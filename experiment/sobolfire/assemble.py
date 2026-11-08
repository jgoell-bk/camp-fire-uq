"""Collect run outputs into per-station data files, in Sobol order (handoff section 8).

Rows are always design_row 1..1536 of one station, whatever order the runs
finished in. Stations are kept separate: station is a scenario, not a factor.
This only gathers data; choosing a quantity of interest and computing Sobol
indices is left to the analysis.
"""
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import xy

from .case import URBAN, nonburnable, template_inputs, template_text
from .design import derived_settings, load_matrix
from .runner import final_dir
from .util import ConfigError, read_json, write_text


def done_records(cfg, prefix):
    runs = [r for r in load_matrix(cfg) if r.prefix == prefix]
    if [r.design_row for r in runs] != list(range(1, len(runs) + 1)):
        raise ConfigError(f"{prefix}: matrix rows are not in design_row order")
    dirs = [final_dir(cfg, r) for r in runs]
    missing = [r.run_id for r, d in zip(runs, dirs) if not (d / "DONE").exists()]
    if missing:
        raise ConfigError(f"{prefix}: {len(missing)} of {len(runs)} runs are not validated "
                          f"(first: {missing[:5]}); every response is needed")
    return runs, dirs, [read_json(d / "DONE") for d in dirs]


def assemble(cfg, prefix, mask_path=None):
    """Write results/<prefix>_runs.csv (one row per run: inputs, burned area, path to
    its arrival-time raster). With a mask raster on the fuel grid, also write
    results/<prefix>_arrival.npz: raw arrival times at the masked cells."""
    runs, dirs, done = done_records(cfg, prefix)
    out = cfg.dir("results")
    table = pd.DataFrame({
        "design_row": [r.design_row for r in runs], "run_id": [r.run_id for r in runs],
        "base_row": [r.base_row for r in runs], "sobol_block": [r.sobol_block for r in runs],
        **{f: [r.coeffs[f] for r in runs] for f in cfg.factors},
        **pd.DataFrame([derived_settings(cfg, r.coeffs) for r in runs]).to_dict("list"),
        **{k: [d[k] for d in done] for k in ("burned_area_ac", "crown_area_ac", "burned_cells",
                                             "tstop_h", "wall_s")},
        "stalled": [any("STALLED" in m for m in d["stop_messages"]) for d in done],
        "toa_path": [str(dir_ / d["toa_file"]) for dir_, d in zip(dirs, done)]})
    write_text(out / f"{prefix}_runs.csv", table.to_csv(index=False, float_format="%.17g"))
    print(f"wrote {out / f'{prefix}_runs.csv'} ({len(table)} rows in design_row order)")
    if mask_path is None:
        return table

    fuel = read_json(cfg.shared_dir / "manifest.json")["fuel_grid"]
    with rasterio.open(mask_path) as m:
        if (m.height, m.width) != (fuel["nrows"], fuel["ncols"]) or abs(m.transform.c - fuel["x0"]) > 1 \
                or abs(m.transform.f - fuel["top"]) > 1:
            raise ConfigError(f"{mask_path} is not on the fuel grid")
        rows, cols = np.nonzero(m.read(1, masked=True).filled(0) > 0)
        x, y = xy(m.transform, rows, cols)
    fbfm_stem = template_inputs(template_text(cfg))[0]["FBFM"]
    with rasterio.open(cfg.shared_dir / "fuels" / f"{fbfm_stem}.tif") as f:
        fbfm = f.read(1)[rows, cols]
    arrival = np.empty((len(runs), len(rows)), np.float32)
    for i, path in enumerate(table.toa_path):
        with rasterio.open(path) as f:
            arrival[i] = f.read(1)[rows, cols]
    arrival[arrival <= 0] = np.nan          # ELMFIRE's -9999 / -1: not burned by the end
    path = out / f"{prefix}_arrival.npz"
    np.savez_compressed(
        path, arrival_s=arrival, design_row=table.design_row.to_numpy(),
        run_id=table.run_id.to_numpy(str), row=rows, col=cols, x=np.asarray(x), y=np.asarray(y),
        fbfm=fbfm, unburnable=nonburnable(fbfm) & (fbfm != URBAN),
        stop_s=float(cfg["simulation"]["stop_s"]), mask=str(mask_path),
        note="arrival_s: seconds after t0 (2018-11-08T14:00Z); NaN = not burned by stop_s. "
             "unburnable: nonburnable fuel other than urban (91), which burns via buildings.")
    print(f"wrote {path}: {len(runs)} runs x {len(rows)} cells")
    return table
