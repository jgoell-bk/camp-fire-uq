"""Collect one station's run outputs into results/, one row per run in offset order
(-50, -49, ..., +50). This only gathers data; the analysis is done elsewhere."""
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import xy

from .case import URBAN, nonburnable, template_inputs, template_text
from .common import StudyError, read_json, run_dir, run_id, write_text


def collect(cfg, prefix, mask_path=None):
    """Write results/<S>_runs.csv (offset, wind direction, burned area, path to the
    arrival-time raster; rotation_deg is empty for fixed-direction runs). With a mask raster on the fuel grid, also write
    results/<S>_arrival.npz: arrival times at the masked cells, one row per run."""
    if prefix not in cfg.stations:
        raise StudyError(f"unknown station {prefix!r}; one of {list(cfg.stations)}")
    ids = [run_id(prefix, d) for d in cfg.offsets]
    dirs = [run_dir(cfg, i) for i in ids]
    missing = [i for i, d in zip(ids, dirs) if not (d / "DONE").exists()]
    if missing:
        raise StudyError(f"{prefix}: {len(missing)} of {len(ids)} runs are not validated (first: {missing[:5]})")
    done = [read_json(d / "DONE") for d in dirs]
    rec = [read_json(d / "run.json") for d in dirs]
    # Runs built before study.direction existed (the first 303, 2026-10-26) are hourly.
    directions = {r.get("direction", "hourly") for r in rec}
    if directions != {cfg["study"]["direction"]}:
        raise StudyError(f"{prefix}: runs have direction {sorted(directions)}, the config says "
                         f"{cfg['study']['direction']!r}")
    table = pd.DataFrame({
        "run_id": ids, "offset_deg": cfg.offsets, "direction": cfg["study"]["direction"],
        # the run's mean wind direction; for fixed runs, the direction of every hour
        "wd_from_deg": [(r["target_from_deg"] + d) % 360 for r, d in zip(rec, cfg.offsets)],
        **{k: [r.get(k) for r in rec] for k in ("aim_deg", "rotation_deg", "mean_from_deg", "target_from_deg")},
        **{k: [d[k] for d in done] for k in ("burned_area_ac", "crown_area_ac", "burned_cells",
                                             "first_arrival_s", "tstop_h", "wall_s")},
        "stalled": [any("STALLED" in m for m in d["stop_messages"]) for d in done],
        "toa_path": [str(dir_ / d["toa_file"]) for dir_, d in zip(dirs, done)]})
    out = cfg.dir("results")
    write_text(out / f"{prefix}_runs.csv", table.to_csv(index=False, float_format="%.17g"))
    print(f"wrote {out / f'{prefix}_runs.csv'} ({len(table)} runs, offsets {ids[0]} .. {ids[-1]})")
    if mask_path is None:
        return table

    fuel = read_json(cfg.shared_dir / "manifest.json")["fuel_grid"]
    with rasterio.open(mask_path) as m:
        if (m.height, m.width) != (fuel["nrows"], fuel["ncols"]) or abs(m.transform.c - fuel["x0"]) > 1 \
                or abs(m.transform.f - fuel["top"]) > 1:
            raise StudyError(f"{mask_path} is not on the fuel grid")
        rows, cols = np.nonzero(m.read(1, masked=True).filled(0) > 0)
        x, y = xy(m.transform, rows, cols)
    fbfm_stem = template_inputs(template_text(cfg))[0]["FBFM"]
    with rasterio.open(cfg.shared_dir / "fuels" / f"{fbfm_stem}.tif") as f:
        fbfm = f.read(1)[rows, cols]
    arrival = np.empty((len(ids), len(rows)), np.float32)
    for i, path in enumerate(table.toa_path):
        with rasterio.open(path) as f:
            arrival[i] = f.read(1)[rows, cols]
    arrival[arrival <= 0] = np.nan          # ELMFIRE's -9999 / -1: not burned by the end
    path = out / f"{prefix}_arrival.npz"
    np.savez_compressed(
        path, arrival_s=arrival, offset_deg=np.asarray(cfg.offsets), run_id=np.asarray(ids),
        row=rows, col=cols, x=np.asarray(x), y=np.asarray(y),
        fbfm=fbfm, unburnable=nonburnable(fbfm) & (fbfm != URBAN),
        stop_s=float(cfg["simulation"]["stop_s"]), mask=str(mask_path),
        note="arrival_s: seconds after t0 (2018-11-08T14:00Z); NaN = not burned by stop_s. "
             "Rows follow offset_deg. unburnable: nonburnable fuel other than urban (91).")
    print(f"wrote {path}: {len(ids)} runs x {len(rows)} cells")
    return table
