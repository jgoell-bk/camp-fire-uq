"""The group's QoI for each run: arrival time in the Paradise wildland border band.

Writes results/border_qoi.csv (one row per run) and, next to the run data,
<scratch_root>/border_toa.csv: every border cell's arrival time in each run, with 1-based
row/col in MATLAB's find(mask) order (column-major), as in
../experiment/matlab/extract_border_toa.m. Times are seconds (CSV) or hours (summary)
after t0 = 06:00 PST 8 Nov 2018; the ignition is at 900 s (0.25 h).
"""
import numpy as np
import pandas as pd
import rasterio

from .common import StudyError, read_json, sha256, write_text


def border_cells(cfg):
    with rasterio.open(cfg.path(cfg["qoi"]["border_mask"])) as f:
        mask = f.read(1) == 1
    cols, rows = np.nonzero(mask.T)            # column-major, like MATLAB's find(mask)
    return mask, rows, cols


def same_as(cfg, rid, toa_sha):
    """'identical' / 'different' against qoi.same_as, or '' if none is configured."""
    ref = (cfg["qoi"].get("same_as") or {}).get(rid)
    if not ref:
        return ""
    done = cfg.path(ref) / "DONE"
    if not done.exists():
        return f"reference not validated ({done})"
    return "identical" if read_json(done)["toa_sha256"] == toa_sha else "different"


def summarize(t, horizon_s):
    """Summary of one run's border-band arrival times t (s, NaN = unburned)."""
    burned = np.isfinite(t)
    h = t[burned] / 3600
    cap = np.where(burned, t, horizon_s) / 3600
    return {"cells": t.size, "burned_cells": int(burned.sum()), "burned_share": burned.mean(),
            "earliest_h": h.min() if h.size else np.nan, "median_burned_h": np.median(h) if h.size else np.nan,
            "mean_burned_h": h.mean() if h.size else np.nan, "mean_capped_h": cap.mean()}


def qoi(cfg):
    mask, rows, cols = border_cells(cfg)
    horizon = float(cfg["simulation"]["stop_s"])
    table, cells = [], {"row": rows + 1, "col": cols + 1}
    for rid in cfg.run_ids:
        done_path = cfg.run_dir(rid) / "DONE"
        if not done_path.exists():
            raise StudyError(f"{rid} is not validated yet (python -m baselines status)")
        done = read_json(done_path)
        toa_path = cfg.run_dir(rid) / done["toa_file"]
        if sha256(toa_path) != done["toa_sha256"]:
            raise StudyError(f"{toa_path} changed after validation")
        with rasterio.open(toa_path) as f:
            a = f.read(1)
        if a.shape != mask.shape:
            raise StudyError(f"{rid}: arrival raster {a.shape} and border mask {mask.shape} differ")
        v = a[rows, cols].astype(float)
        if (v == -1).any():
            raise StudyError(f"{rid}: the ignition marker lies in the border band")
        t = np.where(v > 0, v, np.nan)
        cells[rid] = t
        table.append({"run_id": rid, "wind_from": str(cfg["study"]["runs"][rid]), **summarize(t, horizon),
                      "burned_area_ac": done["burned_area_ac"], "same_inputs_run": same_as(cfg, rid, done["toa_sha256"])})
    out = pd.DataFrame(table)
    path = cfg.dir("results") / "border_qoi.csv"
    write_text(path, out.to_csv(index=False, float_format="%.6g"))
    cell_path = cfg.scratch_root / "border_toa.csv"
    write_text(cell_path, pd.DataFrame(cells).to_csv(index=False, float_format="%.1f", na_rep="NaN"))
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(out.drop(columns=["wind_from"]).to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print(f"hours after 06:00 PST 8 Nov; ignition at 0.25 h; unburned cells count as "
          f"{horizon / 3600:.0f} h in mean_capped_h")
    print(f"wrote {path}\nwrote {cell_path} (per-cell seconds, NaN = unburned)")
    return out
