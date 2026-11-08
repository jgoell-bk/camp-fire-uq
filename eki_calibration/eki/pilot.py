"""Pilot diagnostics for an extracted batch, on the fit constraints only.

`python -m eki pilot batches/I00.txt` reads results/I00_G.npz and writes:
- results/I00_pilot.json: per-run reachability and deadline violations; cells reached by
  all/some/no members; breakdowns; parameter-response relationships and conditioning;
  the distance response;
- results/I00_pilot_cells.csv: one row per fit constraint;
- results/I00_distance.npz: distance to the burned area at each deadline.

It never builds an update and never calls fit_arrays, so it works when the update is
refused. The holdout is not touched. Lateness is T - u, the modelled arrival minus the
report deadline (positive = violation). It is a deadline violation, not an arrival-time
error.

Distance response (evaluated, not adopted): d_i is the distance from the constraint's
sampling-cell centre to the nearest cell burned by u_i (0 < T <= u_i, or the ignition
cell). It is 0 when the cell itself burned by its deadline. It is defined for every
member, including where the cell is never reached. Its limits:
- it is in metres on a 30 m grid;
- it ignores wind, slope and fuel between the cell and the fire;
- it jumps when a spot fire lands nearby;
- it needs a declared length scale before it can sit beside lateness in one objective.
"""
import numpy as np
import pandas as pd
import rasterio

from .common import BASELINE, NOT_REACHED, REACHED, StudyError, read_json, run_dir, utc_now, write_json, write_text
from .holdout import load_split
from .observations import load_table

RADII_CELLS = (10, 34, 100, 334, 1500)     # search half-widths: 0.3, 1, 3, 10, 45 km
HOUR_BINS = ((0, 3), (3, 6), (6, 9), (9, 12), (12, 18), (18, 24.01))
KM_BINS = ((0, 10), (10, 13), (13, 16), (16, 20), (20, 50))


def distance_to_burned(toa, rows, cols, u, cell_m):
    """Distance (m) from each cell to the nearest cell burned by its deadline u."""
    H, W = toa.shape
    out = np.full(len(rows), np.nan)
    for i, (r, c, ui) in enumerate(zip(rows, cols, u)):
        if 0 < toa[r, c] <= ui:
            out[i] = 0.0
            continue
        for h in RADII_CELLS:
            r0, r1, c0, c1 = max(r - h, 0), min(r + h + 1, H), max(c - h, 0), min(c + h + 1, W)
            w = toa[r0:r1, c0:c1]
            rr, cc = np.nonzero(((w > 0) & (w <= ui)) | (w == -1.0))
            if rr.size:
                d = float(np.hypot(rr + r0 - r, cc + c0 - c).min())
                if d <= h or h == RADII_CELLS[-1]:     # nothing outside the square can be closer
                    out[i] = d * cell_m
                    break
    return out


def _q(x, qs=(0.1, 0.5, 0.9)):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return [round(float(v), 2) for v in np.quantile(x, qs)] if x.size else None


def _corr(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.std() == 0 or b.std() == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def pilot(cfg, batch, lateness_scale_s=3600.0):
    out = cfg.dir("results")
    path = out / f"{batch.stem}_G.npz"
    if not path.exists():
        raise StudyError(f"{path} not found; run `python -m eki extract {batch}` first")
    z = np.load(path)
    obs, split = load_table(cfg), load_split(cfg)
    fit = z["candidate_fit"]
    if list(z["obs_id"]) != list(obs.obs_id):
        raise StudyError("extract file and observation table list different observations")
    ids, ok = z["run_id"], z["run_validated"]
    if not ok.all():
        raise StudyError(f"runs failed: {list(ids[~ok])}; retry them before the pilot review")
    mem = ids != BASELINE
    st, T, u = z["status_primary"][:, fit], z["arrival_s_primary"][:, fit], z["deadline_s"][fit]
    o, loc = obs[fit].reset_index(drop=True), split[fit].reset_index(drop=True)
    N = int(fit.sum())
    reached = st == REACHED
    late_h = (T - u) / 3600.0                               # NaN where not reached
    viol = reached & (late_h > 0)
    theta, names = z["theta"], list(z["theta_names"])

    # ---- per run
    runs = []
    for j, rid in enumerate(ids):
        lv = late_h[j][viol[j]]
        runs.append({"run_id": str(rid), "reached": int(reached[j].sum()), "not_reached": int((st[j] == NOT_REACHED).sum()),
                     "violated": int(viol[j].sum()), "met": int((reached[j] & ~viol[j]).sum()),
                     "lateness_h_p10_p50_p90_reached": _q(late_h[j][reached[j]]),
                     "lateness_h_p50_p90_violated": _q(lv, (0.5, 0.9))})
    nmem = reached[mem].sum(axis=0)
    cells = {"members": int(mem.sum()), "reached_by_all": int((nmem == mem.sum()).sum()),
             "reached_by_some": int(((nmem > 0) & (nmem < mem.sum())).sum()), "reached_by_none": int((nmem == 0).sum()),
             "baseline_reached": int(reached[~mem].sum())}

    # ---- breakdowns (members pooled)
    grid = read_json(cfg.shared_dir / "manifest.json")
    ign = grid["ignition"]
    km = np.hypot(loc.cell_x - ign["cell_x"], loc.cell_y - ign["cell_y"]) / 1000
    u_h = u / 3600
    labels = {
        "deadline_h": pd.Series([next((f"{a}-{b:g}" for a, b in HOUR_BINS if a <= v < b), "other") for v in u_h]),
        "km_from_ignition": pd.Series([next((f"{a}-{b}" for a, b in KM_BINS if a <= v < b), "other") for v in km]),
        "nist_region": o.nist_sheet, "fire_type": o.nist_fire_type,
        "mapping": o["r{:g}_status".format(float(cfg["observations"]["proxy_radius_m"]))]}
    breakdown = {}
    for name, lab in labels.items():
        rows = {}
        for g in sorted(lab.unique()):
            k = (lab == g).to_numpy()
            r_, v_ = reached[np.ix_(mem, k)], viol[np.ix_(mem, k)]
            lv = late_h[np.ix_(mem, k)][v_]
            rows[str(g)] = {"cells": int(k.sum()), "reached_share": round(float(r_.mean()), 3),
                            "violated_share_of_reached": round(float(v_.sum() / max(r_.sum(), 1)), 3),
                            "lateness_h_p50_violated": _q(lv, (0.5,))[0] if lv.size else None,
                            "reached_by_no_member": int((nmem[k] == 0).sum())}
        breakdown[name] = rows

    # ---- distance response
    cell_m = float(grid["fuel_grid"]["cellsize"])
    p = "r{:g}_".format(float(cfg["observations"]["proxy_radius_m"]))
    rr, cc = o[p + "row"].astype(int).to_numpy(), o[p + "col"].astype(int).to_numpy()
    dist = np.empty((len(ids), N))
    for j, rid in enumerate(ids):
        d = run_dir(cfg, str(rid))
        with rasterio.open(d / read_json(d / "DONE")["toa_file"]) as f:
            dist[j] = distance_to_burned(f.read(1), rr, cc, u, cell_m)
    np.savez_compressed(out / f"{batch.stem}_distance.npz", run_id=ids, obs_id=o.obs_id.to_numpy(str),
                        distance_m=dist, deadline_s=u,
                        note="distance from the sampling-cell centre to the nearest cell burned by the deadline "
                             "(0 < T <= u, or the ignition cell); 0 = cell burned by its deadline")
    both = reached & np.isfinite(dist)
    distance = {"zero_share_by_run": {str(r): round(float((dist[j] == 0).mean()), 3) for j, r in enumerate(ids)},
                "positive_m_p10_p50_p90_members": _q(dist[mem][dist[mem] > 0]),
                "unreached_cells_m_p10_p50_p90_members": _q(dist[mem][st[mem] == NOT_REACHED]),
                "corr_with_lateness_where_violated": _corr(dist[viol], late_h[viol]) if viol.sum() > 2 else None,
                "met_implies_zero": bool((dist[both & ~viol] == 0).all())}

    # ---- parameter-response (members only) and conditioning
    th = theta[mem]
    zt = (th - th.mean(0)) / th.std(0, ddof=1)
    common = reached[mem].all(axis=0)
    resp = {"reached": reached[mem].sum(axis=1).astype(float),
            "violated_share_of_reached": viol[mem].sum(axis=1) / np.maximum(reached[mem].sum(axis=1), 1),
            "mean_lateness_h_common_cells": np.nanmean(late_h[mem][:, common], axis=1) if common.any() else None,
            "mean_distance_km": dist[mem].mean(axis=1) / 1000,
            "burned_by_deadline_share": (dist[mem] == 0).mean(axis=1)}
    pr = {k: {n: round(_corr(zt[:, i], v), 2) for i, n in enumerate(names)} for k, v in resp.items() if v is not None}
    X = np.hstack([np.ones((len(zt), 1)), zt])
    r2 = {}
    for k, v in resp.items():
        if v is not None:
            beta, *_ = np.linalg.lstsq(X, v, rcond=None)
            r2[k] = round(float(1 - ((v - X @ beta) ** 2).sum() / ((v - v.mean()) ** 2).sum()), 3)
    per_cell = {}
    for label, M in (("arrival_common_cells", T[mem][:, common]), ("distance_all_cells", dist[mem])):
        if M.shape[1] == 0:
            continue
        sd = M.std(axis=0)
        keep = sd > 0
        C = np.array([[_corr(zt[:, i], M[:, k]) for k in np.where(keep)[0]] for i in range(len(names))])
        per_cell[label] = {"cells_varying": int(keep.sum()),
                           **{n: {"median_abs_corr": round(float(np.nanmedian(np.abs(C[i]))), 2),
                                  "share_abs_corr_over_0.35": round(float(np.nanmean(np.abs(C[i]) > 0.35)), 2)}
                              for i, n in enumerate(names)}}
    s_th = np.linalg.svd(zt, compute_uv=False)
    cond = {"theta_singular_values": [round(float(v), 3) for v in s_th],
            "theta_condition": round(float(s_th[0] / s_th[-1]), 2)}
    for label, M in (("arrival_common_cells", T[mem][:, common] / 3600), ("distance_km", dist[mem] / 1000)):
        if M.shape[1] > 1:
            s = np.linalg.svd(M - M.mean(0), compute_uv=False)
            cond[f"{label}_singular_values_top"] = [round(float(v), 3) for v in s[:10]]
            cond[f"{label}_effective_rank_1pct"] = int((s > 0.01 * s[0]).sum())
    if common.any():
        r = np.maximum(0, (T[mem][:, common] - u[common]) / lateness_scale_s) / np.sqrt(common.sum())
        dF = r - r.mean(0)
        ev = np.linalg.eigvalsh(dF.T @ dF / (mem.sum() - 1) + np.eye(common.sum()))
        cond["deadline_CFF_plus_I_condition_common_cells"] = round(float(ev[-1] / ev[0]), 2)
        cond["deadline_CFF_trace_common_cells"] = round(float(np.trace(dF.T @ dF / (mem.sum() - 1))), 4)
        cond["deadline_scale_used_s"] = lateness_scale_s

    # ---- per-cell table and json
    tab = pd.DataFrame({"obs_id": o.obs_id, "nist_region": o.nist_sheet, "fire_type": o.nist_fire_type,
                        "mapping": labels["mapping"], "deadline_h": u_h, "km_from_ignition": km,
                        "members_reached": nmem, "baseline_reached": reached[~mem].sum(axis=0),
                        "members_violated": viol[mem].sum(axis=0),
                        "lateness_h_median_members": np.nanmedian(np.where(reached[mem], late_h[mem], np.nan), axis=0),
                        "baseline_lateness_h": late_h[~mem][0] if (~mem).any() else np.nan,
                        "distance_m_median_members": np.median(dist[mem], axis=0),
                        "baseline_distance_m": dist[~mem][0] if (~mem).any() else np.nan})
    write_text(out / f"{batch.stem}_pilot_cells.csv", tab.to_csv(index=False, float_format="%.6g"))
    result = {"created_utc": utc_now(), "batch": str(batch), "fit_constraints": N, "runs": runs, "cells": cells,
              "breakdown": breakdown, "distance_response": distance, "parameter_response_corr": pr,
              "parameter_response_r2": r2, "per_cell_response": per_cell, "conditioning": cond,
              "notes": ["lateness = modelled arrival - report deadline (h); positive = deadline violated",
                        "correlations over 32 members: |r| > 0.35 is about p < 0.05",
                        "common cells = reached by every member; used for diagnostics only, never as the target"]}
    write_json(out / f"{batch.stem}_pilot.json", result)
    return result, tab
