"""Observation operator: arrival times at the frozen observation cells, run by run.

`python -m eki extract batches/I00.txt` writes results/I00_G.npz and results/I00_reach.csv.

For every run of the batch and every row of prepared/observations.csv, three rules
are sampled:
- primary: the frozen cell of the primary radius (columns r60_*), direct or proxy;
- sensitivity: the frozen cell of the sensitivity radius (r90_*);
- nearby3x3: 'first nearby arrival', the earliest arrival among the eligible cells of
  the 3 x 3 block around the containing cell (a sensitivity only).

Status codes are in eki.common:
- 0 reached;
- 1 eligible cell not reached by SIMULATION_TSTOP, kept as non-arrival;
- 2 ignition marker;
- 3 not mapped (statically nonburnable, excluded or outside);
- 4 run failed.
arrival_s is NaN unless the status is 0: nothing is dropped, capped, filled or moved to
another cell. Reachability (status) and finite arrival times (arrival_s) are separate
arrays. A batch with runs still running or not started is refused.
"""
import numpy as np
import pandas as pd
import rasterio

from .case import coefficients, template_inputs, template_text
from .common import (IGNITION, NOT_MAPPED, NOT_REACHED, REACHED, RUN_FAILED, STATUS_LEGEND, StudyError,
                     read_json, run_dir, sha256, utc_now, write_text)
from .corrections import theta_names
from .holdout import load_split
from .landscape import eligible
from .observations import load_table
from .run import TOA_NODATA, read_batch, state
from .target import selection


def classify(values):
    """Status and arrival (NaN unless reached) for raster values at assigned cells."""
    values = np.asarray(values, float)
    status = np.full(values.shape, NOT_REACHED, np.int8)
    status[values > 0] = REACHED
    status[values == -1.0] = IGNITION
    bad = ~((values > 0) | (values == -1.0) | (values == TOA_NODATA))
    if bad.any():
        raise StudyError(f"unexpected arrival-time values {np.unique(values[bad])[:5]}")
    return status, np.where(status == REACHED, values, np.nan)


def frozen_cells(obs, prefix):
    """Row/col of the frozen cell per observation, and which observations have one."""
    mapped = obs[prefix + "status"].isin(["direct", "proxy"]).to_numpy()
    rows = np.where(mapped, obs[prefix + "row"].fillna(-1), -1).astype(int)
    cols = np.where(mapped, obs[prefix + "col"].fillna(-1), -1).astype(int)
    return mapped, rows, cols


def nearby_blocks(obs, ok):
    """For each observation, the eligible cells of the 3 x 3 block around its containing cell."""
    blocks = []
    for r, c, inside in zip(obs.row, obs.col, obs.in_domain):
        if not inside:
            blocks.append((np.array([], int), np.array([], int)))
            continue
        r0, r1, c0, c1 = max(r - 1, 0), min(r + 2, ok.shape[0]), max(c - 1, 0), min(c + 2, ok.shape[1])
        rr, cc = np.nonzero(ok[r0:r1, c0:c1])
        blocks.append((rr + r0, cc + c0))
    return blocks


def sample_nearby(toa, blocks):
    status = np.full(len(blocks), NOT_MAPPED, np.int8)
    arrival = np.full(len(blocks), np.nan)
    for i, (rr, cc) in enumerate(blocks):
        if rr.size == 0:
            continue
        v = toa[rr, cc]
        if (v > 0).any():
            status[i], arrival[i] = REACHED, float(v[v > 0].min())
        else:
            status[i] = IGNITION if (v == -1.0).any() else NOT_REACHED
    return status, arrival


def extract(cfg, batch):
    obs = load_table(cfg)
    split = load_split(cfg)
    if list(split.obs_id) != list(obs.obs_id):
        raise StudyError("holdout file and observation table list different observations")
    meta = read_json(cfg.path(cfg["paths"]["prepared"]) / "observations.json")
    prefixes = meta["column_prefixes"]
    ids = read_batch(batch)
    dirs = [run_dir(cfg, i) for i in ids]
    states = [state(d) for d in dirs]
    pending = [i for i, s in zip(ids, states) if s not in ("validated", "failed")]
    if pending:
        raise StudyError(f"{len(pending)} runs are still running or not started (first: {pending[:5]})")
    ok_run = np.array([s == "validated" for s in states])
    names = theta_names(cfg)
    theta = np.array([[coefficients(cfg, r)[0][n] for n in names] for r in ids])
    for i, d in enumerate(dirs):
        if ok_run[i] and not np.allclose(read_json(d / "run.json")["theta"], theta[i], rtol=0, atol=1e-12):
            raise StudyError(f"{ids[i]}: run.json theta differs from its ensemble file")

    stems = template_inputs(template_text(cfg))[0]
    fuels = cfg.shared_dir / "fuels"
    with rasterio.open(fuels / f"{stems['FBFM']}.tif") as f:
        fbfm = f.read(1)
    with rasterio.open(fuels / f"{stems['BLDG_AREA']}.tif") as f:
        ok = eligible(fbfm, f.read(1))
    cells = {rule: frozen_cells(obs, prefixes[rule]) for rule in ("primary", "sensitivity")}
    blocks = nearby_blocks(obs, ok)

    rules = ("primary", "sensitivity", "nearby3x3")
    n, m = len(ids), len(obs)
    status = {r: np.full((n, m), NOT_MAPPED, np.int8) for r in rules}
    arrival = {r: np.full((n, m), np.nan, np.float64) for r in rules}
    done = {}
    for i, d in enumerate(dirs):
        if not ok_run[i]:
            for r in rules:
                status[r][i] = RUN_FAILED
            continue
        done[i] = read_json(d / "DONE")
        with rasterio.open(d / done[i]["toa_file"]) as f:
            toa = f.read(1)
        for rule, (mapped, rows, cols) in cells.items():
            s, a = classify(toa[rows[mapped], cols[mapped]])
            status[rule][i, mapped], arrival[rule][i, mapped] = s, a
        status["nearby3x3"][i], arrival["nearby3x3"][i] = sample_nearby(toa, blocks)

    part = split.partition
    sel = selection(cfg, obs)
    primary = (sel.selection == "primary").to_numpy()
    cand = {k: primary & part.eq(k).to_numpy() for k in ("fit", "holdout")}
    out = cfg.dir("results")
    stem = batch.stem
    np.savez_compressed(
        out / f"{stem}_G.npz", run_id=np.asarray(ids), run_validated=ok_run, obs_id=obs.obs_id.to_numpy(str),
        theta=theta, theta_names=np.asarray(names),
        **{f"arrival_s_{r}": arrival[r] for r in rules}, **{f"status_{r}": status[r] for r in rules},
        model_t_s=obs.model_t_s.to_numpy(float), window_end_s=obs.window_end_s.to_numpy(float),
        in_window=obs.in_window.to_numpy(bool), window_end_after_stop=obs.window_end_after_stop.to_numpy(bool),
        event_type=obs.event_type.to_numpy(str), partition=part.to_numpy(str),
        deadline_s=sel.deadline_s.to_numpy(float), selection=sel.selection.to_numpy(str),
        selection_reason=sel.selection_reason.to_numpy(str),
        candidate_fit=cand["fit"], candidate_holdout=cand["holdout"],
        primary_mapping=obs[prefixes["primary"] + "status"].to_numpy(str),
        sensitivity_mapping=obs[prefixes["sensitivity"] + "status"].to_numpy(str),
        fuel_class=obs.fuel_class.to_numpy(str),
        status_legend=np.asarray([f"{k}: {v}" for k, v in STATUS_LEGEND.items()]),
        observations_sha256=sha256(cfg.path(cfg["paths"]["prepared"]) / "observations.csv"),
        created_utc=utc_now(),
        note="Rows follow run_id, columns follow obs_id (prepared/observations.csv order). arrival_s is "
             "seconds after t0 = 2018-11-08T14:00Z, NaN unless status is 0. deadline_s is the report's "
             "upper arrival bound u (eki.target); candidate_fit marks the primary constraints of the fit "
             "partition. The holdout is stored but not summarized: it is scored only once settings are fixed.")

    prim = obs[prefixes["primary"] + "status"].to_numpy(str)
    later = (~obs.in_window | obs.window_end_after_stop).to_numpy() & np.isin(prim, ["direct", "proxy"]) \
        & part.eq("fit").to_numpy()
    groups = {"fit_direct": cand["fit"] & (prim == "direct"), "fit_proxy": cand["fit"] & (prim == "proxy"),
              "fit_after_window": later}
    rows = []
    for i, rid in enumerate(ids):
        row = {"run_id": rid, "validated": bool(ok_run[i]),
               "burned_area_ac": done[i]["burned_area_ac"] if i in done else np.nan,
               "wall_s": done[i]["wall_s"] if i in done else np.nan}
        for g, sel in groups.items():
            row[f"{g}_n"] = int(sel.sum())
            row[f"{g}_reached"] = int((status["primary"][i, sel] == REACHED).sum())
        for rule in ("sensitivity", "nearby3x3"):
            row[f"fit_{rule}_reached"] = int((status[rule][i, cand["fit"]] == REACHED).sum())
        rows.append(row)
    reach = pd.DataFrame(rows)
    write_text(out / f"{stem}_reach.csv", reach.to_csv(index=False, float_format="%.6g"))
    print(f"wrote {out / f'{stem}_G.npz'} ({n} runs x {m} observations) and {stem}_reach.csv")
    if not ok_run.all():
        print(f"runs that failed (status 4 everywhere): {[r for r, v in zip(ids, ok_run) if not v]}")
    cols = ["run_id", "burned_area_ac", "fit_direct_reached", "fit_direct_n", "fit_proxy_reached", "fit_proxy_n",
            "fit_after_window_reached", "fit_after_window_n"]
    print(reach[cols].to_string(index=False, float_format=lambda v: f"{v:.0f}"))
    members = np.array([r != "BASELINE" for r in ids])
    s = status["primary"][np.ix_(members & ok_run, cand["fit"])]
    print(f"fit constraints ({int(cand['fit'].sum())}): reached by every validated member: "
          f"{int((s == REACHED).all(axis=0).sum())}; by none: {int((s != REACHED).all(axis=0).sum())}")
    return reach
