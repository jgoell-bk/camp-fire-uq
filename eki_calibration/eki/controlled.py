"""The first controlled update: propose -> check runs at the means -> verify -> promote.

propose I00 (from results/I00_distance.npz, written by `pilot`):
- Builds the step for the fit constraints: G = sqrt(w) d / L, B from
  eki.objective, alpha by eki.update.controlled_step.
- Writes proposals/I01.csv and .json: alpha, the displacements, the objective parts of
  every I00 member, and early/late response diagnostics.
- Writes ensembles/CHK01.csv with three check runs:
  - CHK01_I00MEAN: the current ensemble mean;
  - CHK01_I00MEAN_REP: the same inputs and seed, a repeat that measures run-to-run
    variability;
  - CHK01_I01MEAN: the proposed ensemble mean.
  It also writes batches/CHK01.txt and prints its sbatch line.
verify CHK01: actual Phi of the runs at the two means (never the ensemble average of
  the member responses), the repeat difference, and diagnostics by deadline window,
  fire type and direct/proxy. It writes results/CHK01_verify.json.
promote I01: copies the proposal to ensembles/I01.csv only if verify found it no worse
  than the measured variability. Then make-batch I01. One step at a time: nothing
  iterates on its own.
The holdout is never read for any of this beyond its block labels for fit constraints.
"""
import shutil

import numpy as np
import pandas as pd
import rasterio

from .common import BASELINE, REACHED, StudyError, read_json, run_dir, sha256, utc_now, write_json, write_text
from .corrections import corrected_weather, theta_names
from .forcing import load_baseline
from .holdout import load_split
from .objective import block_weights, distance_outputs, objective, overspread_outputs, regularization_matrix
from .observations import load_table
from .pilot import distance_to_burned
from .run import state
from .target import candidates
from .update import controlled_step, write_ensemble

WINDOWS = ((0, 6), (6, 12), (12, 24.01))


def fit_constraints(cfg):
    """The fixed fit constraints (obs order), their blocks and weights."""
    obs, split = load_table(cfg), load_split(cfg)
    c = candidates(cfg, obs, split.partition)
    blocks = split.block[c.index].to_numpy(str)
    return c, blocks, block_weights(blocks)


def admissible(cfg):
    baseline, names = load_baseline(cfg), theta_names(cfg)

    def ok(row):
        try:
            corrected_weather(cfg, baseline, dict(zip(names, row)))
            return True
        except StudyError:
            return False
    return ok


def run_responses(cfg, rid, c):
    """Distance (m), arrival (s, NaN unless reached) and status at the fit constraints."""
    d = run_dir(cfg, rid)
    if state(d) != "validated":
        raise StudyError(f"{rid} is {state(d)}, not validated")
    with rasterio.open(d / read_json(d / "DONE")["toa_file"]) as f:
        toa = f.read(1)
    p = "r{:g}_".format(float(cfg["observations"]["proxy_radius_m"]))
    rr, cc = c[p + "row"].astype(int).to_numpy(), c[p + "col"].astype(int).to_numpy()
    cell = float(read_json(cfg.shared_dir / "manifest.json")["fuel_grid"]["cellsize"])
    v = toa[rr, cc]
    return distance_to_burned(toa, rr, cc, c.deadline_s.to_numpy(), cell), np.where(v > 0, v, np.nan), \
        np.where(v > 0, REACHED, 1)


def breakdown(cfg, c, w, d, T):
    """Per group: cells, deadline met share (d = 0), mean and median distance (km,
    unweighted), weighted distance term, median standardized lateness among reached."""
    L, s = float(cfg["fit"]["length_scale_m"]), float(cfg["target"]["lateness_scale_s"])
    u_h = c.deadline_s.to_numpy() / 3600
    p = "r{:g}_".format(float(cfg["observations"]["proxy_radius_m"]))
    groups = {"deadline_h": np.array([next(f"{a}-{b:g}" for a, b in WINDOWS if a <= x < b) for x in u_h]),
              "fire_type": c.nist_fire_type.to_numpy(str), "mapping": c[p + "status"].to_numpy(str)}
    out = {}
    for name, lab in groups.items():
        out[name] = {}
        for g in sorted(set(lab)):
            k = lab == g
            late = (T[k] - c.deadline_s.to_numpy()[k]) / s
            out[name][g] = {"cells": int(k.sum()), "met_share": round(float((d[k] == 0).mean()), 3),
                            "mean_km": round(float(d[k].mean() / 1000), 2),
                            "median_km": round(float(np.median(d[k]) / 1000), 2),
                            "weighted_distance_term": round(float(0.5 * (w[k] * (d[k] / L) ** 2).sum()), 3),
                            "median_std_lateness_reached": round(float(np.nanmedian(late)), 2)
                            if np.isfinite(late).any() else None}
    return out


def early_late_response(cfg, theta, d, c):
    """Standardized regression of the mean distance (km) in each deadline window on all
    nine coefficients over the members. A knot cannot affect weather before its basis
    turns on, so an early-window weight on a later knot reflects correlated
    initialization or chance, not control."""
    names = theta_names(cfg)
    zt = (theta - theta.mean(0)) / theta.std(0, ddof=1)
    X = np.hstack([np.ones((len(zt), 1)), zt])
    u_h = c.deadline_s.to_numpy() / 3600
    out = {}
    for a, b in WINDOWS:
        k = (u_h >= a) & (u_h < b)
        y = d[:, k].mean(axis=1) / 1000
        beta, *_ = np.linalg.lstsq(X, (y - y.mean()) / y.std(ddof=1), rcond=None)
        out[f"{a}-{b:g}h"] = {"cells": int(k.sum()), **{n: round(float(v), 2) for n, v in zip(names, beta[1:])}}
    knots = cfg["corrections"]["knots_h"]
    support = {f"k{r + 1}": [knots[r - 1] if r > 0 else 0.0, knots[r + 1] if r + 1 < len(knots) else
                             cfg["simulation"]["stop_s"] / 3600] for r in range(len(knots))}
    return {"standardized_betas_mean_distance": out, "basis_support_h": support,
            "note": "k1 acts on 0-6 h, k2 on 0-12 h, k3 from 6 h on (model hours after 06:00 PST); "
                    "k3 cannot change any deadline before 6 h"}


def propose(cfg, iteration=0):
    src, nxt = f"I{iteration:02d}", f"I{iteration + 1:02d}"
    chk = f"CHK{iteration + 1:02d}"
    res = cfg.dir("results")
    path = res / f"{src}_distance.npz"
    if not path.exists():
        raise StudyError(f"{path} not found; run `python -m eki pilot batches/{src}.txt` first")
    c, blocks, w = fit_constraints(cfg)
    D = np.load(path)
    if list(D["obs_id"]) != list(c.obs_id) or not np.array_equal(D["deadline_s"], c.deadline_s.to_numpy()):
        raise StudyError(f"{path.name} was made for other constraints or deadlines; rerun pilot")
    ids = list(D["run_id"])
    ens = pd.read_csv(cfg.path(cfg["paths"]["ensembles"]) / f"{src}.csv")
    rows = [ids.index(r) for r in ens.run_id]
    names = theta_names(cfg)
    theta, d = ens[names].to_numpy(float), D["distance_m"][rows]
    L = float(cfg["fit"]["length_scale_m"])
    G, B = distance_outputs(d, w, L), regularization_matrix(cfg)
    over = cfg["fit"].get("overspread")
    acres = base_acres = None
    if over:                                  # v3: one more output, A_out / reference_acres
        if "overspread_acres" not in D.files or float(D["overspread_at_s"]) != float(over["at_s"]):
            raise StudyError(f"{path.name} has no over-spread areas for at_s = {over['at_s']}; rerun pilot")
        acres = D["overspread_acres"][rows]
        G = np.hstack([G, overspread_outputs(cfg, acres)])
        base_acres = D["overspread_acres"][[ids.index(BASELINE)]] if BASELINE in ids else None
    ok = admissible(cfg)
    new, step = controlled_step(cfg, theta, G, B, ok)
    obj = objective(cfg, theta, d, w, acres)
    base = objective(cfg, np.zeros((1, len(names))), D["distance_m"][ids.index(BASELINE)][None], w, base_acres) \
        if BASELINE in ids else None

    out = cfg.path("proposals") / f"{nxt}.csv"
    m0, m1 = theta.mean(axis=0), new.mean(axis=0)
    for label, m in (("current mean", m0), ("proposed mean", m1)):
        if not ok(m):
            raise StudyError(f"the {label} is not admissible")
    _, nb = np.unique(blocks, return_counts=True)
    parent = {"from": src, "distance_file": str(path), "distance_sha256": sha256(path),
              "constraints": int(len(c)), "blocks": int(len(nb)), "constraints_per_block_min_median_max":
                  [int(nb.min()), float(np.median(nb)), int(nb.max())],
              "weights_sum": float(w.sum()), "length_scale_m": L, "fit_settings": cfg["fit"],
              "B_rows": int(len(B)), "step": step,
              "current_objective": {k: {"mean": float(v.mean()), "min": float(v.min()), "max": float(v.max())}
                                    for k, v in obj.items()},
              "current_objective_per_member": {r: {k: float(v[j]) for k, v in obj.items()}
                                               for j, r in enumerate(ens.run_id)},
              "baseline_objective": {k: float(v[0]) for k, v in base.items()} if base else None,
              "current_mean": dict(zip(names, map(float, m0))), "proposed_mean": dict(zip(names, map(float, m1))),
              "early_late_response": early_late_response(cfg, theta, d, c),
              "note": "proposal only: not an ensemble file until verify accepts it and promote copies it"}
    write_ensemble(cfg, out, iteration + 1, new, parent)
    checks = pd.DataFrame({"run_id": [f"{chk}_{src}MEAN", f"{chk}_{src}MEAN_REP", f"{chk}_{nxt}MEAN"],
                           **{n: [m0[i], m0[i], m1[i]] for i, n in enumerate(names)}})
    cpath = cfg.path(cfg["paths"]["ensembles"]) / f"{chk}.csv"
    text = checks.to_csv(index=False, float_format="%.17g")
    if cpath.exists() and cpath.read_text() != text:
        raise StudyError(f"{cpath} exists with different contents")
    write_text(cpath, text)
    write_json(cpath.with_suffix(".json"), {"created_utc": utc_now(), "sha256": sha256(cpath),
                                            "proposal": str(out), "proposal_sha256": sha256(out),
                                            "rows": {f"{chk}_{src}MEAN": f"mean of {src}",
                                                     f"{chk}_{src}MEAN_REP": "identical inputs and seed (variability)",
                                                     f"{chk}_{nxt}MEAN": f"mean of the proposed {nxt}"}})
    return out, cpath, parent


def verify(cfg, chk):
    cpath = cfg.path(cfg["paths"]["ensembles"]) / f"{chk}.csv"
    meta = read_json(cpath.with_suffix(".json"))
    if sha256(meta["proposal"]) != meta["proposal_sha256"]:
        raise StudyError("the proposal changed after the check runs were defined")
    rows = pd.read_csv(cpath)
    names = theta_names(cfg)
    c, blocks, w = fit_constraints(cfg)
    resp = {r: run_responses(cfg, r, c) for r in rows.run_id}
    theta = rows[names].to_numpy(float)
    d = np.array([resp[r][0] for r in rows.run_id])
    acres = None
    if cfg["fit"].get("overspread"):
        from .overspread import perimeter_mask, run_acres
        mask = perimeter_mask(cfg)
        acres = np.array([run_acres(cfg, r, mask) for r in rows.run_id])
    obj = objective(cfg, theta, d, w, acres)
    phi = dict(zip(rows.run_id, obj["total"]))
    cur, rep, prop = rows.run_id
    toa_same = read_json(run_dir(cfg, cur) / "DONE")["toa_sha256"] == read_json(run_dir(cfg, rep) / "DONE")["toa_sha256"]
    variability = abs(phi[rep] - phi[cur])
    change = phi[prop] - phi[cur]
    verdict = "improves" if change < -variability else ("worsens" if change > variability else "no measurable change")
    result = {"created_utc": utc_now(), "check": chk, "proposal": meta["proposal"],
              "proposal_sha256": meta["proposal_sha256"],
              "objective": {r: {k: float(v[j]) for k, v in obj.items()} for j, r in enumerate(rows.run_id)},
              **({"overspread_acres": dict(zip(rows.run_id, map(float, acres)))} if acres is not None else {}),
              "repeat_bit_identical": bool(toa_same), "variability": float(variability),
              "change_proposed_minus_current": float(change), "verdict": verdict,
              "accepted": bool(change <= variability),
              "diagnostics": {r: breakdown(cfg, c, w, resp[r][0], resp[r][1]) for r in rows.run_id},
              "note": "Phi from the runs at the mean parameters; the ensemble average of member responses is "
                      "not used. Step limits do not promise a decrease."}
    write_json(cfg.dir("results") / f"{chk}_verify.json", result)
    return result


def promote(cfg, iteration):
    nxt, chk = f"I{iteration:02d}", f"CHK{iteration:02d}"
    v = read_json(cfg.dir("results") / f"{chk}_verify.json")
    prop = cfg.path("proposals") / f"{nxt}.csv"
    if sha256(prop) != v["proposal_sha256"]:
        raise StudyError(f"{prop} is not the proposal that {chk} verified")
    if not v["accepted"]:
        raise StudyError(f"{chk}: the proposed mean {v['verdict']} the objective "
                         f"({v['change_proposed_minus_current']:+.4f}); increase damping and reconsider")
    dest = cfg.path(cfg["paths"]["ensembles"]) / f"{nxt}.csv"
    if dest.exists():
        raise StudyError(f"{dest} exists; ensemble files are never edited")
    shutil.copyfile(prop, dest)
    write_json(dest.with_suffix(".json"), {"created_utc": utc_now(), "file": dest.name, "sha256": sha256(dest),
                                           "promoted_from": str(prop), "verified_by": f"results/{chk}_verify.json",
                                           "proposal_record": read_json(prop.with_suffix(".json"))})
    return dest
