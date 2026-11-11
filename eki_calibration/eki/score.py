"""Score runs on the fit and holdout constraints: python -m eki score NAME RUN ... [--members I08].

For each run and each partition (fit, holdout) this gives:
- the weighted distance term 1/2 sum w_i (d_i / L)^2, where w_i = 1/(K n_b) over that
  partition's blocks, so the weights sum to 1 in both partitions;
- the unweighted share of deadlines met and the mean distance;
- the breakdown by deadline window, fire type and direct/proxy, as in verify.
Distances come from each run's arrival raster (eki.pilot.distance_to_burned). Nothing
is fitted or changed. Output: results/score_<NAME>.json.

The holdout was kept unscored during the calibration. It is scored once, after the
calibration has finished, and never used to choose settings (Jacob, 2026-11-10).
"""
import numpy as np
import pandas as pd

from .common import StudyError, utc_now, write_json
from .controlled import breakdown, run_responses
from .holdout import load_split
from .objective import block_weights
from .observations import load_table
from .target import candidates


def partition_constraints(cfg, obs, split, which):
    c = candidates(cfg, obs, split.partition, which=which)
    blocks = split.block[c.index].to_numpy(str)
    return c, blocks, block_weights(blocks)


def one(cfg, rid, c, w):
    L = float(cfg["fit"]["length_scale_m"])
    d, T, _ = run_responses(cfg, rid, c)
    return {"distance_term": float(0.5 * (w * (d / L) ** 2).sum()), "met_share": float((d == 0).mean()),
            "mean_km": float(d.mean() / 1000)}, d, T


def score(cfg, name, run_ids, members=None):
    obs, split = load_table(cfg), load_split(cfg)
    out = {"created_utc": utc_now(), "name": name, "runs": list(run_ids), "members_of": members,
           "length_scale_m": float(cfg["fit"]["length_scale_m"]), "partitions": {}}
    member_ids = []
    if members:
        path = cfg.path(cfg["paths"]["ensembles"]) / f"{members}.csv"
        if not path.exists():
            raise StudyError(f"{path} not found")
        member_ids = list(pd.read_csv(path).run_id)
    if cfg["fit"].get("overspread"):                 # v3: area outside the final perimeter, per run
        from .overspread import perimeter_mask, run_acres
        mask = perimeter_mask(cfg)
        out["overspread_acres"] = {r: run_acres(cfg, r, mask) for r in run_ids}
    for which in ("fit", "holdout"):
        c, blocks, w = partition_constraints(cfg, obs, split, which)
        runs = {}
        for rid in run_ids:
            s, d, T = one(cfg, rid, c, w)
            runs[rid] = {**s, "breakdown": breakdown(cfg, c, w, d, T)}
        part = {"constraints": int(len(c)), "blocks": int(len(set(blocks))), "runs": runs}
        if member_ids:
            terms = np.array([one(cfg, r, c, w)[0]["distance_term"] for r in member_ids])
            part["members"] = {"distance_term": dict(zip(member_ids, map(float, terms))),
                               "median": float(np.median(terms)), "min": float(terms.min()),
                               "max": float(terms.max())}
        out["partitions"][which] = part
    write_json(cfg.dir("results") / f"score_{name}.json", out)
    return out


def report(out):
    for rid, a in out.get("overspread_acres", {}).items():
        print(f"{rid:28s}burned outside the final perimeter: {a:,.0f} ac")
    for which, part in out["partitions"].items():
        print(f"\n{which}: {part['constraints']} constraints in {part['blocks']} blocks")
        print(f"  {'run':26s}{'distance term':>14s}{'met':>6s}{'mean km':>9s}   by deadline window: met / mean km")
        for rid, r in part["runs"].items():
            win = "  ".join(f"{k} h {g['met_share']:.2f}/{g['mean_km']:.1f}"
                            for k, g in r["breakdown"]["deadline_h"].items())
            print(f"  {rid:26s}{r['distance_term']:14.3f}{r['met_share']:6.2f}{r['mean_km']:9.2f}   {win}")
        if "members" in part:
            m = part["members"]
            print(f"  {out['members_of']} members' distance term: median {m['median']:.3f}, "
                  f"range {m['min']:.3f}-{m['max']:.3f}")
