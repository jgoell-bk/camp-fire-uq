"""The initial (iteration 0) ensemble: ensembles/I00.csv and I00.json.

For each family (speed, direction, moisture) and member, the K knot values are
drawn together: z ~ N(0, C) with C_rs = rho^|r-s| (Cholesky factor of C times
independent standard normals), u = Phi(z), value = lo + (hi - lo) * u. So every
knot is uniform on its range and neighbouring knots move together. C is the
correlation of the normal scores z; the uniform values' own correlation is
slightly lower (6/pi * arcsin(rho/2), 0.68 for rho = 0.7). The families are drawn
independently, in the order speed, direction, moisture, from one generator.

Speed is drawn as a multiplier on its range and stored as a = ln(multiplier), the
coefficient the parameterization uses.

These ranges only spread the starting ensemble. They are not regularization scales
and do not constrain later updates.
"""
import math

import numpy as np
import pandas as pd

from .common import StudyError, member_id, sha256, utc_now, write_json, write_text
from .corrections import FAMILIES, corrected_weather, theta_names
from .forcing import load_baseline

RANGE_KEYS = {"a": "speed_multiplier", "d": "direction_deg", "m": "moisture_pp"}
_erf = np.vectorize(math.erf)


def normal_cdf(z):
    return 0.5 * (1.0 + _erf(np.asarray(z, float) / math.sqrt(2.0)))


def knot_correlation(rho, k):
    r = np.arange(k)
    return float(rho) ** np.abs(r[:, None] - r[None, :])


def sample(cfg):
    """(theta table, normal scores) for the configured initial ensemble."""
    e = cfg["initial_ensemble"]
    n, k = int(e["size"]), len(cfg["corrections"]["knots_h"])
    L = np.linalg.cholesky(knot_correlation(e["knot_correlation"], k))
    rng = np.random.default_rng(int(e["seed"]))
    cols, scores = {}, {}
    for f in FAMILIES:
        lo, hi = (float(v) for v in e[RANGE_KEYS[f]])
        if not lo < hi or (f == "a" and lo <= 0):
            raise StudyError(f"initial_ensemble.{RANGE_KEYS[f]} must be [lo, hi] with lo < hi"
                             + (" and lo > 0" if f == "a" else ""))
        z = rng.standard_normal((n, k)) @ L.T
        values = lo + (hi - lo) * normal_cdf(z)
        for r in range(k):
            name = f"{f}_k{r + 1}"
            scores[name] = z[:, r]
            if f == "a":
                cols[f"speed_mult_k{r + 1}"] = values[:, r]
                cols[name] = np.log(values[:, r])
            else:
                cols[name] = values[:, r]
    ids = [member_id(0, j) for j in range(1, n + 1)]
    names = theta_names(cfg)
    table = pd.DataFrame({"run_id": ids, "iteration": 0, "member": np.arange(1, n + 1),
                          **{c: cols[c] for c in names},
                          **{c: cols[c] for c in cols if c.startswith("speed_mult")}})
    return table, pd.DataFrame(scores, index=ids)[names]


def rank_report(theta):
    """Singular values of the centred ensemble, each column scaled by its standard
    deviation (so degrees, log multipliers and percentage points compare)."""
    x = np.asarray(theta, float)
    c = x - x.mean(axis=0)
    s = np.linalg.svd(c / c.std(axis=0, ddof=1), compute_uv=False)
    rank = int((s > s[0] * 1e-10).sum())
    return {"members": x.shape[0], "parameters": x.shape[1], "rank": rank,
            "full_rank": rank == x.shape[1], "singular_values": [float(v) for v in s],
            "condition_number": float(s[0] / s[-1]) if s[-1] > 0 else float("inf")}


def correlations(cfg, table, scores):
    """Sample correlation of neighbouring and second-neighbour knots, per family."""
    k = len(cfg["corrections"]["knots_h"])
    out = {}
    for f in FAMILIES:
        for lag in range(1, k):
            pairs = [(f"{f}_k{r}", f"{f}_k{r + lag}") for r in range(1, k - lag + 1)]
            out[f"{f}_lag{lag}"] = {
                "target_normal_scores": float(cfg["initial_ensemble"]["knot_correlation"]) ** lag,
                "normal_scores": float(np.mean([np.corrcoef(scores[p], scores[q])[0, 1] for p, q in pairs])),
                "coefficients": float(np.mean([np.corrcoef(table[p], table[q])[0, 1] for p, q in pairs]))}
    return out


def make_initial(cfg):
    """Sample, check rank and admissibility of every member, write I00.csv and I00.json.
    An existing I00.csv with different contents is never replaced."""
    table, scores = sample(cfg)
    names = theta_names(cfg)
    rank = rank_report(table[names])
    if not rank["full_rank"]:
        raise StudyError(f"initial ensemble has rank {rank['rank']} < {len(names)}")
    baseline = load_baseline(cfg)
    low = {}
    for _, row in table.iterrows():             # raises on the first inadmissible member
        series, _ = corrected_weather(cfg, baseline, {n: row[n] for n in names})
        for v in ("ws", "m1", "m10", "m100"):
            low[v] = min(low.get(v, np.inf), float(series[v].min()))
    out = cfg.dir("ensembles")
    csv_path, meta_path = out / "I00.csv", out / "I00.json"
    text = table.to_csv(index=False, float_format="%.17g")
    if csv_path.exists() and csv_path.read_text() != text:
        raise StudyError(f"{csv_path} exists with different contents; ensemble files are never edited")
    write_text(csv_path, text)
    meta = {"created_utc": utc_now(), "file": csv_path.name, "sha256": sha256(csv_path),
            "config_sha256": cfg.sha256, "settings": cfg["initial_ensemble"],
            "knots_h": cfg["corrections"]["knots_h"], "theta_names": names,
            "method": "per family: z = L e, e ~ N(0, I), L L^T = C, C_rs = rho^|r-s|; "
                      "value = lo + (hi - lo) Phi(z); a = ln(speed multiplier)",
            "rank": rank, "correlations": correlations(cfg, table, scores),
            "minimum_written": {"ws_mph": low["ws"], "m1_pct": low["m1"], "m10_pct": low["m10"],
                                "m100_pct": low["m100"]},
            "baseline_sha256": sha256(cfg.path(cfg["paths"]["prepared"]) / "baseline.csv")}
    write_json(meta_path, meta)
    return table, meta
