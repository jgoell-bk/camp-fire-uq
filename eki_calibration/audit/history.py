"""Audit 2026-11-10: the v3 (and v2) iteration history, and the alpha = 4 step from I08.

    python audit/history.py          (from eki_calibration/, elmfire-vnv active)

Reads saved files only (results/CHK<nn>_verify.json, proposals/I<nn>.json,
ensembles/I<nn>.csv, results/I<nn>_distance.npz). For every check:
- Phi and its parts from the runs at the current and the proposed ensemble means;
- alpha and the largest standardized moves;
- a linearized prediction of the change at the mean: the members' outputs F = [sqrt(w) d / L;
  A_out / A_ref; B theta] regressed on theta (least squares over the 32 members), evaluated
  at the two means. It is a diagnostic only; nothing in the calibration uses it.
Per ensemble: each coefficient's sd in units of its penalty scale, and P_EPS quantiles.
Then, from the unchanged I08 ensemble and its responses: the alpha = 2 step (must reproduce
v3/proposals/I09.csv, the rejected proposal) and the alpha = 4 step through the same EKI
solve, with the existing step limits and admissibility checks.
Writes audit/out/history_*.{csv,json,png}.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib                                                    # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                      # noqa: E402
import numpy as np                                                   # noqa: E402
import pandas as pd                                                  # noqa: E402

from eki.common import Config, read_json                             # noqa: E402
from eki.controlled import admissible, fit_constraints               # noqa: E402
from eki.corrections import theta_names                              # noqa: E402
from eki.objective import (distance_outputs, overspread_outputs, penalty_scale_vector,   # noqa: E402
                           regularization_matrix)
from eki.update import eki_update, step_limits                       # noqa: E402

OUT = ROOT / "audit" / "out"
OUT.mkdir(parents=True, exist_ok=True)
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SER = {"distance": "#2a78d6", "over-spread": "#eb6834", "amplitude": "#1baf7a", "smoothness": "#eda100"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.titlesize": 9, "axes.titlecolor": INK, "figure.dpi": 130})


def outputs(cfg, theta, D, acres, w, B):
    F = [distance_outputs(D, w, cfg["fit"]["length_scale_m"])]
    if cfg["fit"].get("overspread"):
        F.append(overspread_outputs(cfg, acres))
    F.append(theta @ B.T)
    return np.hstack(F)


def members(cfg, n):
    ens = pd.read_csv(cfg.path(cfg["paths"]["ensembles"]) / f"I{n:02d}.csv")
    z = np.load(cfg.dir("results") / f"I{n:02d}_distance.npz")
    ids = list(z["run_id"])
    rows = [ids.index(r) for r in ens.run_id]
    acres = z["overspread_acres"][rows] if "overspread_acres" in z.files else None
    return ens, ens[theta_names(cfg)].to_numpy(float), z["distance_m"][rows], acres


def linear_phi(theta, F, points):
    """1/2 |F_lin(x)|^2 at each point, F_lin = F_bar + (x - theta_bar) A from least squares."""
    tb, fb = theta.mean(0), F.mean(0)
    A, *_ = np.linalg.lstsq(theta - tb, F - fb, rcond=None)
    return [float(0.5 * ((fb + (x - tb) @ A) ** 2).sum()) for x in points]


def history(version):
    cfg = Config(ROOT / version / "config.yaml")
    c, blocks, w = fit_constraints(cfg)
    B = regularization_matrix(cfg)
    names = theta_names(cfg)
    scale = penalty_scale_vector(cfg)
    rows, spread = [], []
    n = 1
    while (cfg.dir("results") / f"CHK{n:02d}_verify.json").exists():
        v = read_json(cfg.dir("results") / f"CHK{n:02d}_verify.json")
        prop = read_json(ROOT / version / f"proposals/I{n:02d}.json")["parent"]
        cur, new = f"CHK{n:02d}_I{n - 1:02d}MEAN", f"CHK{n:02d}_I{n:02d}MEAN"
        oc, on = v["objective"][cur], v["objective"][new]
        ens, th, D, acres = members(cfg, n - 1)
        F = outputs(cfg, th, D, acres, w, B)
        m0, m1 = np.array([prop["current_mean"][k] for k in names]), np.array([prop["proposed_mean"][k] for k in names])
        lin0, lin1 = linear_phi(th, F, [m0, m1])
        row = {"version": version, "check": f"CHK{n:02d}", "from": f"I{n - 1:02d}", "alpha": prop["step"]["alpha"],
               "max_mean_move_sd": max(prop["step"]["mean_move_sd"].values()),
               "max_member_move_sd": max(prop["step"]["max_member_move_sd"].values()),
               "repeat_bit_identical": v["repeat_bit_identical"], "variability": v["variability"],
               "accepted": v["accepted"], "verdict": v["verdict"],
               "members_phi_mean": prop["current_objective"]["total"]["mean"],
               "members_phi_min": prop["current_objective"]["total"]["min"]}
        for part in ("distance", "overspread", "amplitude", "smoothness", "total"):
            if part in oc:
                row[f"cur_{part}"], row[f"new_{part}"] = oc[part], on[part]
        if "overspread_acres" in v:
            row["cur_outside_acres"], row["new_outside_acres"] = v["overspread_acres"][cur], v["overspread_acres"][new]
        row["actual_change"] = on["total"] - oc["total"]
        row["linear_predicted_change"] = lin1 - lin0
        row["linear_phi_at_current_mean"] = lin0
        row["ratio_actual_to_predicted"] = row["actual_change"] / row["linear_predicted_change"]
        rows.append(row)
        n += 1
    for k in range(n):
        p = cfg.path(cfg["paths"]["ensembles"]) / f"I{k:02d}.csv"
        if not p.exists():
            break
        th = pd.read_csv(p)[names].to_numpy(float)
        sd = th.std(0, ddof=1) / scale
        rec = {"version": version, "ensemble": f"I{k:02d}", **{f"sd/scale {nm}": float(s) for nm, s in zip(names, sd)}}
        if "p" in names:
            pe = 0.01 * np.exp(th[:, names.index("p")])
            rec.update({"P_EPS_min": float(pe.min()), "P_EPS_median": float(np.median(pe)), "P_EPS_max": float(pe.max()),
                        "P_EPS_at_mean_p": float(0.01 * np.exp(th[:, names.index("p")].mean()))})
        spread.append(rec)
    return cfg, pd.DataFrame(rows), pd.DataFrame(spread)


cfg3, h3, s3 = history("v3")
_, h2, s2 = history("v2")
hist = pd.concat([h2, h3], ignore_index=True)
hist.to_csv(OUT / "history_checks.csv", index=False, float_format="%.6g")
pd.concat([s2, s3], ignore_index=True).to_csv(OUT / "history_spread.csv", index=False, float_format="%.5g")

# ---------------------------------------------------------------- the step from I08: alpha 2 and 4
c, blocks, w = fit_constraints(cfg3)
B = regularization_matrix(cfg3)
names = theta_names(cfg3)
ens, th, D, acres = members(cfg3, 8)
G = np.hstack([distance_outputs(D, w, cfg3["fit"]["length_scale_m"]), overspread_outputs(cfg3, acres)])
y, R = np.zeros(G.shape[1]), np.ones(G.shape[1])
floor = float(cfg3["fit"]["step"]["sd_floor"]) * penalty_scale_vector(cfg3)
ok = admissible(cfg3)
rejected = pd.read_csv(ROOT / "v3/proposals/I09.csv")[names].to_numpy(float)
F = outputs(cfg3, th, D, acres, w, B)
steps = {}
for alpha in (1.0, 2.0, 4.0, 8.0):
    new = eki_update(th, G, y, R, B, alpha)
    mm, mx, _ = step_limits(th, new, floor)
    lin0, lin1 = linear_phi(th, F, [th.mean(0), new.mean(0)])
    s = cfg3["fit"]["step"]
    steps[alpha] = {"max_mean_move_sd": float(mm.max()), "max_member_move_sd": float(mx.max()),
                    "inadmissible_members": [j for j, r in enumerate(new) if not ok(r)],
                    "within_limits": bool(mm.max() <= s["max_mean_sd"] and mx.max() <= s["max_member_sd"]),
                    "mean_move_sd": dict(zip(names, map(float, mm))),
                    "proposed_mean": dict(zip(names, map(float, new.mean(0)))),
                    "linear_predicted_change": lin1 - lin0}
    if alpha == 2.0:
        steps[alpha]["max_abs_diff_from_v3_proposals_I09"] = float(np.abs(new - rejected).max())
    if alpha == 4.0:
        d2 = rejected.mean(0) - th.mean(0)
        d4 = new.mean(0) - th.mean(0)
        steps[alpha]["mean_step_ratio_to_alpha2_per_coefficient"] = dict(zip(names, map(float, d4 / d2)))
        steps[alpha]["mean_step_norm_ratio_to_alpha2_in_sd_units"] = float(
            np.linalg.norm(d4 / th.std(0, ddof=1)) / np.linalg.norm(d2 / th.std(0, ddof=1)))
        steps[alpha]["cosine_with_alpha2_step_in_sd_units"] = float(
            (d4 / th.std(0, ddof=1)) @ (d2 / th.std(0, ddof=1))
            / np.linalg.norm(d4 / th.std(0, ddof=1)) / np.linalg.norm(d2 / th.std(0, ddof=1)))
(OUT / "history_I08_steps.json").write_text(json.dumps({str(k): v for k, v in steps.items()}, indent=1))

# ---------------------------------------------------------------- figure: Phi at the mean, by part (v3)
fig, ax = plt.subplots(figsize=(7.4, 4.0))
it = np.arange(len(h3) + 1)
parts = {"distance": "distance", "over-spread": "overspread", "amplitude": "amplitude", "smoothness": "smoothness"}
for label, key in parts.items():
    vals = list(h3[f"cur_{key}"]) + [h3[f"new_{key}"].iloc[-1]]
    ax.plot(it[:-1], vals[:-1], color=SER[label], lw=2, marker="o", ms=4, label=label)
    ax.plot(it[-2:], vals[-2:], color=SER[label], lw=1.2, ls=":", marker="o", ms=4, mfc="#fcfcfb")
    ax.annotate(f"{vals[-2]:.2f}", (it[-2], vals[-2]), xytext=(-30, 4), textcoords="offset points", fontsize=7,
                color=INK2)
ax.set_xticks(it, [f"I{k:02d}" for k in it[:-1]] + ["I09\n(rejected)"])
ax.set_ylabel("part of Phi, run at the ensemble mean")
ax.set_yscale("log")
ax.grid(True, color=GRID, lw=0.6), ax.set_axisbelow(True)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.legend(frameon=False, ncol=4, loc="upper right", fontsize=8)
ax.set_title("v3: objective parts at each ensemble mean (dotted: the rejected alpha = 2 proposal)", loc="left")
fig.tight_layout()
fig.savefig(OUT / "history_v3_parts.png")
plt.close(fig)

pd.set_option("display.width", 250)
print(hist[["version", "check", "alpha", "max_mean_move_sd", "cur_total", "new_total", "actual_change",
            "linear_predicted_change", "ratio_actual_to_predicted", "accepted"]].round(4).to_string(index=False))
print(h3[["check", "cur_distance", "cur_overspread", "cur_amplitude", "cur_smoothness", "cur_total",
          "new_distance", "new_overspread", "new_amplitude", "new_smoothness", "new_total",
          "cur_outside_acres", "new_outside_acres"]].round(4).to_string(index=False))
print(s3.round(3).to_string(index=False))
print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk not in ("mean_move_sd", "proposed_mean")}
                  for k, v in steps.items()}, indent=1))
