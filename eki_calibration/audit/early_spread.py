"""Audit 2026-11-10: why the 0-6 h deadlines stay missed. Reads saved outputs only.

    python audit/early_spread.py          (from eki_calibration/, elmfire-vnv active)

Fit constraints only (v2 = v3, 1,062); the holdout is never read. Runs:
- the final runs at the I08 ensemble-mean parameters, CHK08_I08MEAN of v2 and v3,
  and the uncorrected run (BASELINE);
- every ensemble member of v2 and v3 (I00-I08), from results/I<nn>_distance.npz;
- SEN01 (v1's I08 mean at P_EPS 0.01 / 0.003 / 0.001), re-scored on these constraints.
Each observation is scored at its own deadline u_i with the calibration's operator
(eki.pilot.distance_to_burned); this script also finds the nearest burned cell, for
the maps and the residual directions, and checks that its distances equal the operator's.
Contribution of constraint i to Phi: 1/2 w_i (d_i / L)^2.
Writes audit/out/early_*.{csv,json,png}.
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
import rasterio                                                      # noqa: E402

from eki.common import Config, read_json, run_dir                    # noqa: E402
from eki.controlled import WINDOWS, fit_constraints                  # noqa: E402
from eki.overspread import burned_outside_acres, perimeter_mask      # noqa: E402
from eki.pilot import RADII_CELLS, distance_to_burned                # noqa: E402

OUT = ROOT / "audit" / "out"
OUT.mkdir(parents=True, exist_ok=True)
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
COL = {"uncorrected": "#1baf7a", "v2 final": "#2a78d6", "v3 final": "#eb6834"}
SEQ = ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]       # distance classes, light -> dark
KM_EDGES = [0, 1, 3, 6, 10, 1e9]
KM_LABELS = ["< 1 km", "1-3 km", "3-6 km", "6-10 km", "> 10 km"]
plt.rcParams.update({"font.size": 9, "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.titlesize": 9, "axes.titlecolor": INK, "figure.dpi": 130})

CFG = {"v1": Config(ROOT / "config.yaml"), "v2": Config(ROOT / "v2/config.yaml"), "v3": Config(ROOT / "v3/config.yaml")}
REF = CFG["v3"]
c, blocks, w = fit_constraints(REF)
c = c.reset_index(drop=True)
L = float(REF["fit"]["length_scale_m"])
ROW, COLS = c.r60_row.astype(int).to_numpy(), c.r60_col.astype(int).to_numpy()
U = c.deadline_s.to_numpy(float)
MAN = read_json(REF.shared_dir / "manifest.json")
G = MAN["fuel_grid"]
CS, X0, TOP = float(G["cellsize"]), float(G["x0"]), float(G["top"])
IGN = MAN["ignition"]
INSIDE, CELL_AC = perimeter_mask(REF)
AT_S = float(REF["fit"]["overspread"]["at_s"])
EARLY = U < 6 * 3600
T0 = pd.Timestamp(REF["simulation"]["t0_utc"]).tz_convert("America/Los_Angeles")
RUNS = {"uncorrected": ("v3", "BASELINE"), "v2 final": ("v2", "CHK08_I08MEAN"), "v3 final": ("v3", "CHK08_I08MEAN")}


def xy(row, col):
    return X0 + (np.asarray(col) + 0.5) * CS, TOP - (np.asarray(row) + 0.5) * CS


def pst(seconds):
    return (T0 + pd.to_timedelta(seconds, unit="s")).strftime("%H:%M")


def toa_of(version, rid):
    d = run_dir(CFG[version], rid)
    with rasterio.open(d / read_json(d / "DONE")["toa_file"]) as f:
        return f.read(1)


def nearest(toa, rows, cols, u):
    """distance_to_burned, plus the nearest burned cell (row, col)."""
    H, W = toa.shape
    d, nr, nc = np.full(len(rows), np.nan), np.full(len(rows), -1), np.full(len(rows), -1)
    for i, (r, cl, ui) in enumerate(zip(rows, cols, u)):
        if 0 < toa[r, cl] <= ui:
            d[i], nr[i], nc[i] = 0.0, r, cl
            continue
        for h in RADII_CELLS:
            r0, r1, c0, c1 = max(r - h, 0), min(r + h + 1, H), max(cl - h, 0), min(cl + h + 1, W)
            win = toa[r0:r1, c0:c1]
            a, b = np.nonzero(((win > 0) & (win <= ui)) | (win == -1.0))
            if a.size:
                dist = np.hypot(a + r0 - r, b + c0 - cl)
                k = int(dist.argmin())
                if dist[k] <= h or h == RADII_CELLS[-1]:
                    d[i], nr[i], nc[i] = dist[k] * CS, a[k] + r0, b[k] + c0
                    break
    return d, nr, nc


def contrib(d):
    return 0.5 * w * (np.asarray(d) / L) ** 2


def burned_at(toa, t):
    return ((toa > 0) & (toa <= t)) | (toa == -1.0)


# ---------------------------------------------------------------- the three runs
res, toas = {}, {}
for name, (v, rid) in RUNS.items():
    toa = toa_of(v, rid)
    d, nr, nc = nearest(toa, ROW, COLS, U)
    ref_d = distance_to_burned(toa, ROW, COLS, U, CS)
    if not np.array_equal(d, ref_d):
        raise SystemExit(f"{name}: nearest() differs from eki.pilot.distance_to_burned")
    toas[name] = toa
    res[name] = {"d": d, "nr": nr, "nc": nc, "T": toa[ROW, COLS], "q": contrib(d)}
check = {}
for name in ("v2 final", "v3 final"):
    v = RUNS[name][0]
    diag = read_json(ROOT / v / "results/CHK08_verify.json")["diagnostics"]["CHK08_I08MEAN"]["deadline_h"]
    for a, b in WINDOWS:
        k = (U / 3600 >= a) & (U / 3600 < b)
        mine = round(float(res[name]["q"][k].sum()), 3)
        check[f"{name} {a}-{b:g}"] = {"this_script": mine, "CHK08_verify": diag[f"{a}-{b:g}"]["weighted_distance_term"]}
        if mine != diag[f"{a}-{b:g}"]["weighted_distance_term"]:
            raise SystemExit(f"window term differs from CHK08_verify.json: {name} {a}-{b}")

# ---------------------------------------------------------------- windows
win_rows = []
for name in RUNS:
    q, d = res[name]["q"], res[name]["d"]
    for a, b in WINDOWS:
        k = (U / 3600 >= a) & (U / 3600 < b)
        win_rows.append({"run": name, "window_h": f"{a}-{b:g}", "constraints": int(k.sum()),
                         "weight_share": float(w[k].sum()), "distance_term": float(q[k].sum()),
                         "share_of_run_term": float(q[k].sum() / q.sum()), "met_share": float((d[k] == 0).mean()),
                         "mean_km": float(d[k].mean() / 1000)})
windows = pd.DataFrame(win_rows)
windows.to_csv(OUT / "early_windows.csv", index=False, float_format="%.4f")

# ---------------------------------------------------------------- early: per constraint
e = c[EARLY].copy()
e["block"] = blocks[EARLY]
nb = pd.Series(blocks).value_counts()
e["block_n"] = e.block.map(nb).to_numpy()
e["w"] = w[EARLY]
e["deadline_pst"] = [pst(x) for x in U[EARLY]]
ix, iy = IGN["cell_x"], IGN["cell_y"]
ox, oy = xy(ROW[EARLY], COLS[EARLY])
e["km_from_ignition"] = np.hypot(ox - ix, oy - iy) / 1000
for name in RUNS:
    r = res[name]
    e[f"d_km {name}"] = r["d"][EARLY] / 1000
    e[f"q {name}"] = r["q"][EARLY]
    t = r["T"][EARLY]
    e[f"arrival_pst {name}"] = [pst(x) if x > 0 else ("ignition" if x == -1 else "never") for x in t]
# residual direction: from the nearest burned cell to the observation, against the
# direction from the ignition to the observation (0 deg: the observation lies straight
# ahead of the burned area along the ignition -> observation line)
r3 = res["v3 final"]
nx, ny = xy(r3["nr"][EARLY], r3["nc"][EARLY])
gap = np.degrees(np.arctan2(ox - nx, oy - ny))
radial = np.degrees(np.arctan2(ox - ix, oy - iy))
ang = np.abs((gap - radial + 180) % 360 - 180)
e["v3_gap_bearing_deg"] = np.where(e["d_km v3 final"] > 0, gap % 360, np.nan)
e["v3_gap_vs_radial_deg"] = np.where(e["d_km v3 final"] > 0, ang, np.nan)
cols_out = ["obs_id", "deadline_pst", "deadline_s", "nist_sheet", "nist_fire_type", "r60_status", "block", "block_n",
            "w", "km_from_ignition"] + [f"{p} {n}" for n in RUNS for p in ("d_km", "q", "arrival_pst")] + \
           ["v3_gap_bearing_deg", "v3_gap_vs_radial_deg", "nist_description"]
e[cols_out].sort_values("q v3 final", ascending=False).to_csv(OUT / "early_constraints.csv", index=False,
                                                             float_format="%.5g")

# ---------------------------------------------------------------- by region and block
def group_table(key):
    rows = []
    for g, s in e.groupby(key):
        row = {key: g, "constraints": len(s), "blocks": s.block.nunique(), "weight_share_of_all": float(s.w.sum())}
        for name in RUNS:
            row[f"term {name}"] = float(s[f"q {name}"].sum())
            row[f"share {name}"] = float(s[f"q {name}"].sum() / e[f"q {name}"].sum())
            row[f"met {name}"] = float((s[f"d_km {name}"] == 0).mean())
            row[f"mean_km {name}"] = float(s[f"d_km {name}"].mean())
        rows.append(row)
    return pd.DataFrame(rows).sort_values("term v3 final", ascending=False)


regions, blocks_t = group_table("nist_sheet"), group_table("block")
blocks_t["sheets"] = blocks_t.block.map(e.groupby("block").nist_sheet.agg(lambda s: "/".join(sorted(set(s)))))
regions.to_csv(OUT / "early_by_region.csv", index=False, float_format="%.4f")
blocks_t.to_csv(OUT / "early_by_block.csv", index=False, float_format="%.4f")


def concentration(q):
    q = np.sort(np.asarray(q))[::-1]
    cs = np.cumsum(q) / q.sum()
    return {"top1": float(cs[0]), "top5": float(cs[4]), "top10": float(cs[9]), "top20": float(cs[19]),
            "top50": float(cs[49]), "n_for_50pct": int(np.searchsorted(cs, 0.5) + 1),
            "n_for_80pct": int(np.searchsorted(cs, 0.8) + 1), "n": int(len(q))}


conc = {}
for name in RUNS:
    conc[name] = {"constraints": concentration(e[f"q {name}"]),
                  "blocks": concentration(blocks_t[f"term {name}"]),
                  "share_d_over_1km": float((e[f"d_km {name}"] > 1).mean()),
                  "share_d_over_3km": float((e[f"d_km {name}"] > 3).mean()),
                  "share_d_over_5km": float((e[f"d_km {name}"] > 5).mean()),
                  "term_from_blocks_with_1_or_2_constraints": float(e.loc[e.block_n <= 2, f"q {name}"].sum()
                                                                    / e[f"q {name}"].sum())}
miss = e["d_km v3 final"] > 0
qv = e.loc[miss, "q v3 final"]
a3 = e.loc[miss, "v3_gap_vs_radial_deg"]
direction = {"unmet_constraints": int(miss.sum()),
             "share_of_v3_early_term_with_gap_within_45deg_of_radial": float(qv[a3 <= 45].sum() / qv.sum()),
             "within_45_to_90": float(qv[(a3 > 45) & (a3 <= 90)].sum() / qv.sum()),
             "more_than_90": float(qv[a3 > 90].sum() / qv.sum()),
             "median_angle_weighted_by_count": float(np.median(a3))}

# ---------------------------------------------------------------- named points
special = {}
for oid in ("OBS0001", "OBS0005", "OBS0009"):
    s = e[e.obs_id == oid]
    if len(s):
        s = s.iloc[0]
        special[oid] = {"deadline_pst": s.deadline_pst, "km_from_ignition_cell": float(s.km_from_ignition),
                        "description": s.nist_description,
                        **{name: {"d_km": float(s[f"d_km {name}"]), "arrival_pst": s[f"arrival_pst {name}"]}
                           for name in RUNS}}

# ---------------------------------------------------------------- flexibility: every member
def early_term(D):
    return 0.5 * (w[EARLY] * (D[:, EARLY] / L) ** 2).sum(axis=1)


members = []
per_obs_min = np.full(EARLY.sum(), np.inf)
for v in ("v2", "v3"):
    for it in range(9):
        z = np.load(ROOT / v / f"results/I{it:02d}_distance.npz")
        ids, D = list(z["run_id"]), z["distance_m"]
        acres = z["overspread_acres"] if "overspread_acres" in z.files else None
        if acres is None:                      # v2 has no over-spread output: compute it from the rasters
            acres = np.array([burned_outside_acres(toa_of(v, r), INSIDE, CELL_AC, AT_S) for r in ids])
        et = early_term(D)
        per_obs_min = np.minimum(per_obs_min, D[:, EARLY].min(axis=0))
        for j, r in enumerate(ids):
            members.append({"version": v, "ensemble": f"I{it:02d}", "run_id": r, "early_term": float(et[j]),
                            "early_met_share": float((D[j, EARLY] == 0).mean()),
                            "outside_acres": float(acres[j]), "total_term": float(0.5 * (w * (D[j] / L) ** 2).sum())})
sen = []
for rid in pd.read_csv(ROOT / "ensembles/SEN01.csv").run_id:
    toa = toa_of("v1", rid)
    d = distance_to_burned(toa, ROW, COLS, U, CS)
    per_obs_min = np.minimum(per_obs_min, d[EARLY])
    sen.append({"version": "v1 SEN01", "ensemble": "SEN01", "run_id": rid, "early_term": float(contrib(d)[EARLY].sum()),
                "early_met_share": float((d[EARLY] == 0).mean()),
                "outside_acres": burned_outside_acres(toa, INSIDE, CELL_AC, AT_S), "total_term": float(contrib(d).sum())})
mem = pd.DataFrame(members + sen)
mem.to_csv(OUT / "early_members.csv", index=False, float_format="%.5g")
flex = {"runs_scored": int(len(mem)),
        "min_early_term_calibration_members": mem[mem.version != "v1 SEN01"].nsmallest(3, "early_term")[
            ["version", "run_id", "early_term", "early_met_share", "outside_acres"]].to_dict("records"),
        "SEN01": mem[mem.version == "v1 SEN01"][["run_id", "early_term", "early_met_share", "outside_acres"]].to_dict("records"),
        "early_constraints_never_met_by_any_run": int((per_obs_min > 0).sum()),
        "early_constraints": int(EARLY.sum()),
        "oracle_early_term_best_run_per_constraint": float(0.5 * (w[EARLY] * (per_obs_min / L) ** 2).sum()),
        "note": "oracle: each constraint scored by whichever run came closest; a lower bound no single run attains"}

# ---------------------------------------------------------------- timing table
tim = e.nsmallest(8, "deadline_s")[["obs_id", "deadline_pst", "deadline_s", "nist_sheet", "km_from_ignition",
                                    "d_km v2 final", "arrival_pst v2 final", "d_km v3 final", "arrival_pst v3 final"]]
summary = {"window_check_against_CHK08_verify": check, "concentration": conc, "residual_direction_v3": direction,
           "named_points": special, "flexibility": flex, "earliest_constraints": tim.to_dict("records"),
           "ignition_cell": IGN, "statistic": "runs at the I08 ensemble-mean parameters (CHK08_I08MEAN) and BASELINE; "
                                              "members are individual ensemble runs"}
(OUT / "early_summary.json").write_text(json.dumps(summary, indent=1, default=float))

# ================================================================= figures
mask_rows, mask_cols = np.nonzero(burned_at(toas["v3 final"], 6 * 3600) | burned_at(toas["v2 final"], 6 * 3600))
ex, ey = xy(ROW[EARLY], COLS[EARLY])
fx, fy = xy(mask_rows, mask_cols)
xmin, xmax = min(ex.min(), fx.min(), ix) - 1500, max(ex.max(), fx.max(), ix) + 1500
ymin, ymax = min(ey.min(), fy.min(), iy) - 1500, max(ey.max(), fy.max(), iy) + 1500
c0, c1 = int((xmin - X0) // CS), int((xmax - X0) // CS) + 1
r0, r1 = int((TOP - ymax) // CS), int((TOP - ymin) // CS) + 1
extent_km = [(X0 + c0 * CS) / 1000, (X0 + c1 * CS) / 1000, (TOP - r1 * CS) / 1000, (TOP - r0 * CS) / 1000]
o5 = special.get("OBS0005")
o5xy = xy(ROW[c.obs_id == "OBS0005"], COLS[c.obs_id == "OBS0005"])


def km_class(dkm):
    return np.digitize(dkm, KM_EDGES[1:-1], right=False)


def base_map(ax, toa, t, title):
    ax.set_facecolor("#fcfcfb")
    burned = burned_at(toa, t)[r0:r1, c0:c1]
    ax.imshow(np.where(burned, 1.0, np.nan), extent=extent_km, cmap=matplotlib.colors.ListedColormap([AXIS]),
              vmin=0, vmax=1, interpolation="nearest", zorder=1)
    ax.contour(np.flipud(INSIDE[r0:r1, c0:c1].astype(float)), levels=[0.5], colors=INK2, linewidths=0.8,
               extent=extent_km, origin="lower", zorder=2)
    ax.plot(ix / 1000, iy / 1000, marker="*", ms=11, color=INK, mec="#fcfcfb", mew=0.8, zorder=6)
    ax.plot(o5xy[0] / 1000, o5xy[1] / 1000, marker="P", ms=8, color=INK, mec="#fcfcfb", mew=0.8, zorder=6)
    ax.set_title(title, loc="left")
    ax.set_xticks([]), ax.set_yticks([])
    ax.set_aspect("equal")


def obs_layer(ax, k, dkm):
    k = np.asarray(k)
    met = k & (dkm == 0)
    miss_ = k & (dkm > 0)
    x, y = xy(ROW, COLS)
    cls = km_class(dkm)
    order = np.argsort(dkm[miss_])
    ax.scatter(x[miss_][order] / 1000, y[miss_][order] / 1000, c=[SEQ[j] for j in cls[miss_][order]], s=22,
               edgecolors="#fcfcfb", linewidths=0.6, zorder=4)
    ax.scatter(x[met] / 1000, y[met] / 1000, s=24, facecolors="none", edgecolors=INK, linewidths=0.9, zorder=5)


for name in ("v3 final", "v2 final"):
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.9))
    for j, hh in enumerate(range(1, 7)):
        ax = axes.flat[j]
        lo, hi = (hh - 1) * 3600, hh * 3600
        k = EARLY & (U > lo) & (U <= hi) if hh > 1 else EARLY & (U <= hi)
        q = res[name]["q"][k].sum()
        base_map(ax, toas[name], hi, f"deadlines {pst(lo)}-{pst(hi)} PST: {int(k.sum())} reports, "
                                     f"{q:.2f} of the 0-6 h term\nburned area at {pst(hi)} (grey)")
        obs_layer(ax, k, res[name]["d"] / 1000)
    handles = [matplotlib.lines.Line2D([], [], marker="o", ls="", ms=6, mfc=SEQ[j], mec="#fcfcfb", label=lab)
               for j, lab in enumerate(KM_LABELS)]
    handles += [matplotlib.lines.Line2D([], [], marker="o", ls="", ms=6, mfc="none", mec=INK, label="met (0 km)"),
                matplotlib.patches.Patch(color=AXIS, label="burned by the panel's end time"),
                matplotlib.lines.Line2D([], [], color=INK2, lw=0.8, label="CAL FIRE final perimeter"),
                matplotlib.lines.Line2D([], [], marker="*", ls="", ms=9, color=INK, label="model ignition (Pulga)"),
                matplotlib.lines.Line2D([], [], marker="P", ls="", ms=7, color=INK,
                                        label="OBS0005: reported second ignition (Rim Rd / Concow Rd)")]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=8)
    fig.suptitle(f"{name} (run at the I08 ensemble-mean parameters): 0-6 h reports coloured by distance to the area "
                 f"burned by their own deadline", x=0.01, ha="left", color=INK, fontsize=10)
    fig.tight_layout(rect=(0, 0.07, 1, 0.96))
    fig.savefig(OUT / f"early_map_{name.split()[0]}.png")
    plt.close(fig)

# top contributors, each at its exact deadline
top = e.sort_values("q v3 final", ascending=False).head(8)
fig, axes = plt.subplots(2, 4, figsize=(13, 7))
for ax, (_, s) in zip(axes.flat, top.iterrows()):
    i = int(np.flatnonzero(c.obs_id.to_numpy() == s.obs_id)[0])
    r, cl = ROW[i], COLS[i]
    half = int(max(s["d_km v3 final"], 1.0) * 1000 / CS * 1.35) + 40
    rr0, rr1, cc0, cc1 = max(r - half, 0), r + half, max(cl - half, 0), cl + half
    ext = [(X0 + cc0 * CS) / 1000, (X0 + cc1 * CS) / 1000, (TOP - rr1 * CS) / 1000, (TOP - rr0 * CS) / 1000]
    b3 = burned_at(toas["v3 final"], U[i])[rr0:rr1, cc0:cc1]
    b2 = burned_at(toas["v2 final"], U[i])[rr0:rr1, cc0:cc1].astype(float)
    ax.imshow(np.where(b3, 1.0, np.nan), extent=ext, cmap=matplotlib.colors.ListedColormap([AXIS]), vmin=0, vmax=1,
              interpolation="nearest")
    if b2.any():
        ax.contour(np.flipud(b2), levels=[0.5], colors=COL["v2 final"], linewidths=1.2, extent=ext, origin="lower")
    ax.contour(np.flipud(INSIDE[rr0:rr1, cc0:cc1].astype(float)), levels=[0.5], colors=INK2, linewidths=0.7,
               linestyles="dashed", extent=ext, origin="lower")
    px, py = xy(r, cl)
    qx, qy = xy(res["v3 final"]["nr"][i], res["v3 final"]["nc"][i])
    ax.plot([px / 1000, qx / 1000], [py / 1000, qy / 1000], color=COL["v3 final"], lw=1.6)
    ax.plot(px / 1000, py / 1000, "o", ms=7, color=COL["v3 final"], mec=INK, mew=0.8)
    ax.set_xlim(ext[0], ext[1]), ax.set_ylim(ext[2], ext[3]), ax.set_aspect("equal")
    ax.set_xticks([]), ax.set_yticks([])
    ax.set_title(f"{s.obs_id} ({s.nist_sheet.split(' and ')[0]}), due {s.deadline_pst}\nv3 {s['d_km v3 final']:.1f} km "
                 f"(v2 {s['d_km v2 final']:.1f} km); {s['q v3 final']:.3f} of Phi; block of {s.block_n}", loc="left",
                 fontsize=8)
handles = [matplotlib.patches.Patch(color=AXIS, label="v3 final: burned by the report's deadline"),
           matplotlib.lines.Line2D([], [], color=COL["v2 final"], lw=1.2, label="v2 final: edge of area burned by the deadline"),
           matplotlib.lines.Line2D([], [], color=COL["v3 final"], lw=1.6, marker="o", label="report and its v3 distance"),
           matplotlib.lines.Line2D([], [], color=INK2, lw=0.7, ls="--", label="final perimeter")]
fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8)
fig.suptitle("The eight largest 0-6 h contributions to v3's objective, each at its own deadline "
             "(runs at the I08 mean parameters)", x=0.01, ha="left", color=INK, fontsize=10)
fig.tight_layout(rect=(0, 0.05, 1, 0.96))
fig.savefig(OUT / "early_top_contributors.png")
plt.close(fig)

# by region: grouped horizontal bars (v2, v3), direct labels
reg = regions.sort_values("term v3 final")
fig, ax = plt.subplots(figsize=(7.5, 4.2))
yk = np.arange(len(reg))
for off, name in ((0.2, "v2 final"), (-0.2, "v3 final")):
    vals = reg[f"term {name}"].to_numpy()
    ax.barh(yk + off, vals, height=0.38, color=COL[name], label=name)
    for yy, vv in zip(yk + off, vals):
        ax.text(vv + 0.02, yy, f"{vv:.2f}", va="center", fontsize=7, color=INK2)
ax.set_yticks(yk, [f"{s} ({n})" for s, n in zip(reg.nist_sheet, reg.constraints)])
ax.set_xlabel("contribution to the distance term, 0-6 h deadlines (NIST sheet, reports)")
ax.xaxis.grid(True, color=GRID, lw=0.6), ax.set_axisbelow(True)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.legend(frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(OUT / "early_by_region.png")
plt.close(fig)

# concentration: cumulative share of the 0-6 h term against the number of reports
fig, ax = plt.subplots(figsize=(6.5, 3.8))
for name in ("v2 final", "v3 final"):
    q = np.sort(e[f"q {name}"].to_numpy())[::-1]
    cs = np.cumsum(q) / q.sum()
    ax.plot(np.arange(1, len(q) + 1), cs, color=COL[name], lw=2, label=name)
    n50 = conc[name]["constraints"]["n_for_50pct"]
    ax.annotate(f"{n50} reports give 50 %", (n50, 0.5), xytext=(n50 + 25, 0.5 - (0.12 if name == "v3 final" else 0.05)),
                fontsize=7, color=INK2, arrowprops={"arrowstyle": "-", "color": MUTED, "lw": 0.6})
ax.set_xlabel("0-6 h reports, largest contribution first")
ax.set_ylabel("cumulative share of the 0-6 h term")
ax.set_ylim(0, 1.02), ax.set_xlim(0, len(e))
ax.grid(True, color=GRID, lw=0.6), ax.set_axisbelow(True)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.legend(frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(OUT / "early_concentration.png")
plt.close(fig)

# trade-off: every member's 0-6 h term against its area outside the perimeter
fig, ax = plt.subplots(figsize=(7.2, 4.6))
for v, name in (("v2", "v2 final"), ("v3", "v3 final")):
    s = mem[mem.version == v]
    ax.scatter(s.outside_acres / 1000, s.early_term, s=14, color=COL[name], alpha=0.55, edgecolors="none",
               label=f"{v} members, I00-I08")
for name in RUNS:
    acres = burned_outside_acres(toas[name], INSIDE, CELL_AC, AT_S)
    ax.scatter([acres / 1000], [res[name]["q"][EARLY].sum()], s=70, color=COL[name], edgecolors=INK, linewidths=1,
               zorder=5)
    ax.annotate(name, (acres / 1000, res[name]["q"][EARLY].sum()), xytext=(6, 6), textcoords="offset points",
                fontsize=8, color=INK)
s = mem[mem.version == "v1 SEN01"]
ax.scatter(s.outside_acres / 1000, s.early_term, s=40, marker="D", facecolors="none", edgecolors=INK, linewidths=1,
           label="SEN01: v1 I08 mean at P_EPS 0.01 / 0.003 / 0.001", zorder=5)
for _, r in s.iterrows():
    ax.annotate(r.run_id.split("_")[-1], (r.outside_acres / 1000, r.early_term), xytext=(6, -10),
                textcoords="offset points", fontsize=7, color=INK2)
ax.set_xlabel("area burned outside the final perimeter by 24 h (thousand acres)")
ax.set_ylabel("0-6 h part of the distance term")
ax.grid(True, color=GRID, lw=0.6), ax.set_axisbelow(True)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.legend(frameon=False, fontsize=7, loc="upper right")
fig.tight_layout()
fig.savefig(OUT / "early_tradeoff.png")
plt.close(fig)

print(windows.to_string(index=False))
print("\nby region (0-6 h):")
print(regions[["nist_sheet", "constraints", "blocks", "weight_share_of_all", "term v2 final", "term v3 final",
               "share v3 final", "met v3 final", "mean_km v2 final", "mean_km v3 final"]].to_string(index=False))
print("\ntop blocks:")
print(blocks_t[["block", "sheets", "constraints", "term v2 final", "term v3 final", "share v3 final"]].head(10).to_string(index=False))
print("\ntop constraints:")
print(e.sort_values("q v3 final", ascending=False)[["obs_id", "deadline_pst", "nist_sheet", "block_n", "km_from_ignition",
                                                    "d_km v2 final", "d_km v3 final", "q v2 final", "q v3 final",
                                                    "arrival_pst v3 final"]].head(12).to_string(index=False))
print(json.dumps({k: summary[k] for k in ("concentration", "residual_direction_v3", "named_points", "flexibility")},
                 indent=1, default=float))
