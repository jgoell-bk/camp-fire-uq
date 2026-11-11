"""Compare the final runs of v1, v2 and v3 on common yardsticks (fit partition only).

    python -m tools.compare_finals        (from eki_calibration/, elmfire-vnv active)

Each version's Phi uses its own terms (v1: original observations; v3: + over-spread and
the P_EPS penalty), so Phi is not comparable across versions. Here every run is measured
on the same things:
- the cleaned fit deadlines (v2 = v3: 1,062 constraints, v3/config.yaml): distance term,
  weighted RMS km = L sqrt(2 x term), share met, mean km, by deadline window;
- area burned outside the CAL FIRE final perimeter by 24 h (v3/prepared/final_perimeter.tif),
  and the total area burned by 24 h.
The holdout is not touched. Output: v3/results/compare_finals.json.
"""
import numpy as np
import rasterio

from eki.common import Config, read_json, run_dir, utc_now, write_json
from eki.controlled import WINDOWS, fit_constraints
from eki.overspread import burned_outside_acres, perimeter_mask
from eki.pilot import distance_to_burned

RUNS = [("baseline", "config.yaml", "BASELINE"),
        ("v1 final", "config.yaml", "CHK08_I08MEAN"),
        ("v2 final", "v2/config.yaml", "CHK08_I08MEAN"),
        ("v3 final", "v3/config.yaml", "CHK08_I08MEAN")]

ref = Config("v3/config.yaml")
c, blocks, w = fit_constraints(ref)
L = float(ref["fit"]["length_scale_m"])
p = "r{:g}_".format(float(ref["observations"]["proxy_radius_m"]))
rr, cc = c[p + "row"].astype(int).to_numpy(), c[p + "col"].astype(int).to_numpy()
u = c.deadline_s.to_numpy()
cell = float(read_json(ref.shared_dir / "manifest.json")["fuel_grid"]["cellsize"])
inside, cell_acres = perimeter_mask(ref)
at_s = float(ref["fit"]["overspread"]["at_s"])
u_h = u / 3600

out = {"created_utc": utc_now(), "constraints": int(len(c)), "length_scale_m": L, "runs": {}}
print(f"cleaned fit deadlines: {len(c)} in {len(set(blocks))} blocks\n")
print(f"{'run':10s}{'RMS km':>8s}{'met':>6s}{'mean km':>9s}{'0-6 h':>12s}{'6-12 h':>12s}{'12-24 h':>12s}"
      f"{'burned 24 h':>13s}{'outside':>10s}")
for label, cfg_file, rid in RUNS:
    cfg = Config(cfg_file)
    d_run = run_dir(cfg, rid)
    with rasterio.open(d_run / read_json(d_run / "DONE")["toa_file"]) as f:
        toa = f.read(1)
    d = distance_to_burned(toa, rr, cc, u, cell)
    term = float(0.5 * (w * (d / L) ** 2).sum())
    win = {}
    for a, b in WINDOWS:
        k = (u_h >= a) & (u_h < b)
        win[f"{a}-{b:g}"] = {"n": int(k.sum()), "met_share": float((d[k] == 0).mean()),
                             "mean_km": float(d[k].mean() / 1000)}
    burned = ((toa > 0) & (toa <= at_s)) | (toa == -1.0)
    r = {"config": cfg_file, "run_id": rid, "distance_term": term, "rms_km": L * np.sqrt(2 * term) / 1000,
         "met_share": float((d == 0).mean()), "mean_km": float(d.mean() / 1000), "by_window": win,
         "burned_acres_24h": float(burned.sum() * cell_acres),
         "outside_acres_24h": burned_outside_acres(toa, inside, cell_acres, at_s)}
    out["runs"][label] = r
    ws = "".join(f"{g['met_share']:6.2f}/{g['mean_km']:4.1f}" for g in win.values())
    print(f"{label:10s}{r['rms_km']:8.2f}{r['met_share']:6.2f}{r['mean_km']:9.2f}{ws}"
          f"{r['burned_acres_24h']:13,.0f}{r['outside_acres_24h']:10,.0f}")
print("\nwindows: share met / mean km; areas in acres")
write_json(ref.dir("results") / "compare_finals.json", out)
print("wrote v3/results/compare_finals.json")
