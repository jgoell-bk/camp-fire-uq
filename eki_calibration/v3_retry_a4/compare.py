"""After the retry's check run and `verify`: the I08 control, v3's rejected alpha = 2
proposal and the alpha = 4 retry, side by side (defined before the run, audit 2026-11-10).

    python -m eki --config v3_retry_a4/config.yaml verify CHK09
    python v3_retry_a4/compare.py

All three are single runs at ensemble-mean parameters (never averages of member scores),
with SEED 2023, on v3's objective. Acceptance is the existing rule (eki.controlled.verify):
accepted iff Phi(proposed mean) - Phi(I08 mean) <= variability = |Phi(repeat) - Phi(I08 mean)|,
which is 0 when the repeat is bit-identical. Acceptance means this one step did not raise
Phi at the mean for this seed; it is not evidence of convergence.
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

v3 = json.loads((ROOT / "v3/results/CHK09_verify.json").read_text())
path = HERE / "results/CHK09_verify.json"
if not path.exists():
    sys.exit(f"{path} not found: run `python -m eki --config v3_retry_a4/config.yaml verify CHK09` first")
rt = json.loads(path.read_text())
cols = {"I08 mean (control)": (v3, "CHK09_I08MEAN"), "rejected I09, alpha 2": (v3, "CHK09_I09MEAN"),
        "retry I09, alpha 4": (rt, "CHK09_I09MEAN"), "retry control": (rt, "CHK09_I08MEAN")}
parts = ("distance", "overspread", "amplitude", "smoothness", "total")
print(f"{'':22s}" + "".join(f"{k:>24s}" for k in cols))
for p in parts:
    print(f"{p:22s}" + "".join(f"{v['objective'][r][p]:24.4f}" for v, r in cols.values()))
print(f"{'outside acres':22s}" + "".join(f"{v['overspread_acres'][r]:24,.0f}" for v, r in cols.values()))
for win in ("0-6", "6-12", "12-24.01"):
    for key, lab in (("weighted_distance_term", "term"), ("met_share", "met"), ("mean_km", "mean km")):
        print(f"{win + ' h ' + lab:22s}" + "".join(
            f"{v['diagnostics'][r]['deadline_h'][win][key]:24.3f}" for v, r in cols.values()))
ctrl = rt["objective"]["CHK09_I08MEAN"]["total"]
if abs(ctrl - v3["objective"]["CHK09_I08MEAN"]["total"]) > 0:
    print("WARNING: the retry's control differs from v3's I08-mean run")
print(f"\nrepeat bit-identical: {rt['repeat_bit_identical']}; variability {rt['variability']:g}")
print(f"retry minus control: {rt['change_proposed_minus_current']:+.4f} -> {rt['verdict']}; "
      f"accepted under the existing rule: {rt['accepted']}")
print(f"(v3's alpha = 2 proposal: {v3['change_proposed_minus_current']:+.4f}, {v3['verdict']})")
print("Accepted or not, this is one step at one seed: it is not a convergence test.")
