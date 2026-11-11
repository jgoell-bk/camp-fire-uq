"""Reuse v3's validated I08-mean control runs in the retry, only if their inputs match.

    python v3_retry_a4/reuse_controls.py        (from eki_calibration/, after `propose I08`)

The retry's check batch CHK09 is CHK09_I08MEAN, CHK09_I08MEAN_REP (the I08 mean, twice:
the repeat measures run-to-run variability) and CHK09_I09MEAN (the alpha = 4 proposal).
v3 already ran the first two (validated, bit-identical to each other). For each of them
this script:
1. checks that the retry's CHK09.csv row equals v3's CHK09.csv row as text;
2. builds the run here with the retry config into a temporary folder (eki.case.build_case)
   and compares it with v3's run: theta, every weather raster checksum, weather.json bands,
   P_EPS, SEED, SIMULATION_TSTART/TSTOP, template, baseline, ELMFIRE binary, shared input
   files, and the namelist with its four folder paths replaced by placeholders;
3. only if all of that is equal, copies v3's run folder into this retry's scratch root
   (v3's copy is not touched) and adds REUSED.json with the comparison.
Expected differences, recorded but not required to match: config_sha256 (the retry config
differs in fit.step.alpha0, the scratch root and auto.max_ensemble), coefficients_from (the
CHK09.csv path and checksum: its third row is the new proposal), shared_manifest_sha256
(creation time), built_utc and host.
It also checks proposals/I09.csv against an independent alpha = 4 EKI step from the I08
files and records how it differs from v3's rejected (alpha = 2) proposal.
Then only CHK09_I09MEAN is left to run.
"""
import json
import re
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

import numpy as np                                                   # noqa: E402
import pandas as pd                                                  # noqa: E402

from eki.case import build_case                                      # noqa: E402
from eki.common import Config, read_json, run_dir, sha256, utc_now   # noqa: E402
from eki.controlled import fit_constraints                           # noqa: E402
from eki.corrections import theta_names                              # noqa: E402
from eki.objective import distance_outputs, overspread_outputs, regularization_matrix   # noqa: E402
from eki.update import eki_update                                    # noqa: E402

CONTROLS = ("CHK09_I08MEAN", "CHK09_I08MEAN_REP")
PATH_KEYS = ("FUELS_AND_TOPOGRAPHY_DIRECTORY", "WEATHER_DIRECTORY", "MISCELLANEOUS_INPUTS_DIRECTORY", "SCRATCH")


def normalized_namelist(text):
    for k in PATH_KEYS:
        text, n = re.subn(rf"(?m)^{k} = '[^']*'", f"{k} = '<{k}>'", text)
        if n != 1:
            raise SystemExit(f"{k} not found exactly once in a namelist")
    return text


def check_proposal(cfg):
    c, blocks, w = fit_constraints(cfg)
    names = theta_names(cfg)
    ens = pd.read_csv(HERE / "ensembles/I08.csv")
    z = np.load(HERE / "results/I08_distance.npz")
    ids = list(z["run_id"])
    rows = [ids.index(r) for r in ens.run_id]
    th = ens[names].to_numpy(float)
    G = np.hstack([distance_outputs(z["distance_m"][rows], w, cfg["fit"]["length_scale_m"]),
                   overspread_outputs(cfg, z["overspread_acres"][rows])])
    new = eki_update(th, G, np.zeros(G.shape[1]), np.ones(G.shape[1]), regularization_matrix(cfg), 4.0)
    prop = pd.read_csv(HERE / "proposals/I09.csv")[names].to_numpy(float)
    rej = pd.read_csv(ROOT / "v3/proposals/I09.csv")[names].to_numpy(float)
    rec = read_json(HERE / "proposals/I09.json")["parent"]["step"]
    out = {"alpha_chosen": rec["alpha"], "alphas_tried": [t["alpha"] for t in rec["tried"]],
           "max_abs_diff_from_independent_alpha4_step": float(np.abs(prop - new).max()),
           "max_abs_diff_from_v3_rejected_alpha2": float(np.abs(prop - rej).max()),
           "proposed_mean": dict(zip(names, map(float, prop.mean(0)))),
           "v3_rejected_mean": dict(zip(names, map(float, rej.mean(0)))),
           "I08_mean": dict(zip(names, map(float, th.mean(0))))}
    if rec["alpha"] != 4.0 or out["max_abs_diff_from_independent_alpha4_step"] > 1e-12:
        raise SystemExit(f"proposals/I09.csv is not the alpha = 4 step: {out}")
    return out


def main():
    cfg, v3 = Config(HERE / "config.yaml"), Config(ROOT / "v3/config.yaml")
    report = {"created_utc": utc_now(), "proposal": check_proposal(cfg), "controls": {}}
    mine = pd.read_csv(HERE / "ensembles/CHK09.csv", dtype=str).set_index("run_id")
    theirs = pd.read_csv(ROOT / "v3/ensembles/CHK09.csv", dtype=str).set_index("run_id")
    for rid in CONTROLS:
        if not mine.loc[rid].equals(theirs.loc[rid]):
            raise SystemExit(f"{rid}: the retry's CHK09.csv row differs from v3's")
        src, dst = run_dir(v3, rid), run_dir(cfg, rid)
        done = read_json(src / "DONE")
        if sha256(src / done["toa_file"]) != done["toa_sha256"]:
            raise SystemExit(f"{src}: arrival raster does not match its DONE checksum")
        old = read_json(src / "run.json")
        with tempfile.TemporaryDirectory() as tmp:
            new = build_case(cfg, rid, Path(tmp) / "case", Path(tmp) / "case" / "scratch")
            new_nl = (Path(tmp) / "case" / "elmfire.data.in").read_text()
            new_wx = read_json(Path(tmp) / "case" / "inputs" / "weather.json")
        old_nl, old_wx = (src / "elmfire.data.in").read_text(), read_json(src / "weather.json")
        must = {"theta": old["theta"] == new["theta"],
                "theta_names": old["theta_names"] == new["theta_names"],
                "weather_raster_sha256": old["weather_raster_sha256"] == new["weather_raster_sha256"],
                "weather_json_bands": old_wx["bands"] == new_wx["bands"],
                "namelist_without_paths": normalized_namelist(old_nl) == normalized_namelist(new_nl),
                **{k: old[k] == new[k] for k in ("P_EPS", "SEED", "SIMULATION_TSTART", "SIMULATION_TSTOP",
                                                 "template_sha256", "baseline_sha256", "elmfire",
                                                 "namelist_overrides", "knots_h", "moisture_weights")}}
        shared_old, shared_new = read_json(v3.shared_dir / "manifest.json"), read_json(cfg.shared_dir / "manifest.json")
        must["shared_input_files"] = {k: v["sha256"] for k, v in shared_old["files"].items()} == \
                                     {k: v["sha256"] for k, v in shared_new["files"].items()}
        expected = {k: {"v3": old[k], "retry": new[k]} for k in
                    ("config_sha256", "coefficients_from", "shared_manifest_sha256", "built_utc", "host")}
        if not all(must.values()):
            raise SystemExit(f"{rid}: inputs differ from v3's run, not reused: "
                             f"{[k for k, v in must.items() if not v]}")
        if dst.exists():
            if read_json(dst / "DONE")["toa_sha256"] != done["toa_sha256"]:
                raise SystemExit(f"{dst} exists and is not v3's run")
        else:
            shutil.copytree(src, dst)
        rec = {"reused_from": str(src), "copied_utc": utc_now(), "toa_sha256": done["toa_sha256"],
               "equal": must, "expected_differences": expected,
               "note": "v3's validated run, reused because the retry config builds identical inputs (paths aside)"}
        (dst / "REUSED.json").write_text(json.dumps(rec, indent=1) + "\n")
        report["controls"][rid] = rec
    same = read_json(run_dir(cfg, CONTROLS[0]) / "DONE")["toa_sha256"] == \
        read_json(run_dir(cfg, CONTROLS[1]) / "DONE")["toa_sha256"]
    report["controls_bit_identical_to_each_other"] = bool(same)
    report["left_to_run"] = [r for r in pd.read_csv(HERE / "ensembles/CHK09.csv").run_id
                             if not (run_dir(cfg, r) / "DONE").exists()]
    (HERE / "REUSE.json").write_text(json.dumps(report, indent=1) + "\n")
    print(json.dumps({"proposal": {k: v for k, v in report["proposal"].items() if not k.endswith("mean")},
                      "controls": {r: all(v["equal"].values()) for r, v in report["controls"].items()},
                      "controls_bit_identical": same, "left_to_run": report["left_to_run"]}, indent=1))


if __name__ == "__main__":
    main()
