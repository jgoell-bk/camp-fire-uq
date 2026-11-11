"""Set up v3_retry_a4: one damped retry of v3's rejected I08 -> I09 step (audit, 2026-11-10).

    python v3_retry_a4/setup_retry.py            (from eki_calibration/, elmfire-vnv active)

What changes from v3: fit.step.alpha0 1 -> 4 (alpha doubles from alpha0 until the step
limits hold, so the smallest alpha allowed is 4; v3's rejected step used alpha 2). The
objective, observations, constraints, parameterization, penalty scales, over-spread
reference, seed and acceptance rule are v3's. auto.max_ensemble is set to 7 so that
`python -m eki auto` refuses to run here (the newest ensemble, I08, is past it): this
folder is for one manual check, not a loop.

Steps (each checked, nothing in v3/ is written):
1. config.yaml = v3/config.yaml with only the lines listed in EDITS changed.
2. Copy, byte for byte: v3's template, prepared/ (incl. the perimeter raster),
   ensembles/I08.csv(.json) and results/I08_* (the I08 members' responses).
3. Build this folder's scratch shared/ inputs (eki.case.prepare_shared) and compare every
   file's checksum with v3's shared manifest.
Writes PROVENANCE.json. Then run `python -m eki --config v3_retry_a4/config.yaml propose I08`
and `python v3_retry_a4/reuse_controls.py` (see README.md).
"""
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from eki.common import Config, read_json, sha256, utc_now   # noqa: E402

V3 = ROOT / "v3"
EDITS = [  # (old line, new line), each must occur exactly once in v3/config.yaml
    ("    alpha0: 1.0\n",
     "    alpha0: 4.0                    # v3_retry_a4: 1.0 in v3. Smallest alpha allowed = 4 (v3's rejected I09 used 2)\n"),
    ("  scratch_root: $SCRATCH/elmfire-class/eki_calibration_v3   # shared/ and runs/ (v3's own)\n",
     "  scratch_root: $SCRATCH/elmfire-class/eki_calibration_v3_retry_a4   # shared/ and runs/ (this retry's own)\n"),
    ("  max_ensemble: 10            # after running and scoring I10\n",
     "  max_ensemble: 7             # v3_retry_a4 guard: auto refuses here (I08 > I07); one manual check only\n"),
]
HEADER = ("# v3_retry_a4 (audit, 2026-11-10): v3 with one change, fit.step.alpha0 = 4, to retry the\n"
          "# rejected I08 -> I09 step with more damping. Everything else, including the objective, the\n"
          "# seed and the acceptance rule, is v3's. Starting state: v3's I08 (copied). See README.md here.\n"
          "# Run from eki_calibration/ with --config v3_retry_a4/config.yaml.\n#\n")
COPIES = ["elmfire.data.template", "ensembles/I08.csv", "ensembles/I08.json",
          "results/I08_G.npz", "results/I08_distance.npz", "results/I08_pilot.json",
          "results/I08_pilot_cells.csv", "results/I08_reach.csv"] + \
         [f"prepared/{p.name}" for p in sorted((V3 / "prepared").iterdir()) if p.is_file()]


def main():
    text = (V3 / "config.yaml").read_text()
    for old, new in EDITS:
        if text.count(old) != 1:
            raise SystemExit(f"expected exactly one {old!r} in v3/config.yaml")
        text = text.replace(old, new)
    cfg_path = HERE / "config.yaml"
    if cfg_path.exists() and cfg_path.read_text() != HEADER + text:
        raise SystemExit(f"{cfg_path} exists with other contents")
    cfg_path.write_text(HEADER + text)

    copies = {}
    for rel in COPIES:
        src, dst = V3 / rel, HERE / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists() and sha256(dst) != sha256(src):
            raise SystemExit(f"{dst} exists and differs from {src}")
        if not dst.exists():
            shutil.copy2(src, dst)
        if sha256(dst) != sha256(src):
            raise SystemExit(f"copy of {src} does not match")
        copies[rel] = {"source": str(src.relative_to(ROOT)), "sha256": sha256(dst)}

    cfg = Config(cfg_path)
    v3 = Config(V3 / "config.yaml")
    diff = {k: (v3[k], cfg[k]) for k in v3 if v3[k] != cfg[k]}
    allowed = {"fit", "paths", "auto"}
    if set(diff) - allowed:
        raise SystemExit(f"unexpected config differences: {sorted(set(diff) - allowed)}")
    f3, fr = dict(v3["fit"]), dict(cfg["fit"])
    s3, sr = f3.pop("step"), fr.pop("step")
    if f3 != fr or {k: v for k, v in s3.items() if k != "alpha0"} != {k: v for k, v in sr.items() if k != "alpha0"}:
        raise SystemExit("fit differs from v3 beyond fit.step.alpha0")

    from eki.case import prepare_shared
    if not (cfg.shared_dir / "manifest.json").exists():
        prepare_shared(cfg)
    mine, theirs = read_json(cfg.shared_dir / "manifest.json"), read_json(v3.shared_dir / "manifest.json")
    same = {k: mine[k] == theirs[k] for k in ("elmfire", "fuel_grid", "weather_grid", "ignition", "ignition_csv")}
    same["files"] = {rel: mine["files"][rel]["sha256"] == theirs["files"][rel]["sha256"] for rel in theirs["files"]}
    if not all(v if isinstance(v, bool) else all(v.values()) for v in same.values()) or \
            set(mine["files"]) != set(theirs["files"]):
        raise SystemExit(f"shared inputs differ from v3's: {same}")

    record = {"created_utc": utc_now(), "purpose": "one damped retry of v3's rejected I08 -> I09 step",
              "config_sha256": cfg.sha256, "v3_config_sha256": v3.sha256,
              "config_changes": {"fit.step.alpha0": [s3["alpha0"], sr["alpha0"]],
                                 "paths.scratch_root": [v3["paths"]["scratch_root"], cfg["paths"]["scratch_root"]],
                                 "auto.max_ensemble": [v3["auto"]["max_ensemble"], cfg["auto"]["max_ensemble"]]},
              "copied_from_v3": copies, "shared_inputs_equal_v3": True,
              "shared_manifest_sha256": sha256(cfg.shared_dir / "manifest.json"),
              "v3_files_written": "none"}
    (HERE / "PROVENANCE.json").write_text(json.dumps(record, indent=1) + "\n")
    print(json.dumps(record["config_changes"], indent=1))
    print(f"copied {len(copies)} files from v3 (checksums equal); shared inputs equal v3's; wrote PROVENANCE.json")


if __name__ == "__main__":
    main()
