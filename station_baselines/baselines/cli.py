"""python -m baselines {prepare, build, run-one, run-all, status, qoi}"""
import argparse
import sys
from pathlib import Path

from .common import DEFAULT_CONFIG, Config, StudyError, write_text


def prepare(cfg):
    from .case import prepare_shared
    from .forcing import check_reference, weather_table
    table = weather_table(cfg)
    ref = check_reference(cfg, table)
    path = cfg.dir("prepared") / "weather.csv"
    write_text(path, table.to_csv(float_format="%.17g"))
    print(f"wrote {path}: {len(table)} rows ({', '.join(cfg.run_ids)}); MEAN matches the "
          f"groupmate's file to {max(ref.values()):.1e}")
    manifest = prepare_shared(cfg)
    print(f"wrote {cfg.shared_dir} ({len(manifest['files'])} files); ignition cell "
          f"row {manifest['ignition']['row']}, col {manifest['ignition']['col']}, fuel {manifest['ignition']['fbfm']}")


def main(argv=None):
    p = argparse.ArgumentParser(prog="python -m baselines", description=__doc__)
    p.add_argument("--config", default=str(DEFAULT_CONFIG))
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("prepare", help="write prepared/weather.csv and <scratch_root>/shared/")
    b = sub.add_parser("build", help="build one case into a folder without running it")
    b.add_argument("run_id")
    b.add_argument("--out", required=True)
    o = sub.add_parser("run-one", help="build, run and validate one run in this session")
    o.add_argument("run_id")
    o.add_argument("--force", action="store_true")
    a = sub.add_parser("run-all", help="every run, WORKERS at a time (what the Slurm script calls)")
    a.add_argument("--workers", type=int, default=3)
    a.add_argument("--force", action="store_true")
    sub.add_parser("status")
    sub.add_parser("qoi", help="border-band arrival summary of the validated runs")
    args = p.parse_args(argv)
    cfg = Config(args.config)
    try:
        if args.cmd == "prepare":
            prepare(cfg)
        elif args.cmd == "build":
            from .case import build_case
            out = Path(args.out).resolve()
            cfg.run_dir(args.run_id)                  # checks the ID
            build_case(cfg, args.run_id, out, out / "scratch")
            print(f"built {args.run_id} in {out}")
        elif args.cmd == "run-one":
            from .run import run_one
            return 0 if run_one(cfg, args.run_id, args.force) != "failed" else 1
        elif args.cmd == "run-all":
            from .run import run_all
            return run_all(cfg, args.workers, args.force)
        elif args.cmd == "status":
            from .run import status
            status(cfg)
        elif args.cmd == "qoi":
            from .qoi import qoi
            qoi(cfg)
    except StudyError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    return 0
