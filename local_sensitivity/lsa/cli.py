"""Command line: python -m lsa <command> (run from the local_sensitivity folder)."""
import argparse
import sys
from pathlib import Path

from .common import DEFAULT_CONFIG, Config, StudyError, all_run_ids, sha256, utc_now, write_json, write_text


def prepare(cfg):
    """Shared inputs on scratch, forcing tables and aim.json in prepared/."""
    from .case import prepare_shared
    from .forcing import aim, forcing_tables

    shared = prepare_shared(cfg)
    ign, wx = shared["ignition"], shared["weather_grid"]
    print(f"baseline case and ELMFIRE found; ignition cell row {ign['row']} col {ign['col']} "
          f"(fuel model {ign['fbfm']})")
    print(f"weather grid {wx['ncols']} x {wx['nrows']} cells of {wx['cellsize']:.0f} m covers the fuel grid")
    print(f"wrote {cfg.shared_dir}")
    tables = forcing_tables(cfg)
    out = cfg.dir("prepared")
    files = {}
    for prefix, t in tables.items():
        path = out / f"forcing_{prefix}.csv"
        write_text(path, t.to_csv(date_format="%Y-%m-%dT%H:%M:%SZ"))
        files[prefix] = sha256(path)
    a = aim(cfg, tables, shared)
    src = cfg["sources"]
    write_json(out / "aim.json", {**a, "created_utc": utc_now(), "forcing_sha256": files,
                                  "sources_sha256": {k: sha256(cfg.path(src[k])) for k in ("wind_csv", "moisture_csv")}})
    print(f"wrote {out}: forcing_<station>.csv ({cfg.nbands} hourly bands from "
          f"{cfg['simulation']['t0_utc']}) and aim.json")
    print(f"target: from {a['target_from_deg']:.2f} deg (blowing toward {a['target_toward_deg']:.2f}), "
          f"{a['target_distance_m'] / 1000:.2f} km from the ignition; mean = {a['mean']}")
    if cfg["study"]["direction"] == "fixed":
        first, last = ((a["target_from_deg"] + d) % 360 for d in (cfg.offsets[0], cfg.offsets[-1]))
        print(f"  direction fixed: every hour of run WD<offset> blows from target + offset, "
              f"{first:.2f} .. {last:.2f} deg; speed and moisture hourly")
        return
    for prefix, s in a["stations"].items():
        print(f"  {prefix}: mean wind from {s['mean_from_deg']:7.2f} deg -> aim {s['aim_deg']:+7.2f} deg; "
              f"rotations {s['aim_deg'] + cfg.offsets[0]:+.2f} .. {s['aim_deg'] + cfg.offsets[-1]:+.2f}")


def main(argv=None):
    p = argparse.ArgumentParser(prog="python -m lsa", description=__doc__)
    p.add_argument("--config", default=DEFAULT_CONFIG, help="default: %(default)s")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("prepare", help="check the case; write shared inputs (scratch), forcing tables and aim.json")
    c = sub.add_parser("build", help="build one case without running it, for inspection")
    c.add_argument("run_id")
    c.add_argument("--out", required=True, help="empty folder to build into")
    c = sub.add_parser("run-one", help="build, run and validate one case here")
    c.add_argument("run_id")
    c.add_argument("--force", action="store_true", help="rerun even if validated (old attempt is archived)")
    c = sub.add_parser("run-batch", help="run one array task's chunk of a batch file (called by Slurm)")
    c.add_argument("batch")
    c.add_argument("--index", type=int, required=True, help="array task id; chunk size is slurm.runs_per_task")
    c.add_argument("--workers", type=int, required=True)
    c.add_argument("--force", action="store_true")
    c = sub.add_parser("make-batch", help="write a batch file and print the sbatch command")
    c.add_argument("which", nargs="?", default="all", help="all (default) | retry")
    c.add_argument("--of", help="for retry: the batch file whose failed or unfinished runs to rerun")
    c = sub.add_parser("status", help="count validated/failed/running/not started runs")
    c.add_argument("batch", nargs="?", help="batch file (default: every run of the study)")
    c = sub.add_parser("collect", help="one station's outputs in offset order -> results/")
    c.add_argument("station", help="JG, OP or CO")
    c.add_argument("--mask", help="raster on the fuel grid; arrival times at nonzero cells are collected")
    args = p.parse_args(argv)

    try:
        cfg = Config(args.config)
        if args.cmd == "prepare":
            prepare(cfg)
        elif args.cmd == "build":
            from .case import build_case
            out = Path(args.out).resolve()
            if out.exists() and any(out.iterdir()):
                raise StudyError(f"{out} is not empty")
            build_case(cfg, args.run_id, out, out / "scratch")
            print(f"built {args.run_id} in {out}")
        elif args.cmd == "run-one":
            from .run import run_one
            return 0 if run_one(cfg, args.run_id, args.force) != "failed" else 1
        elif args.cmd == "run-batch":
            from .run import run_batch
            return run_batch(cfg, args.batch, args.index, args.workers, args.force)
        elif args.cmd == "make-batch":
            from .run import make_batch
            make_batch(cfg, args.which, args.of)
        elif args.cmd == "status":
            from .run import read_batch, status
            status(cfg, read_batch(args.batch) if args.batch else all_run_ids(cfg))
        elif args.cmd == "collect":
            from .collect import collect
            collect(cfg, args.station, args.mask)
    except StudyError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    return 0
