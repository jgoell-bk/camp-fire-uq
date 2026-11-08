"""Command line: python -m sobolfire <command> (run from the experiment folder)."""
import argparse
import sys

from .config import DEFAULT, Config
from .util import ConfigError


def prepare(cfg, write):
    from .case import prepare_shared
    from .design import load_matrix
    from .forcing import prepare_forcing

    runs = load_matrix(cfg)
    print(f"run matrix: {len(runs)} rows; checksum, order, bounds and Sobol structure OK")
    shared = prepare_shared(cfg, write)
    ign, wx = shared["ignition"], shared["weather_grid"]
    print(f"baseline case and ELMFIRE found; ignition cell row {ign['row']} col {ign['col']} "
          f"(fuel model {ign['fbfm']}) at x={ign['cell_x']:.1f} y={ign['cell_y']:.1f}")
    print(f"weather grid {wx['ncols']} x {wx['nrows']} cells of {wx['cellsize']:.0f} m covers the fuel grid")
    tables = prepare_forcing(cfg, write)
    print(f"forcing: {cfg.nbands} hourly bands from {cfg['simulation']['t0_utc']} for "
          f"{', '.join(tables)}; QA, coverage and checkpoint OK")
    if write:
        print(f"wrote {cfg.path(cfg['paths']['prepared'])} and {cfg.shared_dir}")
    else:
        print("nothing written (validate)")


def all_ids(cfg):
    from .design import load_matrix
    return [f"{p}_BASELINE" for p in cfg.stations] + [r.run_id for r in load_matrix(cfg)]


def main(argv=None):
    p = argparse.ArgumentParser(prog="python -m sobolfire", description=__doc__)
    p.add_argument("--config", default=DEFAULT, help="default: %(default)s")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("validate", help="check matrix, sources, case, ignition and grids; write nothing")
    sub.add_parser("prepare", help="validate, then write prepared/ tables and the shared inputs on scratch")
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
    c.add_argument("which", help="pilot | all | controls | spotting | JG | OP | CO | retry")
    c.add_argument("--of", help="for retry: the batch file whose failed or unfinished runs to rerun")
    c.add_argument("--name", help="batch file name (default: the selection)")
    c = sub.add_parser("status", help="count validated/failed/running/not started runs")
    c.add_argument("batch", nargs="?", help="batch file (default: all matrix runs and controls)")
    c = sub.add_parser("assemble", help="collect a station's run outputs in design_row order")
    c.add_argument("station", help="JG, OP or CO")
    c.add_argument("--mask", help="raster on the fuel grid; arrival times at nonzero cells are collected")
    c = sub.add_parser("restore-inputs", help="rebuild validated runs' weather rasters into <run>/inputs/, "
                       "checked against the checksums recorded when they ran")
    c.add_argument("batch", nargs="?", help="batch file (default: all matrix runs and controls)")
    c.add_argument("--workers", type=int, default=1)
    c = sub.add_parser("compare", help="arrival-time differences between two validated runs")
    c.add_argument("run_a")
    c.add_argument("run_b")
    args = p.parse_args(argv)

    try:
        cfg = Config(args.config)
        if args.cmd in ("validate", "prepare"):
            prepare(cfg, write=args.cmd == "prepare")
        elif args.cmd == "build":
            from pathlib import Path
            from .case import build_case
            from .design import load_matrix, resolve
            matrix = {r.run_id: r for r in load_matrix(cfg)}
            out = Path(args.out).resolve()
            if out.exists() and any(out.iterdir()):
                raise ConfigError(f"{out} is not empty")
            build_case(cfg, resolve(cfg, matrix, args.run_id), out, out / "scratch")
            print(f"built {args.run_id} in {out}")
        elif args.cmd == "run-one":
            from .design import load_matrix, resolve
            from .runner import run_one
            matrix = {r.run_id: r for r in load_matrix(cfg)}
            return 0 if run_one(cfg, resolve(cfg, matrix, args.run_id), args.force) != "failed" else 1
        elif args.cmd == "run-batch":
            from .runner import run_batch
            return run_batch(cfg, args.batch, args.index, args.workers, args.force)
        elif args.cmd == "make-batch":
            from .runner import make_batch
            make_batch(cfg, args.which, args.of, args.name)
        elif args.cmd == "status":
            from .runner import read_batch, status
            status(cfg, read_batch(args.batch) if args.batch else all_ids(cfg))
        elif args.cmd == "assemble":
            from .assemble import assemble
            assemble(cfg, args.station, args.mask)
        elif args.cmd == "restore-inputs":
            from .runner import read_batch, restore_all
            return restore_all(cfg, read_batch(args.batch) if args.batch else all_ids(cfg), args.workers)
        elif args.cmd == "compare":
            from .runner import compare
            compare(cfg, args.run_a, args.run_b)
    except ConfigError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    return 0
