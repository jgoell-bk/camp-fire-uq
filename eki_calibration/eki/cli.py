"""Command line: python -m eki <command> (run from the eki_calibration folder)."""
import argparse
import sys
from pathlib import Path

from .common import DEFAULT_CONFIG, Config, StudyError, sha256, utc_now, write_json, write_text


def prepare(cfg):
    """Shared inputs on scratch; prepared/baseline.csv, checked against the reference MEAN."""
    from .case import prepare_shared
    from .corrections import check_knots
    from .forcing import baseline_table, check_reference

    check_knots(cfg)
    shared = prepare_shared(cfg)
    ign, wx = shared["ignition"], shared["weather_grid"]
    print(f"baseline case and ELMFIRE found; ignition cell row {ign['row']} col {ign['col']} "
          f"(fuel model {ign['fbfm']}), ignition at t = {cfg['simulation']['ignition_s']} s")
    print(f"weather grid {wx['ncols']} x {wx['nrows']} cells of {wx['cellsize']:.0f} m covers the fuel grid")
    print(f"wrote {cfg.shared_dir}")
    table = baseline_table(cfg)
    diff = check_reference(cfg, table)
    out = cfg.dir("prepared")
    path = out / "baseline.csv"
    write_text(path, table.to_csv(date_format="%Y-%m-%dT%H:%M:%SZ", float_format="%.17g"))
    src = cfg["sources"]
    write_json(out / "baseline.json", {
        "created_utc": utc_now(), "sha256": sha256(path), "stations": cfg.stations,
        "rule": "equal-weight mean of the stations' (u, v) wind vectors; Jarbo Gap NFDRS dead fuel moisture",
        "reference_check": diff,
        "sources_sha256": {k: sha256(cfg.path(src[k])) for k in ("wind_csv", "moisture_csv", "reference_mean_csv")}})
    print(f"wrote {path} ({cfg.nbands} hourly bands from {cfg['simulation']['t0_utc']}); matches "
          f"reference_mean_csv within {diff['max_abs_ws_diff_mph']:.1e} mph and {diff['max_abs_wd_diff_deg']:.1e} deg")
    print(f"  mean wind {table.ws_mph.min():.1f}-{table.ws_mph.max():.1f} mph; "
          f"1-h moisture {table.m1_pct.min():.2f}-{table.m1_pct.max():.2f} %")


def observations(cfg):
    from .observations import write_table
    df, meta = write_table(cfg)
    c = meta["counts"]
    p, s = meta["column_prefixes"]["primary"], meta["column_prefixes"]["sensitivity"]
    print(f"wrote {cfg.path(cfg['paths']['prepared']) / 'observations.csv'}: {len(df)} reports, "
          f"{c['in_window']} in the model window, {c['exact_duplicates']} exact duplicates")
    print(f"  containing cell: {c['fuel_class']}")
    print(f"  primary ({meta['radii_m']['primary']:g} m): {c[p + 'status']}")
    print(f"  sensitivity ({meta['radii_m']['sensitivity']:g} m): {c[s + 'status']}")
    review = df[df[p + "status"] == "review"]
    for _, r in review.iterrows():
        print(f"  review {r.obs_id}: {r.lon:.6f}, {r.lat:.6f} at {r.observed_time_pst}; {r[p + 'reason']}")


def holdout(cfg):
    from .holdout import write_split
    part, info, cov = write_split(cfg)
    v = cfg["holdout"]["version"]
    print(f"wrote {cfg.path(cfg['paths']['prepared']) / f'holdout_v{v}.csv'}: {info['blocks']} occupied 1 km blocks "
          f"in {info['units']} units; held out {info['holdout_units']} units = {info['holdout_blocks']} blocks "
          f"({info['holdout_blocks'] / info['blocks']:.0%}); links spanning blocks: {info['links_spanning_blocks']}")
    keys = list(cov["fit"])
    print(f"  {'':32s}{'fit':>14s}{'holdout':>14s}")
    for k in keys:
        f, h = cov["fit"][k], cov["holdout"][k]
        print(f"  {k:32s}{str(f):>14s}{str(h):>14s}")


def constraints(cfg):
    """prepared/constraints.csv: every report's deadline, selection and reason."""
    from .holdout import load_split
    from .observations import load_table
    from .target import selection
    obs = load_table(cfg)
    s = selection(cfg, obs)
    s.insert(1, "partition", load_split(cfg).partition)
    path = cfg.dir("prepared") / "constraints.csv"
    write_text(path, s.to_csv(index=False, float_format="%.10g"))
    counts = s.groupby(["selection_rule", "partition"]).size().unstack(fill_value=0)
    write_json(path.with_suffix(".json"), {"created_utc": utc_now(), "sha256": sha256(path),
                                           "config_target": cfg["target"], "counts": counts.to_dict()})
    print(f"wrote {path}")
    print(counts.to_string())


def pilot(cfg, batch):
    from .pilot import pilot as run_pilot
    res, tab = run_pilot(cfg, batch)
    print(f"fit constraints: {res['fit_constraints']}; cells reached by all / some / no members: "
          f"{res['cells']['reached_by_all']} / {res['cells']['reached_by_some']} / {res['cells']['reached_by_none']}; "
          f"baseline reaches {res['cells']['baseline_reached']}")
    print(f"{'run':10s}{'reached':>9s}{'not':>6s}{'met':>6s}{'late':>6s}  lateness h p10/p50/p90 (reached)")
    for r in res["runs"]:
        print(f"{r['run_id']:10s}{r['reached']:9d}{r['not_reached']:6d}{r['met']:6d}{r['violated']:6d}  "
              f"{r['lateness_h_p10_p50_p90_reached']}")
    print(f"wrote results/{batch.stem}_pilot.json, _pilot_cells.csv, _distance.npz")


def propose(cfg, src):
    from .controlled import propose as run_propose
    from .run import make_batch
    it = int(src[1:])
    out, cpath, rec = run_propose(cfg, it)
    s = rec["step"]
    print(f"proposal {out} (not an ensemble file yet): {rec['constraints']} fit constraints in {rec['blocks']} "
          f"blocks (per block min/median/max {rec['constraints_per_block_min_median_max']}), weights sum "
          f"{rec['weights_sum']:.6f}; L = {rec['length_scale_m']:g} m")
    for t in s["tried"]:
        print(f"  alpha {t['alpha']:>7g}: max mean move {t['max_mean_move_sd']:.3f} sd, max member move "
              f"{t['max_member_move_sd']:.3f} sd, inadmissible members {len(t['inadmissible_members'])}")
    print(f"  chosen alpha {s['alpha']:g}; collapsed coefficients {s['collapsed_coefficients'] or 'none'}; "
          f"boundary projections {s['boundary_projections']}")
    o = rec["current_objective"]
    print(f"  current members' Phi: mean {o['total']['mean']:.3f} (distance {o['distance']['mean']:.3f}, amplitude "
          f"{o['amplitude']['mean']:.3f}, smoothness {o['smoothness']['mean']:.3f}); min {o['total']['min']:.3f}")
    if rec["baseline_objective"]:
        print(f"  BASELINE Phi {rec['baseline_objective']['total']:.3f}")
    print(f"  {'coefficient':12s}{'current mean':>14s}{'proposed mean':>15s}{'mean move sd':>14s}{'max member sd':>15s}")
    for n in rec["current_mean"]:
        print(f"  {n:12s}{rec['current_mean'][n]:14.4f}{rec['proposed_mean'][n]:15.4f}"
              f"{s['mean_move_sd'][n]:14.3f}{s['max_member_move_sd'][n]:15.3f}")
    print(f"wrote {cpath} (check runs at the means); the check batch:")
    make_batch(cfg, cpath.stem)


def verify(cfg, chk):
    from .controlled import verify as run_verify
    v = run_verify(cfg, chk)
    for r, o in v["objective"].items():
        print(f"  {r:22s} Phi {o['total']:.4f} = distance {o['distance']:.4f} + amplitude {o['amplitude']:.4f} "
              f"+ smoothness {o['smoothness']:.4f}")
    print(f"  repeat bit-identical: {v['repeat_bit_identical']}; variability {v['variability']:.4g}")
    print(f"  proposed - current: {v['change_proposed_minus_current']:+.4f} -> {v['verdict']}; "
          f"accepted for promotion: {v['accepted']}")
    print(f"wrote results/{chk}_verify.json")


def promote(cfg, nxt):
    from .controlled import promote as run_promote
    from .run import make_batch
    dest = run_promote(cfg, int(nxt[1:]))
    print(f"wrote {dest}")
    make_batch(cfg, nxt)


def ensemble(cfg):
    from .ensemble import make_initial
    table, meta = make_initial(cfg)
    r = meta["rank"]
    print(f"wrote {cfg.path(cfg['paths']['ensembles']) / 'I00.csv'}: {r['members']} members x "
          f"{r['parameters']} coefficients, rank {r['rank']}, condition number {r['condition_number']:.2f} "
          "(centred, each column scaled by its standard deviation)")
    for k, v in meta["correlations"].items():
        print(f"  {k}: normal scores {v['normal_scores']:+.2f} (target {v['target_normal_scores']:.2f}), "
              f"coefficients {v['coefficients']:+.2f}")
    w = meta["minimum_written"]
    print(f"  lowest written values: wind {w['ws_mph']:.2f} mph, moisture 1-h {w['m1_pct']:.2f} %, "
          f"10-h {w['m10_pct']:.2f} %, 100-h {w['m100_pct']:.2f} % (floor {cfg['corrections']['moisture_min_pct']} %)")


def main(argv=None):
    p = argparse.ArgumentParser(prog="python -m eki", description=__doc__)
    p.add_argument("--config", default=DEFAULT_CONFIG, help="default: %(default)s")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("prepare", help="check the case; write shared inputs (scratch) and prepared/baseline.csv")
    sub.add_parser("observations", help="write prepared/observations.csv (times, cells, frozen mapping)")
    sub.add_parser("holdout", help="write prepared/holdout_v<version>.csv, the frozen spatial-block split")
    sub.add_parser("constraints", help="write prepared/constraints.csv (deadlines, primary constraint per cell, reasons)")
    c = sub.add_parser("pilot", help="pilot diagnostics of an extracted batch (fit constraints only)")
    c.add_argument("batch")
    c = sub.add_parser("propose", help="controlled update from an ensemble: proposals/I<n+1>.csv and CHK check runs")
    c.add_argument("ensemble", help="the source ensemble, e.g. I00 (needs results/I00_distance.npz)")
    c = sub.add_parser("verify", help="actual objective at the current and proposed means (CHK runs)")
    c.add_argument("check", help="e.g. CHK01")
    c = sub.add_parser("promote", help="copy a verified proposal to ensembles/ and write its batch")
    c.add_argument("ensemble", help="e.g. I01")
    sub.add_parser("ensemble", help="write ensembles/I00.csv, the initial ensemble (checks rank, admissibility)")
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
    c.add_argument("which", nargs="?", default="I00", help="I00 (default: baseline + 32 members) | I01 ... | CHK01 ... | retry")
    c.add_argument("--of", help="for retry: the batch file whose failed or unfinished runs to rerun")
    c = sub.add_parser("status", help="count validated/failed/running/not started runs of a batch")
    c.add_argument("batch")
    c = sub.add_parser("extract", help="arrival times at the observation cells -> results/<batch>_G.npz")
    c.add_argument("batch")
    args = p.parse_args(argv)

    try:
        cfg = Config(args.config)
        if args.cmd == "prepare":
            prepare(cfg)
        elif args.cmd == "observations":
            observations(cfg)
        elif args.cmd == "holdout":
            holdout(cfg)
        elif args.cmd == "constraints":
            constraints(cfg)
        elif args.cmd == "pilot":
            pilot(cfg, Path(args.batch))
        elif args.cmd == "propose":
            propose(cfg, args.ensemble)
        elif args.cmd == "verify":
            verify(cfg, args.check)
        elif args.cmd == "promote":
            promote(cfg, args.ensemble)
        elif args.cmd == "ensemble":
            ensemble(cfg)
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
            status(cfg, read_batch(args.batch))
        elif args.cmd == "extract":
            from .extract import extract
            extract(cfg, Path(args.batch))
    except StudyError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    return 0
