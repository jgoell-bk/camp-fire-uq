"""The whole EKI loop in one job: python -m eki auto (run by slurm/eki_auto.slurm).

Starting from the newest ensemble file (ensembles/I<n>.csv), each pass does what the
manual commands do, in the same order and with the same files:
1. run I<n>'s members here, `workers` at a time (validated runs are skipped);
2. `extract` and `pilot` it, and log the members' objective;
3. stop if I<n> is auto.max_ensemble, or if the update that made I<n> lowered Phi at
   the mean by less than auto.min_relative_gain (read from results/CHK<n>_verify.json);
4. `propose` I<n+1> with its check runs CHK<n+1>, run them, `verify`;
5. stop if verify does not accept the update; otherwise `promote` I<n+1> and go to 1.

It also stops at the first failed run or error and never retries on its own. Every
step is appended to results/auto_log.json. Resubmitting the job resumes: validated
runs are skipped, and an existing proposal or promotion is reused, not redone.
"""
import os
import re
import traceback
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

from . import cli
from .common import StudyError, read_json, run_dir, sha256, utc_now, write_json
from .controlled import WINDOWS, fit_constraints
from .corrections import theta_names
from .extract import extract
from .objective import objective
from .run import check_ids, make_batch, read_batch, run_one, selection, state


def latest_ensemble(cfg):
    names = [p.stem for p in cfg.path(cfg["paths"]["ensembles"]).glob("I*.csv")]
    ns = [int(s[1:]) for s in names if re.fullmatch(r"I\d{2}", s)]
    if not ns:
        raise StudyError("no ensembles/I00.csv; run `python -m eki ensemble` first")
    return max(ns)


def update_gain(cfg, n):
    """Relative decrease of Phi at the mean from the update that made I<n>, measured by
    CHK<n> (None for I00 or if CHK<n> was never verified)."""
    path = cfg.dir("results") / f"CHK{n:02d}_verify.json"
    if n == 0 or not path.exists():
        return None
    v = read_json(path)
    current = v["objective"][f"CHK{n:02d}_I{n - 1:02d}MEAN"]["total"]
    return -v["change_proposed_minus_current"] / current


def decide(n, top, gain, min_gain):
    """After I<n> has been run and scored: None to continue, or the reason to stop."""
    if n >= top:
        return f"reached I{top:02d}, the last ensemble allowed (auto.max_ensemble)"
    if gain is not None and gain < min_gain:
        return (f"the update that made I{n:02d} lowered Phi at the mean by {gain:.1%}, "
                f"less than {min_gain:.0%} (auto.min_relative_gain)")
    return None


def run_all(cfg, name, workers):
    """Run every run of ensemble or check batch `name` here. Returns the failed IDs."""
    ids, _ = selection(cfg, name)
    check_ids(cfg, ids)
    batch = cfg.dir("batches") / f"{name}.txt"
    if not batch.exists() or read_batch(batch) != ids:
        make_batch(cfg, name, quiet=True)
    todo = [r for r in ids if state(run_dir(cfg, r)) != "validated"]
    print(f"{name}: {len(ids)} runs, {len(ids) - len(todo)} already validated, {workers} at a time", flush=True)
    with ThreadPoolExecutor(workers) as pool:
        statuses = list(pool.map(lambda r: run_one(cfg, r), ids))
    return batch, [r for r, s in zip(ids, statuses) if s == "failed"]


def member_summary(cfg, n, c, w):
    """Members' Phi (from results/I<n>_distance.npz) and, per deadline window, the
    members' median share of deadlines met and median mean distance."""
    name = f"I{n:02d}"
    D = np.load(cfg.dir("results") / f"{name}_distance.npz")
    ids = list(D["run_id"])
    ens = pd.read_csv(cfg.path(cfg["paths"]["ensembles"]) / f"{name}.csv")
    rows = [ids.index(r) for r in ens.run_id]
    d = D["distance_m"][rows]
    acres = D["overspread_acres"][rows] if cfg["fit"].get("overspread") else None
    phi = objective(cfg, ens[theta_names(cfg)].to_numpy(float), d, w, acres)["total"]
    u_h = c.deadline_s.to_numpy() / 3600
    windows = {}
    for a, b in WINDOWS:
        k = (u_h >= a) & (u_h < b)
        windows[f"{a}-{b:g}h"] = {"met_share": round(float(np.median((d[:, k] == 0).mean(1))), 3),
                                  "mean_km": round(float(np.median(d[:, k].mean(1)) / 1000), 2)}
    out = {"members_phi": {"mean": float(phi.mean()), "median": float(np.median(phi)),
                           "min": float(phi.min()), "max": float(phi.max())},
           "windows_members_median": windows}
    if acres is not None:
        out["overspread_acres_members_median"] = float(np.median(acres))
    if "p" in ens.columns:
        out["p_eps_members_median"] = float(np.median(ens.p_eps_value))
    return out


def auto(cfg, workers):
    a = cfg["auto"]
    top, min_gain = int(a["max_ensemble"]), float(a["min_relative_gain"])
    log_path = cfg.dir("results") / "auto_log.json"
    log = read_json(log_path) if log_path.exists() else {"jobs": []}
    job = {"started_utc": utc_now(), "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
           "host": os.uname().nodename, "settings": a, "steps": []}
    log["jobs"].append(job)

    def note(**step):
        job["steps"].append({"utc": utc_now(), **step})
        write_json(log_path, log)

    def stop(reason, code):
        job.update(stopped=reason, exit_code=code, ended_utc=utc_now())
        write_json(log_path, log)
        print(f"\nauto stopped: {reason}\nlog: {log_path}", flush=True)
        summary(log_path)
        return code

    try:
        n = latest_ensemble(cfg)
        if n > top:
            return stop(f"the newest ensemble, I{n:02d}, is past auto.max_ensemble (I{top:02d})", 0)
        c, _, w = fit_constraints(cfg)
        while True:
            name = f"I{n:02d}"
            print(f"\n===== {name}: run, extract, pilot ({utc_now()})", flush=True)
            batch, failed = run_all(cfg, name, workers)
            if failed:
                return stop(f"{len(failed)} runs of {name} failed (first: {failed[:5]}); see their attempt.json. "
                            "Find the cause before resubmitting", 1)
            extract(cfg, batch)
            cli.pilot(cfg, batch)
            gain = update_gain(cfg, n)
            note(step="ensemble", ensemble=name, update_gain=gain, **member_summary(cfg, n, c, w))
            reason = decide(n, top, gain, min_gain)
            if reason:
                return stop(reason, 0)

            nxt, chk = f"I{n + 1:02d}", f"CHK{n + 1:02d}"
            print(f"\n===== {name} -> {nxt}: propose, check runs, verify ({utc_now()})", flush=True)
            prop = cfg.path("proposals") / f"{nxt}.csv"
            cpath = cfg.path(cfg["paths"]["ensembles"]) / f"{chk}.csv"
            if prop.exists() and cpath.exists():
                print(f"{prop.name} and {cpath.name} exist (an earlier job proposed them): reused", flush=True)
            elif prop.exists() or cpath.exists():
                raise StudyError(f"only one of {prop} and {cpath} exists; check them by hand")
            else:
                cli.propose(cfg, name, quiet=True)
            _, failed = run_all(cfg, chk, workers)
            if failed:
                return stop(f"check runs failed: {failed}; see their attempt.json", 1)
            cli.verify(cfg, chk)
            v = read_json(cfg.dir("results") / f"{chk}_verify.json")
            current = v["objective"][f"{chk}_{name}MEAN"]["total"]
            proposed = v["objective"][f"{chk}_{nxt}MEAN"]["total"]
            note(step="check", check=chk, phi_current_mean=current, phi_proposed_mean=proposed,
                 relative_gain=(current - proposed) / current, verdict=v["verdict"], accepted=v["accepted"],
                 repeat_bit_identical=v["repeat_bit_identical"],
                 alpha=read_json(prop.with_suffix(".json"))["parent"]["step"]["alpha"])
            if not v["accepted"]:
                return stop(f"{chk}: the proposed mean {v['verdict']} Phi ({proposed - current:+.4f}); "
                            f"{name} is the final ensemble", 0)
            dest = cfg.path(cfg["paths"]["ensembles"]) / f"{nxt}.csv"
            if dest.exists():
                if sha256(dest) != sha256(prop):
                    raise StudyError(f"{dest} exists and differs from {prop}")
                print(f"{dest.name} was already promoted", flush=True)
            else:
                cli.promote(cfg, nxt, quiet=True)
            n += 1
    except Exception as e:                       # record why, then stop: nothing retries on its own
        job["traceback"] = traceback.format_exc()
        return stop(f"error: {type(e).__name__}: {e}", 1)


def summary(log_path):
    """Every logged ensemble and check, oldest first."""
    log = read_json(log_path)
    print(f"{'step':9s}{'name':7s}{'Phi':>26s}  details")
    for job in log["jobs"]:
        for s in job["steps"]:
            if s["step"] == "ensemble":
                p = s["members_phi"]
                phi = f"members mean {p['mean']:.2f} min {p['min']:.2f}"
                win = "  ".join(f"{k} {x['met_share']:.2f} met {x['mean_km']:.1f} km"
                                for k, x in s["windows_members_median"].items())
                more = "".join(f"  {label} {s[key]:{fmt}}" for key, label, fmt in
                               (("overspread_acres_members_median", "outside", ",.0f"),
                                ("p_eps_members_median", "P_EPS", ".4g")) if key in s)
                print(f"{'ensemble':9s}{s['ensemble']:7s}{phi:>26s}  {win}{more}")
            else:
                phi = f"mean {s['phi_current_mean']:.2f} -> {s['phi_proposed_mean']:.2f}"
                print(f"{'check':9s}{s['check']:7s}{phi:>26s}  gain {s['relative_gain']:.1%}, "
                      f"alpha {s['alpha']:g}, {s['verdict']}")
