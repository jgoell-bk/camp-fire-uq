"""Running cases, as in ../local_sensitivity/lsa/run.py.

Each run is built and executed in node-local scratch ($L_SCRATCH_JOB). Its records,
outputs and weather rasters are copied back to <scratch_root>/runs/<group>/<run_id>/
(group BASELINE or I00, I01, ...). A folder holds attempt.json (status
running/failed/validated) and, only after the outputs pass validation, DONE. Earlier
attempts are moved to attempts/<k>/.
"""
import csv
import math
import os
import re
import shutil
import socket
import subprocess
import time
import traceback
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

from .case import build_case, coefficients
from .common import (BASELINE, CODE_ROOT, StudyError, read_json, run_dir, run_group, sha256, utc_now,
                     write_json)

TOA_NODATA = -9999.0      # ELMFIRE writes unburned cells as NoData
KEEP = ["elmfire.data.in", "run.json", "inputs/weather.json", "outputs/fire_size_stats.csv",
        "logs/elmfire.stdout", "logs/elmfire.stderr"]


def state(directory):
    if (directory / "DONE").exists():
        return "validated"
    if (directory / "attempt.json").exists():
        return read_json(directory / "attempt.json")["status"]   # running or failed
    return "not_started"


def archive_previous(directory):
    """Move an earlier attempt's files into attempts/<k>/ so its evidence is kept."""
    entries = [p for p in directory.iterdir() if p.name != "attempts"]
    if not entries:
        return
    k = 1
    while (directory / "attempts" / str(k)).exists():
        k += 1
    dest = directory / "attempts" / str(k)
    dest.mkdir(parents=True)
    for p in entries:
        shutil.move(str(p), str(dest / p.name))


def describe(record):
    """One log line for a built run; every field is always a number."""
    k = len(record["knots_h"])
    t = record["theta"]
    mult = "/".join(f"{v:.2f}" for v in record["speed_mult_knots"])
    d = "/".join(f"{v:+.1f}" for v in t[k:2 * k])
    m = "/".join(f"{v:+.2f}" for v in t[2 * k:3 * k])
    return f"speed x{mult}, direction {d} deg, moisture {m} pp at the knots"


def run_one(cfg, rid, force=False):
    """Build, run and validate one case. Returns skipped, validated or failed."""
    final = run_dir(cfg, rid)
    if (final / "DONE").exists() and not force:
        print(f"{rid}: already validated, skipped", flush=True)
        return "skipped"
    if final.exists():
        archive_previous(final)
    final.mkdir(parents=True, exist_ok=True)
    attempt = {"run_id": rid, "status": "running", "started_utc": utc_now(), "host": socket.gethostname(),
               **{k.lower(): os.environ.get(k) for k in ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID")}}
    write_json(final / "attempt.json", attempt)
    base = os.environ.get("L_SCRATCH_JOB") or str(cfg.scratch_root / "work")
    work = Path(base) / "eki" / rid
    shutil.rmtree(work, ignore_errors=True)
    print(f"{rid}: start work={work}", flush=True)
    try:
        record = build_case(cfg, rid, work, work / "scratch")
        print(f"{rid}: {describe(record)}", flush=True)
        t0 = time.time()
        with open(work / "logs/elmfire.stdout", "w") as out, open(work / "logs/elmfire.stderr", "w") as err:
            proc = subprocess.Popen([str(cfg.path(cfg["case"]["elmfire"])), "elmfire.data.in"],
                                    cwd=work, stdout=out, stderr=err)
            _, status, usage = os.wait4(proc.pid, 0)
        code = os.waitstatus_to_exitcode(status)
        proc.returncode = code                      # reaped by wait4; tell Popen
        attempt.update(exit_code=code, wall_s=round(time.time() - t0, 1),
                       cpu_s=round(usage.ru_utime + usage.ru_stime, 1),
                       max_rss_mb=round(usage.ru_maxrss / 1024, 1))
        shutil.rmtree(work / "scratch", ignore_errors=True)
        if code != 0:
            raise RuntimeError(f"ELMFIRE exited with code {code}")
        attempt["validation"] = validate_outputs(cfg, work)
        attempt["status"] = "validated"
    except Exception as e:
        attempt.update(status="failed", error=f"{type(e).__name__}: {e}", traceback=traceback.format_exc())
    attempt["ended_utc"] = utc_now()
    try:
        copy_back(work, final)
    except Exception as e:
        attempt.update(status="failed", error=f"copy back failed: {e}")
    write_json(final / "attempt.json", attempt)
    if attempt["status"] == "validated":
        write_json(final / "DONE", {"run_id": rid, "validated_utc": attempt["ended_utc"],
                                    "wall_s": attempt["wall_s"], "max_rss_mb": attempt["max_rss_mb"],
                                    **attempt["validation"]})
    shutil.rmtree(work, ignore_errors=True)
    print(f"{rid}: {attempt['status']}" + (f" ({attempt.get('error')})" if attempt["status"] == "failed" else
          f" in {attempt['wall_s']} s, {attempt['validation']['burned_area_ac']:.0f} ac"), flush=True)
    return attempt["status"]


def copy_back(work, final):
    outputs = work / "outputs"
    toa = [str(p.relative_to(work)) for p in outputs.glob("time_of_arrival_*.tif")] if outputs.is_dir() else []
    for rel in KEEP + toa:
        if (work / rel).exists():
            shutil.copy2(work / rel, final / Path(rel).name)
    weather = sorted((work / "inputs").glob("*.tif")) if (work / "inputs").is_dir() else []
    if weather:                                   # the run's weather rasters, kept for reference
        (final / "inputs").mkdir(exist_ok=True)
        for p in weather:
            shutil.copy2(p, final / "inputs" / p.name)


def read_fire_size_stats(path):
    with open(path, newline="") as f:
        rows = [[c.strip() for c in r] for r in csv.reader(f) if r]
    if len(rows) != 2:
        raise ValueError(f"{path.name}: expected one fire, found {len(rows) - 1} rows")
    stats = dict(zip(*rows))
    return {k: float(stats[k]) for k in
            ("tstop (h)", "Wall clock time (s)", "Total fire area (ac)", "Crown fire area (ac)")}


def validate_outputs(cfg, work):
    """A zero exit code is not enough: check the final arrival-time raster (grid,
    values, coverage of SIMULATION_TSTOP, ignition at the right cell and time) and
    the fire size table."""
    s = cfg["simulation"]
    shared = read_json(cfg.shared_dir / "manifest.json")
    fuel, ign = shared["fuel_grid"], shared["ignition"]
    toa_files = sorted((work / "outputs").glob("time_of_arrival_0000001_*.tif"))
    if len(toa_files) != 1:
        raise ValueError(f"expected one final time_of_arrival raster, found {[p.name for p in toa_files]}")
    toa_path = toa_files[0]
    t_end = int(toa_path.stem.rsplit("_", 1)[1])
    if t_end < s["stop_s"]:
        raise ValueError(f"final dump at {t_end} s, before SIMULATION_TSTOP = {s['stop_s']} s")
    with rasterio.open(toa_path) as f:
        a = f.read(1)
        t = f.transform
        grid_ok = ((f.height, f.width) == (fuel["nrows"], fuel["ncols"])
                   and f.crs.to_string() == fuel["crs"] and math.isclose(t.a, fuel["cellsize"])
                   # ELMFIRE stores the corner in single precision: allow 1 m
                   and abs(t.c - fuel["x0"]) < 1 and abs(t.f - fuel["top"]) < 1)
    if not grid_ok:
        raise ValueError(f"{toa_path.name} is not on the fuel grid")
    burned = a > 0
    if not np.isfinite(a).all() or not np.isin(a[~burned], [-1.0, TOA_NODATA]).all():
        raise ValueError(f"unexpected arrival-time values: {np.unique(a[~burned])[:5]}")
    if a[burned].max(initial=0) > t_end + 1:
        raise ValueError(f"arrival time {a.max()} after the final dump {t_end}")
    # ELMFIRE marks the ignition cell itself with -1; its neighbours carry arrival times.
    if a[ign["row"], ign["col"]] != -1.0:
        raise ValueError(f"ignition cell holds {a[ign['row'], ign['col']]}, not ELMFIRE's ignition marker -1")
    first = float(a[burned].min()) if burned.any() else None
    if first is not None and first < s["ignition_s"]:
        raise ValueError(f"cells burned at {first} s, before the ignition at {s['ignition_s']} s")
    stats = read_fire_size_stats(work / "outputs/fire_size_stats.csv")
    stdout = (work / "logs/elmfire.stdout").read_text(errors="replace")
    stops = sorted({line.split("]", 1)[-1].strip() for line in stdout.splitlines() if "STOPP" in line})
    return {"toa_file": toa_path.name, "toa_sha256": sha256(toa_path), "t_end_s": t_end,
            "first_arrival_s": first, "burned_cells": int(burned.sum()),
            "burned_area_ac": stats["Total fire area (ac)"],
            "crown_area_ac": stats["Crown fire area (ac)"], "tstop_h": stats["tstop (h)"],
            "elmfire_wall_s": stats["Wall clock time (s)"], "stop_messages": stops}


# ---------------------------------------------------------------- batches

def read_batch(path):
    with open(path) as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]


def chunk(ids, size, index):
    """Array task `index` runs batch lines index*size+1 .. (index+1)*size."""
    return ids[index * size:(index + 1) * size]


def check_ids(cfg, ids):
    """Every ID must parse and have coefficients before anything runs."""
    if len(set(ids)) != len(ids):
        raise StudyError("duplicate run IDs in the batch")
    for rid in ids:
        coefficients(cfg, rid)


def run_batch(cfg, batch, index, workers, force=False):
    size = cfg["slurm"]["runs_per_task"]
    ids = chunk(read_batch(batch), size, index)
    if not ids:
        raise StudyError(f"array task {index} has no runs in {batch} (chunk size {size})")
    if not (cfg.shared_dir / "manifest.json").exists():
        raise StudyError("shared inputs not found; run `python -m eki prepare` first")
    check_ids(cfg, ids)
    print(f"array task {index}: lines {index * size + 1}-{index * size + len(ids)} of {batch}, "
          f"{len(ids)} runs, {workers} at a time", flush=True)
    with ThreadPoolExecutor(workers) as pool:
        statuses = list(pool.map(lambda rid: run_one(cfg, rid, force), ids))
    counts = Counter(statuses)
    print(f"array task {index} finished: {dict(counts)}", flush=True)
    return 1 if counts["failed"] else 0


def selection(cfg, which, of=None):
    """(run IDs, batch name).
    I00: the baseline and every iteration-0 member; I01, I02 ...: that ensemble file's
    members; CHK01 ...: the check runs in ensembles/CHK01.csv; retry: the failed or
    unfinished runs of the batch file `of`."""
    if which == "retry":
        if not of:
            raise StudyError("make-batch retry needs --of <batch file>")
        return [i for i in read_batch(of) if state(run_dir(cfg, i)) != "validated"], \
            f"retry_{utc_now().replace(':', '')}"
    if re.fullmatch(r"I\d{2}|CHK\d{2}", which):
        path = cfg.path(cfg["paths"]["ensembles"]) / f"{which}.csv"
        if not path.exists():
            raise StudyError(f"{path} not found")
        ids = list(pd.read_csv(path).run_id)
        return ([BASELINE] + ids if which == "I00" else ids), which
    raise StudyError(f"unknown selection {which!r} (I00, I01 ..., CHK01 ..., or retry)")


def make_batch(cfg, which="I00", of=None):
    """Write batches/<name>.txt and print the sbatch line."""
    ids, name = selection(cfg, which, of)
    if not ids:
        print("nothing to run")
        return None
    check_ids(cfg, ids)
    path = cfg.dir("batches") / f"{name}.txt"
    sl = cfg["slurm"]
    size = sl["runs_per_task"]
    ntasks = math.ceil(len(ids) / size)
    path.write_text(f"# {len(ids)} runs, selection {which!r}, created {utc_now()}\n" + "\n".join(ids) + "\n")
    print(f"wrote {path} ({len(ids)} runs)")
    print(f"{size} runs per array task -> tasks 0-{ntasks - 1}")
    for i in range(ntasks):
        part = chunk(ids, size, i)
        print(f"  task {i}: {part[0]} .. {part[-1]} ({len(part)} runs)")
    script = sl["check_script"] if name.startswith("CHK") else sl["script"]
    print("Submit from the eki_calibration folder (CPUs, time and memory are set in the Slurm file):")
    print(f"  cd {CODE_ROOT}")
    print(f"  sbatch --array=0-{ntasks - 1}%{sl['max_concurrent_tasks']} {script} "
          f"{os.path.relpath(path, CODE_ROOT)}")
    return path


def status(cfg, ids):
    counts, stalled, unfinished = Counter(), [], []
    groups = []
    for rid in ids:
        group = run_group(rid)
        if group not in groups:
            groups.append(group)
        d = run_dir(cfg, rid)
        st = state(d)
        counts[(group, st)] += 1
        if st in ("failed", "running"):
            unfinished.append(f"{rid} ({st})")
        if st == "validated" and any("STALLED" in m for m in read_json(d / "DONE")["stop_messages"]):
            stalled.append(rid)
    states = ["validated", "failed", "running", "not_started"]
    print(f"{'group':10s}" + "".join(f"{s:>13s}" for s in states))
    for g in groups:
        print(f"{g:10s}" + "".join(f"{counts[(g, s)]:13d}" for s in states))
    print("'running' includes runs whose job ended before they finished.")
    if stalled:
        print(f"{len(stalled)} validated runs stopped early (fire front stalled), e.g. {stalled[:5]}")
    if unfinished:
        print(f"failed or unfinished: {unfinished[:10]}{' ...' if len(unfinished) > 10 else ''}")
    return counts
