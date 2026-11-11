"""Running the runs, as in ../eki_calibration/eki/run.py.

Each run is built and executed in node-local scratch ($L_SCRATCH_JOB). Its records,
outputs and weather rasters are copied back to <scratch_root>/runs/<run_id>/. A folder
holds attempt.json (status running/failed/validated) and, only after the outputs pass
validation, DONE. Earlier attempts are moved to attempts/<k>/.
"""
import csv
import math
import os
import shutil
import socket
import subprocess
import time
import traceback
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import rasterio

from .case import build_case
from .common import StudyError, read_json, sha256, utc_now, write_json

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


def run_one(cfg, rid, force=False):
    """Build, run and validate one case. Returns skipped, validated or failed."""
    final = cfg.run_dir(rid)
    if (final / "DONE").exists() and not force:
        print(f"{rid}: already validated, skipped", flush=True)
        return "skipped"
    if final.exists():
        archive_previous(final)
    final.mkdir(parents=True, exist_ok=True)
    attempt = {"run_id": rid, "status": "running", "started_utc": utc_now(), "host": socket.gethostname(),
               "slurm_job_id": os.environ.get("SLURM_JOB_ID")}
    write_json(final / "attempt.json", attempt)
    base = os.environ.get("L_SCRATCH_JOB") or str(cfg.scratch_root / "work")
    work = Path(base) / "station_baselines" / rid
    shutil.rmtree(work, ignore_errors=True)
    print(f"{rid}: start work={work}, wind from {cfg['study']['runs'][rid]}", flush=True)
    try:
        build_case(cfg, rid, work, work / "scratch")
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


def run_all(cfg, workers, force=False):
    """Every run in the config, `workers` at a time. Exit code 1 if any failed."""
    if not (cfg.shared_dir / "manifest.json").exists():
        raise StudyError("shared inputs not found; run `python -m baselines prepare` first")
    ids = cfg.run_ids
    print(f"{len(ids)} runs ({', '.join(ids)}), {workers} at a time", flush=True)
    with ThreadPoolExecutor(workers) as pool:
        statuses = list(pool.map(lambda rid: run_one(cfg, rid, force), ids))
    counts = Counter(statuses)
    print(f"finished: {dict(counts)}", flush=True)
    return 1 if counts["failed"] else 0


def status(cfg):
    for rid in cfg.run_ids:
        d = cfg.run_dir(rid)
        st = state(d)
        extra = ""
        if st == "validated":
            done = read_json(d / "DONE")
            extra = f"  {done['wall_s']:.0f} s, {done['burned_area_ac']:.0f} ac"
        elif st == "failed":
            extra = f"  {read_json(d / 'attempt.json').get('error')}"
        print(f"{rid:6s} {st}{extra}")
