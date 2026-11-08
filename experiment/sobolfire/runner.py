"""Running cases (handoff section 7).

Each run is built and executed in node-local scratch ($L_SCRATCH_JOB). Its records,
outputs and weather rasters (inputs/) are copied back to its folder on $SCRATCH:

  runs/<prefix>/b<base_row>/<run_id>/   matrix rows, 12 per base-row folder
  controls/<run_id>/                    nominal controls and repeats of them

A folder holds attempt.json (status running/failed/validated) and, only after the
outputs pass validation, DONE. Earlier attempts are moved to attempts/<k>/.
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
from .design import REP, derived_settings, level_text, load_matrix, resolve
from .util import ConfigError, read_json, sha256, utc_now, write_json

TOA_NODATA = -9999.0      # ELMFIRE writes unburned cells as NoData (and rarely -1)
KEEP = ["elmfire.data.in", "run.json", "inputs/weather.json", "outputs/fire_size_stats.csv",
        "logs/elmfire.stdout", "logs/elmfire.stderr"]


def final_dir(cfg, run):
    if run.is_control:
        return cfg.scratch_root / "controls" / run.run_id
    return cfg.scratch_root / "runs" / run.prefix / f"b{run.base_row:03d}" / run.run_id


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


def slurm_ids():
    keys = ["SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID"]
    return {k.lower(): os.environ.get(k) for k in keys}


def run_one(cfg, run, force=False):
    """Build, run and validate one case. Returns skipped, validated or failed."""
    final = final_dir(cfg, run)
    if (final / "DONE").exists() and not force:
        print(f"{run.run_id}: already validated, skipped", flush=True)
        return "skipped"
    if final.exists():
        archive_previous(final)
    final.mkdir(parents=True, exist_ok=True)
    attempt = {"run_id": run.run_id, "status": "running", "started_utc": utc_now(),
               "host": socket.gethostname(), **slurm_ids()}
    write_json(final / "attempt.json", attempt)
    base = os.environ.get("L_SCRATCH_JOB") or str(cfg.scratch_root / "work")
    work = Path(base) / "camp_sobol" / run.run_id
    shutil.rmtree(work, ignore_errors=True)
    print(f"{run.run_id}: start  station={cfg.stations[run.prefix][0]} design_row={run.design_row} "
          f"block={run.sobol_block} coeffs={run.coeffs} derived={derived_settings(cfg, run.coeffs)} "
          f"overrides={run.overrides} work={work}", flush=True)
    try:
        build_case(cfg, run, work, work / "scratch")
        exe = cfg.path(cfg["case"]["elmfire"])
        t0 = time.time()
        with open(work / "logs/elmfire.stdout", "w") as out, open(work / "logs/elmfire.stderr", "w") as err:
            proc = subprocess.Popen([str(exe), "elmfire.data.in"], cwd=work, stdout=out, stderr=err)
            _, status, usage = os.wait4(proc.pid, 0)
            proc.returncode = os.waitstatus_to_exitcode(status)
        attempt.update(exit_code=proc.returncode, wall_s=round(time.time() - t0, 1),
                       cpu_s=round(usage.ru_utime + usage.ru_stime, 1),
                       max_rss_mb=round(usage.ru_maxrss / 1024, 1))
        shutil.rmtree(work / "scratch", ignore_errors=True)
        if proc.returncode != 0:
            raise RuntimeError(f"ELMFIRE exited with code {proc.returncode}")
        attempt["validation"] = validate_outputs(cfg, work)
        attempt["status"] = "validated"
    except Exception as e:
        attempt.update(status="failed", error=f"{type(e).__name__}: {e}",
                       traceback=traceback.format_exc())
    attempt["ended_utc"] = utc_now()
    try:
        copy_back(work, final)
    except Exception as e:
        attempt.update(status="failed", error=f"copy back failed: {e}")
    write_json(final / "attempt.json", attempt)
    if attempt["status"] == "validated":
        write_json(final / "DONE", {"run_id": run.run_id, "validated_utc": attempt["ended_utc"],
                                    "wall_s": attempt["wall_s"], "max_rss_mb": attempt["max_rss_mb"],
                                    **attempt["validation"]})
    shutil.rmtree(work, ignore_errors=True)
    print(f"{run.run_id}: {attempt['status']}" + (f" ({attempt.get('error')})" if attempt["status"] == "failed" else
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


def restore_inputs(cfg, run):
    """Rebuild a validated run's weather rasters into <run folder>/inputs/ with the same
    code that built them, and keep them only if they match the checksums recorded in
    its weather.json. Returns present, restored, not_validated or mismatch: <files>."""
    d = final_dir(cfg, run)
    if not (d / "DONE").exists():
        return "not_validated"
    want = read_json(d / "weather.json")["raster_sha256"]
    inputs = d / "inputs"
    if all((inputs / n).exists() and sha256(inputs / n) == h for n, h in want.items()):
        return "present"
    tmp = d / ".restore"
    shutil.rmtree(tmp, ignore_errors=True)
    try:
        build_case(cfg, run, tmp, tmp / "scratch")
        bad = [n for n, h in want.items() if sha256(tmp / "inputs" / n) != h]
        if bad:
            return f"mismatch: {bad}"
        inputs.mkdir(exist_ok=True)
        for n in want:
            os.replace(tmp / "inputs" / n, inputs / n)
        return "restored"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def restore_all(cfg, ids, workers):
    matrix = {r.run_id: r for r in load_matrix(cfg)}
    runs = [resolve(cfg, matrix, i) for i in ids]
    with ThreadPoolExecutor(workers) as pool:
        results = list(pool.map(lambda r: restore_inputs(cfg, r), runs))
    counts = Counter(s if not s.startswith("mismatch") else "mismatch" for s in results)
    print(f"{len(runs)} runs: {dict(counts)}", flush=True)
    for r, s in zip(runs, results):
        if s.startswith("mismatch"):
            print(f"{r.run_id}: {s}", flush=True)
    return 1 if counts["mismatch"] else 0


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
    # ELMFIRE marks the ignition cell itself with -1 (seen in the baseline and the
    # pilot); its neighbours carry arrival times. So -1 here confirms the location.
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


def run_batch(cfg, batch, index, workers, force=False):
    size = cfg["slurm"]["runs_per_task"]
    ids = chunk(read_batch(batch), size, index)
    if not ids:
        raise ConfigError(f"array task {index} has no runs in {batch} (chunk size {size})")
    for needed in (cfg.shared_dir / "manifest.json", cfg.path(cfg["paths"]["prepared"]) / "provenance.json"):
        if not needed.exists():
            raise ConfigError(f"{needed} not found; run `python -m sobolfire prepare` first")
    matrix = {r.run_id: r for r in load_matrix(cfg)}
    runs = [resolve(cfg, matrix, i) for i in ids]     # every ID resolves before anything runs
    print(f"array task {index}: lines {index * size + 1}-{index * size + len(ids)} of {batch}, "
          f"{len(runs)} runs, {workers} at a time", flush=True)
    with ThreadPoolExecutor(workers) as pool:
        statuses = list(pool.map(lambda r: run_one(cfg, r, force), runs))
    counts = Counter(statuses)
    print(f"array task {index} finished: {dict(counts)}", flush=True)
    return 1 if counts["failed"] else 0


def selection(cfg, which, of=None):
    """Run IDs for make-batch: pilot, all, controls, spotting (the spotting_study),
    a station prefix, or retry (failed or unfinished runs in the batch file `of`)."""
    matrix = load_matrix(cfg)
    stations = list(cfg.stations)
    if which == "all":
        return [r.run_id for r in matrix]
    if which in stations:
        return [r.run_id for r in matrix if r.prefix == which]
    if which == "controls":
        return [f"{p}_BASELINE" for p in stations]
    if which == "spotting":      # nominal controls plus one spotting setting changed at a time
        study = cfg.get("spotting_study") or {}
        return [f"{p}_BASELINE" for p in stations] + [
            f"{p}_BASELINE__{key}-{level_text(v)}" for p in stations
            for key, levels in study.items() for v in levels]
    if which == "pilot":
        p = cfg["pilot"]
        return ([f"{s}_BASELINE" for s in stations] + [f"{r}{REP}2" for r in p["repeats"]]
                + [r.run_id for r in matrix if r.base_row in p["base_rows"]])
    if which == "retry":
        if not of:
            raise ConfigError("make-batch retry needs --of <batch file>")
        by_id = {r.run_id: r for r in matrix}
        return [i for i in read_batch(of)
                if state(final_dir(cfg, resolve(cfg, by_id, i))) in ("failed", "running", "not_started")]
    raise ConfigError(f"unknown selection {which!r}")


def make_batch(cfg, which, of=None, name=None):
    ids = selection(cfg, which, of)
    if not ids:
        print("nothing to run")
        return None
    name = name or (f"retry_{utc_now().replace(':', '')}" if which == "retry" else which)
    path = cfg.dir("batches") / f"{name}.txt"
    sl = cfg["slurm"]
    size = sl["runs_per_task"]
    ntasks = math.ceil(len(ids) / size)
    path.write_text(f"# {len(ids)} runs, selection {which!r}, created {utc_now()}\n" + "\n".join(ids) + "\n")
    print(f"wrote {path} ({len(ids)} runs)")
    print(f"{size} runs per array task -> tasks 0-{ntasks - 1}; task i runs lines {size}*i+1 .. {size}*(i+1)")
    for i in sorted({0, ntasks - 1}):
        part = chunk(ids, size, i)
        print(f"  task {i}: {part[0]} .. {part[-1]} ({len(part)} runs)")
    print("Submit from the experiment folder (CPUs, time and memory are set in the Slurm file):")
    print(f"  cd {cfg.root}")
    print(f"  sbatch --array=0-{ntasks - 1}%{sl['max_concurrent_tasks']} slurm/sobol_batch.slurm "
          f"{path.relative_to(cfg.root)}")
    return path


def status(cfg, ids):
    matrix = {r.run_id: r for r in load_matrix(cfg)}
    counts, stalled, unfinished = Counter(), [], []
    for i in ids:
        run = resolve(cfg, matrix, i)
        d = final_dir(cfg, run)
        st = state(d)
        counts[(run.prefix, st)] += 1
        if st in ("failed", "running"):
            unfinished.append(f"{i} ({st})")
        if st == "validated" and any("STALLED" in m for m in read_json(d / "DONE")["stop_messages"]):
            stalled.append(i)
    states = ["validated", "failed", "running", "not_started"]
    print(f"{'station':8s}" + "".join(f"{s:>13s}" for s in states))
    for p in cfg.stations:
        print(f"{p:8s}" + "".join(f"{counts[(p, s)]:13d}" for s in states))
    print("'running' includes runs whose job ended before they finished.")
    if stalled:
        print(f"{len(stalled)} validated runs stopped early (fire front stalled), e.g. {stalled[:5]}")
    if unfinished:
        print(f"failed or unfinished: {unfinished[:10]}{' ...' if len(unfinished) > 10 else ''}")
    return counts


def compare(cfg, id_a, id_b):
    """Arrival-time differences between two validated runs (seed-policy check)."""
    matrix = {r.run_id: r for r in load_matrix(cfg)}
    arrays = []
    for i in (id_a, id_b):
        d = final_dir(cfg, resolve(cfg, matrix, i))
        if not (d / "DONE").exists():
            raise ConfigError(f"{i} is not validated")
        with rasterio.open(d / read_json(d / "DONE")["toa_file"]) as f:
            arrays.append(f.read(1))
    a, b = arrays
    both = (a > 0) & (b > 0)
    print(f"identical: {np.array_equal(a, b)}; burned in only one: {int(((a > 0) != (b > 0)).sum())} cells; "
          f"max |difference| where both burned: {float(np.abs(a - b)[both].max(initial=0)):.3f} s")
