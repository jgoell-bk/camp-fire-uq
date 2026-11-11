"""Over-spread (v3, Jacob 2026-11-10): the area a run burned outside the final fire perimeter.

The deadline term only penalizes a fire that is late, so on its own a faster fire
always scores better. The final perimeter bounds the fire at every time: a cell outside
it never burned. So

    A_out = (cell area) x #{cells outside the final perimeter with 0 < T <= at_s, or T = -1}

(T = -1 marks the ignition cell, as in eki.pilot.distance_to_burned), and the objective
gains 1/2 (A_out / reference_acres)^2 (eki.objective).

The perimeter is the CAL FIRE final perimeter of the Camp Fire (153,336 ac, contained
26 Nov 2018), rasterized once on the fuel grid by `python -m eki perimeter`: a cell is
inside when its centre is inside the polygon (gdal_rasterize's default rule). It is a
lower bound on the over-spread at 24 h, since by then the real fire filled only part of
the final perimeter. It never uses model output.
"""
import shutil
import subprocess
import tempfile

import numpy as np
import rasterio

from .common import StudyError, read_json, run_dir, sha256, utc_now, write_json
from .weather import Grid, grid_of

ACRE_M2 = 4046.8564224


def settings(cfg):
    o = cfg["fit"].get("overspread")
    if not o:
        raise StudyError("fit.overspread is not configured")
    return o


def raster_path(cfg):
    return cfg.path(cfg["paths"]["prepared"]) / "final_perimeter.tif"


def fuel_grid(cfg):
    return Grid(**read_json(cfg.shared_dir / "manifest.json")["fuel_grid"])


def prepare_perimeter(cfg):
    """Rasterize sources.final_perimeter on the fuel grid -> prepared/final_perimeter.tif
    (1 inside, 0 outside) and its .json. Needs the GDAL command-line tools
    (module load physics gdal/3.10.2)."""
    for tool in ("ogr2ogr", "gdal_rasterize"):
        if shutil.which(tool) is None:
            raise StudyError(f"{tool} not found; module load physics gdal/3.10.2 first")
    src = cfg.path(cfg["sources"]["final_perimeter"])
    if not src.exists():
        raise StudyError(f"{src} not found")
    g = fuel_grid(cfg)
    out = raster_path(cfg)
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        vec, ras = f"{tmp}/perimeter.gpkg", f"{tmp}/perimeter.tif"
        cmds = [["ogr2ogr", "-t_srs", g.crs, "-nlt", "MULTIPOLYGON", vec, str(src)],
                ["gdal_rasterize", "-q", "-burn", "1", "-init", "0", "-ot", "Byte", "-a_srs", g.crs,
                 "-te", repr(g.x0), repr(g.bottom), repr(g.x0 + g.ncols * g.cellsize), repr(g.top),
                 "-ts", str(g.ncols), str(g.nrows), "-co", "COMPRESS=DEFLATE", vec, ras]]
        for c in cmds:
            r = subprocess.run(c, capture_output=True, text=True)
            if r.returncode:
                raise StudyError(f"{c[0]} failed: {r.stderr.strip()}")
        if grid_of(ras) != g:
            raise StudyError(f"rasterized perimeter is not on the fuel grid: {grid_of(ras)} != {g}")
        shutil.copyfile(ras, out)
    with rasterio.open(out) as f:
        inside = f.read(1).astype(bool)
    info = {"created_utc": utc_now(), "file": out.name, "sha256": sha256(out), "source": str(src),
            "source_sha256": sha256(src), "commands": cmds, "fuel_grid": g.as_dict(),
            "inside_cells": int(inside.sum()), "inside_acres": float(inside.sum() * g.cellsize ** 2 / ACRE_M2),
            "rule": "a cell is inside when its centre is inside the polygon (gdal_rasterize default)"}
    write_json(out.with_suffix(".json"), info)
    return info


def perimeter_mask(cfg):
    """(inside mask, cell area in acres) of prepared/final_perimeter.tif, checked to be
    on the fuel grid."""
    path = raster_path(cfg)
    if not path.exists():
        raise StudyError(f"{path} not found; run `python -m eki --config ... perimeter` first")
    g = fuel_grid(cfg)
    if grid_of(path) != g:
        raise StudyError(f"{path} is not on the fuel grid")
    with rasterio.open(path) as f:
        return f.read(1).astype(bool), g.cellsize ** 2 / ACRE_M2


def burned_outside_acres(toa, inside, cell_acres, at_s):
    """A_out: area burned by at_s outside the perimeter (T = -1 is the ignition cell)."""
    if toa.shape != inside.shape:
        raise StudyError(f"arrival raster {toa.shape} and perimeter {inside.shape} differ in shape")
    burned = ((toa > 0) & (toa <= at_s)) | (toa == -1.0)
    return float((burned & ~inside).sum() * cell_acres)


def run_acres(cfg, rid, mask=None):
    """A_out of a validated run."""
    from .run import state
    inside, cell_acres = mask if mask is not None else perimeter_mask(cfg)
    d = run_dir(cfg, rid)
    if state(d) != "validated":
        raise StudyError(f"{rid} is {state(d)}, not validated")
    with rasterio.open(d / read_json(d / "DONE")["toa_file"]) as f:
        toa = f.read(1)
    return burned_outside_acres(toa, inside, cell_acres, float(settings(cfg)["at_s"]))
