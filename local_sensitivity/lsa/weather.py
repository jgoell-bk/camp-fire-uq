"""Spatially uniform weather rasters on a coarse grid in the case CRS.

ELMFIRE only requires the weather rasters to match each other; it maps each fuel
cell to a weather cell by coordinates (MAP_FINE_TO_COARSE in elmfire_subs.f90) and
does not check the CRS. A small grid that covers the fuel grid with a margin is
enough for uniform weather. Same code as ../experiment/sobolfire/weather.py.
"""
import math
import os
from dataclasses import asdict, dataclass

import numpy as np
import rasterio
from rasterio.transform import from_origin

from .common import StudyError, sha256

VARIABLES = ["ws", "wd", "m1", "m10", "m100"]   # file stems, as in the namelist template
NODATA = -9999.0


@dataclass(frozen=True)
class Grid:
    crs: str
    x0: float          # left edge
    top: float         # top edge
    cellsize: float
    ncols: int
    nrows: int

    @property
    def transform(self):
        return from_origin(self.x0, self.top, self.cellsize, self.cellsize)

    @property
    def bottom(self):
        return self.top - self.nrows * self.cellsize

    def as_dict(self):
        return asdict(self)


def grid_of(path):
    with rasterio.open(path) as f:
        t = f.transform
        if t.b or t.d or not math.isclose(t.a, -t.e):
            raise StudyError(f"{path}: not a north-up grid with square cells")
        return Grid(f.crs.to_string(), t.c, t.f, t.a, f.width, f.height)


def weather_grid(fuel, cellsize, margin):
    """Coarse grid covering the fuel grid, with `margin` extra cells on every side."""
    ncols = math.ceil(fuel.ncols * fuel.cellsize / cellsize) + 2 * margin
    nrows = math.ceil(fuel.nrows * fuel.cellsize / cellsize) + 2 * margin
    return Grid(fuel.crs, fuel.x0 - margin * cellsize, fuel.top + margin * cellsize,
                float(cellsize), ncols, nrows)


def _nint(x):
    """Fortran NINT: round half away from zero."""
    return np.where(x >= 0, np.floor(x + 0.5), np.ceil(x - 0.5)).astype(int)


def check_coverage(fuel, wx):
    """Repeat ELMFIRE's fuel-to-weather index arithmetic for every fuel cell centre
    (MAP_FINE_TO_COARSE and GET_BILINEAR_INTERPOLATE_COEFFS, elmfire_subs.f90).
    Each centre must lie inside the weather grid, and the bilinear stencil
    I1 = 1 + NINT(d), I2 = I1 + 1 must not reach past the last weather cell: where
    I1 is clamped in both directions ELMFIRE divides 0 by 0."""
    problems = []
    for axis, n_fuel, start_fuel, n_wx, start_wx in [
            ("column", fuel.ncols, fuel.x0, wx.ncols, wx.x0),
            ("row", fuel.nrows, fuel.bottom, wx.nrows, wx.bottom)]:   # ELMFIRE rows count up from the bottom
        centres = start_fuel + (np.arange(1, n_fuel + 1) - 0.5) * fuel.cellsize
        d = (centres - start_wx) / wx.cellsize
        containing = np.ceil(d).astype(int)
        if containing.min() < 1 or containing.max() > n_wx:
            problems.append(f"fuel {axis}s extend outside the weather grid")
        if (1 + _nint(d)).max() >= n_wx:
            problems.append(f"bilinear stencil reaches the last weather {axis}")
    if fuel.crs != wx.crs:
        problems.append(f"CRS {wx.crs} differs from the fuel CRS {fuel.crs}")
    if problems:
        raise StudyError("weather grid: " + "; ".join(problems) + "; increase weather_grid.margin_cells")


def write_weather(directory, grid, series):
    """One Float32 multiband GeoTIFF per variable; band k holds series[var][k-1] in
    every cell. Each file is read back and checked under a temporary name, then
    renamed. Returns {file name: sha256}."""
    hashes = {}
    for var in VARIABLES:
        values = np.asarray(series[var], dtype=np.float32)
        cube = np.broadcast_to(values[:, None, None], (len(values), grid.nrows, grid.ncols))
        final = directory / f"{var}.tif"
        tmp = directory / f".{var}.tif.partial"
        with rasterio.open(tmp, "w", driver="GTiff", dtype="float32", count=len(values),
                           width=grid.ncols, height=grid.nrows, crs=grid.crs,
                           transform=grid.transform, nodata=NODATA) as dst:
            dst.write(cube)
        with rasterio.open(tmp) as f:
            ok = (f.count == len(values) and (f.height, f.width) == (grid.nrows, grid.ncols)
                  and f.crs.to_string() == grid.crs and f.dtypes[0] == "float32"
                  and np.allclose(tuple(f.transform)[:6], tuple(grid.transform)[:6])
                  and np.array_equal(f.read(), cube))
        if not ok:
            raise RuntimeError(f"{tmp}: read-back check failed")
        os.replace(tmp, final)
        hashes[final.name] = sha256(final)
    return hashes
