"""Local sensitivity of fire arrival time to wind direction (Camp Fire, ELMFIRE).

Two run sets (RUN_SET below; ../README.md):
- hourly: each run <S>_WD<theta> turned the station's recorded hourly wind
  directions by aim + theta degrees, so that theta = 0 points their mean from the
  ignition straight at Paradise;
- fixed: every hour of run <S>_WD<theta> blows from one direction, the direction
  that points straight at Paradise, plus theta.
In both, positive theta turns the wind clockwise, and wind speed and moisture are
the recorded hourly values.

For every 30 m cell, this script estimates the derivative of the arrival time T
with respect to theta at one point theta0, by a central finite difference,

    dT/dtheta (theta0)  ~=  [ T(theta0 + h) - T(theta0 - h) ] / (2 h),

It draws the result as a filled contour map, one panel per station, with the
arrival-time isochrones of the theta0 run on top (figures/). It also writes it as
GeoTIFFs for QGIS (qgis/, section 5).

Before the first use, collect each station's runs into a table (seconds each):
    python -m lsa collect JG
    python -m lsa collect OP
    python -m lsa collect CO
and, for the fixed-direction runs,
    python -m lsa --config fixed_direction/config.yaml collect JG     (and OP, CO)
Then, from the local_sensitivity folder:
    python analysis/lsa_gradient.py
"""
import json
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from matplotlib.colors import LinearSegmentedColormap

# %% 0. Settings --------------------------------------------------------------
RUN_SET = "fixed"  # "hourly" (recorded hourly directions, turned) or "fixed" (one direction)
THETA0 = 0        # degrees: where the derivative is taken (0 = mean wind aimed at Paradise)
H = 1             # degrees: finite-difference step. The runs are 1 degree apart, so any
                  # whole number works. H = 1 is the most local; a larger H averages the
                  # slope over theta0 - H .. theta0 + H and is less noisy (see section 2).
STATIONS = {"JG": "Jarbo Gap", "OP": "Openshaw", "CO": "Cohasset"}
ISOCHRONES_H = np.arange(3, 25, 3)   # arrival-time contour lines of the theta0 run, hours
CLIP_PERCENTILE = 95  # colour scale spans +/- this percentile of |dT/dtheta| (all stations)
ZOOM_KM = None        # None = whole domain, or (xmin, xmax, ymin, ymax) in km from the ignition

STUDY = Path(__file__).resolve().parents[1]            # the local_sensitivity folder
WILDFIRE = STUDY.parent
DATA = {"hourly": STUDY, "fixed": STUDY / "fixed_direction"}[RUN_SET]   # the run set's folder
RESULTS = DATA / "results"                             # <S>_runs.csv from `python -m lsa collect`
FIGURES = DATA / "figures"
QGIS_DIR = DATA / "qgis"                               # GeoTIFFs and styles, to download for QGIS
MASK_TIF = WILDFIRE / "experiment/qgis_setup/paradise_border_mask.tif"   # the QoI cells
TOWN_SHP = WILDFIRE / "experiment/external/tl_2018_06_place/tl_2018_06_place.shp"

# Colours: a diverging scale (two hues, neutral grey at zero) for the signed
# derivative, and neutral inks for everything else. The QGIS style in section 5
# uses the same five colours.
DIVERGING_COLOURS = ["#104281", "#2a78d6", "#f0efec", "#e34948", "#8c1d1d"]
DIVERGING = LinearSegmentedColormap.from_list("blue_grey_red", DIVERGING_COLOURS)
INK, INK_2 = "#0b0b0b", "#52514e"
plt.rcParams.update({"axes.edgecolor": INK_2, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": INK_2, "ytick.color": INK_2, "font.size": 9})


# %% 1. Read the arrival-time rasters ------------------------------------------
def read_toa_hours(path):
    """Arrival time in hours after t0 (06:00 PST 8 Nov 2018). ELMFIRE writes -9999
    where a cell did not burn within 24 h and -1 at the ignition cell itself; both
    become NaN here, because neither is an arrival time."""
    with rasterio.open(path) as f:
        t = f.read(1).astype(np.float64)
    t[t <= 0] = np.nan
    return t / 3600.0


aim = json.loads((DATA / "prepared/aim.json").read_text())
ix, iy = aim["ignition_cell_x_y"]                      # ignition cell centre, EPSG:26910 metres

toa0, t_minus, t_plus = {}, {}, {}
for s in STATIONS:
    runs = pd.read_csv(RESULTS / f"{s}_runs.csv").set_index("offset_deg")   # one row per theta
    for theta in (THETA0 - H, THETA0, THETA0 + H):
        if theta not in runs.index:
            raise SystemExit(f"{s}: no run at theta = {theta}; the runs span "
                             f"{runs.index.min()} .. {runs.index.max()} degrees")
    t_minus[s] = read_toa_hours(runs.loc[THETA0 - H, "toa_path"])
    toa0[s] = read_toa_hours(runs.loc[THETA0, "toa_path"])
    t_plus[s] = read_toa_hours(runs.loc[THETA0 + H, "toa_path"])

# Cell-centre coordinates, from the raster's geotransform, in km from the ignition.
# Row 0 is the northern edge (the transform's e is negative), so y decreases down the rows.
with rasterio.open(runs.loc[THETA0, "toa_path"]) as f:
    tr, nrows, ncols = f.transform, f.height, f.width
x_km = (tr.c + (np.arange(ncols) + 0.5) * tr.a - ix) / 1000.0
y_km = (tr.f + (np.arange(nrows) + 0.5) * tr.e - iy) / 1000.0


# %% 2. Central finite difference ----------------------------------------------
# The derivative only exists where the cell burned in both neighbouring runs:
# NaN - number = NaN, so cells unburned in either run drop out by themselves.
# Units: hours per degree x 60 = minutes of arrival time per degree of turn.
#
# Two things to keep in mind when reading the result:
# - Spot fires start in jumps. Whether embers reach a cell can change between
#   two runs 1 degree apart, and that changes the arrival times downwind of it
#   abruptly, not smoothly. With H = 1, such jumps show up as speckle or sharp
#   patches. Compare H = 1 with H = 3 or 5 to see which features are robust.
# - Near the edge of the burned area the derivative is large and then undefined,
#   because a small turn decides whether the cell burns at all within 24 h.
grad = {s: (t_plus[s] - t_minus[s]) / (2 * H) * 60.0 for s in STATIONS}


# %% 3. Numbers for the Paradise border band (the QoI cells) ------------------
with rasterio.open(MASK_TIF) as f:
    band = f.read(1) == 1
    assert band.shape == (nrows, ncols), "mask is not on the arrival-time grid"

print(f"dT/dtheta at theta0 = {THETA0:+d} deg, h = {H} deg, over the {band.sum()} border-band cells")
for s, name in STATIONS.items():
    g = grad[s][band]
    ok = np.isfinite(g)
    only_one = np.isfinite(t_minus[s][band]) != np.isfinite(t_plus[s][band])
    print(f"  {s} ({name}): burned at theta0 {np.isfinite(toa0[s][band]).mean():5.1%}, "
          f"derivative defined {ok.mean():5.1%}, burned in only one of the two runs {only_one.mean():5.1%}")
    if ok.any():
        p5, p50, p95 = np.percentile(g[ok], [5, 50, 95])
        print(f"      min/deg: median {p50:+.2f}, mean {g[ok].mean():+.2f}, 5-95% {p5:+.2f} .. {p95:+.2f}")


# %% 4. Contour maps ------------------------------------------------------------
# One colour scale for all three stations, so the panels can be compared; it is
# symmetric about zero and clipped at a percentile so a few extreme cells don't
# wash out the rest (they get the end colours: extend="both").
everything = np.concatenate([np.abs(g[np.isfinite(g)]) for g in grad.values()])
lim = np.percentile(everything, CLIP_PERCENTILE)
levels = np.linspace(-lim, lim, 18)     # 17 bins: 8 per sign and a grey one around zero

town = gpd.read_file(TOWN_SHP, where="GEOID='0655520'").to_crs("EPSG:26910")
town_km = town.geometry.translate(-ix, -iy).scale(1e-3, 1e-3, origin=(0, 0))

fig, axes = plt.subplots(1, len(STATIONS), figsize=(5 * len(STATIONS), 5.6),
                         sharex=True, sharey=True, layout="constrained")
for ax, (s, name) in zip(np.atleast_1d(axes), STATIONS.items()):
    filled = ax.contourf(x_km, y_km, grad[s], levels=levels, cmap=DIVERGING, extend="both")
    iso = ax.contour(x_km, y_km, toa0[s], levels=ISOCHRONES_H, colors=INK_2, linewidths=0.6)
    ax.clabel(iso, fmt="%d h", fontsize=7)
    ax.contour(x_km, y_km, band.astype(float), levels=[0.5], colors=INK, linewidths=0.8,
               linestyles="dashed")                                     # the QoI band
    town_km.boundary.plot(ax=ax, color=INK, linewidth=1.0)              # Paradise town limits
    ax.plot(0, 0, marker="*", markersize=11, color=INK, linestyle="none")   # ignition
    # Direction the mean wind blows toward in the theta0 run, as a 5 km arrow.
    toward = np.radians(aim["target_toward_deg"] + THETA0)
    ax.annotate("", xy=(5 * np.sin(toward), 5 * np.cos(toward)), xytext=(0, 0),
                arrowprops=dict(arrowstyle="->", color=INK, linewidth=1.2))
    ax.set_title(f"{name} ({s})")
    ax.set_aspect("equal")
    ax.set_xlabel("km east of the ignition")
    if ZOOM_KM:
        ax.set_xlim(ZOOM_KM[:2])
        ax.set_ylim(ZOOM_KM[2:])
np.atleast_1d(axes)[0].set_ylabel("km north of the ignition")

cbar = fig.colorbar(filled, ax=axes, orientation="horizontal", shrink=0.6, aspect=40)
cbar.set_label("dT/dθ, minutes per degree   (> 0: the cell burns later as the wind turns clockwise)")
fig.suptitle(f"Sensitivity of arrival time to wind direction ({RUN_SET} direction runs) at θ0 = {THETA0:+d}° "
             f"(central difference, h = {H}°)\n"
             "Lines: arrival time of the θ0 run. Dashed: Paradise border band; solid: town limits. "
             "Star: ignition; arrow: mean wind. Blank: unburned in a neighbouring run.", fontsize=9)

FIGURES.mkdir(exist_ok=True)
out = FIGURES / f"lsa_gradient_theta{THETA0:+03d}_h{H}.png"
fig.savefig(out, dpi=200)
print(f"wrote {out}")


# %% 5. GeoTIFFs for QGIS -------------------------------------------------------
# Per station, two single-band GeoTIFFs in EPSG:26910:
#   <S>_dTdtheta_theta+00_h1.tif   dT/dtheta, minutes per degree (the map above)
#   <S>_toa_theta+00_hours.tif     arrival time of the theta0 run, hours after 06:00 PST
# Cells without a value (unburned, or the ignition cell) hold -9999, declared as
# NoData, so QGIS leaves them empty.
#
# The grid (corner, cell size, CRS) is taken from the border-band mask, which sits
# exactly on the fuel grid like every layer in ../experiment/qgis_setup/. ELMFIRE's
# own rasters store the corner in single precision, about 0.2 m off; the values
# are the same cells either way.
with rasterio.open(MASK_TIF) as f:
    grid = {"crs": f.crs, "transform": f.transform, "width": f.width, "height": f.height}


def write_geotiff(path, values, description):
    """values: a 2-D array on the fuel grid, NaN where there is no value."""
    data = np.where(np.isfinite(values), values, -9999.0).astype(np.float32)
    with rasterio.open(path, "w", driver="GTiff", count=1, dtype="float32", nodata=-9999.0,
                       compress="deflate", **grid) as dst:
        dst.write(data, 1)
        dst.set_band_description(1, description)       # shown in QGIS layer properties
    print(f"wrote {path}")


# A QGIS style (.qml) for the gradient: the same diverging colours as the PNG, from
# -lim to +lim, with grey at zero. Values beyond +/- lim get the end colours. QGIS
# should apply a .qml with the same name as the .tif when the layer is added; if it
# doesn't, use Layer Properties > Symbology > Style > Load Style.
stops = np.linspace(-lim, lim, len(DIVERGING_COLOURS))
items = "\n".join(f'          <item alpha="255" value="{v:.6g}" label="{v:+.2f} min/deg" color="{c}"/>'
                  for v, c in zip(stops, DIVERGING_COLOURS))
qml = f"""<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.28" styleCategories="Symbology">
  <pipe>
    <rasterrenderer type="singlebandpseudocolor" band="1" opacity="1" alphaBand="-1"
                    classificationMin="{-lim:.6g}" classificationMax="{lim:.6g}" nodataColor="">
      <rastershader>
        <colorrampshader colorRampType="INTERPOLATED" classificationMode="1" clip="0"
                         minimumValue="{-lim:.6g}" maximumValue="{lim:.6g}">
{items}
        </colorrampshader>
      </rastershader>
    </rasterrenderer>
  </pipe>
</qgis>
"""

QGIS_DIR.mkdir(exist_ok=True)
tag = f"theta{THETA0:+03d}"
for s in STATIONS:
    grad_tif = QGIS_DIR / f"{s}_dTdtheta_{tag}_h{H}.tif"
    write_geotiff(grad_tif, grad[s], f"{s} dT/dtheta at theta0 = {THETA0:+d} deg (h = {H} deg), minutes per degree")
    grad_tif.with_suffix(".qml").write_text(qml)
    write_geotiff(QGIS_DIR / f"{s}_toa_{tag}_hours.tif", toa0[s],
                  f"{s} arrival time of the theta0 = {THETA0:+d} deg run, hours after 06:00 PST 8 Nov 2018")
print(f"colour range in the .qml files: +/- {lim:.2f} min/deg ({CLIP_PERCENTILE}th percentile of |dT/dtheta|)")
