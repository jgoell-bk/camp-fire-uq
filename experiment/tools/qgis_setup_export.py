"""Write qgis_setup/: layers for a slide showing the simulation setup.

  cd /oak/stanford/groups/gorle/jgoell/Wildfire/experiment
  source /oak/stanford/groups/gorle/jgoell/venvs/elmfire-vnv/bin/activate
  python -I tools/qgis_setup_export.py external/tl_2018_06_place/tl_2018_06_place.shp

The argument is the Census TIGER/Line 2018 California places shapefile
(https://www2.census.gov/geo/tiger/TIGER2018/PLACE/tl_2018_06_place.zip).
See qgis_setup/README.txt for the layers and their sources.
"""
import shutil
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer
from rasterio.features import rasterize, shapes
from shapely.geometry import Point, box, shape

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sobolfire.case import nonburnable       # noqa: E402
from sobolfire.config import Config          # noqa: E402
from sobolfire.util import read_json         # noqa: E402

CRS = "EPSG:26910"
# NIFC "Public View - Interagency RAWS" (retrieved 2026-10-20); IDs match the forcing data.
STATIONS = [  # prefix, name, WIMS id, lat, lon, elevation ft, role in the experiment
    ("JG", "Jarbo Gap", "041214", 39.73591, -121.48898, 2535, "wind scenario; fuel moisture for all runs"),
    ("OP", "Openshaw", "041215", 39.58983, -121.63516, 268, "wind scenario"),
    ("CO", "Cohasset", "041211", 39.87184, -121.76898, 1733, "wind scenario"),
]
PARADISE_CENTRE = (-121.6219, 39.7596)   # stand-in point used in the pilot summaries
BORDER_M = 1000                          # width of the wildland border band outside the town limits


def hillshade(dem, cellsize, azimuth=315.0, altitude=45.0):
    z = np.where(dem > -9000, dem, np.nan)
    gy, gx = np.gradient(z, cellsize)            # gy: dz per metre southward; gx: eastward
    slope = np.arctan(np.hypot(gx, gy))
    aspect = np.arctan2(gy, -gx)                 # ESRI hillshade convention
    az, alt = np.radians(360 - azimuth + 90), np.radians(altitude)
    hs = np.sin(alt) * np.cos(slope) + np.cos(alt) * np.sin(slope) * np.cos(az - aspect)
    return np.where(np.isfinite(hs), np.clip(hs, 0, 1) * 255, 0).astype("uint8")


def main(places_shp):
    cfg = Config()
    out = cfg.root / "qgis_setup"
    out.mkdir(exist_ok=True)
    gpkg = out / "setup.gpkg"
    gpkg.unlink(missing_ok=True)
    shared = read_json(cfg.shared_dir / "manifest.json")
    fuel, ign = shared["fuel_grid"], shared["ignition"]
    with rasterio.open(cfg.shared_dir / "fuels" / "fbfm40b.tif") as f:
        fb, profile = f.read(1), f.profile
    T = profile["transform"]
    to_utm = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)

    def layer(name, records, geoms):
        gpd.GeoDataFrame(records, geometry=geoms, crs=CRS).to_file(gpkg, layer=name)

    domain = box(fuel["x0"], fuel["top"] - fuel["nrows"] * fuel["cellsize"],
                 fuel["x0"] + fuel["ncols"] * fuel["cellsize"], fuel["top"])
    layer("domain", {"name": ["ELMFIRE domain, 30 m cells"]}, [domain])
    layer("ignition", {"name": ["Ignition, 06:30 PST 8 Nov 2018"], "lat": [ign["lat"]], "lon": [ign["lon"]]},
          [Point(ign["cell_x"], ign["cell_y"])])

    # Stations, with the observed wind over the 24 h window (vector mean, blowing-from)
    recs, pts = [], []
    for prefix, name, wims, lat, lon, elev, role in STATIONS:
        w = pd.read_csv(cfg.path(cfg["paths"]["prepared"]) / f"forcing_{prefix}.csv")
        u = (w.ws_mph * np.sin(np.radians(w.wd_deg))).mean()
        v = (w.ws_mph * np.cos(np.radians(w.wd_deg))).mean()
        recs.append({"name": name, "label": f"{name} ({prefix})", "wims_id": wims, "elev_ft": elev, "role": role,
                     "lat": lat, "lon": lon, "mean_ws_mph": round(w.ws_mph.mean(), 1),
                     "mean_wind_from_deg": round(float(np.degrees(np.arctan2(u, v)) % 360), 1),
                     "inside_domain": domain.contains(Point(to_utm.transform(lon, lat)))})
        pts.append(Point(to_utm.transform(lon, lat)))
    layer("stations", pd.DataFrame(recs).to_dict("list"), pts)

    places = gpd.read_file(places_shp).to_crs(CRS)
    town = places[places.GEOID == "0655520"]            # Paradise town, Butte County
    magalia = places[places.GEOID == "0645120"]         # Magalia CDP, just north
    cols = ["NAMELSAD", "GEOID", "ALAND", "geometry"]
    town[cols].to_file(gpkg, layer="paradise_town_2018")
    magalia[cols].to_file(gpkg, layer="magalia_cdp_2018")
    px, py = to_utm.transform(*PARADISE_CENTRE)
    layer("paradise_2km_standin", {"name": ["2 km zone used in the pilot summaries"]}, [Point(px, py).buffer(2000)])

    # Model cells whose centres fall inside the town limits: a mask for `assemble --mask`
    mask = rasterize([(g, 1) for g in town.geometry], out_shape=fb.shape, transform=T, dtype="uint8")
    p8 = {**profile, "dtype": "uint8", "nodata": 0, "compress": "deflate"}
    with rasterio.open(out / "paradise_town_mask.tif", "w", **p8) as f:
        f.write(mask, 1)
    cells = [shape(g) for g, val in shapes(mask, mask=mask.astype(bool), transform=T) if val == 1]
    layer("paradise_town_cells", {"name": ["model cells inside Paradise town limits"] * len(cells)}, cells)

    # Border band: wildland cells (burnable, not urban) within BORDER_M outside the limits
    limits = town.geometry.iloc[0]
    ring = limits.buffer(BORDER_M).difference(limits)
    border = rasterize([(ring, 1)], out_shape=fb.shape, transform=T, dtype="uint8") & (~nonburnable(fb)).astype("uint8")
    with rasterio.open(out / "paradise_border_mask.tif", "w", **p8) as f:
        f.write(border, 1)
    cells = [shape(g) for g, val in shapes(border, mask=border.astype(bool), transform=T) if val == 1]
    layer("paradise_border_cells", {"name": [f"wildland cells within {BORDER_M} m outside Paradise limits"] * len(cells)},
          cells)

    urban = (fb == 91).astype("uint8")                   # cells with buildings (urban fuel)
    with rasterio.open(out / "building_cells.tif", "w", **p8) as f:
        f.write(urban, 1)
    with rasterio.open(cfg.shared_dir / "fuels" / "dem.tif") as f:
        dem, dprof = f.read(1), f.profile
    with rasterio.open(out / "hillshade.tif", "w", **{**dprof, "dtype": "uint8", "nodata": None, "compress": "deflate"}) as f:
        f.write(hillshade(dem, dprof["transform"].a), 1)
    shutil.copy(cfg.shared_dir / "fuels" / "dem.tif", out / "dem.tif")

    print(f"wrote {out}")
    print(f"  Paradise town: {int(mask.sum())} model cells inside the limits, "
          f"{int((mask.astype(bool) & (fb == 91)).sum())} of them with buildings")
    print(f"  border band: {int(border.sum())} wildland cells within {BORDER_M} m outside the limits, "
          f"{len(cells)} polygons")
    for r in recs:
        print(f"  {r['label']:16s} inside domain: {r['inside_domain']}, 24 h mean wind {r['mean_ws_mph']} mph "
              f"from {r['mean_wind_from_deg']} deg")


if __name__ == "__main__":
    main(sys.argv[1])
