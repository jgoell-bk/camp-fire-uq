"""Write QGIS-ready files to qgis/ (small; copy the folder to a laptop with QGIS).

  <run_id>_arrival_h.tif   arrival time in hours after t0 (06:00 PST 8 Nov 2018);
                           unburned cells and ELMFIRE's ignition marker are NoData
  arrival_hours.qml        0-24 h colour style for those rasters
  context.gpkg             ignition, central Paradise, 2 km zone, VIIRS detections
                           in the first 24 h (hours after t0), model domain
  fuel_model.tif, dem.tif  copies of the case's fuel models and elevation

  cd /oak/stanford/groups/gorle/jgoell/Wildfire/experiment
  source /oak/stanford/groups/gorle/jgoell/venvs/elmfire-vnv/bin/activate
  python tools/qgis_export.py                      # default runs (below)
  python tools/qgis_export.py JG_0123 CO_0456      # any validated runs
"""
import shutil
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer
from shapely.geometry import Point, box

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sobolfire.config import Config                       # noqa: E402
from sobolfire.design import load_matrix, resolve         # noqa: E402
from sobolfire.runner import final_dir                    # noqa: E402
from sobolfire.util import read_json                      # noqa: E402

# A row of base row 1, and the same row with only q_spot taken from B (P_EPS 0.0125 vs 0.0037)
DEFAULT = ["JG_0001", "JG_0007", "OP_0001", "OP_0007", "CO_0001", "CO_0007"]
PARADISE = (-121.6219, 39.7596)
RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]


def arrival_hours(cfg, by_id, run_id, out):
    d = final_dir(cfg, resolve(cfg, by_id, run_id))
    with rasterio.open(d / read_json(d / "DONE")["toa_file"]) as f:
        a, profile = f.read(1), f.profile
    h = np.where(a > 0, a / 3600.0, -9999.0).astype("float32")
    profile.update(dtype="float32", nodata=-9999.0, compress="deflate")
    with rasterio.open(out / f"{run_id}_arrival_h.tif", "w", **profile) as f:
        f.write(h, 1)


def style(out):
    items = "\n".join(f'          <item alpha="255" value="{24 * i / (len(RAMP) - 1):g}" '
                      f'label="{24 * i / (len(RAMP) - 1):g} h" color="{c}"/>' for i, c in enumerate(RAMP))
    (out / "arrival_hours.qml").write_text(f"""<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.28" styleCategories="Symbology">
  <pipe>
    <rasterrenderer type="singlebandpseudocolor" band="1" opacity="1" alphaBand="-1"
                    classificationMin="0" classificationMax="24" nodataColor="">
      <rastershader>
        <colorrampshader colorRampType="INTERPOLATED" classificationMode="1" clip="0"
                         minimumValue="0" maximumValue="24">
{items}
        </colorrampshader>
      </rastershader>
    </rasterrenderer>
  </pipe>
</qgis>
""")


def context(cfg, out):
    crs = "EPSG:26910"
    shared = read_json(cfg.shared_dir / "manifest.json")
    ign, fuel = shared["ignition"], shared["fuel_grid"]
    px, py = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(*PARADISE)
    gpkg = out / "context.gpkg"
    gpkg.unlink(missing_ok=True)
    gpd.GeoDataFrame({"name": ["ignition 06:30 PST"]}, geometry=[Point(ign["cell_x"], ign["cell_y"])],
                     crs=crs).to_file(gpkg, layer="ignition")
    gpd.GeoDataFrame({"name": ["central Paradise"]}, geometry=[Point(px, py)], crs=crs).to_file(gpkg, layer="paradise_centre")
    gpd.GeoDataFrame({"name": ["2 km stand-in zone"]}, geometry=[Point(px, py).buffer(2000)],
                     crs=crs).to_file(gpkg, layer="paradise_2km")
    gpd.GeoDataFrame({"name": ["model domain"]}, geometry=[box(fuel["x0"], fuel["top"] - fuel["nrows"] * fuel["cellsize"],
                     fuel["x0"] + fuel["ncols"] * fuel["cellsize"], fuel["top"])], crs=crs).to_file(gpkg, layer="domain")
    case = cfg.path(cfg["case"]["dir"])
    v = gpd.read_file(case / "data/viirs_observation/viirs_camp.shp").to_crs(crs)
    t = pd.to_datetime(v.ACQ_DATE.astype(str) + v.ACQ_TIME.astype(str).str.zfill(4), format="%Y-%m-%d%H%M", utc=True)
    v["hours_after_t0"] = (t - pd.Timestamp(cfg["simulation"]["t0_utc"])).dt.total_seconds() / 3600
    v[v.hours_after_t0 <= 24].to_file(gpkg, layer="viirs_first_24h")
    for name, stem in (("fuel_model.tif", "fbfm40b"), ("dem.tif", "dem")):
        shutil.copy(cfg.shared_dir / "fuels" / f"{stem}.tif", out / name)


if __name__ == "__main__":
    cfg = Config()
    out = cfg.root / "qgis"
    out.mkdir(exist_ok=True)
    by_id = {r.run_id: r for r in load_matrix(cfg)}
    for run_id in sys.argv[1:] or DEFAULT:
        arrival_hours(cfg, by_id, run_id, out)
        print(f"wrote {out / (run_id + '_arrival_h.tif')}")
    style(out)
    context(cfg, out)
    print(f"wrote {out}/arrival_hours.qml, context.gpkg, fuel_model.tif, dem.tif")
