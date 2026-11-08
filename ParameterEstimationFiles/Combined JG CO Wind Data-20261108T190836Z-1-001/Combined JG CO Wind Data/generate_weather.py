"""
Generate Weather

Blend the hourly event wind of two RAWS stations (Cohasset, Jarbo Gap) into two new
wind series for the camp_sobol_d6 runs:
  IDW  : inverse-distance weighted by distance from the ignition point
  MEAN : equal-weight (50/50) mean
Wind is averaged as vectors: convert to u/v, average u and v separately, recombine.
Moisture is identical across stations (all runs use Jarbo Gap), so it is not handled here.

Input  : 'observed' (unperturbed) block of a reference run's weather.json per station,
         ignition point from the shared manifest.json.
Output : one CSV per blend with columns band, time_utc, time_pst, ws_mph, wd_deg
         (ws in mph at 20 ft; wd in degrees the wind blows FROM, true north).

To run from shell:

ml physics py-pyproj/3.7.2_py312
ml py-pandas/2.2.1_py312
cd /scratch/users/rwdavies/elmfire
python3 generate_weather.py

@author: ryandavies19, Claude (Opus 5.5)
"""

import json
import numpy as np
import pandas as pd
from pyproj import Transformer

#############
# Functions #
#############

def wind_to_uv(ws, wd):
    """Speed and meteorological direction (blowing from, deg) -> u (east), v (north) components."""
    wd_rad = np.radians(wd)
    return -ws * np.sin(wd_rad), -ws * np.cos(wd_rad)

def uv_to_wind(u, v):
    """u (east), v (north) components -> speed and meteorological direction (blowing from, deg)."""
    return np.hypot(u, v), np.degrees(np.arctan2(-u, -v)) % 360

def load_observed_wind(path):
    """Read the unperturbed hourly wind (band, time_utc, time_pst, ws, wd) from a run's weather.json."""
    with open(path) as f:
        bands = json.load(f)["bands"]
    return pd.DataFrame([{"band": b["band"], "time_utc": b["time_utc"], "time_pst": b["time_pst"],
                          "ws": b["observed"]["ws"], "wd": b["observed"]["wd"]} for b in bands])

# Executable
if __name__ == "__main__":

    ##################
    ##### Inputs #####
    ##################

    run_root = "/scratch/users/jgoell/elmfire-class/camp_sobol_d6"
    manifest_path = f"{run_root}/shared/manifest.json" # ignition x, y in EPSG:26910
    station_weather = { # station: [lat, lon] (forecast.weather.gov), weather.json of a reference run
        "Cohasset": {"latlon": [39.87184, -121.76898],
                     "weather_json": f"{run_root}/runs/CO/b001/CO_0001/weather.json"},
        "JarboGap": {"latlon": [39.73591, -121.48898],
                     "weather_json": f"{run_root}/runs/JG/b001/JG_0001/weather.json"},
    }
    output_paths = { # blend name: output csv path
        "IDW":  "/scratch/users/rwdavies/elmfire/outputs/weather_IDW.csv",
        "MEAN": "/scratch/users/rwdavies/elmfire/outputs/weather_MEAN.csv",
    }
    target_crs = "EPSG:26910" # NAD83 / UTM zone 10N [m], same CRS as manifest.json

    ################
    ##### Code #####
    ################

    # ignition point (same cell the runs ignite from)
    with open(manifest_path) as f:
        ignition = json.load(f)["ignition"]
    ignition_xy = np.array([ignition["x"], ignition["y"]])

    # station distances from ignition [m] and inverse-distance weights
    to_target = Transformer.from_crs("EPSG:4326", target_crs, always_xy=True)
    dist = {name: np.hypot(*(np.array(to_target.transform(lon, lat)) - ignition_xy))
            for name, (lat, lon) in ((n, s["latlon"]) for n, s in station_weather.items())}
    inv = {name: 1 / d for name, d in dist.items()}
    weight_sets = {
        "IDW":  {name: w / sum(inv.values()) for name, w in inv.items()},
        "MEAN": {name: 1 / len(station_weather) for name in station_weather},
    }

    # load station wind and check the hourly bands line up
    wind = {name: load_observed_wind(s["weather_json"]) for name, s in station_weather.items()}
    ref = next(iter(wind.values()))
    for name, df in wind.items():
        assert df["time_utc"].equals(ref["time_utc"]), f"{name} time bands do not match"

    # blend in u/v space and write
    for blend, weights in weight_sets.items():
        u = sum(weights[n] * wind_to_uv(wind[n]["ws"], wind[n]["wd"])[0] for n in wind)
        v = sum(weights[n] * wind_to_uv(wind[n]["ws"], wind[n]["wd"])[1] for n in wind)
        ws, wd = uv_to_wind(u, v)
        out = ref[["band", "time_utc", "time_pst"]].assign(ws_mph=ws.round(4), wd_deg=wd.round(4))
        out.to_csv(output_paths[blend], index=False)

    #######################
    ##### Diagnostics #####
    #######################

    for name in station_weather:
        print(f"{name:<10} {dist[name]/1000:6.2f} km  IDW weight {weight_sets['IDW'][name]:.4f}")
    for blend, path in output_paths.items():
        print(f"{blend:<5} -> {path}")
