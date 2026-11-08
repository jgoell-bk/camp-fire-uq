"""Hourly station weather for the 24 h window, and each station's aim.

`prepare` writes, per wind station, prepared/forcing_<S>.csv: the observed wind
speed and direction at that station and Jarbo Gap NFDRS dead fuel moisture, one
row per hourly weather band. It also writes prepared/aim.json: each station's
mean wind direction and the aim, the rotation that turns that mean into the
direction blowing from the ignition straight at the target.

Directions are as ELMFIRE reads them: degrees the wind blows from, measured from
the grid's north. The case sets no GRID_DECLINATION and the station azimuths are
written unchanged, as in the Sobol runs. At the ignition, grid north is 1.0
degree from true north, so the target bearing is computed on the grid.
"""
import numpy as np
import pandas as pd
from pyproj import Transformer

from .common import StudyError, read_json

MOISTURE = {"m1_pct": "oneHR_TL_FuelMoisture",
            "m10_pct": "tenHR_TL_FuelMoisture",
            "m100_pct": "hundredHR_TL_FuelMoisture"}
MEANS = ("vector", "unit")


def band_times(cfg):
    """UTC time of each weather band: band k (from 1) is t0 + (k - 1) * dt."""
    s = cfg["simulation"]
    t0 = pd.Timestamp(s["t0_utc"])
    times = pd.DatetimeIndex([t0 + pd.Timedelta(seconds=k * s["dt_weather_s"])
                              for k in range(cfg.nbands)], name="time_utc")
    cutoff = pd.Timestamp(cfg["sources"]["jarbo_gap_reliable_until_utc"])
    if times[-1] > cutoff:
        raise StudyError(f"the window ends {times[-1]}, after the Jarbo Gap reliability cutoff {cutoff} "
                         "(Brewer and Clements 2020, sec. 3.3)")
    return times


def exactly(rows, times, what):
    """The rows at `times`, one per time, in band order."""
    if rows.time.duplicated().any() or set(rows.time) != set(times):
        raise StudyError(f"{what}: the window's hourly records are missing or duplicated")
    return rows.set_index("time").loc[times]


def read_wind(cfg, times):
    path = cfg.path(cfg["sources"]["wind_csv"])
    df = pd.read_csv(path, dtype=str)
    df["station_id"] = df.StationId.str.strip().str.zfill(6)    # the file drops the leading zero
    df["time"] = pd.to_datetime(df.DateTime, utc=True, format="ISO8601")
    tables = {}
    for prefix, (name, station_id) in cfg.stations.items():
        w = exactly(df[(df.station_id == station_id) & df.time.isin(times)], times, f"{path.name}, {name}")
        ws, wd = pd.to_numeric(w["WindSpeed(mph)"]), pd.to_numeric(w["WindAzimuth(degrees)"])
        flags = pd.to_numeric(w.WSflag).abs() + pd.to_numeric(w.WAflag).abs()
        if (set(w.StationName) != {name} or (w.ObservationType != "O").any() or (flags != 0).any()
                or not (ws >= 0).all() or not ((wd >= 0) & (wd <= 360)).all()):
            raise StudyError(f"{path.name}, {name}: wrong name, observation type, QA flag or value in the window")
        tables[prefix] = pd.DataFrame({"ws_mph": ws.to_numpy(float), "wd_deg": wd.to_numpy(float)}, index=times)
    return tables


def read_moisture(cfg, times):
    src = cfg["sources"]
    path = cfg.path(src["moisture_csv"])
    name, station_id = src["moisture_station"]
    days = set(times.strftime("%Y-%m-%d"))
    cols = ["stationName", "stationId", "observationTime", "NFDRType", "fuelModelType", *MOISTURE.values()]
    df = pd.concat(c[c.observationTime.str[:10].isin(days)]       # 23 years of hours: read in chunks
                   for c in pd.read_csv(path, usecols=cols, dtype=str, chunksize=250_000))
    df = df[(df.stationId.str.strip().str.zfill(6) == station_id) & (df.NFDRType == "O")
            & (df.fuelModelType == src["nfdrs_fuel_model"])].copy()
    df["time"] = pd.to_datetime(df.observationTime, utc=True, format="ISO8601")
    df = exactly(df[df.time.isin(times)], times, f"{path.name}, {name}")
    table = pd.DataFrame({k: pd.to_numeric(df[v]).to_numpy(float) for k, v in MOISTURE.items()}, index=times)
    if set(df.stationName) != {name} or not (table > 0).all().all():
        raise StudyError(f"{path.name}: wrong station name, or missing or nonpositive moisture")
    return table


def forcing_tables(cfg):
    """One table per wind station: band, PST time, observed wind, moisture."""
    times = band_times(cfg)
    moisture = read_moisture(cfg, times)
    tables = {}
    for prefix, wind in read_wind(cfg, times).items():
        t = wind.join(moisture)
        t.insert(0, "band", np.arange(1, len(times) + 1))
        t.insert(1, "time_pst", times.tz_convert("America/Los_Angeles").strftime("%Y-%m-%d %H:%M %Z"))
        tables[prefix] = t
    return tables


def wrap360(deg):
    """deg in [0, 360). Plain `% 360` turns -1e-15 into 360.0."""
    deg = float(deg) % 360
    return 0.0 if deg >= 360 else deg


def mean_from(wd_deg, ws_mph, how):
    """Mean of hourly wind-from directions, in degrees [0, 360).
    vector: direction of the mean wind vector (hours weighted by speed);
    unit: direction of the mean unit vector (every hour counts the same).
    Either way 350 and 10 average to 0, where the arithmetic mean gives 180."""
    if how not in MEANS:
        raise StudyError(f"study.mean must be one of {MEANS}, not {how!r}")
    theta = np.radians(np.asarray(wd_deg, float))
    w = np.asarray(ws_mph, float) if how == "vector" else np.ones_like(theta)
    return wrap360(np.degrees(np.arctan2((w * np.sin(theta)).sum(), (w * np.cos(theta)).sum())))


def aim(cfg, tables, shared):
    """Target direction from the ignition cell, and each station's mean direction
    and aim (target - mean, in [-180, 180))."""
    ign, crs = shared["ignition"], shared["fuel_grid"]["crs"]
    x, y = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(*cfg["study"]["target_lonlat"])
    dx, dy = x - ign["cell_x"], y - ign["cell_y"]
    toward = wrap360(np.degrees(np.arctan2(dx, dy)))
    target = wrap360(toward + 180)
    how = cfg["study"]["mean"]
    stations = {}
    for prefix, t in tables.items():
        mean = mean_from(t.wd_deg, t.ws_mph, how)
        stations[prefix] = {"mean_from_deg": mean, "aim_deg": (target - mean + 180) % 360 - 180}
    return {"mean": how, "target_lonlat": cfg["study"]["target_lonlat"], "target_x_y": [x, y],
            "ignition_cell_x_y": [ign["cell_x"], ign["cell_y"]], "target_distance_m": float(np.hypot(dx, dy)),
            "target_toward_deg": toward, "target_from_deg": target, "stations": stations}


def load_prepared(cfg, prefix):
    """A station's forcing table and aim record, as written by prepare."""
    prepared = cfg.path(cfg["paths"]["prepared"])
    path = prepared / f"forcing_{prefix}.csv"
    if not path.exists() or not (prepared / "aim.json").exists():
        raise StudyError(f"{path} or aim.json not found; run `python -m lsa prepare` first")
    table = pd.read_csv(path, index_col="time_utc")
    if list(table.index) != list(band_times(cfg).strftime("%Y-%m-%dT%H:%M:%SZ")):
        raise StudyError(f"{path.name} does not match the configured window; rerun prepare")
    aims = read_json(prepared / "aim.json")
    if aims["mean"] != cfg["study"]["mean"] or aims["target_lonlat"] != cfg["study"]["target_lonlat"]:
        raise StudyError("aim.json was made with another study.mean or target; rerun prepare")
    return table, aims
