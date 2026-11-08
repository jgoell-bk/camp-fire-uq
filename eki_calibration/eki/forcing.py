"""The baseline weather: hourly JG+CO mean wind and Jarbo Gap dead fuel moisture.

`prepare` writes prepared/baseline.csv, one row per hourly weather band. Wind is the
equal-weight mean of the two stations' wind vectors: each hour's speed and
direction become u (east) and v (north), u and v are averaged, and the mean vector
gives speed and direction again. This is the MEAN blend of the groupmate's
generate_weather.py, and prepare checks that it reproduces that file. Directions are
degrees the wind blows from, as the stations report them and as ELMFIRE reads them
(the case sets no GRID_DECLINATION, as in the earlier studies).

Reading the station files is the same code as ../local_sensitivity/lsa/forcing.py.
"""
import numpy as np
import pandas as pd

from .common import StudyError

MOISTURE = {"m1_pct": "oneHR_TL_FuelMoisture",
            "m10_pct": "tenHR_TL_FuelMoisture",
            "m100_pct": "hundredHR_TL_FuelMoisture"}
COLUMNS = ["band", "time_pst", "ws_mph", "wd_deg", "m1_pct", "m10_pct", "m100_pct"]


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


def to_uv(ws, wd):
    """Speed and direction blowing from (deg) -> u (east), v (north) of the moving air."""
    r = np.radians(wd)
    return -ws * np.sin(r), -ws * np.cos(r)


def from_uv(u, v):
    """u, v -> speed and direction blowing from, in [0, 360)."""
    wd = np.degrees(np.arctan2(-u, -v)) % 360.0
    return np.hypot(u, v), np.where(wd >= 360.0, 0.0, wd)


def vector_mean(winds):
    """Equal-weight mean of several stations' wind vectors, hour by hour."""
    uv = [to_uv(w.ws_mph.to_numpy(float), w.wd_deg.to_numpy(float)) for w in winds]
    u = np.mean([p[0] for p in uv], axis=0)
    v = np.mean([p[1] for p in uv], axis=0)
    return from_uv(u, v)


def baseline_table(cfg):
    """band, PST time, mean wind, moisture, and each station's own wind (for reference)."""
    times = band_times(cfg)
    winds = read_wind(cfg, times)
    ws, wd = vector_mean(list(winds.values()))
    t = pd.DataFrame({"band": np.arange(1, len(times) + 1),
                      "time_pst": times.tz_convert("America/Los_Angeles").strftime("%Y-%m-%d %H:%M %Z"),
                      "ws_mph": ws, "wd_deg": wd}, index=times)
    t = t.join(read_moisture(cfg, times))
    for prefix, w in winds.items():
        t[f"{prefix}_ws_mph"], t[f"{prefix}_wd_deg"] = w.ws_mph, w.wd_deg
    return t


def check_reference(cfg, table):
    """Compare with the groupmate's MEAN csv (rounded to 4 decimals). Returns the largest
    differences; raises if any exceeds the rounding."""
    ref = pd.read_csv(cfg.path(cfg["sources"]["reference_mean_csv"]))
    if list(ref.time_utc) != list(pd.to_datetime(table.index, utc=True).strftime("%Y-%m-%dT%H:%M:%SZ")):
        raise StudyError("reference_mean_csv has different hours from the baseline window")
    dws = float(np.abs(ref.ws_mph.to_numpy() - table.ws_mph.to_numpy()).max())
    dwd = float(np.abs((ref.wd_deg.to_numpy() - table.wd_deg.to_numpy() + 180) % 360 - 180).max())
    if dws > 1e-4 or dwd > 1e-4:
        raise StudyError(f"baseline differs from reference_mean_csv: speed by {dws}, direction by {dwd} deg")
    return {"max_abs_ws_diff_mph": dws, "max_abs_wd_diff_deg": dwd}


def load_baseline(cfg):
    path = cfg.path(cfg["paths"]["prepared"]) / "baseline.csv"
    if not path.exists():
        raise StudyError(f"{path} not found; run `python -m eki prepare` first")
    table = pd.read_csv(path, index_col="time_utc")
    if list(table.index) != list(band_times(cfg).strftime("%Y-%m-%dT%H:%M:%SZ")):
        raise StudyError(f"{path.name} does not match the configured window; rerun prepare")
    return table
