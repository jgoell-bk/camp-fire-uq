"""Each run's hourly weather: its station's wind (or the JG+CO vector mean) and Jarbo
Gap NFDRS dead fuel moisture.

The station-file reading and the vector mean are copied from
../eki_calibration/eki/forcing.py. Directions are degrees the wind blows from, as the
stations report them and as ELMFIRE reads them.
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
    """{station ID: DataFrame(ws_mph, wd_deg)} for every station in sources.stations."""
    path = cfg.path(cfg["sources"]["wind_csv"])
    df = pd.read_csv(path, dtype=str)
    df["station_id"] = df.StationId.str.strip().str.zfill(6)    # the file drops the leading zero
    df["time"] = pd.to_datetime(df.DateTime, utc=True, format="ISO8601")
    tables = {}
    for station_id, name in cfg["sources"]["stations"].items():
        w = exactly(df[(df.station_id == station_id) & df.time.isin(times)], times, f"{path.name}, {name}")
        ws, wd = pd.to_numeric(w["WindSpeed(mph)"]), pd.to_numeric(w["WindAzimuth(degrees)"])
        flags = pd.to_numeric(w.WSflag).abs() + pd.to_numeric(w.WAflag).abs()
        if (set(w.StationName) != {name} or (w.ObservationType != "O").any() or (flags != 0).any()
                or not (ws >= 0).all() or not ((wd >= 0) & (wd <= 360)).all()):
            raise StudyError(f"{path.name}, {name}: wrong name, observation type, QA flag or value in the window")
        tables[station_id] = pd.DataFrame({"ws_mph": ws.to_numpy(float), "wd_deg": wd.to_numpy(float)}, index=times)
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


def weather_table(cfg):
    """One table for all runs: run_id, time_utc and the COLUMNS, one row per run and band."""
    times = band_times(cfg)
    winds = read_wind(cfg, times)
    moisture = read_moisture(cfg, times)
    parts = []
    for rid, source in cfg["study"]["runs"].items():
        ids = [source] if isinstance(source, str) else list(source)
        missing = [i for i in ids if i not in winds]
        if missing:
            raise StudyError(f"run {rid}: station {missing} is not in sources.stations")
        if len(ids) == 1:
            ws, wd = winds[ids[0]].ws_mph.to_numpy(), winds[ids[0]].wd_deg.to_numpy()
        else:
            ws, wd = vector_mean([winds[i] for i in ids])
        t = pd.DataFrame({"run_id": rid, "band": np.arange(1, len(times) + 1),
                          "time_pst": times.tz_convert("America/Los_Angeles").strftime("%Y-%m-%d %H:%M %Z"),
                          "ws_mph": ws, "wd_deg": wd}, index=times).join(moisture)
        parts.append(t)
    table = pd.concat(parts)
    table.index = table.index.strftime("%Y-%m-%dT%H:%M:%SZ")
    table.index.name = "time_utc"
    return table


def check_reference(cfg, table):
    """Compare the MEAN run with the groupmate's MEAN csv (rounded to 4 decimals)."""
    ref = pd.read_csv(cfg.path(cfg["sources"]["reference_mean_csv"]))
    mean = table[table.run_id == "MEAN"]
    if list(ref.time_utc) != list(mean.index):
        raise StudyError("reference_mean_csv has different hours from the window")
    dws = float(np.abs(ref.ws_mph.to_numpy() - mean.ws_mph.to_numpy()).max())
    dwd = float(np.abs((ref.wd_deg.to_numpy() - mean.wd_deg.to_numpy() + 180) % 360 - 180).max())
    if dws > 1e-4 or dwd > 1e-4:
        raise StudyError(f"MEAN differs from reference_mean_csv: speed by {dws}, direction by {dwd} deg")
    return {"max_abs_ws_diff_mph": dws, "max_abs_wd_diff_deg": dwd}


def load_weather(cfg, rid):
    """The prepared hourly table of one run, in band order."""
    path = cfg.path(cfg["paths"]["prepared"]) / "weather.csv"
    if not path.exists():
        raise StudyError(f"{path} not found; run `python -m baselines prepare` first")
    table = pd.read_csv(path, index_col="time_utc")
    t = table[table.run_id == rid]
    if list(t.index) != list(band_times(cfg).strftime("%Y-%m-%dT%H:%M:%SZ")):
        raise StudyError(f"{path.name} has no complete window for {rid}; rerun prepare")
    return t
