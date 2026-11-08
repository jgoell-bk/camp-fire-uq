"""Hourly station forcing for the simulation window (handoff section 4).

One table per wind station: observed wind speed and direction at that station and
Jarbo Gap NFDRS dead fuel moisture, one row per weather band, before perturbation.
"""
import platform

import numpy as np
import pandas as pd

from .util import ConfigError, read_json, sha256, utc_now, write_json, write_text

WIND_COLUMNS = ["StationName", "StationId", "DateTime", "ObservationType",
                "WindSpeed(mph)", "WindAzimuth(degrees)", "WSflag", "WAflag"]
MOISTURE = {"m1_pct": "oneHR_TL_FuelMoisture",
            "m10_pct": "tenHR_TL_FuelMoisture",
            "m100_pct": "hundredHR_TL_FuelMoisture"}
TABLE_COLUMNS = ["band", "time_pst", "ws_mph", "wd_deg", "ws_flag", "wa_flag",
                 *MOISTURE, "wind_source_time", "moisture_source_time"]
UNITS = {
    "ws_mph": "mph at 20 ft (RAWS standard height; ELMFIRE WS_AT_10M = .FALSE.)",
    "wd_deg": "degrees, direction the wind blows from, true north (no grid-declination "
              "correction, as in the baseline case)",
    "moisture": "percent (ELMFIRE DEAD_MC_IN_PERCENT = .TRUE.)",
}


def band_times(cfg):
    """UTC time of each weather band: band k (from 1) is t0 + (k - 1) * dt."""
    s = cfg["simulation"]
    t0 = pd.Timestamp(s["t0_utc"])
    if t0.tzinfo is None:
        raise ConfigError("simulation.t0_utc must include a time zone (end it with Z)")
    return pd.DatetimeIndex([t0 + pd.Timedelta(seconds=k * s["dt_weather_s"])
                             for k in range(cfg.nbands)], name="time_utc")


def window_problems(times_found, times):
    problems = []
    duplicated = times_found[times_found.duplicated()]
    if len(duplicated):
        problems.append(f"duplicate times, e.g. {duplicated.iloc[0]}")
    missing = times.difference(pd.DatetimeIndex(times_found))
    if len(missing):
        problems.append(f"{len(missing)} band times missing, first {missing[0]}")
    return problems


def read_wind(cfg, times):
    path = cfg.path(cfg["sources"]["wind_csv"])
    df = pd.read_csv(path, dtype=str)
    missing = set(WIND_COLUMNS) - set(df.columns)
    if missing:
        raise ConfigError(f"{path.name}: missing columns {sorted(missing)}")
    df["station_id"] = df["StationId"].str.strip().str.zfill(6)
    df["time"] = pd.to_datetime(df["DateTime"], utc=True, format="ISO8601")
    for col in ["WindSpeed(mph)", "WindAzimuth(degrees)", "WSflag", "WAflag"]:
        df[col] = pd.to_numeric(df[col])
    tables, summary = {}, {}
    for prefix, (name, station_id) in cfg.stations.items():
        station = df[df.station_id == station_id]
        if set(station.StationName) != {name}:
            raise ConfigError(f"{path.name}: station {station_id} is named "
                              f"{sorted(set(station.StationName))}, expected {name!r}")
        win = station[station.time.isin(times)]
        ws, wa = win["WindSpeed(mph)"], win["WindAzimuth(degrees)"]
        problems = window_problems(win.time, times)
        if (win.ObservationType != "O").any():
            problems.append("observation type other than 'O'")
        if not np.isfinite(ws).all() or (ws < 0).any():
            problems.append("missing or negative wind speed")
        if not np.isfinite(wa).all() or ((wa < 0) | (wa > 360)).any():
            problems.append("missing or out-of-range wind azimuth")
        if (win.WSflag != 0).any() or (win.WAflag != 0).any():
            problems.append("nonzero WSflag/WAflag")
        if problems:
            raise ConfigError(f"{path.name}, {name}: " + "; ".join(problems))
        win = win.set_index("time").loc[times]
        tables[prefix] = pd.DataFrame({
            "ws_mph": win["WindSpeed(mph)"].to_numpy(float),
            "wd_deg": win["WindAzimuth(degrees)"].to_numpy(float),
            "ws_flag": win.WSflag.to_numpy(int), "wa_flag": win.WAflag.to_numpy(int),
            "wind_source_time": win.DateTime.to_numpy()}, index=times)
        summary[prefix] = {"station": name, "station_id": station_id,
                           "records": len(station), "first": str(station.time.min()),
                           "last": str(station.time.max())}
    return tables, summary


def read_moisture(cfg, times):
    src = cfg["sources"]
    path = cfg.path(src["moisture_csv"])
    name, station_id = src["moisture_station"]
    days = set(times.strftime("%Y-%m-%d"))
    cols = ["stationName", "stationId", "observationTime", "NFDRType", "fuelModelType",
            *MOISTURE.values()]
    parts, records = [], 0
    for chunk in pd.read_csv(path, usecols=cols, dtype=str, chunksize=250_000):
        records += len(chunk)
        parts.append(chunk[chunk.observationTime.str[:10].isin(days)])
    df = pd.concat(parts)
    df = df[(df.stationId.str.strip().str.zfill(6) == station_id) & (df.NFDRType == "O")].copy()
    df["time"] = pd.to_datetime(df.observationTime, utc=True, format="ISO8601")
    for col in MOISTURE.values():
        df[col] = pd.to_numeric(df[col])
    win = df[df.time.isin(times)]
    if set(win.stationName) != {name}:
        raise ConfigError(f"{path.name}: station names {sorted(set(win.stationName))}, expected {name!r}")
    # The NFDRS fuel-model rows of one hour carry the same dead fuel moistures;
    # confirm that before keeping one of them.
    spread = win.groupby("time")[list(MOISTURE.values())].agg(lambda s: s.max() - s.min())
    if (spread.to_numpy() != 0).any():
        worst = spread.max().to_dict()
        raise ConfigError(f"{path.name}: NFDRS fuel models disagree in the window, max spread {worst}")
    sel = win[win.fuelModelType == src["nfdrs_fuel_model"]]
    values = sel[list(MOISTURE.values())].to_numpy(float)
    problems = window_problems(sel.time, times)
    if not np.isfinite(values).all() or (values <= 0).any():
        problems.append("missing or nonpositive moisture")
    if problems:
        raise ConfigError(f"{path.name}: " + "; ".join(problems))
    sel = sel.set_index("time").loc[times]
    table = pd.DataFrame({k: sel[v].to_numpy(float) for k, v in MOISTURE.items()}, index=times)
    table["moisture_source_time"] = sel.observationTime.to_numpy()
    summary = {"station": name, "station_id": station_id, "rows_read": records,
               "fuel_models_in_window": sorted(set(win.fuelModelType)),
               "kept_fuel_model": src["nfdrs_fuel_model"]}
    return table, summary


def check_cutoff(cfg, times):
    src = cfg["sources"]
    cutoff = pd.Timestamp(src["jarbo_gap_reliable_until_utc"])
    if times[-1] > cutoff and not src["allow_after_jarbo_gap_cutoff"]:
        raise ConfigError(
            f"the window ends {times[-1]}, after the Jarbo Gap reliability cutoff {cutoff} "
            "(Brewer and Clements 2020, sec. 3.3). Set sources.allow_after_jarbo_gap_cutoff "
            "to use it anyway.")


def check_checkpoint(cfg, tables):
    cp = cfg["checkpoint"]
    t = pd.Timestamp(cp["time_utc"])
    problems = []
    for prefix, table in tables.items():
        if t not in table.index:
            raise ConfigError(f"checkpoint time {t} is outside the simulation window")
        wind = table.loc[t, ["ws_mph", "wd_deg"]].to_numpy(float)
        if not np.allclose(wind, cp["wind"][prefix]):
            problems.append(f"{prefix} wind {wind.tolist()} != {cp['wind'][prefix]}")
        moisture = table.loc[t, list(MOISTURE)].to_numpy(float)
        if not np.allclose(moisture, cp["moisture"]):
            problems.append(f"{prefix} moisture {moisture.tolist()} != {cp['moisture']}")
    if problems:
        raise ConfigError(f"checkpoint at {t}: " + "; ".join(problems))


def prepare_forcing(cfg, write=True):
    """Build, check and (if write) save one forcing table per wind station."""
    times = band_times(cfg)
    check_cutoff(cfg, times)
    wind, wind_summary = read_wind(cfg, times)
    moisture, moisture_summary = read_moisture(cfg, times)
    tables = {}
    for prefix, w in wind.items():
        table = w.join(moisture)
        table.insert(0, "band", np.arange(1, len(times) + 1))
        table.insert(1, "time_pst", times.tz_convert("America/Los_Angeles").strftime("%Y-%m-%d %H:%M %Z"))
        tables[prefix] = table[TABLE_COLUMNS]
    check_checkpoint(cfg, tables)
    if not write:
        return tables
    out = cfg.dir("prepared")
    outputs = {}
    for prefix, table in tables.items():
        path = out / f"forcing_{prefix}.csv"
        write_text(path, table.to_csv(date_format="%Y-%m-%dT%H:%M:%SZ"))
        outputs[prefix] = {"path": str(path), "sha256": sha256(path)}
    src = cfg["sources"]
    write_json(out / "provenance.json", {
        "created_utc": utc_now(),
        "config": {"file": str(cfg.file), "sha256": cfg.sha256},
        "window": {"t0_utc": cfg["simulation"]["t0_utc"], "bands": cfg.nbands,
                   "dt_weather_s": cfg["simulation"]["dt_weather_s"],
                   "band_times_utc": list(times.strftime("%Y-%m-%dT%H:%M:%SZ"))},
        "sources": {
            "wind_csv": {"path": str(cfg.path(src["wind_csv"])),
                         "sha256": sha256(cfg.path(src["wind_csv"])),
                         "used": "WindSpeed(mph), WindAzimuth(degrees); not gusts",
                         "stations": wind_summary},
            "moisture_csv": {"path": str(cfg.path(src["moisture_csv"])),
                             "sha256": sha256(cfg.path(src["moisture_csv"])),
                             "used": list(MOISTURE.values()), **moisture_summary}},
        "checks": ["unique, complete hourly records for every band", "ObservationType O",
                   "WS/WA QA flags zero", "NFDRS fuel-model rows identical",
                   "positive moisture, nonnegative speed, azimuth in [0, 360]",
                   f"checkpoint values at {cfg['checkpoint']['time_utc']}"],
        "units": UNITS,
        "warnings": [
            "Brewer and Clements (2020), sec. 3.3: Jarbo Gap observations are questionable "
            "after about 06:00 PST 9 Nov 2018 (2018-11-09T14:00Z), possibly affected by the "
            "fire. Zero QA flags do not establish that they are reliable. The NFDRS moistures "
            "are model estimates driven by the same station's weather, not direct sensors."],
        "outputs": outputs,
        "software": {"python": platform.python_version(), "pandas": pd.__version__,
                     "numpy": np.__version__},
    })
    return tables


def load_forcing(cfg, prefix):
    """A station's prepared table; checks it still matches the configured window."""
    path = cfg.path(cfg["paths"]["prepared"]) / f"forcing_{prefix}.csv"
    if not path.exists():
        raise ConfigError(f"{path} not found; run `python -m sobolfire prepare` first")
    table = pd.read_csv(path, index_col="time_utc")
    expected = list(band_times(cfg).strftime("%Y-%m-%dT%H:%M:%SZ"))
    if list(table.index) != expected or list(table.band) != list(range(1, cfg.nbands + 1)):
        raise ConfigError(f"{path.name} does not match the configured window; rerun prepare")
    return table, sha256(path)


def provenance(cfg):
    return read_json(cfg.path(cfg["paths"]["prepared"]) / "provenance.json")
