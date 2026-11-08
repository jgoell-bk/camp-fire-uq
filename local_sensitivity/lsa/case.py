"""Inputs shared by every run, and one run's ELMFIRE case."""
import re
import shutil
import socket

import numpy as np
import rasterio
from pyproj import Transformer

from .common import StudyError, parse_run_id, read_json, sha256, utc_now, write_json, write_text
from .forcing import load_prepared
from .weather import VARIABLES, Grid, check_coverage, grid_of, weather_grid, write_weather

URBAN = 91   # nonburnable for surface spread, but burns through the building model
# study.direction: hourly = the recorded hourly directions turned by aim + offset;
# fixed = one direction, target + offset, in every hour. Speed and moisture are the
# recorded hourly values in both.
DIRECTIONS = ("hourly", "fixed")


def nonburnable(code):
    """ELMFIRE's ISNONBURNABLE for the Rothermel model (elmfire_init.f90)."""
    code = np.asarray(code)
    return ((code >= 90) & (code <= 100)) | (code == 256) | (code <= 0) | (code > 303)


def template_text(cfg):
    return cfg.path(cfg["case"]["template"]).read_text()


def template_inputs(text):
    """Fuel/building raster stems and misc files named in the namelist template."""
    named = dict(re.findall(r"^\s*(\w+)_FILENAME\s*=\s*'([^']*)'", text, re.M))
    weather = {k: named.pop(k) for k in ("WS", "WD", "M1", "M10", "M100")}
    if sorted(weather.values()) != sorted(VARIABLES):
        raise StudyError(f"template weather files {weather} != {VARIABLES}")
    named.pop("IGNITIONS_CSV")
    misc = re.findall(r"^\s*(?:BUILDING_)?FUEL_MODEL_FILE\s*=\s*'([^']*)'", text, re.M)
    return named, misc


def baseline_dirs(cfg):
    """Fuel and misc folders of the baseline case, read from its own namelist."""
    case = cfg.path(cfg["case"]["dir"])
    text = (case / "elmfire.data.in").read_text()
    get = lambda key: re.search(rf"^\s*{key}\s*=\s*'([^']*)'", text, re.M).group(1)
    return case / get("FUELS_AND_TOPOGRAPHY_DIRECTORY"), case / get("MISCELLANEOUS_INPUTS_DIRECTORY")


def ignition_cell(cfg, fbfm_path):
    """The fuel cell containing the configured ignition point; must be burnable."""
    lon, lat = cfg["simulation"]["ignition_lonlat"]
    with rasterio.open(fbfm_path) as f:
        x, y = Transformer.from_crs("EPSG:4326", f.crs, always_xy=True).transform(lon, lat)
        row, col = (int(i) for i in f.index(x, y))
        if not (0 <= row < f.height and 0 <= col < f.width):
            raise StudyError(f"ignition point ({lon}, {lat}) is outside the fuel grid")
        code = int(f.read(1)[row, col])
        cell_x, cell_y = (float(v) for v in f.xy(row, col))
    if nonburnable(code):
        raise StudyError(f"ignition cell (row {row}, col {col}) has nonburnable fuel model {code}")
    return {"lon": lon, "lat": lat, "x": x, "y": y, "row": row, "col": col,
            "cell_x": cell_x, "cell_y": cell_y, "fbfm": code}


def prepare_shared(cfg):
    """Check the baseline case, the ignition cell and the weather grid; copy the
    case's fuel and misc inputs to <scratch_root>/shared/ (read by every run) and
    add the fixed-ignition CSV. Returns the shared manifest."""
    fuels, misc_dir = baseline_dirs(cfg)
    rasters, misc = template_inputs(template_text(cfg))
    sources = {f"fuels/{stem}.tif": fuels / f"{stem}.tif" for stem in rasters.values()}
    sources.update({f"misc/{name}": misc_dir / name for name in misc})
    missing = [str(p) for p in sources.values() if not p.exists()]
    if missing:
        raise StudyError(f"baseline case files missing: {missing}")
    exe = cfg.path(cfg["case"]["elmfire"])
    if not exe.is_file():
        raise StudyError(f"ELMFIRE executable not found: {exe}")
    fbfm = fuels / f"{rasters['FBFM']}.tif"
    fuel = grid_of(fbfm)
    for rel, path in sources.items():
        if rel.endswith(".tif") and grid_of(path) != fuel:
            raise StudyError(f"{path.name} is not on the fuel grid")
    wx = weather_grid(fuel, cfg["weather_grid"]["cellsize_m"], cfg["weather_grid"]["margin_cells"])
    check_coverage(fuel, wx)
    ignition = ignition_cell(cfg, fbfm)
    s = cfg["simulation"]
    # Columns read by ELMFIRE: case, starting weather band, x, y, area stop (acres),
    # time stop (hours). The fire starts at SIMULATION_TSTART in band 1.
    ignition_csv = ("icase,band,x,y,astop_ac,tstop_h\n"
                    f"1,1,{ignition['cell_x']:.3f},{ignition['cell_y']:.3f},9.0E9,"
                    f"{s['stop_s'] / 3600:.6f}\n")
    manifest = {"created_utc": utc_now(), "baseline_case": str(cfg.path(cfg["case"]["dir"])),
                "elmfire": {"path": str(exe), "sha256": sha256(exe)},
                "fuel_grid": fuel.as_dict(), "weather_grid": wx.as_dict(),
                "ignition": ignition, "ignition_csv": ignition_csv, "files": {}}
    dest = cfg.shared_dir
    for rel, src in sources.items():
        target = dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        tmp = target.with_name(f".{target.name}.partial")
        shutil.copy(src, tmp)
        if sha256(tmp) != sha256(src):
            raise RuntimeError(f"copy of {src} does not match the source")
        tmp.replace(target)
        manifest["files"][rel] = {"source": str(src), "sha256": sha256(target)}
    write_text(dest / "fuels" / "ignition.csv", ignition_csv)
    manifest["files"]["fuels/ignition.csv"] = {"sha256": sha256(dest / "fuels" / "ignition.csv")}
    write_json(dest / "manifest.json", manifest)
    return manifest


def render_namelist(cfg, weather_dir, scratch_dir):
    """The template filled in for one run (only paths differ between runs)."""
    s = cfg["simulation"]
    return template_text(cfg).format(
        fuels_dir=cfg.shared_dir / "fuels", misc_dir=cfg.shared_dir / "misc",
        weather_dir=weather_dir, scratch_dir=scratch_dir,
        dt_weather_s=f"{float(s['dt_weather_s']):.1f}", nbands=cfg.nbands,
        ignition_s=f"{float(s['ignition_s']):.1f}", stop_s=f"{float(s['stop_s']):.1f}",
        seed=int(s["seed"]))


def build_case(cfg, rid, workdir, scratch_dir):
    """Write a run's case into workdir: inputs/ (uniform weather rasters and
    weather.json), outputs/, logs/, elmfire.data.in and run.json. Returns run.json."""
    prefix, offset = parse_run_id(cfg, rid)
    direction = cfg["study"]["direction"]
    if direction not in DIRECTIONS:
        raise StudyError(f"study.direction must be one of {DIRECTIONS}, not {direction!r}")
    shared_manifest = cfg.shared_dir / "manifest.json"
    if not shared_manifest.exists():
        raise StudyError(f"{shared_manifest} not found; run `python -m lsa prepare` first")
    shared = read_json(shared_manifest)
    forcing, aims = load_prepared(cfg, prefix)
    aim = aims["stations"][prefix]["aim_deg"]
    inputs = workdir / "inputs"
    for d in (inputs, workdir / "outputs", workdir / "logs", scratch_dir):
        d.mkdir(parents=True, exist_ok=True)
    # The run's mean wind blows from target + offset either way (turning every hour
    # by one angle turns their mean by the same angle).
    wd_from = float(np.mod(aims["target_from_deg"] + offset, 360.0))
    if direction == "hourly":     # the recorded hourly directions, all turned by aim + offset
        wd = np.mod(forcing.wd_deg.to_numpy(float) + aim + offset, 360.0)
        wd_equation = "wd = (wd_obs + aim_deg + offset_deg) mod 360"
    else:                         # fixed: the same direction in every hour
        wd = np.full(len(forcing), wd_from)
        wd_equation = "wd = (target_from_deg + offset_deg) mod 360 in every hour"
    series = {"ws": forcing.ws_mph.to_numpy(float), "wd": wd,
              "m1": forcing.m1_pct.to_numpy(float), "m10": forcing.m10_pct.to_numpy(float),
              "m100": forcing.m100_pct.to_numpy(float)}
    raster_sha = write_weather(inputs, Grid(**shared["weather_grid"]), series)
    observed = {"ws": "ws_mph", "wd": "wd_deg", "m1": "m1_pct", "m10": "m10_pct", "m100": "m100_pct"}
    bands = [{"band": int(row.band), "time_utc": t, "time_pst": row.time_pst,
              "observed": {v: float(getattr(row, col)) for v, col in observed.items()},
              "written": {v: float(np.float32(series[v][i])) for v in VARIABLES}}
             for i, (t, row) in enumerate(forcing.iterrows())]
    write_json(inputs / "weather.json", {
        "run_id": rid, "direction": direction, "aim_deg": aim, "offset_deg": offset,
        "equations": [wd_equation, "ws, m1, m10, m100 = observed values"],
        "units": {"ws": "mph at 20 ft", "wd": "degrees the wind blows from, grid north",
                  "m1/m10/m100": "percent"},
        "grid": shared["weather_grid"], "uniform": True, "bands": bands, "raster_sha256": raster_sha})
    write_text(workdir / "elmfire.data.in", render_namelist(cfg, inputs, scratch_dir))
    name, station_id = cfg.stations[prefix]
    a = aims["stations"][prefix]
    record = {
        "run_id": rid, "station_prefix": prefix, "wind_station": name, "wind_station_id": station_id,
        "moisture_station": cfg["sources"]["moisture_station"],
        "direction": direction, "offset_deg": offset, "wd_from_deg": wd_from, "aim_deg": aim,
        "rotation_deg": aim + offset if direction == "hourly" else None,
        "mean_from_deg": a["mean_from_deg"], "mean": aims["mean"],
        "target_from_deg": aims["target_from_deg"], "target_lonlat": aims["target_lonlat"],
        "config_sha256": cfg.sha256, "template_sha256": sha256(cfg.path(cfg["case"]["template"])),
        "namelist_sha256": sha256(workdir / "elmfire.data.in"),
        "shared_manifest_sha256": sha256(shared_manifest),
        "forcing_sha256": sha256(cfg.path(cfg["paths"]["prepared"]) / f"forcing_{prefix}.csv"),
        "weather_raster_sha256": raster_sha, "elmfire": shared["elmfire"],
        "built_utc": utc_now(), "host": socket.gethostname(),
    }
    write_json(workdir / "run.json", record)
    return record
