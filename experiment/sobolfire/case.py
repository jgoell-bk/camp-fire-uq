"""Shared case inputs and per-run ELMFIRE cases (handoff sections 2, 5 and 6)."""
import functools
import re
import shutil
import socket

import numpy as np
import rasterio
from pyproj import Transformer

from .design import derived_settings
from .forcing import UNITS, load_forcing
from .util import ConfigError, read_json, sha256, utc_now, write_json, write_text
from .weather import VARIABLES, Grid, check_coverage, grid_of, weather_grid, write_weather

URBAN = 91   # nonburnable for surface spread, but burns through the building model


def nonburnable(code):
    """ELMFIRE's ISNONBURNABLE for the Rothermel model (elmfire_init.f90)."""
    code = np.asarray(code)
    return ((code >= 90) & (code <= 100)) | (code == 256) | (code <= 0) | (code > 303)


def perturb(forcing, coeffs):
    """The handoff's equations, applied once to the hourly series. Moisture stays in
    percent: 7 % plus 1 percentage point is 8 %."""
    out = {
        "ws": coeffs["a_U"] * forcing["ws_mph"].to_numpy(float),
        "wd": np.mod(forcing["wd_deg"].to_numpy(float) + coeffs["b_theta_deg"], 360.0),
        "m1": forcing["m1_pct"].to_numpy(float) + coeffs["b_1_pp"],
        "m10": forcing["m10_pct"].to_numpy(float) + coeffs["b_10_pp"],
        "m100": forcing["m100_pct"].to_numpy(float) + coeffs["b_100_pp"],
    }
    problems = [f"nonfinite {k}" for k, v in out.items() if not np.isfinite(v).all()]
    if (out["ws"] < 0).any():
        problems.append("negative wind speed")
    problems += [f"nonpositive {k}" for k in ("m1", "m10", "m100") if (out[k] <= 0).any()]
    if problems:
        raise ValueError("perturbed forcing rejected: " + ", ".join(problems))
    return out


def template_text(cfg):
    return cfg.path(cfg["case"]["template"]).read_text()


def template_inputs(text):
    """Fuel/building raster stems and misc files named in the namelist template."""
    named = dict(re.findall(r"^\s*(\w+)_FILENAME\s*=\s*'([^']*)'", text, re.M))
    weather = {k: named.pop(k) for k in ("WS", "WD", "M1", "M10", "M100")}
    if sorted(weather.values()) != sorted(VARIABLES):
        raise ConfigError(f"template weather files {weather} != {VARIABLES}")
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
            raise ConfigError(f"ignition point ({lon}, {lat}) is outside the fuel grid")
        code = int(f.read(1)[row, col])
        cell_x, cell_y = (float(v) for v in f.xy(row, col))
    if nonburnable(code):
        raise ConfigError(f"ignition cell (row {row}, col {col}) has nonburnable fuel model {code}")
    return {"lon": lon, "lat": lat, "x": x, "y": y, "row": row, "col": col,
            "cell_x": cell_x, "cell_y": cell_y, "fbfm": code}


def prepare_shared(cfg, write=True):
    """Check the baseline case, the ignition cell and the weather grid; if write,
    copy the case's fuel and misc inputs to scratch (read by every run) and add
    the fixed-ignition CSV. Returns the shared manifest."""
    fuels, misc_dir = baseline_dirs(cfg)
    rasters, misc = template_inputs(template_text(cfg))
    sources = {f"fuels/{stem}.tif": fuels / f"{stem}.tif" for stem in rasters.values()}
    sources.update({f"misc/{name}": misc_dir / name for name in misc})
    missing = [str(p) for p in sources.values() if not p.exists()]
    if missing:
        raise ConfigError(f"baseline case files missing: {missing}")
    exe = cfg.path(cfg["case"]["elmfire"])
    if not exe.is_file():
        raise ConfigError(f"ELMFIRE executable not found: {exe}")
    fbfm = fuels / f"{rasters['FBFM']}.tif"
    fuel = grid_of(fbfm)
    for rel, path in sources.items():
        if rel.endswith(".tif") and grid_of(path) != fuel:
            raise ConfigError(f"{path.name} is not on the fuel grid")
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
                "elmfire": {"path": str(exe), "sha256": exe_sha256(str(exe))},
                "fuel_grid": fuel.as_dict(), "weather_grid": wx.as_dict(),
                "ignition": ignition, "ignition_csv": ignition_csv, "files": {}}
    if not write:
        return manifest
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


@functools.lru_cache
def exe_sha256(path):
    return sha256(path)


def render_namelist(cfg, weather_dir, scratch_dir, overrides=None):
    """The template filled in for one run; `overrides` replaces the value of
    settings that appear on exactly one template line (e.g. P_EPS = 0.01)."""
    s = cfg["simulation"]
    text = template_text(cfg).format(
        fuels_dir=cfg.shared_dir / "fuels", misc_dir=cfg.shared_dir / "misc",
        weather_dir=weather_dir, scratch_dir=scratch_dir,
        dt_weather_s=f"{float(s['dt_weather_s']):.1f}", nbands=cfg.nbands,
        ignition_s=f"{float(s['ignition_s']):.1f}", stop_s=f"{float(s['stop_s']):.1f}",
        seed=int(s["seed"]), spotting=".TRUE." if s["spotting"] else ".FALSE.")
    for key, value in (overrides or {}).items():
        text, n = re.subn(rf"^{key} = .*$", f"{key} = {namelist_literal(value)}", text, flags=re.M)
        if n != 1:
            raise ConfigError(f"{key} is set on {n} template lines; overrides need exactly one")
    return text


def namelist_literal(value):
    """Fortran text for a setting. Nine significant digits are more than ELMFIRE's
    single-precision REAL holds, so 0.010000000000000002 is written as 0.01."""
    return (".TRUE." if value else ".FALSE.") if isinstance(value, bool) else f"{float(value):.9g}"


def build_case(cfg, run, workdir, scratch_dir):
    """Write a run's case into workdir: inputs/ (uniform weather rasters and their
    sidecar), outputs/, logs/, elmfire.data.in and run.json. Returns the run record."""
    shared_manifest = cfg.shared_dir / "manifest.json"
    if not shared_manifest.exists():
        raise ConfigError(f"{shared_manifest} not found; run `python -m sobolfire prepare` first")
    shared = read_json(shared_manifest)
    inputs = workdir / "inputs"
    for d in (inputs, workdir / "outputs", workdir / "logs", scratch_dir):
        d.mkdir(parents=True, exist_ok=True)
    forcing, forcing_sha = load_forcing(cfg, run.prefix)
    series = perturb(forcing, run.coeffs)
    raster_sha = write_weather(inputs, Grid(**shared["weather_grid"]), series)
    observed = {"ws": "ws_mph", "wd": "wd_deg", "m1": "m1_pct", "m10": "m10_pct", "m100": "m100_pct"}
    bands = [{"band": int(row.band), "time_utc": t, "time_pst": row.time_pst,
              "observed": {v: float(getattr(row, col)) for v, col in observed.items()},
              "written": {v: float(np.float32(series[v][i])) for v in VARIABLES}}
             for i, (t, row) in enumerate(forcing.iterrows())]
    write_json(inputs / "weather.json", {
        "run_id": run.run_id, "coefficients": run.coeffs, "units": UNITS,
        "equations": ["ws = a_U * ws_obs", "wd = (wd_obs + b_theta_deg) mod 360",
                      "m1 = m1_obs + b_1_pp", "m10 = m10_obs + b_10_pp",
                      "m100 = m100_obs + b_100_pp"],
        "grid": shared["weather_grid"], "uniform": True, "bands": bands,
        "raster_sha256": raster_sha, "forcing_sha256": forcing_sha})
    derived = derived_settings(cfg, run.coeffs)
    clash = sorted(set(derived) & set(run.overrides))
    if clash:
        raise ConfigError(f"{run.run_id}: {clash} set both by a factor and by a spotting_study override")
    settings = {**derived, **run.overrides}
    namelist = render_namelist(cfg, inputs, scratch_dir, settings)
    write_text(workdir / "elmfire.data.in", namelist)
    name, station_id = cfg.stations[run.prefix]
    s = cfg["simulation"]
    record = {
        "run_id": run.run_id, "station_prefix": run.prefix,
        "wind_station": name, "wind_station_id": station_id,
        "moisture_station": cfg["sources"]["moisture_station"],
        "design_row": run.design_row, "base_row": run.base_row, "sobol_block": run.sobol_block,
        "coefficients": run.coeffs, "derived_settings": derived,
        "namelist_overrides": run.overrides,
        "namelist_values": {k: namelist_literal(v) for k, v in settings.items()},
        "coefficients_csv_text": {f: run.text[f] for f in cfg.factors} if run.text else None,
        "sampling_seed": run.text.get("sampling_seed"),
        "matrix_sha256": cfg["design"]["sha256"], "config_sha256": cfg.sha256,
        "template_sha256": sha256(cfg.path(cfg["case"]["template"])),
        "namelist_sha256": sha256(workdir / "elmfire.data.in"),
        "shared_manifest_sha256": sha256(shared_manifest), "forcing_sha256": forcing_sha,
        "weather_raster_sha256": raster_sha, "elmfire": shared["elmfire"],
        "model_randomness": {"RANDOMIZE_RANDOM_SEED": False, "SEED": s["seed"],
                             "WIND_FLUCTUATIONS": False, "ENABLE_SPOTTING": s["spotting"],
                             "ignition": "fixed cell from ignition.csv"},
        "built_utc": utc_now(), "host": socket.gethostname(),
    }
    write_json(workdir / "run.json", record)
    return record
