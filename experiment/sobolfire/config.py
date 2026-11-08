"""Experiment configuration: one YAML file, paths resolved against its folder."""
import math
import os
from pathlib import Path

import yaml

from .util import ConfigError, sha256

DEFAULT = Path(__file__).resolve().parents[1] / "config.yaml"

REQUIRED = {
    "design": ["matrix", "sha256", "base_samples", "factors", "stations"],
    "sources": ["wind_csv", "moisture_csv", "moisture_station", "nfdrs_fuel_model",
                "jarbo_gap_reliable_until_utc", "allow_after_jarbo_gap_cutoff"],
    "simulation": ["t0_utc", "ignition_s", "stop_s", "dt_weather_s", "ignition_lonlat",
                   "seed", "spotting"],
    "weather_grid": ["cellsize_m", "margin_cells"],
    "case": ["dir", "template", "elmfire"],
    "paths": ["prepared", "batches", "results", "scratch_root"],
    "checkpoint": ["time_utc", "wind", "moisture"],
    "pilot": ["base_rows", "repeats"],
    "slurm": ["runs_per_task", "max_concurrent_tasks"],
}


class Config(dict):
    """The parsed YAML, plus the file's location and checksum."""

    def __init__(self, path=DEFAULT):
        path = Path(path).resolve()
        with open(path) as f:
            super().__init__(yaml.safe_load(f))
        self.file = path
        self.root = path.parent
        self.sha256 = sha256(path)
        missing = [f"{section}.{key}" for section, keys in REQUIRED.items() for key in keys
                   if (self.get(section) or {}).get(key) is None]
        if missing:
            raise ConfigError(f"{path.name}: missing settings {missing}")

    def path(self, value):
        """Expand $VARS and ~; relative paths are relative to the config folder."""
        expanded = os.path.expandvars(os.path.expanduser(str(value)))
        if "$" in expanded:
            raise ConfigError(f"unset environment variable in path {value!r}")
        p = Path(expanded)
        return p if p.is_absolute() else self.root / p

    @property
    def factors(self):
        return list(self["design"]["factors"])

    @property
    def stations(self):
        return self["design"]["stations"]

    @property
    def nbands(self):
        """Hourly bands from t0 through SIMULATION_TSTOP, both ends included."""
        s = self["simulation"]
        return math.ceil(s["stop_s"] / s["dt_weather_s"]) + 1

    @property
    def scratch_root(self):
        return self.path(self["paths"]["scratch_root"])

    @property
    def shared_dir(self):
        return self.scratch_root / "shared"

    def dir(self, key):
        """A small-file output folder under paths (prepared, batches, results)."""
        d = self.path(self["paths"][key])
        d.mkdir(parents=True, exist_ok=True)
        return d
