"""Configuration and small helpers (trimmed from ../eki_calibration/eki/common.py)."""
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path

import yaml

CODE_ROOT = Path(__file__).resolve().parents[1]     # station_baselines/: `python -m baselines` runs here
DEFAULT_CONFIG = CODE_ROOT / "config.yaml"


class StudyError(ValueError):
    """A problem with inputs or settings; reported without a traceback."""


class Config(dict):
    """The parsed config.yaml, plus its location and checksum."""

    def __init__(self, path=DEFAULT_CONFIG):
        path = Path(path).resolve()
        with open(path) as f:
            super().__init__(yaml.safe_load(f))
        self.file, self.root, self.sha256 = path, path.parent, sha256(path)

    def path(self, value):
        """Expand $VARS and ~; relative paths are relative to the config folder."""
        expanded = os.path.expandvars(os.path.expanduser(str(value)))
        if "$" in expanded:
            raise StudyError(f"unset environment variable in path {value!r}")
        p = Path(expanded)
        return p if p.is_absolute() else self.root / p

    def dir(self, key):
        d = self.path(self["paths"][key])
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def run_ids(self):
        return list(self["study"]["runs"])

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

    def run_dir(self, rid):
        if rid not in self["study"]["runs"]:
            raise StudyError(f"{rid!r} is not a run of this study ({', '.join(self.run_ids)})")
        return self.scratch_root / "runs" / rid


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def write_text(path, text):
    """Write via a temporary file and rename, so readers never see a partial file."""
    path = Path(path)
    tmp = path.with_name(f".{path.name}.partial")
    tmp.write_text(text)
    os.replace(tmp, path)


def write_json(path, obj):
    write_text(path, json.dumps(obj, indent=1, default=str) + "\n")


def read_json(path):
    with open(path) as f:
        return json.load(f)
