"""Configuration, run IDs and small helpers used by the other modules."""
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path

import yaml

CODE_ROOT = Path(__file__).resolve().parents[1]     # local_sensitivity/: `python -m lsa` runs here
DEFAULT_CONFIG = CODE_ROOT / "config.yaml"          # the hourly-direction runs
# The fixed-direction runs use fixed_direction/config.yaml (python -m lsa --config ...).


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
        """A small-file folder under paths (prepared, batches, results), created if needed."""
        d = self.path(self["paths"][key])
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def stations(self):
        return self["study"]["stations"]

    @property
    def offsets(self):
        first, last, step = self["study"]["offsets_deg"]
        return list(range(first, last + 1, step))

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


def run_id(prefix, offset):
    return f"{prefix}_WD{offset:+03d}"          # JG_WD-50, JG_WD+00, JG_WD+07


def all_run_ids(cfg):
    """Every run of the study, station by station, offsets in increasing order."""
    return [run_id(p, d) for p in cfg.stations for d in cfg.offsets]


def parse_run_id(cfg, rid):
    """(station prefix, offset in degrees) of a run ID of this study."""
    prefix = rid.partition("_")[0]
    offsets = {run_id(prefix, d): d for d in cfg.offsets}
    if prefix not in cfg.stations or rid not in offsets:
        raise StudyError(f"{rid!r} is not a run of this study (IDs look like "
                         f"{run_id(next(iter(cfg.stations)), cfg.offsets[0])})")
    return prefix, offsets[rid]


def run_dir(cfg, rid):
    prefix, _ = parse_run_id(cfg, rid)
    return cfg.scratch_root / "runs" / prefix / rid


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
