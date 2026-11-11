"""Configuration, run IDs and small helpers used by the other modules."""
import hashlib
import json
import math
import os
import re
from datetime import datetime, timezone
from pathlib import Path

import yaml

CODE_ROOT = Path(__file__).resolve().parents[1]     # eki_calibration/: `python -m eki` runs here
DEFAULT_CONFIG = CODE_ROOT / "config.yaml"
BASELINE = "BASELINE"                                # zero correction: the uncorrected JG+CO mean
MEMBER_RE = re.compile(r"^I(\d{2})_M(\d{2})$")       # I00_M01: iteration 0, member 1

# Observation-operator status, one per run and observation (eki.extract)
REACHED = 0       # arrival > 0, seconds after t0
NOT_REACHED = 1   # eligible cell, not reached by SIMULATION_TSTOP (ELMFIRE -9999); kept as non-arrival
IGNITION = 2      # ELMFIRE's ignition marker (-1)
NOT_MAPPED = 3    # statically nonburnable / excluded / outside the domain under this rule
RUN_FAILED = 4    # the run did not validate (job or model failure)
STATUS_LEGEND = {REACHED: "reached", NOT_REACHED: "eligible cell not reached by SIMULATION_TSTOP",
                 IGNITION: "ignition marker", NOT_MAPPED: "not mapped (statically nonburnable, excluded or outside)",
                 RUN_FAILED: "run failed (not validated)"}


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
        """A small-file folder under paths (prepared, ensembles, batches, results), created if needed."""
        d = self.path(self["paths"][key])
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def stations(self):
        return self["study"]["stations"]

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


# Check runs at ensemble means, e.g. CHK01_I00MEAN, CHK01_I01MEAN, CHK01_I00MEAN_REP
# (the repeat has identical inputs and seed, in its own folder).
CHECK_RE = re.compile(r"^(CHK\d{2})_I\d{2}MEAN(_REP)?$")
# Sensitivity runs (not part of the calibration), e.g. SEN01_I08MEAN_PEPS0.001: a row of
# ensembles/SEN01.csv, whose extra columns may override namelist settings (case.OVERRIDABLE).
SENS_RE = re.compile(r"^(SEN\d{2})_[A-Za-z0-9._]+$")


def member_id(iteration, member):
    return f"I{iteration:02d}_M{member:02d}"


def parse_run_id(rid):
    """(iteration, member) of a member run; (None, None) for BASELINE and check runs.
    Use run_group, not this, to tell the kinds apart."""
    if rid == BASELINE or CHECK_RE.match(rid) or SENS_RE.match(rid):
        return None, None
    m = MEMBER_RE.match(rid)
    if not m or int(m.group(2)) < 1:
        raise StudyError(f"{rid!r} is not a run ID of this study "
                         "(BASELINE, I00_M01 ... I00_M32, or a check run such as CHK01_I00MEAN)")
    return int(m.group(1)), int(m.group(2))


def run_group(rid):
    """BASELINE, I00, I01, ... or CHK01, ...: the runs folder and the coefficient file."""
    if rid == BASELINE:
        return BASELINE
    c = CHECK_RE.match(rid) or SENS_RE.match(rid)
    if c:
        return c.group(1)
    iteration, _ = parse_run_id(rid)
    return f"I{iteration:02d}"


def run_dir(cfg, rid):
    return cfg.scratch_root / "runs" / run_group(rid) / rid


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
