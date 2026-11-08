"""The immutable Sobol run matrix (handoff section 3): loading, checks, run lookup."""
import csv
from dataclasses import dataclass, field, replace

from .util import ConfigError, sha256

BASE_COLUMNS = ["run_id", "wind_station", "wind_station_id", "moisture_station",
                "moisture_station_id", "design_row", "base_row", "sobol_block", "sampling_seed"]
REP = "__rep"   # JG_0005__rep2: same inputs as JG_0005, run again in its own folder


def columns(cfg):
    """Matrix CSV columns: identifiers, the factors, then derived settings."""
    return BASE_COLUMNS + cfg.factors + list(cfg["design"].get("derived") or {})


def derived_settings(cfg, coeffs):
    """Namelist settings computed from a factor as low * base**factor (design.derived)."""
    return {name: d["low"] * d["base"] ** coeffs[d["factor"]]
            for name, d in (cfg["design"].get("derived") or {}).items()}


@dataclass(frozen=True)
class Run:
    run_id: str
    prefix: str                 # JG, OP or CO
    design_row: int | None      # None for the nominal controls
    base_row: int | None
    sobol_block: str            # A, AB_<factor>, BA_<factor>, B, or CONTROL
    coeffs: dict                # factor -> value, CSV column order
    text: dict = field(default_factory=dict)   # the CSV row as read (empty for controls)
    overrides: dict = field(default_factory=dict)   # namelist settings changed for this run

    @property
    def is_control(self):
        return self.design_row is None


def block_order(factors):
    return ["A"] + [f"AB_{f}" for f in factors] + [f"BA_{f}" for f in factors] + ["B"]


def load_matrix(cfg):
    """Read the run matrix, check it, and return its Runs in file order."""
    path = cfg.path(cfg["design"]["matrix"])
    digest = sha256(path)
    if digest != cfg["design"]["sha256"]:
        raise ConfigError(f"{path.name}: sha256 {digest} differs from the configured one")
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != columns(cfg):
            raise ConfigError(f"{path.name}: columns {reader.fieldnames} != {columns(cfg)}")
        runs = [Run(r["run_id"], r["run_id"].split("_")[0], int(r["design_row"]),
                    int(r["base_row"]), r["sobol_block"],
                    {k: float(r[k]) for k in cfg.factors}, r) for r in reader]
    problems = check_matrix(cfg, runs)
    if problems:
        raise ConfigError(f"{path.name}: {len(problems)} problems: " + "; ".join(problems[:5]))
    return runs


def check_matrix(cfg, runs):
    """Station-major blocks of N*(2D+2) rows, IDs and blocks in the documented order,
    coefficients within bounds, derived columns equal to their formula, Saltelli
    structure, identical vectors in every station."""
    factors, stations = cfg.factors, cfg.stations
    blocks = block_order(factors)
    n = cfg["design"]["base_samples"] * len(blocks)
    moisture_name, moisture_id = cfg["sources"]["moisture_station"]
    problems = []
    if len(runs) != n * len(stations):
        problems.append(f"{len(runs)} rows, expected {n * len(stations)}")
    if len({r.run_id for r in runs}) != len(runs):
        problems.append("duplicate run_id")
    if len({r.text["sampling_seed"] for r in runs}) != 1:
        problems.append("more than one sampling_seed")
    reference = None
    for k, (prefix, (name, station_id)) in enumerate(stations.items()):
        block = runs[k * n:(k + 1) * n]
        for i, r in enumerate(block, start=1):
            expected = dict(run_id=f"{prefix}_{i:04d}", wind_station=name,
                            wind_station_id=station_id, moisture_station=moisture_name,
                            moisture_station_id=moisture_id, design_row=str(i),
                            base_row=str((i - 1) // len(blocks) + 1),
                            sobol_block=blocks[(i - 1) % len(blocks)])
            problems += [f"{r.run_id}: {col}={r.text[col]!r}, expected {want!r}"
                         for col, want in expected.items() if r.text[col] != want]
            for f in factors:
                lo, hi = cfg["design"]["factors"][f][:2]
                if not lo <= r.coeffs[f] <= hi:
                    problems.append(f"{r.run_id}: {f}={r.coeffs[f]} outside [{lo}, {hi}]")
            problems += [f"{r.run_id}: {k}={r.text[k]} but {v!r} from its factor"
                         for k, v in derived_settings(cfg, r.coeffs).items() if float(r.text[k]) != v]
        vectors = [r.coeffs for r in block]
        if reference is None:
            reference = vectors
            problems += check_saltelli(block, factors)
        elif vectors != reference:
            problems.append(f"{prefix}: factor vectors differ from the first station's")
    return problems


def check_saltelli(block, factors):
    """AB_k is A with factor k taken from B; BA_k is B with factor k taken from A."""
    problems = []
    size = 2 * len(factors) + 2
    for start in range(0, len(block), size):
        rows = {r.sobol_block: r.coeffs for r in block[start:start + size]}
        a, b = rows.get("A"), rows.get("B")
        if a is None or b is None:
            problems.append(f"base row {start // size + 1}: missing A or B")
            continue
        for f in factors:
            if rows.get(f"AB_{f}") != {**a, f: b[f]}:
                problems.append(f"base row {start // size + 1}: AB_{f} is not A with {f} from B")
            if rows.get(f"BA_{f}") != {**b, f: a[f]}:
                problems.append(f"base row {start // size + 1}: BA_{f} is not B with {f} from A")
    return problems


def level_text(value):
    """How a spotting_study level appears in a run ID: 0.001, 300, TRUE."""
    return str(value).upper() if isinstance(value, bool) else f"{value:g}"


def study_override(cfg, run_id, suffix):
    """{SETTING: level} for a suffix like P_EPS-0.001; the level must be listed
    under spotting_study in the config."""
    key, _, text = suffix.rpartition("-")
    levels = (cfg.get("spotting_study") or {}).get(key, [])
    matches = [v for v in levels if level_text(v) == text]
    if not matches:
        raise ConfigError(f"{run_id}: {suffix!r} is not a level listed under spotting_study")
    return {key: matches[0]}


def resolve(cfg, matrix, run_id):
    """The Run for an ID. `matrix` maps run_id -> Run. Accepted forms:
      JG_0005                          a matrix row
      JG_BASELINE                      the nominal control
      JG_0005__rep2                    a repeat, in its own folder
      JG_BASELINE__P_EPS-0.001         one spotting_study setting changed"""
    base_id, sep, suffix = run_id.partition("__")
    overrides = {}
    if sep and not (suffix.startswith(REP[2:]) and suffix[len(REP) - 2:].isdigit()):
        overrides = study_override(cfg, run_id, suffix)
    prefix, _, rest = base_id.partition("_")
    if prefix not in cfg.stations:
        raise ConfigError(f"{run_id}: unknown station prefix {prefix!r}")
    if rest == "BASELINE":
        nominal = {f: float(v[2]) for f, v in cfg["design"]["factors"].items()}
        return Run(run_id, prefix, None, None, "CONTROL", nominal, overrides=overrides)
    if base_id not in matrix:
        raise ConfigError(f"{run_id}: not in the run matrix")
    return replace(matrix[base_id], run_id=run_id, overrides=overrides)
