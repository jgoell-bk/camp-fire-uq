"""Provisional spatial-block holdout split: prepared/holdout_v<version>.csv and .json.

Each report is located by its frozen sampling cell (primary mapping), or by its
containing cell if it has none. Blocks are holdout.block_m squares in the case CRS,
counted from the fuel grid's top-left corner, and a cell belongs to the block holding
its centre. Reports that must stay together are linked:
- reports sharing a sampling cell, which is automatic, since a cell is in one block;
- repeats of one NIST record;
- the groups in observations.same_event_groups.
Blocks joined by a link form one unit. The units are shuffled with holdout.seed, and
the first round(fraction x units) are held out. Model output is never used.

holdout.blocks_from (optional, v2): instead of a new draw, every block keeps its
partition in that earlier split file, so the fit/holdout blocks stay as they were
(Jacob, 2026-11-10). Every block must be in that file with one partition, and no unit
may span both partitions.

This tests spatial generalization within this fire; it is not validation on an
independent fire. The file is frozen: a different split must get a new version.
"""
import numpy as np
import pandas as pd

from .common import StudyError, read_json, sha256, utc_now, write_json, write_text
from .observations import load_table, radius_prefix
from .target import candidates


class Union:
    def __init__(self):
        self.parent = {}

    def find(self, a):
        self.parent.setdefault(a, a)
        while self.parent[a] != a:
            self.parent[a] = self.parent[self.parent[a]]
            a = self.parent[a]
        return a

    def join(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


def locate(cfg, obs):
    """Cell (row, col), its centre and its block for every report."""
    grid = read_json(cfg.shared_dir / "manifest.json")["fuel_grid"]
    p = radius_prefix(float(cfg["observations"]["proxy_radius_m"]))
    mapped = obs[p + "status"].isin(["direct", "proxy"])
    row = np.where(mapped, obs[p + "row"], np.where(obs.in_domain, obs.row, np.nan))
    col = np.where(mapped, obs[p + "col"], np.where(obs.in_domain, obs.col, np.nan))
    basis = np.where(mapped, "sampling cell", np.where(obs.in_domain, "containing cell", "outside"))
    cx = grid["x0"] + (col + 0.5) * grid["cellsize"]
    cy = grid["top"] - (row + 0.5) * grid["cellsize"]
    b = float(cfg["holdout"]["block_m"])
    bx, by = np.floor((cx - grid["x0"]) / b), np.floor((grid["top"] - cy) / b)
    block = [f"B{int(j):02d}_{int(i):02d}" if np.isfinite(i) else "" for i, j in zip(bx, by)]
    return pd.DataFrame({"obs_id": obs.obs_id, "location_basis": basis, "cell_row": row, "cell_col": col,
                         "cell_x": cx, "cell_y": cy, "block": block})


def split(cfg, obs):
    h = cfg["holdout"]
    loc = locate(cfg, obs)
    if (loc.block == "").any():
        raise StudyError(f"{int((loc.block == '').sum())} reports have no location in the domain")
    u = Union()
    for blk in loc.block:
        u.find(blk)
    block_of = dict(zip(loc.obs_id, loc.block))
    links = [[o, r] for o, r in zip(obs.obs_id, obs.nist_record_repeat_of) if r]
    groups = cfg["observations"].get("same_event_groups") or []
    for g in groups:
        if any(o not in block_of for o in g):
            raise StudyError(f"same_event_groups entry {g} names an unknown obs_id")
    spanning = 0
    for g in links + groups:
        blocks = {block_of[o] for o in g}
        spanning += len(blocks) > 1
        first = next(iter(blocks))
        for b in blocks:
            u.join(first, b)
    loc["unit"] = [u.find(b) for b in loc.block]
    units = sorted(set(loc.unit))
    if h.get("blocks_from"):
        src = cfg.path(h["blocks_from"])
        prev = pd.read_csv(src, keep_default_na=False, na_values=[""])
        roles = prev.groupby("block").partition.agg(lambda x: set(x))
        if (roles.map(len) > 1).any():
            raise StudyError(f"{src.name}: some blocks have two partitions")
        roles = roles.map(lambda x: next(iter(x)))
        unknown = sorted(set(loc.block) - set(roles.index))
        if unknown:
            raise StudyError(f"{len(unknown)} blocks are not in {src.name} (first: {unknown[0]})")
        loc["partition"] = loc.block.map(roles)
        if (loc.groupby("unit").partition.nunique() > 1).any():
            raise StudyError(f"a linked unit spans fit and holdout blocks of {src.name}")
        n_hold = int(loc[loc.partition == "holdout"].unit.nunique())
    else:
        rng = np.random.default_rng(int(h["seed"]))
        order = rng.permutation(len(units))
        n_hold = int(round(float(h["fraction"]) * len(units)))
        held = {units[i] for i in order[:n_hold]}
        loc["partition"] = np.where(loc.unit.isin(held), "holdout", "fit")
    return loc, {"blocks": int(loc.block.nunique()), "units": len(units), "holdout_units": n_hold,
                 "holdout_blocks": int(loc[loc.partition == "holdout"].block.nunique()),
                 "links_spanning_blocks": int(spanning)}


def coverage(cfg, obs, part):
    """Counts per partition from observation metadata only (no model output)."""
    ign = read_json(cfg.shared_dir / "manifest.json")["ignition"]
    rows = {}
    for name in ("fit", "holdout"):
        sel = part.partition == name
        o, loc = obs[sel], part[sel]
        cand = candidates(cfg, obs, part.partition, which=name)
        hours = o.model_t_s[o.in_window] / 3600
        dist = np.hypot(loc.cell_x - ign["cell_x"], loc.cell_y - ign["cell_y"]) / 1000
        p = radius_prefix(float(cfg["observations"]["proxy_radius_m"]))
        rows[name] = {
            "blocks": int(loc.block.nunique()), "reports": int(sel.sum()), "in_window": int(o.in_window.sum()),
            "after_window": int((~o.in_window).sum()),
            "primary_direct": int((o[p + "status"] == "direct").sum()),
            "primary_proxy": int((o[p + "status"] == "proxy").sum()),
            "primary_excluded": int((o[p + "status"] == "excluded").sum()),
            "residual": int((o.event_type == "residual").sum()),
            "candidate_target": len(cand),
            "candidate_target_proxy": int((cand[p + "status"] == "proxy").sum()),
            **{f"in_window_{a}-{b}h": int(((hours >= a) & (hours < b)).sum())
               for a, b in ((0, 3), (3, 6), (6, 9), (9, 12), (12, 18), (18, 24.01))},
            **{f"fire_type_{k}": int((o.nist_fire_type == k).sum()) for k in "VRCSO"},
            "km_from_ignition_p10_p50_p90": [round(float(v), 1) for v in np.quantile(dist, [0.1, 0.5, 0.9])],
            "x_km_range": [round(float(loc.cell_x.min() / 1000), 1), round(float(loc.cell_x.max() / 1000), 1)],
            "y_km_range": [round(float(loc.cell_y.min() / 1000), 1), round(float(loc.cell_y.max() / 1000), 1)]}
    return rows


def write_split(cfg):
    obs = load_table(cfg)
    part, info = split(cfg, obs)
    out = cfg.dir("prepared")
    v = int(cfg["holdout"]["version"])
    path = out / f"holdout_v{v}.csv"
    text = part.to_csv(index=False, float_format="%.10g")
    if path.exists() and path.read_text() != text:
        raise StudyError(f"{path} exists with different contents; the split is frozen (use a new holdout.version)")
    write_text(path, text)
    cov = coverage(cfg, obs, part)
    write_json(out / f"holdout_v{v}.json", {
        "created_utc": utc_now(), "file": path.name, "sha256": sha256(path), "settings": cfg["holdout"],
        "observations_sha256": sha256(out / "observations.csv"), **info,
        "same_event_groups": cfg["observations"].get("same_event_groups") or [],
        "rule": "1 km blocks from the fuel grid's top-left corner, located by the frozen sampling cell "
                "(else the containing cell); blocks linked by a shared NIST record or a same_event_groups "
                "entry form one unit; " + (f"each block keeps its partition in {cfg['holdout']['blocks_from']}"
                                           if cfg["holdout"].get("blocks_from") else
                                           "units shuffled with the seed; the first round(fraction x units) held out"),
        "scope": "spatial generalization within this fire, not validation on an independent fire",
        "coverage": cov})
    return part, info, cov


def load_split(cfg):
    v = int(cfg["holdout"]["version"])
    path = cfg.path(cfg["paths"]["prepared"]) / f"holdout_v{v}.csv"
    if not path.exists():
        raise StudyError(f"{path} not found; run `python -m eki holdout` first")
    return pd.read_csv(path, keep_default_na=False, na_values=[""])
