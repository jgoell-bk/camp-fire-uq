"""The canonical observation table: prepared/observations.csv and observations.json.

One row per GeoJSON feature, in file order, with a stable obs_id (OBS0001 ...). The
source files are only read. For every report the table records:

- the NIST TN 2135 Appendix F record it came from (eki.nist). Every field is kept, plus
  the record's sheet and row and how the match was made;
- its time on the model clock and its observation window [start, start + window];
- its location on the fuel grid, and the frozen cell mapping:
  - direct: the containing cell is eligible (landscape.eligible); that cell is used.
  - proxy: the containing cell is nonburnable; the eligible cell whose centre is
    nearest the reported point is used, if within the radius.
  - excluded: no eligible cell centre within the radius, outside the domain, or
    excluded on review.
  - review: the containing cell is water/barren/snow/agriculture and no decision is
    recorded in observations.review_decisions.

The mapping is made for the primary radius (r60_*) and for the sensitivity radius
(r90_*). It uses the fixed landscape and the recorded review decisions only, never a
run's output, and is never changed between runs.

Matching: the GeoJSON is the workbook's located rows in sheet order. Feature i is
matched to located row i, and that is accepted only if every pair agrees exactly on
longitude, latitude (6 decimals), day and minute; otherwise nothing is matched.
Where several located rows share a key, nist_match_basis says so. The GeoJSON
features in such a group are identical, so the assignment within the group changes no
observation.

Subset (sources.observations_subset_geojson, optional): a cleaned copy of the GeoJSON
that keeps only some of its points (v2, 2026-11-10). The table is built from the full
GeoJSON as above, so every NIST match and obs_id stays the same, and is then
restricted to the kept points (select_subset). Each kept point must equal an original
point exactly in longitude, latitude, day, minute and elapsed minutes. Where k of m
identical original points are kept, the first k in workbook order are used (Jacob,
2026-11-10).

event_type comes from NIST fields only:
- residual: NIST 'Residual Fire?', fire lingering after the front passed;
- fire_present: any other report; fire was observed at the place during the window.
NIST does not label arrivals. A report bounds arrival from above (arrival <= end of
window); it is not an exact arrival time.
"""
import json
import math

import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer
from rasterio.transform import xy

from .case import template_inputs, template_text
from .common import StudyError, sha256, utc_now, write_json, write_text
from .landscape import FUEL_NAMES, eligible, fuel_class
from .nist import COLUMNS, FIRE_TYPES, load_records

CRS84 = ("urn:ogc:def:crs:OGC:1.3:CRS84", "urn:ogc:def:crs:OGC::CRS84")
DOC_SEARCH_M = 300.0     # nearest_eligible_* is looked up this far, for the record
NIST_FIELDS = {"Source #": "nist_source_id", "Date": "nist_date_serial", "Time": "nist_time_serial",
               "Time Source": "nist_time_source", "Obs Window (min)": "nist_obs_window_min",
               "Location Source": "nist_location_source", "Info Source": "nist_info_source",
               "Fire Behavior Observations": "nist_description", "Type of Fire": "nist_fire_type",
               "Residual Fire?": "nist_residual_fire", "SSI": "nist_ssi"}


def read_features(path):
    with open(path) as f:
        fc = json.load(f)
    crs = fc.get("crs", {}).get("properties", {}).get("name")
    if crs not in CRS84:
        raise StudyError(f"{path}: CRS {crs!r}, expected CRS84 (lon, lat)")
    rows = []
    for i, feat in enumerate(fc["features"]):
        g, p = feat["geometry"], feat["properties"]
        if g["type"] != "Point":
            raise StudyError(f"{path}: feature {i} is a {g['type']}, not a Point")
        rows.append({"feature_index": i, "lon": float(g["coordinates"][0]), "lat": float(g["coordinates"][1]),
                     "date_original": p["Date"], "time_original": p["Time"], "elapsed_original": p["elapsed"],
                     "elapsed_minutes_original": float(p["elapsed_minutes"])})
    return pd.DataFrame(rows)


def model_times(cfg, df):
    """Local time, UTC time and model time of every report, with consistency checks."""
    tz = "America/Los_Angeles"
    local = pd.to_datetime(df.date_original.str[:10] + " " + df.time_original).dt.tz_localize(tz)
    if (local.map(lambda t: t.utcoffset()) != pd.Timedelta(hours=-8)).any():
        raise StudyError("some reports are not in PST (UTC-8)")
    origin = pd.Timestamp(cfg["sources"]["observations_origin_pst"]).tz_localize(tz)
    s = cfg["simulation"]
    t0 = pd.Timestamp(s["t0_utc"])
    if origin.tz_convert("UTC") != t0 + pd.Timedelta(seconds=s["ignition_s"]):
        raise StudyError(f"observations_origin_pst {origin} is not the model ignition time t0 + ignition_s")
    elapsed = (local - origin).dt.total_seconds()
    if not (np.isclose(elapsed / 60, df.elapsed_minutes_original, atol=1e-6)
            & np.isclose(elapsed, pd.to_timedelta(df.elapsed_original).dt.total_seconds(), atol=1e-6)).all():
        raise StudyError("elapsed times in the GeoJSON do not equal report time minus observations_origin_pst")
    utc = local.dt.tz_convert("UTC")
    return pd.DataFrame({"observed_time_pst": local.dt.strftime("%Y-%m-%dT%H:%M:%S%z"),
                         "observed_time_utc": utc.dt.strftime("%Y-%m-%dT%H:%M:%SZ"),
                         "elapsed_s": elapsed, "model_t_s": (utc - t0).dt.total_seconds()})


def match_nist(df, records):
    """NIST columns for every GeoJSON row (see the module docstring), and the records
    that are not in the GeoJSON."""
    placed = records[records.placed].reset_index(drop=True)
    if len(placed) != len(df):
        raise StudyError(f"{len(df)} GeoJSON features but {len(placed)} located workbook rows")
    geo = pd.DataFrame({"lon6": df.lon.round(6), "lat6": df.lat.round(6),
                        "day": df.date_original.str[8:10].astype(int),
                        "minute": df.time_original.str[:2].astype(int) * 60 + df.time_original.str[3:5].astype(int)})
    nist = pd.DataFrame({"lon6": placed.lon.round(6), "lat6": placed.lat.round(6),
                         "day": placed.day.astype(int), "minute": placed.minute_of_day.astype(int)})
    bad = ~(geo.to_numpy() == nist.to_numpy()).all(axis=1)
    if bad.any():
        i = int(np.argmax(bad))
        raise StudyError(f"feature {i} does not match workbook row {placed.sheet[i]}:{placed.sheet_row[i]} "
                         f"({geo.iloc[i].tolist()} vs {nist.iloc[i].tolist()}); no match is made")
    key = list(nist.columns)
    size = nist.groupby(key)["lon6"].transform("size")
    content = placed[COLUMNS].astype(str).agg("\x1f".join, axis=1)     # the whole record, all 13 fields
    differ = content.groupby([nist[k] for k in key]).transform("nunique") > 1
    out = pd.DataFrame({"nist_sheet": placed.sheet, "nist_sheet_row": placed.sheet_row,
                        **{new: placed[old] for old, new in NIST_FIELDS.items()},
                        "nist_fire_type_name": placed.fire_type.map(FIRE_TYPES),
                        "nist_window_min": placed.window_min, "nist_residual": placed.residual,
                        "nist_ssi_flag": placed.ssi, "nist_key_group_size": size.to_numpy()})
    out["nist_match_basis"] = np.select(
        [size == 1, ~differ], ["unique key", "key and order; candidates identical"],
        "key and order; candidates differ, but the features are identical")
    codes = pd.factorize(content)[0]                                   # same record -> same code
    first_pos = pd.Series(np.arange(len(content))).groupby(codes).transform("min").to_numpy()
    out["nist_record_repeat_of"] = [f"OBS{j + 1:04d}" if j != i else "" for i, j in enumerate(first_pos)]
    unplaced = records[~records.placed][["sheet", "sheet_row", "Source #", "Date", "Time",
                                         "Latitude", "Longitude", "Fire Behavior Observations"]]
    return out, unplaced.to_dict("records")


def nearest_eligible(x, y, row, col, ok, transform, search_m):
    """(row, col, distance m) of the eligible cell whose centre is nearest (x, y),
    within search_m; ties go to the lower row, then column. None if there is none."""
    half = int(math.ceil(search_m / abs(transform.a))) + 1
    r0, r1 = max(row - half, 0), min(row + half + 1, ok.shape[0])
    c0, c1 = max(col - half, 0), min(col + half + 1, ok.shape[1])
    rr, cc = np.nonzero(ok[r0:r1, c0:c1])
    if rr.size == 0:
        return None
    rr, cc = rr + r0, cc + c0
    cx, cy = xy(transform, rr, cc)
    d = np.hypot(np.asarray(cx) - x, np.asarray(cy) - y)
    i = np.lexsort((cc, rr, d))[0]
    return (int(rr[i]), int(cc[i]), float(d[i])) if d[i] <= search_m else None


def map_point(x, y, row, col, in_domain, fbfm, ok, transform, radius, review_fuels):
    """Status, assigned cell, its fuel, displacement and reason for one point."""
    none = dict(row=None, col=None, fuel=None, displacement_m=None)
    if not in_domain:
        return dict(status="excluded", reason="outside the model domain", **none)
    code = int(fbfm[row, col])
    if ok[row, col]:
        cx, cy = xy(transform, row, col)
        return dict(status="direct", row=row, col=col, fuel=code, displacement_m=float(np.hypot(cx - x, cy - y)),
                    reason=f"containing cell is eligible (fuel {code})")
    near = nearest_eligible(x, y, row, col, ok, transform, DOC_SEARCH_M)
    where = f"nearest eligible cell centre {near[2]:.1f} m" if near else f"no eligible cell within {DOC_SEARCH_M:.0f} m"
    if code in review_fuels:
        return dict(status="review", reason=f"containing cell fuel {code} ({FUEL_NAMES.get(code, 'nonburnable')}): "
                                            f"held for manual review; {where}", **none)
    if near is None or near[2] > radius:
        return dict(status="excluded", reason=f"containing cell fuel {code} is nonburnable and {where}, "
                                              f"beyond {radius:g} m", **none)
    r, c, d = near
    return dict(status="proxy", row=r, col=c, fuel=int(fbfm[r, c]), displacement_m=d,
                reason=f"containing cell fuel {code} ({FUEL_NAMES.get(code, 'nonburnable')}) is nonburnable; {where}")


def radius_prefix(radius):
    return f"r{radius:g}_"


def select_subset(df, path):
    """Positions in df (the full GeoJSON's rows, file = workbook order) of the points the
    subset GeoJSON keeps, in df order, and a record of what was kept. See the module
    docstring for the rule."""
    with open(path) as f:
        fc = json.load(f)
    crs = fc.get("crs", {}).get("properties", {}).get("name")
    if crs not in CRS84:
        raise StudyError(f"{path}: CRS {crs!r}, expected CRS84 (lon, lat)")
    key = lambda lon, lat, date, time: (round(float(lon), 6), round(float(lat), 6), str(date)[:10], str(time))
    want, minutes = {}, {}
    for i, feat in enumerate(fc["features"]):
        g, p = feat["geometry"], feat["properties"]
        if g["type"] != "Point":
            raise StudyError(f"{path}: feature {i} is a {g['type']}, not a Point")
        k = key(*g["coordinates"][:2], p["Date"], p["Time"])
        want[k] = want.get(k, 0) + 1
        minutes.setdefault(k, set()).add(round(float(p["elapsed_minutes"]), 6))
    have = {}
    for i, k in enumerate(map(key, df.lon, df.lat, df.date_original, df.time_original)):
        have.setdefault(k, []).append(i)
    missing = [k for k in want if k not in have]
    if missing:
        raise StudyError(f"{len(missing)} subset points match no original point exactly (first: {missing[0]})")
    keep, partial = [], []
    for k, n in want.items():
        rows = have[k]
        if n > len(rows):
            raise StudyError(f"the subset has {n} copies of {k}, the original only {len(rows)}")
        if minutes[k] != {round(float(df.elapsed_minutes_original[rows[0]]), 6)}:
            raise StudyError(f"subset point {k}: elapsed minutes {minutes[k]} differ from the original")
        keep += rows[:n]                                   # first n in workbook order
        if n < len(rows):
            partial.append({"key": list(k), "kept": int(n), "of": len(rows),
                            "kept_obs_ids": [df.obs_id[i] for i in rows[:n]],
                            "dropped_obs_ids": [df.obs_id[i] for i in rows[n:]]})
    keep = sorted(keep)
    dropped = sorted(set(range(len(df))) - set(keep))
    info = {"subset_source": str(path), "subset_sha256": sha256(path), "subset_features": len(fc["features"]),
            "original_features": len(df), "kept": len(keep), "dropped": len(dropped),
            "dropped_obs_ids": [df.obs_id[i] for i in dropped],
            "partly_kept_groups": partial,
            "rule": "kept = original rows equal to a subset point in lon, lat (6 decimals), day, minute and elapsed "
                    "minutes; where k of m identical rows are kept, the first k in workbook order (Jacob, 2026-11-10)"}
    return keep, info


def apply_review(df, cfg, radius, fbfm, ok, transform):
    """Mapping of the review points that have a recorded decision, for one radius."""
    out = {}
    for oid, d in (cfg["observations"].get("review_decisions") or {}).items():
        hit = df.index[df.obs_id == oid]
        if len(hit) != 1:
            raise StudyError(f"review decision for unknown {oid}")
        i = hit[0]
        if not np.allclose([df.lon[i], df.lat[i]], d["lonlat"], atol=1e-6):
            raise StudyError(f"review decision {oid}: coordinates {d['lonlat']} do not match the table")
        if d["decision"] == "proxy":
            m = map_point(df.x[i], df.y[i], df.row[i], df.col[i], df.in_domain[i], fbfm, ok, transform, radius, set())
            m["reason"] = f"reviewed: {d['note']}; {m['reason']}"
        elif d["decision"] == "excluded":
            m = dict(status="excluded", row=None, col=None, fuel=None, displacement_m=None,
                     reason=f"reviewed: {d['note']}")
        else:
            raise StudyError(f"review decision {oid}: {d['decision']!r} is not proxy or excluded")
        out[i] = m
    return out


def build_table(cfg):
    src = cfg.path(cfg["sources"]["observations_geojson"])
    wb = cfg.path(cfg["sources"]["nist_workbook"])
    df = read_features(src)
    df = pd.concat([df, model_times(cfg, df)], axis=1)
    df.insert(0, "obs_id", [f"OBS{i + 1:04d}" for i in range(len(df))])
    nist, unplaced = match_nist(df, load_records(wb))
    df = pd.concat([df, nist], axis=1)
    subset = cfg["sources"].get("observations_subset_geojson")
    subset_info = None
    if subset:
        keep, subset_info = select_subset(df, cfg.path(subset))
        df = df.iloc[keep].reset_index(drop=True)          # feature_index stays the full GeoJSON's
    s = cfg["simulation"]
    df["window_start_s"] = df.model_t_s
    df["window_end_s"] = df.model_t_s + 60.0 * df.nist_window_min
    df["in_window"] = (df.model_t_s > s["ignition_s"]) & (df.model_t_s <= s["stop_s"])
    df["window_end_after_stop"] = df.window_end_s > s["stop_s"]
    df["event_type"] = np.where(df.nist_residual, "residual", "fire_present")

    stems = template_inputs(template_text(cfg))[0]
    fuels = cfg.shared_dir / "fuels"
    if not (fuels / f"{stems['FBFM']}.tif").exists():
        raise StudyError(f"{fuels} not found; run `python -m eki prepare` first")
    with rasterio.open(fuels / f"{stems['FBFM']}.tif") as f:
        fbfm, transform, crs = f.read(1), f.transform, f.crs
    with rasterio.open(fuels / f"{stems['BLDG_AREA']}.tif") as f:
        bldg = f.read(1)
    ok = eligible(fbfm, bldg)
    x, y = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(df.lon.to_numpy(), df.lat.to_numpy())
    rows, cols = rasterio.transform.rowcol(transform, x, y)
    df["x"], df["y"], df["row"], df["col"] = x, y, np.asarray(rows), np.asarray(cols)
    df["in_domain"] = (df.row >= 0) & (df.row < fbfm.shape[0]) & (df.col >= 0) & (df.col < fbfm.shape[1])
    review = set(int(v) for v in cfg["observations"]["review_fuels"])
    df["fuel"] = [int(fbfm[r, c]) if d else None for r, c, d in zip(df.row, df.col, df.in_domain)]
    df["fuel_class"] = [fuel_class(v, review) if v is not None else "outside" for v in df.fuel]
    df["eligible_cell"] = [bool(ok[r, c]) if d else False for r, c, d in zip(df.row, df.col, df.in_domain)]

    near = [nearest_eligible(px, py, r, c, ok, transform, DOC_SEARCH_M) if (d and not e) else None
            for px, py, r, c, d, e in zip(df.x, df.y, df.row, df.col, df.in_domain, df.eligible_cell)]
    df["nearest_eligible_row"] = [n[0] if n else None for n in near]
    df["nearest_eligible_col"] = [n[1] if n else None for n in near]
    df["nearest_eligible_fuel"] = [int(fbfm[n[0], n[1]]) if n else None for n in near]
    df["nearest_eligible_m"] = [n[2] if n else None for n in near]

    o = cfg["observations"]
    radii = [float(o["proxy_radius_m"]), float(o["sensitivity_radius_m"])]
    decisions = o.get("review_decisions") or {}
    df["reviewed"] = df.obs_id.isin(list(decisions))
    df["review_decision"] = df.obs_id.map({k: v["decision"] for k, v in decisions.items()}).fillna("")
    for radius in radii:
        p = radius_prefix(radius)
        m = [map_point(px, py, r, c, d, fbfm, ok, transform, radius, review)
             for px, py, r, c, d in zip(df.x, df.y, df.row, df.col, df.in_domain)]
        for i, v in apply_review(df, cfg, radius, fbfm, ok, transform).items():
            m[i] = v
        m = pd.DataFrame(m)
        for k in ("status", "row", "col", "fuel", "displacement_m", "reason"):
            df[p + k] = m[k]

    key = ["lon", "lat", "date_original", "time_original"]
    first = df.groupby(key, sort=False).obs_id.transform("first")
    df["exact_duplicate_of"] = np.where(first != df.obs_id, first, "")
    df["reports_at_point"] = df.groupby(["lon", "lat"]).obs_id.transform("size")
    p = radius_prefix(radii[0])
    mapped = df[p + "status"].isin(["direct", "proxy"])
    rank = (df[mapped].sort_values(["model_t_s", "feature_index"])
            .groupby([p + "row", p + "col"]).cumcount() + 1)
    df[p + "cell_report_rank"] = rank.reindex(df.index)

    meta = {"created_utc": utc_now(), "source": str(src), "source_sha256": sha256(src),
            "nist_workbook": str(wb), "nist_workbook_sha256": sha256(wb),
            "nist_rows_not_in_geojson": unplaced, "subset": subset_info,
            "config_sha256": cfg.sha256, "features": len(df),
            "fuel_sha256": sha256(fuels / f"{stems['FBFM']}.tif"),
            "bldg_area_sha256": sha256(fuels / f"{stems['BLDG_AREA']}.tif"),
            "time_rule": "model_t_s = report start time (PST) - t0; equals elapsed + ignition_s "
                         f"(origin {cfg['sources']['observations_origin_pst']} PST = ignition); "
                         "window_end_s = model_t_s + 60 * NIST Obs Window (min)",
            "in_window": "ignition_s < model_t_s <= stop_s; window_end_after_stop flags windows that end later",
            "event_type": "residual = NIST Residual Fire? x (after the front); fire_present = any other report "
                          "(fire observed at the place during the window: an upper bound on arrival)",
            "match": "feature i = located workbook row i, accepted only because every pair agrees exactly on "
                     "lon, lat (6 decimals), day and minute",
            "eligible": "surface-burnable fuel model, or fuel 91 with building area > 0 (fixed landscape)",
            "mapping": "direct: containing cell eligible; proxy: nonburnable containing cell, nearest eligible "
                       "cell centre to the reported point within the radius (ties: lower row, then column); "
                       "excluded: none within the radius, outside the domain, or excluded on review; review: "
                       f"containing fuel in review_fuels {sorted(review)} without a recorded decision",
            "review_decisions": decisions,
            "radii_m": {"primary": radii[0], "sensitivity": radii[1]},
            "column_prefixes": {"primary": radius_prefix(radii[0]), "sensitivity": radius_prefix(radii[1])},
            "counts": {"in_window": int(df.in_window.sum()),
                       "window_end_after_stop": int(df.window_end_after_stop.sum()),
                       "event_type": df.event_type.value_counts().to_dict(),
                       "nist_fire_type": df.nist_fire_type.value_counts().to_dict(),
                       "nist_match_basis": df.nist_match_basis.value_counts().to_dict(),
                       "nist_record_repeats": int((df.nist_record_repeat_of != "").sum()),
                       "fuel_class": df.fuel_class.value_counts().to_dict(),
                       **{radius_prefix(r) + "status": df[radius_prefix(r) + "status"].value_counts().to_dict()
                          for r in radii},
                       "exact_duplicates": int((df.exact_duplicate_of != "").sum())}}
    return df, meta


def write_table(cfg):
    df, meta = build_table(cfg)
    out = cfg.dir("prepared")
    write_text(out / "observations.csv", df.to_csv(index=False, float_format="%.10g"))
    meta["sha256"] = sha256(out / "observations.csv")
    write_json(out / "observations.json", meta)
    return df, meta


TEXT = {c: str for c in ("exact_duplicate_of", "nist_record_repeat_of", "nist_description", "nist_source_id",
                         "nist_residual_fire", "nist_ssi", "review_decision", "nist_time_source",
                         "nist_location_source", "nist_info_source")}


def load_table(cfg):
    path = cfg.path(cfg["paths"]["prepared"]) / "observations.csv"
    if not path.exists():
        raise StudyError(f"{path} not found; run `python -m eki observations` first")
    df = pd.read_csv(path, dtype=TEXT, keep_default_na=False, na_values=[""])
    for c in TEXT:
        df[c] = df[c].fillna("")
    return df
