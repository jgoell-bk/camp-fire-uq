"""The NIST TN 2135 Appendix F workbook (mds2-2954), read with the standard library.

One sheet per focus region, 13 columns each (TN 2135 Table 6):
  Source #                     source reference (TD#, 911 call, photo...); NOT unique
  Date, Time                   Excel serials; Time = start of the observation
  Time Source                  how the time was determined
  Obs Window (min)             duration, or window in which it occurred (end - start)
  Latitude, Longitude          map location (the observation's, else the observer's)
  Location Source, Info Source
  Fire Behavior Observations   the description
  Type of Fire                 S spot, V vegetation, R residential structure,
                               C commercial structure, O other/undetermined
  Residual Fire?               x = fire lingering after the fuels burned (after the front)
  SSI                          x = documented structure-to-structure ignition
A record is identified by its sheet and row. The same record can appear on several
sheets (regions overlap). The workbook's dates are in 2021 (serial 44508 = 2021-11-08):
only the day is used, as day 8 or 9 of November 2018.
"""
import re
import zipfile
import xml.etree.ElementTree as ET

import pandas as pd

from .common import StudyError

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
COLUMNS = ["Source #", "Date", "Time", "Time Source", "Obs Window (min)", "Latitude", "Longitude",
           "Location Source", "Info Source", "Fire Behavior Observations", "Type of Fire", "Residual Fire?", "SSI"]
FIRE_TYPES = {"S": "spot", "V": "vegetation", "R": "residential structure", "C": "commercial structure",
              "O": "other or undetermined"}
SERIAL_NOV8 = 44508        # the workbook's 2021-11-08


def _col(ref):
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group(0):
        n = n * 26 + ord(ch) - 64
    return n - 1


def read_sheets(path):
    """{sheet name: [rows as lists of strings]}, in workbook order."""
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS):
            shared.append("".join(t.text or "" for t in si.iter(f"{{{NS['m']}}}t")))
    rels = {r.get("Id"): r.get("Target") for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
    sheets = {}
    for s in ET.fromstring(z.read("xl/workbook.xml")).find("m:sheets", NS):
        target = rels[s.get(f"{{{NS['r']}}}id")].lstrip("/")
        rows = []
        for row in ET.fromstring(z.read(target if target.startswith("xl/") else "xl/" + target)).iter(f"{{{NS['m']}}}row"):
            cells = {}
            for c in row.findall("m:c", NS):
                v, kind = c.find("m:v", NS), c.get("t")
                if kind == "s" and v is not None:
                    cells[_col(c.get("r"))] = shared[int(v.text)]
                elif kind == "inlineStr":
                    cells[_col(c.get("r"))] = "".join(x.text or "" for x in c.iter(f"{{{NS['m']}}}t"))
                else:
                    cells[_col(c.get("r"))] = v.text if v is not None else ""
            if cells:
                rows.append([cells.get(i, "") for i in range(max(cells) + 1)])
        sheets[s.get("name")] = rows
    return sheets


def load_records(path):
    """Every non-blank data row: sheet, sheet_row (1-based, as in Excel), the 13 fields as
    text, and parsed day, minute of day, lat, lon. `placed` is False where the location or
    time is missing."""
    frames = []
    for name, rows in read_sheets(path).items():
        if [c.strip() for c in rows[0][:13]] != COLUMNS:
            raise StudyError(f"{path}: sheet {name!r} does not have the expected 13 columns")
        for i, r in enumerate(rows[1:], start=2):
            r = (r + [""] * 13)[:13]
            if any(str(c).strip() for c in r):
                frames.append({"sheet": name, "sheet_row": i, **dict(zip(COLUMNS, r))})
    df = pd.DataFrame(frames).astype({c: str for c in COLUMNS})   # fields kept exactly as written
    s = {c: df[c].str.strip() for c in ("Latitude", "Longitude", "Time", "Date", "Obs Window (min)",
                                        "Type of Fire", "Residual Fire?", "SSI")}
    df["placed"] = (s["Latitude"] != "") & (s["Longitude"] != "") & (s["Time"] != "")
    p = df.placed
    df.loc[p, "lat"] = s["Latitude"][p].astype(float)
    df.loc[p, "lon"] = s["Longitude"][p].astype(float)
    seconds = s["Time"][p].astype(float) * 86400
    if ((seconds - seconds.round()).abs() > 1e-3).any() or (seconds.round() % 60 != 0).any():
        raise StudyError(f"{path}: times are not on whole minutes")
    df.loc[p, "minute_of_day"] = (seconds / 60).round()
    df.loc[p, "day"] = s["Date"][p].astype(float).astype(int) - SERIAL_NOV8 + 8
    if not df.loc[p, "day"].isin([8, 9]).all():
        raise StudyError(f"{path}: dates other than 8 or 9 November")
    w = s["Obs Window (min)"]
    df["window_min"] = pd.to_numeric(w.where(w != ""))                 # blank stays missing
    df["fire_type"] = s["Type of Fire"]
    if not df.fire_type.isin(list(FIRE_TYPES)).all():
        raise StudyError(f"{path}: unknown Type of Fire codes {sorted(set(df.fire_type) - set(FIRE_TYPES))}")
    for flag, col in (("residual", "Residual Fire?"), ("ssi", "SSI")):
        if not s[col].isin(["", "x"]).all():
            raise StudyError(f"{path}: {col} values other than blank or x")
        df[flag] = s[col] == "x"
    return df
