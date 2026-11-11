"""Audit 2026-11-10: area burned by 24 h, inside and outside the CAL FIRE final perimeter.

    python audit/acreage.py          (from eki_calibration/, elmfire-vnv active)

Runs at the I08 ensemble-mean parameters (CHK08_I08MEAN) of v1, v2 and v3, the
rejected v3 proposal mean (CHK09_I09MEAN) and the uncorrected run (BASELINE). For each:
- burned = cells with 0 < T <= 86400 s, plus the ignition cell (T = -1); every other
  value is NoData (-9999: not reached, or not burnable);
- inside / outside = burned cells with the perimeter mask (v3/prepared/final_perimeter.tif,
  cell-centre rule) 1 / 0; cell area 30 m x 30 m = 0.2224 ac;
- ELMFIRE's own total (DONE burned_area_ac, from fire_size_stats.csv) as a cross-check;
- burned cells within 30 cells (900 m) of the domain edge, where ELMFIRE stops spread
  (EDGEBUFFER), so a run touching it has a truncated area.
Writes audit/out/acreage.{csv,json} and acreage_map.png.
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib                                                    # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                      # noqa: E402
import numpy as np                                                   # noqa: E402
import pandas as pd                                                  # noqa: E402
import rasterio                                                      # noqa: E402

from eki.common import Config, read_json, run_dir                    # noqa: E402
from eki.overspread import ACRE_M2, perimeter_mask, raster_path      # noqa: E402

OUT = ROOT / "audit" / "out"
OUT.mkdir(parents=True, exist_ok=True)
INK, INK2, AXIS, GRID = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9"
INSIDE_C, OUTSIDE_C = "#2a78d6", "#eb6834"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.titlecolor": INK, "figure.dpi": 130})

CFG = {"v1": Config(ROOT / "config.yaml"), "v2": Config(ROOT / "v2/config.yaml"), "v3": Config(ROOT / "v3/config.yaml")}
REF = CFG["v3"]
inside, cell_ac = perimeter_mask(REF)
grid = read_json(REF.shared_dir / "manifest.json")["fuel_grid"]
EDGE = 30
RUNS = [("uncorrected", "v3", "BASELINE"), ("v1 final", "v1", "CHK08_I08MEAN"),
        ("v2 final", "v2", "CHK08_I08MEAN"), ("v3 final", "v3", "CHK08_I08MEAN"),
        ("v3 rejected I09 mean", "v3", "CHK09_I09MEAN")]

# the perimeter: polygon area in the grid CRS, the attribute, the raster
src = REF.path(REF["sources"]["final_perimeter"])
with tempfile.TemporaryDirectory() as tmp:
    subprocess.run(["ogr2ogr", "-t_srs", grid["crs"], "-nlt", "MULTIPOLYGON", f"{tmp}/p.gpkg", str(src)], check=True)
    layer = subprocess.run(["ogrinfo", "-q", f"{tmp}/p.gpkg"], capture_output=True, text=True).stdout.split()[1]
    q = subprocess.run(["ogrinfo", "-q", "-dialect", "SQLite", "-sql",
                        f'SELECT SUM(ST_Area(geom)) AS a FROM "{layer}"', f"{tmp}/p.gpkg"],
                       capture_output=True, text=True).stdout
    poly_ac = float(q.split("=")[-1]) / ACRE_M2
attr = subprocess.run(["ogrinfo", "-q", "-al", str(src)], capture_output=True, text=True).stdout
gis_acres = float([ln for ln in attr.splitlines() if "GIS_ACRES" in ln][0].split("=")[-1])
with rasterio.open(raster_path(REF)) as f:
    t = f.transform
    on_grid = (abs(t.a - grid["cellsize"]) < 1e-9 and abs(t.c - grid["x0"]) < 1e-6 and abs(t.f - grid["top"]) < 1e-6
               and f.width == grid["ncols"] and f.height == grid["nrows"] and f.crs.to_string() == grid["crs"])
perimeter = {"polygon_acres_in_grid_crs": poly_ac, "GIS_ACRES_attribute": gis_acres,
             "mask_cells": int(inside.sum()), "mask_acres": float(inside.sum() * cell_ac),
             "cell_acres": cell_ac, "mask_on_fuel_grid": bool(on_grid),
             "mask_cells_on_domain_edge_band": int(inside[:EDGE].sum() + inside[-EDGE:].sum()
                                                   + inside[:, :EDGE].sum() + inside[:, -EDGE:].sum())}

rows, burned_maps = [], {}
for label, v, rid in RUNS:
    d = run_dir(CFG[v], rid)
    done = read_json(d / "DONE")
    with rasterio.open(d / done["toa_file"]) as f:
        toa = f.read(1)
        nodata = f.nodata
    vals = np.unique(toa[toa <= 0])
    burned = ((toa > 0) & (toa <= 86400.0)) | (toa == -1.0)
    edge = int(burned[:EDGE].sum() + burned[-EDGE:].sum() + burned[:, :EDGE].sum() + burned[:, -EDGE:].sum())
    sides = {"north": int(burned[:EDGE].sum()), "south": int(burned[-EDGE:].sum()),
             "west": int(burned[:, :EDGE].sum()), "east": int(burned[:, -EDGE:].sum())}
    total, ins = float(burned.sum() * cell_ac), float((burned & inside).sum() * cell_ac)
    rows.append({"run": label, "run_id": f"{v}/{rid}", "statistic": "single run" if rid == "BASELINE" else
                 "single run at ensemble-mean parameters", "total_acres": total, "inside_acres": ins,
                 "outside_acres": total - ins, "inside_share_of_final_perimeter": ins / perimeter["mask_acres"],
                 "outside_share_of_burned": (total - ins) / total,
                 "elmfire_reported_acres": float(done["burned_area_ac"]),
                 "nan_cells": int(np.isnan(toa).sum()), "nonpositive_values": [float(x) for x in vals],
                 "nodata": float(nodata), "burned_cells_within_900m_of_edge": edge, "edge_sides": sides})
    burned_maps[label] = burned
tab = pd.DataFrame(rows)
tab.drop(columns=["edge_sides", "nonpositive_values"]).to_csv(OUT / "acreage.csv", index=False, float_format="%.1f")
(OUT / "acreage.json").write_text(json.dumps({"perimeter": perimeter, "runs": rows}, indent=1))

# map: burned by 24 h, inside (blue) and outside (orange) the perimeter
r_, c_ = np.nonzero(inside | np.any(np.stack(list(burned_maps.values())), axis=0))
r0, r1, c0, c1 = max(r_.min() - 20, 0), r_.max() + 20, max(c_.min() - 20, 0), c_.max() + 20
x0, top, cs = grid["x0"], grid["top"], grid["cellsize"]
ext = [(x0 + c0 * cs) / 1000, (x0 + c1 * cs) / 1000, (top - r1 * cs) / 1000, (top - r0 * cs) / 1000]
show = RUNS[:4]
fig, axes = plt.subplots(1, 4, figsize=(14, 4.6))
for ax, (label, _, _) in zip(axes, show):
    b, s = burned_maps[label][r0:r1, c0:c1], inside[r0:r1, c0:c1]
    img = np.full(b.shape, np.nan)
    img[b & s], img[b & ~s] = 0, 1
    ax.imshow(img, extent=ext, cmap=matplotlib.colors.ListedColormap([INSIDE_C, OUTSIDE_C]), vmin=0, vmax=1,
              interpolation="nearest")
    ax.contour(np.flipud(s.astype(float)), levels=[0.5], colors=INK, linewidths=0.7, extent=ext, origin="lower")
    r = tab[tab.run == label].iloc[0]
    ax.set_title(f"{label}\ninside {r.inside_acres:,.0f} ac, outside {r.outside_acres:,.0f} ac", loc="left")
    ax.set_xticks([]), ax.set_yticks([]), ax.set_aspect("equal")
handles = [matplotlib.patches.Patch(color=INSIDE_C, label="burned by 24 h, inside the final perimeter"),
           matplotlib.patches.Patch(color=OUTSIDE_C, label="burned by 24 h, outside it"),
           matplotlib.lines.Line2D([], [], color=INK, lw=0.7, label="CAL FIRE final perimeter (contained 26 Nov)")]
fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False)
fig.suptitle("Area burned by 06:00 PST 9 Nov (runs at the I08 ensemble-mean parameters; uncorrected = zero corrections)",
             x=0.01, ha="left", color=INK, fontsize=10)
fig.tight_layout(rect=(0, 0.07, 1, 0.95))
fig.savefig(OUT / "acreage_map.png")
plt.close(fig)

print(json.dumps(perimeter, indent=1))
print(tab[["run", "total_acres", "inside_acres", "outside_acres", "inside_share_of_final_perimeter",
           "elmfire_reported_acres", "burned_cells_within_900m_of_edge", "nan_cells", "nonpositive_values"]]
      .to_string(index=False))
for r in rows:
    if r["burned_cells_within_900m_of_edge"]:
        print(r["run"], "edge cells by side:", r["edge_sides"])
