"""Tests for the consequential logic; none of them run ELMFIRE.

    cd <experiment folder> && python -m unittest -v tests.test_pipeline
"""
import csv
import os
import random
import re
import shutil
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

from sobolfire.assemble import assemble
from sobolfire.case import (namelist_literal, nonburnable, perturb, render_namelist,
                            template_inputs, template_text)
from sobolfire.config import Config
from sobolfire.design import check_matrix, derived_settings, load_matrix, resolve
from sobolfire.case import build_case
from sobolfire.runner import (archive_previous, chunk, copy_back, final_dir, restore_inputs, state,
                              validate_outputs)
from sobolfire.util import ConfigError, read_json, write_json
from sobolfire.weather import Grid, check_coverage, weather_grid, write_weather

# The case's fuel grid (gdalinfo of fbfm40b.tif)
FUEL = Grid("EPSG:26910", 598179.544, 4421292.763, 30.0, 1457, 1476)


class Base(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.cfg = Config()
        cls.cfg["paths"]["scratch_root"] = str(Path(cls.tmp.name) / "scratch")
        cls.cfg["paths"]["results"] = str(Path(cls.tmp.name) / "results")
        cls.matrix = load_matrix(cls.cfg)
        cls.by_id = {r.run_id: r for r in cls.matrix}

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()


class TestDesign(Base):
    def test_matrix_shape(self):
        n = 128 * (2 * len(self.cfg.factors) + 2)           # 1,792 rows per station for 6 factors
        self.assertEqual(len(self.matrix), 3 * n)
        self.assertEqual([r.run_id for r in self.matrix[:2]], ["JG_0001", "JG_0002"])
        self.assertEqual(self.matrix[n].run_id, "OP_0001")
        self.assertEqual(self.by_id["CO_0001"].text["wind_station_id"], "041211")

    def test_p_eps_from_q_spot(self):
        """P_EPS = 0.001 * 50**q_spot, log-uniform on [0.001, 0.05]; 0.01 at nominal."""
        p = [derived_settings(self.cfg, r.coeffs)["P_EPS"] for r in self.matrix]
        self.assertTrue(0.001 <= min(p) and max(p) <= 0.05)
        r = self.by_id["JG_0001"]
        self.assertEqual(derived_settings(self.cfg, r.coeffs)["P_EPS"], 0.001 * 50 ** r.coeffs["q_spot"])
        self.assertEqual(float(r.text["P_EPS"]), 0.001 * 50 ** r.coeffs["q_spot"])
        nominal = derived_settings(self.cfg, resolve(self.cfg, self.by_id, "JG_BASELINE").coeffs)["P_EPS"]
        self.assertEqual(namelist_literal(nominal), "0.01")
        tampered = [replace(self.matrix[0], text={**self.matrix[0].text, "P_EPS": "0.02"})] + self.matrix[1:]
        self.assertTrue(any("P_EPS" in p for p in check_matrix(self.cfg, tampered)))

    def test_tampering_is_caught(self):
        changed = list(self.matrix)
        r = changed[3]                                      # JG_0004, block AB_b_1_pp
        changed[3] = replace(r, coeffs={**r.coeffs, "b_1_pp": round(r.coeffs["b_1_pp"], 3)})
        self.assertTrue(any("AB_b_1_pp" in p for p in check_matrix(self.cfg, changed)))
        swapped = list(self.matrix)
        swapped[0], swapped[1] = swapped[1], swapped[0]
        self.assertTrue(check_matrix(self.cfg, swapped))

    def test_resolve(self):
        control = resolve(self.cfg, self.by_id, "OP_BASELINE")
        self.assertTrue(control.is_control)
        self.assertEqual(control.coeffs, {"a_U": 1.0, "b_theta_deg": 0.0, "b_1_pp": 0.0,
                                          "b_10_pp": 0.0, "b_100_pp": 0.0,
                                          "q_spot": 0.588591910067779})
        rep = resolve(self.cfg, self.by_id, "JG_0005__rep2")
        self.assertEqual(rep.coeffs, self.by_id["JG_0005"].coeffs)
        self.assertNotEqual(final_dir(self.cfg, rep), final_dir(self.cfg, self.by_id["JG_0005"]))
        with self.assertRaises(ConfigError):
            resolve(self.cfg, self.by_id, "JG_9999")


class TestCase(Base):
    def forcing(self, ws, wd, m):
        n = len(ws)
        return pd.DataFrame({"ws_mph": ws, "wd_deg": wd, "m1_pct": [m] * n,
                             "m10_pct": [m] * n, "m100_pct": [m] * n})

    def test_equations(self):
        out = perturb(self.forcing([10, 0], [350, 5], 7.0),
                      {"a_U": 1.2, "b_theta_deg": 20, "b_1_pp": 1, "b_10_pp": -1, "b_100_pp": -2})
        np.testing.assert_allclose(out["ws"], [12, 0])
        np.testing.assert_allclose(out["wd"], [10, 25])          # 350 + 20 wraps through north
        np.testing.assert_allclose([out["m1"][0], out["m10"][0], out["m100"][0]], [8, 6, 5])
        out = perturb(self.forcing([5], [10], 7.0),
                      {"a_U": 1, "b_theta_deg": -20, "b_1_pp": 0, "b_10_pp": 0, "b_100_pp": 0})
        np.testing.assert_allclose(out["wd"], [350])

    def test_nonphysical_rejected(self):
        with self.assertRaises(ValueError):
            perturb(self.forcing([5], [10], 1.5),
                    {"a_U": 1, "b_theta_deg": 0, "b_1_pp": 0, "b_10_pp": 0, "b_100_pp": -2})

    def test_weather_grid_coverage(self):
        check_coverage(FUEL, weather_grid(FUEL, 3000, 2))
        with self.assertRaises(ConfigError):                   # the 0/0 edge case
            check_coverage(FUEL, weather_grid(FUEL, 3000, 1))

    def test_rasters_uniform(self):
        grid = weather_grid(FUEL, 3000, 2)
        series = {v: np.arange(25, dtype=float) + k for k, v in enumerate(["ws", "wd", "m1", "m10", "m100"])}
        with tempfile.TemporaryDirectory() as d:
            write_weather(Path(d), grid, series)
            with rasterio.open(Path(d) / "wd.tif") as f:
                a = f.read()
                self.assertEqual((f.count, f.height, f.width), (25, grid.nrows, grid.ncols))
                self.assertEqual(f.crs.to_string(), "EPSG:26910")
            for k in range(25):
                self.assertTrue(np.all(a[k] == np.float32(series["wd"][k])))
            self.assertFalse(list(Path(d).glob(".*partial")))

    def test_namelist(self):
        text = render_namelist(self.cfg, Path("/w/inputs"), Path("/l/scratch"))
        self.assertNotIn("{", text)
        spotting = ".TRUE." if self.cfg["simulation"]["spotting"] else ".FALSE."
        self.assertIn(f"ENABLE_SPOTTING = {spotting}", text)
        self.assertIn("WIND_FLUCTUATIONS = .FALSE.", text)
        self.assertIn("RANDOMIZE_RANDOM_SEED = .FALSE.", text)
        self.assertIn("SIMULATION_TSTOP = 86400.0", text)
        self.assertIn("NUM_METEOROLOGY_TIMES = 25", text)
        self.assertEqual(re.findall(r"RASTER_TO_PERTURB\(\d\) = '(\w+)'", text), ["ADJ", "WAF"])
        # gfortran finds groups by '&': comments must not contain one
        self.assertEqual(len(re.findall("&", text)), 8)
        rasters, misc = template_inputs(template_text(self.cfg))
        self.assertIn("fbfm40b", rasters.values())
        self.assertEqual(misc, ["fuel_models.csv", "building_fuel_models.csv"])

    def test_spotting_study_ids(self):
        run = resolve(self.cfg, self.by_id, "OP_BASELINE__P_EPS-0.001")
        self.assertTrue(run.is_control)
        self.assertEqual(run.overrides, {"P_EPS": 0.001})
        text = render_namelist(self.cfg, Path("/w"), Path("/s"), run.overrides)
        self.assertIn("\nP_EPS = 0.001\n", text)
        self.assertNotIn("P_EPS = 0.01\n", text)
        crosswind = resolve(self.cfg, self.by_id, "JG_BASELINE__USE_CROSSWIND_DISTRIBUTION-TRUE")
        self.assertIn("USE_CROSSWIND_DISTRIBUTION = .TRUE.",
                      render_namelist(self.cfg, Path("/w"), Path("/s"), crosswind.overrides))
        self.assertEqual(resolve(self.cfg, self.by_id, "JG_BASELINE__rep2").overrides, {})
        for bad in ("JG_BASELINE__P_EPS-0.002", "JG_BASELINE__SEED-1"):   # not listed levels
            with self.assertRaises(ConfigError):
                resolve(self.cfg, self.by_id, bad)

    def test_nonburnable(self):
        self.assertEqual(list(nonburnable([91, 98, 102, 182, 256, -9999])),
                         [True, True, False, False, True, True])


class TestRunner(Base):
    def test_chunks_cover_batch_once(self):
        ids = [f"id{i}" for i in range(100)]
        parts = [chunk(ids, 24, i) for i in range(5)]
        self.assertEqual(sum(parts, []), ids)
        self.assertEqual(chunk(ids, 24, 5), [])

    def test_states_and_archive(self):
        d = Path(self.tmp.name) / "run"
        d.mkdir()
        self.assertEqual(state(d), "not_started")
        write_json(d / "attempt.json", {"status": "failed"})
        (d / "elmfire.stdout").write_text("first attempt")
        self.assertEqual(state(d), "failed")
        archive_previous(d)
        self.assertEqual(sorted(p.name for p in d.iterdir()), ["attempts"])
        self.assertEqual((d / "attempts/1/elmfire.stdout").read_text(), "first attempt")
        write_json(d / "attempt.json", {"status": "running"})
        archive_previous(d)
        self.assertTrue((d / "attempts/2/attempt.json").exists())

    def test_validate_outputs(self):
        """Pilot layout: ELMFIRE writes -1 at the ignition cell, arrival times around it."""
        self.cfg.shared_dir.mkdir(parents=True, exist_ok=True)
        write_json(self.cfg.shared_dir / "manifest.json", {
            "fuel_grid": FUEL.as_dict(), "ignition": {"row": 447, "col": 1186}})
        work = Path(self.tmp.name) / "work"
        (work / "outputs").mkdir(parents=True)
        (work / "logs").mkdir()
        (work / "logs/elmfire.stdout").write_text("")
        (work / "outputs/fire_size_stats.csv").write_text(
            "tstop (h),Wall clock time (s),Total fire area (ac),Crown fire area (ac)\n"
            "  24.00,  46.0,  9545.3,  10.0\n")
        toa = np.full((FUEL.nrows, FUEL.ncols), -9999, np.float32)
        toa[447, 1185], toa[447, 1186] = 2195, -1

        def check():
            with rasterio.open(work / "outputs/time_of_arrival_0000001_0086400.tif", "w", driver="GTiff",
                               dtype="float32", count=1, width=FUEL.ncols, height=FUEL.nrows,
                               crs=FUEL.crs, transform=FUEL.transform, nodata=-9999) as f:
                f.write(toa, 1)
            return validate_outputs(self.cfg, work)

        result = check()
        self.assertEqual((result["first_arrival_s"], result["burned_cells"]), (2195.0, 1))
        toa[447, 1186] = 1900                      # no ignition marker at the configured cell
        with self.assertRaises(ValueError):
            check()
        toa[447, 1186], toa[447, 1185] = -1, 1000  # something burned before the ignition
        with self.assertRaises(ValueError):
            check()

    def test_copy_back_keeps_weather(self):
        work, final = Path(self.tmp.name) / "cb_work", Path(self.tmp.name) / "cb_final"
        (work / "inputs").mkdir(parents=True)
        final.mkdir()
        (work / "inputs/ws.tif").write_bytes(b"raster")
        copy_back(work, final)
        self.assertEqual((final / "inputs/ws.tif").read_bytes(), b"raster")

    def test_restore_inputs(self):
        real = Path(os.path.expandvars("$SCRATCH/elmfire-class/camp_sobol_d6/shared/manifest.json"))
        if not real.exists():
            self.skipTest("needs the prepared shared inputs")
        self.cfg.shared_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(real, self.cfg.shared_dir / "manifest.json")
        run = self.by_id["JG_0007"]
        work = Path(self.tmp.name) / "work7"
        build_case(self.cfg, run, work, work / "scratch")           # what the run itself wrote
        d = final_dir(self.cfg, run)
        d.mkdir(parents=True, exist_ok=True)
        shutil.copy(work / "inputs/weather.json", d / "weather.json")
        write_json(d / "DONE", {})
        self.assertEqual(restore_inputs(self.cfg, run), "restored")
        for name in ("ws.tif", "wd.tif", "m1.tif", "m10.tif", "m100.tif"):
            self.assertEqual((d / "inputs" / name).read_bytes(), (work / "inputs" / name).read_bytes())
        self.assertEqual(restore_inputs(self.cfg, run), "present")
        rec = read_json(d / "weather.json")
        rec["raster_sha256"]["ws.tif"] = "0" * 64
        write_json(d / "weather.json", rec)
        self.assertTrue(restore_inputs(self.cfg, run).startswith("mismatch"))
        self.assertFalse((d / ".restore").exists())

    def test_assemble_order_and_completeness(self):
        jg = [r for r in self.matrix if r.prefix == "JG"]
        order = list(jg)
        random.Random(0).shuffle(order)                        # completion order is arbitrary
        for r in order[:-1]:
            d = final_dir(self.cfg, r)
            d.mkdir(parents=True)
            write_json(d / "DONE", {"toa_file": "t.tif", "burned_area_ac": float(r.design_row),
                                    "crown_area_ac": 0.0, "burned_cells": 1, "tstop_h": 24.0,
                                    "wall_s": 1.0, "stop_messages": []})
        with self.assertRaises(ConfigError):                   # one run missing
            assemble(self.cfg, "JG")
        d = final_dir(self.cfg, order[-1])
        d.mkdir(parents=True)
        write_json(d / "DONE", {"toa_file": "t.tif", "burned_area_ac": float(order[-1].design_row),
                                "crown_area_ac": 0.0, "burned_cells": 1, "tstop_h": 24.0,
                                "wall_s": 1.0, "stop_messages": []})
        table = assemble(self.cfg, "JG")
        n = len(jg)
        self.assertEqual(table.design_row.tolist(), list(range(1, n + 1)))
        self.assertEqual(table.burned_area_ac.tolist(), [float(i) for i in range(1, n + 1)])
        with open(Path(self.cfg["paths"]["results"]) / "JG_runs.csv") as f:
            first = next(csv.DictReader(f))
        self.assertEqual(float(first["a_U"]), jg[0].coeffs["a_U"])   # full precision kept
        self.assertEqual(float(first["P_EPS"]), float(jg[0].text["P_EPS"]))   # q_spot and P_EPS both kept


if __name__ == "__main__":
    unittest.main()
