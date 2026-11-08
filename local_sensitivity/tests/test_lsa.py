"""Tests for the consequential logic; none of them run ELMFIRE.

    cd <local_sensitivity folder> && python -m unittest -v tests.test_lsa
"""
import os
import re
import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer

from lsa.case import build_case, render_namelist
from lsa.collect import collect
from lsa.common import (CODE_ROOT, Config, StudyError, all_run_ids, parse_run_id, read_json, run_dir,
                        write_json)
from lsa.forcing import aim, mean_from
from lsa.run import archive_previous, chunk, run_one, state, validate_outputs
from lsa.weather import Grid

FUEL = Grid("EPSG:26910", 598179.544, 4421292.763, 30.0, 1457, 1476)   # the case's fuel grid
REAL_SHARED = Config().shared_dir / "manifest.json"
REAL_AIM = Config().path(Config()["paths"]["prepared"]) / "aim.json"


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.cfg = Config()
        self.cfg["paths"]["scratch_root"] = str(Path(self.tmp.name) / "scratch")
        self.cfg["paths"]["results"] = str(Path(self.tmp.name) / "results")

    def tearDown(self):
        self.tmp.cleanup()


class TestIds(Base):
    def test_run_ids(self):
        ids = all_run_ids(self.cfg)
        self.assertEqual(len(ids), 3 * 101)
        self.assertEqual([ids[0], ids[50], ids[100], ids[101]], ["JG_WD-50", "JG_WD+00", "JG_WD+50", "OP_WD-50"])
        self.assertEqual(len(set(ids)), len(ids))
        self.assertEqual(parse_run_id(self.cfg, "CO_WD-07"), ("CO", -7))
        self.assertEqual(run_dir(self.cfg, "OP_WD+05"), self.cfg.scratch_root / "runs/OP/OP_WD+05")
        for bad in ("JG_WD+51", "JG_WD+7", "JG_WD0", "XX_WD+00", "JG_0001", "JG_BASELINE"):
            with self.assertRaises(StudyError):
                parse_run_id(self.cfg, bad)


class TestAim(Base):
    def test_mean_direction(self):
        self.assertAlmostEqual(mean_from([350, 10], [5, 5], "unit"), 0, places=9)      # not 180
        self.assertAlmostEqual(mean_from([350, 10], [5, 5], "vector"), 0, places=9)
        self.assertAlmostEqual(mean_from([0, 90], [10, 10], "unit"), 45, places=9)
        self.assertAlmostEqual(mean_from([0, 90], [10, 0], "unit"), 45, places=9)     # calm hour counts
        self.assertAlmostEqual(mean_from([0, 90], [10, 0], "vector"), 0, places=9)    # calm hour doesn't
        with self.assertRaises(StudyError):
            mean_from([0], [1], "arithmetic")

    def test_aim_points_at_target(self):
        """Ignition 1 km north of the target: wind must blow from 0 (north)."""
        x, y = Transformer.from_crs("EPSG:4326", "EPSG:26910", always_xy=True).transform(
            *self.cfg["study"]["target_lonlat"])
        shared = {"ignition": {"cell_x": x, "cell_y": y + 1000.0}, "fuel_grid": {"crs": "EPSG:26910"}}
        tables = {"JG": pd.DataFrame({"wd_deg": [30.0, 30.0], "ws_mph": [5.0, 10.0]}),
                  "OP": pd.DataFrame({"wd_deg": [340.0, 20.0], "ws_mph": [5.0, 5.0]}),
                  "CO": pd.DataFrame({"wd_deg": [200.0, 200.0], "ws_mph": [1.0, 1.0]})}
        a = aim(self.cfg, tables, shared)
        turn = lambda deg: (deg + 180) % 360 - 180                                    # angle on the circle
        self.assertLess(abs(turn(a["target_from_deg"])), 1e-6)
        self.assertTrue(0 <= a["target_from_deg"] < 360)
        self.assertAlmostEqual(a["target_distance_m"], 1000, places=6)
        self.assertAlmostEqual(a["stations"]["JG"]["aim_deg"], -30, places=6)
        self.assertLess(abs(turn(a["stations"]["OP"]["aim_deg"])), 1e-6)
        self.assertAlmostEqual(a["stations"]["CO"]["aim_deg"], 160, places=6)        # shortest turn, in [-180, 180)


@unittest.skipUnless(REAL_SHARED.exists() and REAL_AIM.exists(), "needs `python -m lsa prepare` first")
class TestBuild(Base):
    def test_rotated_rasters(self):
        """WD+00: the written directions' mean is the target direction; WD+10 is 10
        degrees more in every band; speed and moisture are the observed values."""
        self.cfg.shared_dir.mkdir(parents=True)
        shutil.copy(REAL_SHARED, self.cfg.shared_dir / "manifest.json")
        work = Path(self.tmp.name) / "build"
        written = {}
        for rid in ("JG_WD+00", "JG_WD+10", "OP_WD-50"):
            rec = build_case(self.cfg, rid, work / rid, work / rid / "scratch")
            with rasterio.open(work / rid / "inputs/wd.tif") as f:
                wd = f.read()
                self.assertTrue((wd == wd[:, :1, :1]).all())                      # uniform in space
            with rasterio.open(work / rid / "inputs/ws.tif") as f:
                ws = f.read()[:, 0, 0]
            written[rid] = wd[:, 0, 0].astype(float)
            obs = pd.read_csv(self.cfg.path(self.cfg["paths"]["prepared"]) / f"forcing_{rec['station_prefix']}.csv")
            np.testing.assert_allclose(ws, obs.ws_mph, rtol=0, atol=1e-5)
            np.testing.assert_allclose(written[rid], np.mod(obs.wd_deg + rec["rotation_deg"], 360).astype(np.float32))
            mean = mean_from(written[rid], ws, self.cfg["study"]["mean"])
            diff = (mean - rec["target_from_deg"] - rec["offset_deg"] + 180) % 360 - 180
            self.assertLess(abs(diff), 1e-3)
            self.assertEqual(rec["rotation_deg"], rec["aim_deg"] + rec["offset_deg"])
        turn = np.mod(written["JG_WD+10"] - written["JG_WD+00"], 360)
        np.testing.assert_allclose(np.where(turn > 180, turn - 360, turn), 10, atol=1e-3)
        with rasterio.open(work / "JG_WD+00/inputs/m100.tif") as f:
            obs = pd.read_csv(self.cfg.path(self.cfg["paths"]["prepared"]) / "forcing_JG.csv")
            np.testing.assert_allclose(f.read()[:, 0, 0], obs.m100_pct, atol=1e-5)

    def test_fixed_direction(self):
        """fixed: every band blows from target + offset; speed and moisture stay hourly."""
        self.cfg.shared_dir.mkdir(parents=True)
        shutil.copy(REAL_SHARED, self.cfg.shared_dir / "manifest.json")
        self.cfg["study"]["direction"] = "fixed"
        work = Path(self.tmp.name) / "fixed"
        for rid in ("JG_WD+10", "CO_WD-50", "OP_WD+00"):
            rec = build_case(self.cfg, rid, work / rid, work / rid / "scratch")
            want = np.float32((rec["target_from_deg"] + rec["offset_deg"]) % 360)
            self.assertAlmostEqual(rec["wd_from_deg"], float(want), places=4)
            self.assertIsNone(rec["rotation_deg"])
            obs = pd.read_csv(self.cfg.path(self.cfg["paths"]["prepared"]) / f"forcing_{rec['station_prefix']}.csv")
            with rasterio.open(work / rid / "inputs/wd.tif") as f:
                self.assertTrue((f.read() == want).all())                        # every band, every cell
            for var, col in (("ws", "ws_mph"), ("m1", "m1_pct"), ("m100", "m100_pct")):
                with rasterio.open(work / rid / f"inputs/{var}.tif") as f:
                    np.testing.assert_allclose(f.read()[:, 0, 0], obs[col], atol=1e-5)
        self.cfg["study"]["direction"] = "constant"
        with self.assertRaises(StudyError):
            build_case(self.cfg, "JG_WD+00", work / "bad", work / "bad" / "scratch")

    def test_run_one_reaches_validation(self):
        """run_one through build, log lines and the ELMFIRE call, for both kinds of run,
        with /bin/true standing in for ELMFIRE. A run must get as far as checking the
        outputs (which then fails: there are none). On 2026-10-28 every fixed-direction
        run crashed before ELMFIRE, on a log line that formatted rotation_deg = None."""
        self.cfg.shared_dir.mkdir(parents=True)
        shutil.copy(REAL_SHARED, self.cfg.shared_dir / "manifest.json")
        self.cfg["case"]["elmfire"] = "/bin/true"
        saved = os.environ.get("L_SCRATCH_JOB")
        os.environ["L_SCRATCH_JOB"] = str(Path(self.tmp.name) / "node")
        try:
            for direction in ("hourly", "fixed"):
                self.cfg["study"]["direction"] = direction
                self.assertEqual(run_one(self.cfg, "CO_WD+05"), "failed")
                attempt = read_json(run_dir(self.cfg, "CO_WD+05") / "attempt.json")
                self.assertEqual(attempt["exit_code"], 0)                      # ELMFIRE was called
                self.assertIn("time_of_arrival", attempt["error"])             # failed at the output check
                self.assertEqual(read_json(run_dir(self.cfg, "CO_WD+05") / "run.json")["direction"], direction)
        finally:
            if saved is None:
                os.environ.pop("L_SCRATCH_JOB")
            else:
                os.environ["L_SCRATCH_JOB"] = saved


class TestFixedConfig(unittest.TestCase):
    def test_same_as_hourly_but_direction(self):
        """fixed_direction/config.yaml differs from config.yaml only where intended."""
        hourly, fixed = Config(), Config(CODE_ROOT / "fixed_direction/config.yaml")
        self.assertEqual((hourly["study"]["direction"], fixed["study"]["direction"]), ("hourly", "fixed"))
        self.assertEqual({k: v for k, v in hourly["study"].items() if k != "direction"},
                         {k: v for k, v in fixed["study"].items() if k != "direction"})
        for section in ("simulation", "weather_grid"):
            self.assertEqual(hourly[section], fixed[section])
        for key in ("wind_csv", "moisture_csv"):
            self.assertEqual(hourly.path(hourly["sources"][key]).resolve(), fixed.path(fixed["sources"][key]).resolve())
            self.assertTrue(fixed.path(fixed["sources"][key]).exists())
        for key in ("moisture_station", "nfdrs_fuel_model", "jarbo_gap_reliable_until_utc"):
            self.assertEqual(hourly["sources"][key], fixed["sources"][key])
        for key in ("dir", "template", "elmfire"):
            self.assertEqual(hourly.path(hourly["case"][key]).resolve(), fixed.path(fixed["case"][key]).resolve())
        self.assertNotEqual(hourly.scratch_root, fixed.scratch_root)
        self.assertNotEqual(hourly.path(hourly["paths"]["results"]), fixed.path(fixed["paths"]["results"]))
        for cfg in (hourly, fixed):
            self.assertTrue((CODE_ROOT / cfg["slurm"]["script"]).exists())
        self.assertEqual(hourly["slurm"]["runs_per_task"], fixed["slurm"]["runs_per_task"])


class TestNamelist(Base):
    def test_namelist(self):
        text = render_namelist(self.cfg, Path("/w/inputs"), Path("/l/scratch"))
        self.assertNotIn("{", text)
        for line in ("ENABLE_SPOTTING = .TRUE.", "P_EPS = 0.01", "WIND_FLUCTUATIONS = .FALSE.",
                     "RANDOMIZE_RANDOM_SEED = .FALSE.", "SEED = 2023", "SIMULATION_TSTART = 1800.0",
                     "SIMULATION_TSTOP = 86400.0", "NUM_METEOROLOGY_TIMES = 25", "NUM_ENSEMBLE_MEMBERS = 1"):
            self.assertIn(f"\n{line}\n", text)
        self.assertEqual(re.findall(r"RASTER_TO_PERTURB\(\d\) = '(\w+)'", text), ["ADJ", "WAF"])
        self.assertEqual(len(re.findall("&", text)), 8)    # gfortran finds groups by '&'


class TestRunner(Base):
    def test_chunks_cover_batch_once(self):
        ids = all_run_ids(self.cfg)
        parts = [chunk(ids, self.cfg["slurm"]["runs_per_task"], i) for i in range(3)]
        self.assertEqual(sum(parts, []), ids)
        self.assertEqual([p[0] for p in parts], ["JG_WD-50", "OP_WD-50", "CO_WD-50"])   # one station per task

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

    def test_validate_outputs(self):
        """ELMFIRE writes -1 at the ignition cell, arrival times around it."""
        self.cfg.shared_dir.mkdir(parents=True)
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


class TestCollect(Base):
    def test_offset_order_and_completeness(self):
        ids = [i for i in all_run_ids(self.cfg) if i.startswith("JG")]
        for rid in reversed(ids):                  # completion order is arbitrary
            d = run_dir(self.cfg, rid)
            d.mkdir(parents=True)
            offset = parse_run_id(self.cfg, rid)[1]
            write_json(d / "run.json", {"aim_deg": 30.5, "rotation_deg": 30.5 + offset,
                                        "mean_from_deg": 35.5, "target_from_deg": 66.0})
            if rid != "JG_WD+13":
                write_json(d / "DONE", {"toa_file": "t.tif", "burned_area_ac": float(offset),
                                        "crown_area_ac": 0.0, "burned_cells": 1, "first_arrival_s": 2000.0,
                                        "tstop_h": 24.0, "wall_s": 1.0, "stop_messages": []})
        with self.assertRaises(StudyError):        # one run missing
            collect(self.cfg, "JG")
        write_json(run_dir(self.cfg, "JG_WD+13") / "DONE", {
            "toa_file": "t.tif", "burned_area_ac": 13.0, "crown_area_ac": 0.0, "burned_cells": 1,
            "first_arrival_s": 2000.0, "tstop_h": 24.0, "wall_s": 1.0, "stop_messages": []})
        table = collect(self.cfg, "JG")
        self.assertEqual(table.offset_deg.tolist(), list(range(-50, 51)))
        self.assertEqual(table.burned_area_ac.tolist(), [float(d) for d in range(-50, 51)])
        saved = pd.read_csv(Path(self.cfg["paths"]["results"]) / "JG_runs.csv")
        self.assertEqual(saved.run_id.tolist(), ids)
        self.assertEqual(saved.rotation_deg.tolist(), [30.5 + d for d in range(-50, 51)])
        self.assertEqual(saved.wd_from_deg.tolist(), [(66.0 + d) % 360 for d in range(-50, 51)])
        self.assertEqual(set(saved.direction), {"hourly"})     # run.json without the key: hourly
        self.cfg["study"]["direction"] = "fixed"               # these runs are not fixed-direction runs
        with self.assertRaises(StudyError):
            collect(self.cfg, "JG")


if __name__ == "__main__":
    unittest.main()
