"""Tests for the consequential logic; none of them run ELMFIRE.

    cd <eki_calibration folder> && python -m unittest -v tests.test_eki

Tests that need prepared/ (baseline.csv, observations.csv), ensembles/I00.csv or the
shared manifest on scratch are skipped when those don't exist yet.
"""
import json
import os
import re
import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin

from eki.auto import auto, decide, latest_ensemble, update_gain
from eki.case import apply_overrides, build_case, coefficients, namelist_overrides
from eki.common import (BASELINE, IGNITION, NOT_MAPPED, NOT_REACHED, REACHED, RUN_FAILED, Config, StudyError,
                        parse_run_id, read_json, run_dir, run_group)
from eki.corrections import basis, band_hours, corrected_weather, curves, theta_names
from eki.ensemble import knot_correlation, normal_cdf, rank_report, sample
from eki.extract import classify, sample_nearby
from eki.forcing import check_reference, from_uv, load_baseline, to_uv, vector_mean
from eki.holdout import Union, load_split, split
from eki.landscape import eligible
from eki.nist import COLUMNS as NIST_COLUMNS
from eki.observations import load_table, map_point, match_nist, model_times, nearest_eligible, select_subset
from eki.run import run_one
from eki.pilot import distance_to_burned
from eki.target import candidates, deadline_residuals, fit_arrays, selection
from eki.objective import (block_weights, distance_outputs, mass_matrix, objective, regularization_matrix,
                           regularization_terms)
from eki.run import selection as selection_runs
from eki.update import controlled_step, eki_update, write_ensemble

REAL = Config()
REAL_SHARED = REAL.shared_dir / "manifest.json"
PREPARED = REAL.path(REAL["paths"]["prepared"])
HAVE_BASELINE = (PREPARED / "baseline.csv").exists()
HAVE_OBS = (PREPARED / "observations.csv").exists()
HAVE_I00 = (REAL.path(REAL["paths"]["ensembles"]) / "I00.csv").exists()


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.cfg = Config()
        self.cfg["paths"]["scratch_root"] = str(Path(self.tmp.name) / "scratch")
        self.cfg["paths"]["results"] = str(Path(self.tmp.name) / "results")

    def tearDown(self):
        self.tmp.cleanup()

    def theta(self, **values):
        t = {n: 0.0 for n in theta_names(self.cfg)}
        t.update(values)
        return t

    def flat_baseline(self, ws=10.0, wd=350.0, m=(5.0, 7.0, 10.0)):
        n = self.cfg.nbands
        return pd.DataFrame({"ws_mph": np.full(n, ws), "wd_deg": np.full(n, wd),
                             "m1_pct": np.full(n, m[0]), "m10_pct": np.full(n, m[1]), "m100_pct": np.full(n, m[2])})


class TestIds(Base):
    def test_run_ids(self):
        self.assertEqual(parse_run_id(BASELINE), (None, None))
        self.assertEqual(parse_run_id("I00_M07"), (0, 7))
        self.assertEqual(run_dir(self.cfg, "I00_M07"), self.cfg.scratch_root / "runs/I00/I00_M07")
        self.assertEqual(run_dir(self.cfg, BASELINE), self.cfg.scratch_root / "runs/BASELINE/BASELINE")
        for bad in ("I00_M00", "I0_M01", "I00_M1", "JG_0001", "baseline", "I00_M01x"):
            with self.assertRaises(StudyError):
                parse_run_id(bad)


class TestCorrections(Base):
    def test_basis(self):
        k = self.cfg["corrections"]["knots_h"]
        self.assertEqual(k, [0, 6, 12])
        B = basis(band_hours(self.cfg), k)
        self.assertEqual(B.shape, (25, 3))
        np.testing.assert_allclose(B.sum(axis=1), 1.0)
        np.testing.assert_allclose(B[[0, 6, 12]], np.eye(3))        # bands 1, 7, 13 are the knots
        np.testing.assert_allclose(B[3], [0.5, 0.5, 0.0])            # 09:00 PST, halfway
        np.testing.assert_allclose(B[9], [0.0, 0.5, 0.5])            # 15:00 PST
        np.testing.assert_allclose(B[13:], np.tile([0, 0, 1.0], (12, 1)))   # constant after 18:00
        with self.assertRaises(StudyError):
            basis([0.0], [0, 6, 6])

    def test_zero_correction_is_the_baseline(self):
        base = self.flat_baseline()
        base["ws_mph"] = np.linspace(2, 16, self.cfg.nbands)
        base["wd_deg"] = np.linspace(0, 359, self.cfg.nbands)
        s, _ = corrected_weather(self.cfg, base, self.theta())
        for v, col in (("ws", "ws_mph"), ("wd", "wd_deg"), ("m1", "m1_pct"), ("m10", "m10_pct"), ("m100", "m100_pct")):
            np.testing.assert_array_equal(s[v], base[col].to_numpy())

    def test_speed_direction_moisture(self):
        base = self.flat_baseline(ws=10.0, wd=350.0)
        th = self.theta(a_k1=np.log(2.0), a_k2=np.log(2.0), a_k3=np.log(0.5),
                        d_k1=30.0, d_k2=-360.0 - 20.0, d_k3=50.0, m_k1=1.0, m_k2=-1.0, m_k3=2.0)
        s, c = corrected_weather(self.cfg, base, th)
        self.assertAlmostEqual(s["ws"][0], 20.0)
        self.assertAlmostEqual(s["ws"][24], 5.0)
        self.assertAlmostEqual(s["wd"][0], 20.0)                    # 350 + 30 wraps to 20
        self.assertAlmostEqual(s["wd"][6], 330.0)                   # 350 - 380 wraps to 330
        self.assertAlmostEqual(c["d"][3], (30.0 - 380.0) / 2)       # interpolated unwrapped
        self.assertTrue(((s["wd"] >= 0) & (s["wd"] < 360)).all())
        for k, m0 in (("m1", 5.0), ("m10", 7.0), ("m100", 10.0)):    # the same offset on all classes
            np.testing.assert_allclose(s[k], m0 + c["m"])
        self.assertAlmostEqual(s["m1"][6], 4.0)
        self.assertAlmostEqual(s["m1"][24], 7.0)

    def test_inadmissible_moisture_is_refused(self):
        base = self.flat_baseline(m=(2.0, 7.0, 10.0))
        corrected_weather(self.cfg, base, self.theta(m_k2=-1.0))     # 1.0 % at band 7: allowed
        with self.assertRaises(StudyError):
            corrected_weather(self.cfg, base, self.theta(m_k2=-1.5))  # 0.5 %: refused, not clipped

    def test_theta_vector_and_dict_agree(self):
        names = theta_names(self.cfg)
        vec = np.arange(1, 10) / 10
        a, b = curves(self.cfg, vec), curves(self.cfg, dict(zip(names, vec)))
        for f in a:
            np.testing.assert_array_equal(a[f], b[f])


class TestBaseline(Base):
    def test_vector_mean(self):
        w = lambda ws, wd: pd.DataFrame({"ws_mph": [float(ws)], "wd_deg": [float(wd)]})
        ws, wd = vector_mean([w(10, 350), w(10, 10)])
        self.assertAlmostEqual(wd[0], 0.0, places=9)                 # not 180
        self.assertAlmostEqual(ws[0], 10 * np.cos(np.radians(10)), places=9)
        ws, _ = vector_mean([w(10, 90), w(10, 270)])
        self.assertAlmostEqual(ws[0], 0.0, places=9)                 # opposing winds cancel
        u, v = to_uv(np.array([5.0]), np.array([45.0]))
        s, d = from_uv(u, v)
        self.assertAlmostEqual(s[0], 5.0)
        self.assertAlmostEqual(d[0], 45.0)

    @unittest.skipUnless(HAVE_BASELINE, "prepare not run")
    def test_matches_reference_mean(self):
        diff = check_reference(self.cfg, load_baseline(self.cfg))
        self.assertLess(diff["max_abs_ws_diff_mph"], 1e-4)


class TestEnsemble(Base):
    def test_sample(self):
        t1, z1 = sample(self.cfg)
        t2, _ = sample(self.cfg)
        pd.testing.assert_frame_equal(t1, t2)                        # reproducible from the seed
        e = self.cfg["initial_ensemble"]
        self.assertEqual(len(t1), e["size"])
        self.assertEqual(list(t1.run_id[:2]), ["I00_M01", "I00_M02"])
        for r in (1, 2, 3):
            mult = t1[f"speed_mult_k{r}"]
            np.testing.assert_allclose(t1[f"a_k{r}"], np.log(mult))
            self.assertTrue(mult.between(*e["speed_multiplier"]).all())
            self.assertTrue(t1[f"d_k{r}"].between(*e["direction_deg"]).all())
            self.assertTrue(t1[f"m_k{r}"].between(*e["moisture_pp"]).all())
        self.assertTrue(rank_report(t1[theta_names(self.cfg)])["full_rank"])

    def test_correlation_and_marginals_in_a_large_sample(self):
        self.cfg["initial_ensemble"]["size"] = 20000
        t, z = sample(self.cfg)
        self.assertAlmostEqual(np.corrcoef(z.d_k1, z.d_k2)[0, 1], 0.7, delta=0.02)
        self.assertAlmostEqual(np.corrcoef(z.d_k1, z.d_k3)[0, 1], 0.49, delta=0.02)
        self.assertAlmostEqual(abs(np.corrcoef(z.d_k1, z.a_k1)[0, 1]), 0.0, delta=0.03)   # families independent
        u = (t.d_k2 + 60) / 120                                       # uniform on [0, 1]
        np.testing.assert_allclose(np.quantile(u, [0.1, 0.5, 0.9]), [0.1, 0.5, 0.9], atol=0.015)

    def test_helpers(self):
        self.assertAlmostEqual(float(normal_cdf(0.0)), 0.5)
        self.assertAlmostEqual(float(normal_cdf(1.959963985)), 0.975, places=6)
        np.testing.assert_allclose(knot_correlation(0.7, 3), [[1, .7, .49], [.7, 1, .7], [.49, .7, 1]])
        self.assertFalse(rank_report(np.tile([[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]], 1))["full_rank"])


class TestObservations(Base):
    def test_model_times(self):
        df = pd.DataFrame({"date_original": ["2018-11-08T00:00:00", "2018-11-09T00:00:00"],
                           "time_original": ["06:25:00", "05:45:00"],
                           "elapsed_original": ["0 days 00:10:00", "0 days 23:30:00"],
                           "elapsed_minutes_original": [10.0, 1410.0]})
        t = model_times(self.cfg, df)
        self.assertEqual(list(t.model_t_s), [1500.0, 85500.0])      # 06:25 PST = 14:25Z = t0 + 1500 s
        self.assertEqual(list(t.observed_time_utc), ["2018-11-08T14:25:00Z", "2018-11-09T13:45:00Z"])
        bad = df.assign(elapsed_minutes_original=[11.0, 1410.0])
        with self.assertRaises(StudyError):
            model_times(self.cfg, bad)
        self.cfg["simulation"]["ignition_s"] = 1800                  # origin no longer = ignition
        with self.assertRaises(StudyError):
            model_times(self.cfg, df)

    def test_mapping_rules(self):
        """A 5 x 5 grid of 30 m cells: nonburnable 256 at the centre, eligible cells at
        known places. Mapping depends on the landscape and the point only."""
        tr = from_origin(0.0, 150.0, 30.0, 30.0)
        fbfm = np.full((5, 5), 256)
        fbfm[2, 3] = 102          # east neighbour of the centre
        fbfm[0, 0] = 91           # urban corner
        fbfm[4, 4] = 98           # water
        bldg = np.zeros((5, 5))
        bldg[0, 0] = 100.0
        ok = eligible(fbfm, bldg)
        self.assertTrue(ok[2, 3] and ok[0, 0] and not ok[2, 2] and not ok[4, 4])
        review = {92, 93, 98, 99}
        m = map_point(75.0, 75.0, 2, 2, True, fbfm, ok, tr, 60.0, review)       # centre of cell (2, 2)
        self.assertEqual((m["status"], m["row"], m["col"], m["fuel"]), ("proxy", 2, 3, 102))
        self.assertAlmostEqual(m["displacement_m"], 30.0)
        m = map_point(15.0, 135.0, 0, 0, True, fbfm, ok, tr, 60.0, review)      # in the urban cell
        self.assertEqual((m["status"], m["row"], m["col"]), ("direct", 0, 0))
        m = map_point(135.0, 15.0, 4, 4, True, fbfm, ok, tr, 60.0, review)      # on water
        self.assertEqual(m["status"], "review")
        self.assertIsNone(m["row"])
        fbfm[2, 3] = 256
        ok = eligible(fbfm, bldg)
        m = map_point(75.0, 75.0, 2, 2, True, fbfm, ok, tr, 60.0, review)       # urban corner is 85 m away
        self.assertEqual(m["status"], "excluded")
        self.assertEqual(map_point(75.0, 75.0, 2, 2, True, fbfm, ok, tr, 90.0, review)["status"], "proxy")
        self.assertEqual(map_point(1e6, 1e6, -1, -1, False, fbfm, ok, tr, 60.0, review)["status"], "excluded")

    def test_nearest_uses_the_reported_point(self):
        tr = from_origin(0.0, 90.0, 30.0, 30.0)
        ok = np.zeros((3, 3), bool)
        ok[1, 0] = ok[1, 2] = True                                    # west and east of the centre
        r, c, d = nearest_eligible(50.0, 45.0, 1, 1, ok, tr, 60.0)    # point east of the centre's middle
        self.assertEqual((r, c), (1, 2))
        r, c, d = nearest_eligible(45.0, 45.0, 1, 1, ok, tr, 60.0)    # exact tie: lower column wins
        self.assertEqual((r, c), (1, 0))

    @unittest.skipUnless(HAVE_OBS, "observations not run")
    def test_real_table(self):
        obs = pd.read_csv(PREPARED / "observations.csv")
        self.assertEqual(len(obs), 2836)
        self.assertEqual(list(obs.obs_id[:2]), ["OBS0001", "OBS0002"])
        self.assertTrue(obs.obs_id.is_unique)
        self.assertEqual((obs.model_t_s - obs.elapsed_s).unique().tolist(), [900.0])
        direct, proxy = obs[obs.r60_status == "direct"], obs[obs.r60_status == "proxy"]
        self.assertTrue((direct.r60_displacement_m <= 15 * np.sqrt(2) + 1e-6).all())   # inside its own cell
        self.assertTrue((proxy.r60_displacement_m <= 60).all())
        self.assertTrue(proxy[~proxy.reviewed].fuel.eq(256).all())
        self.assertTrue(((direct.row == direct.r60_row) & (direct.col == direct.r60_col)).all())
        with rasterio.open(REAL.shared_dir / "fuels/fbfm40b.tif") as f, \
                rasterio.open(REAL.shared_dir / "fuels/baa_m.tif") as b:
            ok = eligible(f.read(1), b.read(1))
        mapped = obs[obs.r60_status.isin(["direct", "proxy"])]
        self.assertTrue(ok[mapped.r60_row.astype(int), mapped.r60_col.astype(int)].all())
        self.assertFalse((obs.r60_status == "review").any())          # every review point has a decision
        status = dict(zip(obs.obs_id, obs.r60_status))
        self.assertEqual([status[o] for o in ("OBS0073", "OBS2440", "OBS2482")], ["proxy", "proxy", "excluded"])
        self.assertAlmostEqual(float(obs.set_index("obs_id").r60_displacement_m["OBS0073"]), 29.8, places=1)
        both = obs[obs.r60_status == "proxy"]                         # 90 m keeps every 60 m proxy cell
        self.assertTrue(((both.r60_row == both.r90_row) & (both.r60_col == both.r90_col)).all())

    @unittest.skipUnless(HAVE_OBS, "observations not run")
    def test_nist_fields(self):
        obs = load_table(REAL)
        self.assertTrue((obs.nist_sheet != "").all() and obs.nist_sheet_row.min() >= 2)
        first = obs.iloc[0]
        self.assertEqual((first.nist_sheet, int(first.nist_sheet_row), first.nist_source_id), ("Concow", 2, "911-001-1"))
        self.assertIn("transmission lines", first.nist_description)
        self.assertEqual(int((obs.event_type == "residual").sum()), int(obs.nist_residual.sum()))
        np.testing.assert_allclose((obs.window_end_s - obs.window_start_s).dropna(),
                                   60 * obs.nist_window_min.dropna())
        rep = obs[obs.nist_record_repeat_of != ""]                     # a repeat has the same record
        src = obs.set_index("obs_id").loc[rep.nist_record_repeat_of]
        self.assertTrue((rep.nist_description.to_numpy() == src.nist_description.to_numpy()).all())
        self.assertTrue((rep.lon.to_numpy() == src.lon.to_numpy()).all())


class TestNistMatch(Base):
    def records(self):
        return pd.DataFrame({"placed": [True, True, False], "lon": [-121.1, -121.2, np.nan],
                             "lat": [39.1, 39.2, np.nan], "day": [8.0, 9.0, np.nan],
                             "minute_of_day": [385.0, 60.0, np.nan], "sheet": ["A", "A", "B"],
                             "sheet_row": [2, 3, 2], "fire_type": ["V", "S", "O"],
                             "window_min": [0.0, 15.0, np.nan], "residual": [False, True, False],
                             "ssi": [False, False, False],
                             **{c: ["a", "b", ""] for c in NIST_COLUMNS}})

    def geo(self):
        return pd.DataFrame({"lon": [-121.1, -121.2], "lat": [39.1, 39.2],
                             "date_original": ["2018-11-08T00:00:00", "2018-11-09T00:00:00"],
                             "time_original": ["06:25:00", "01:00:00"]})

    def test_match_in_order_and_refuse_otherwise(self):
        out, unplaced = match_nist(self.geo(), self.records())
        self.assertEqual(list(out.nist_sheet_row), [2, 3])
        self.assertEqual(list(out.nist_match_basis), ["unique key", "unique key"])
        self.assertEqual(len(unplaced), 1)
        bad = self.geo()
        bad.loc[1, "time_original"] = "01:01:00"                        # one minute off: no match at all
        with self.assertRaises(StudyError):
            match_nist(bad, self.records())


class TestHoldout(Base):
    @unittest.skipUnless(HAVE_OBS and (PREPARED / "holdout_v1.csv").exists(), "holdout not run")
    def test_real_split(self):
        obs, part = load_table(REAL), load_split(REAL)
        self.assertEqual(list(part.obs_id), list(obs.obs_id))
        self.assertEqual(set(part.partition), {"fit", "holdout"})
        self.assertTrue((part.groupby("block").partition.nunique() == 1).all())
        cells = part[part.location_basis == "sampling cell"]
        self.assertTrue((cells.groupby(["cell_row", "cell_col"]).partition.nunique() == 1).all())
        p = dict(zip(part.obs_id, part.partition))
        self.assertEqual(p["OBS2440"], p["OBS2434"])
        rep = obs[obs.nist_record_repeat_of != ""]
        self.assertTrue(all(p[a] == p[b] for a, b in zip(rep.obs_id, rep.nist_record_repeat_of)))
        frac = part[part.partition == "holdout"].block.nunique() / part.block.nunique()
        self.assertAlmostEqual(frac, 0.2, delta=0.01)
        again, _ = split(REAL, obs)                                     # reproducible from the seed
        self.assertEqual(list(again.partition), list(part.partition))

    def test_links_join_blocks(self):
        u = Union()
        u.join("B01_01", "B05_05")
        u.join("B05_05", "B09_09")
        self.assertEqual(u.find("B09_09"), "B01_01")
        self.assertEqual(u.find("B02_02"), "B02_02")


class TestTarget(Base):
    def reports(self):
        """Five reports in one cell and one in another; times in model seconds."""
        return pd.DataFrame({
            "obs_id": ["OBS0001", "OBS0002", "OBS0003", "OBS0004", "OBS0005", "OBS0006"],
            "feature_index": [0, 1, 2, 3, 4, 5],
            "model_t_s": [3000.0, 4000.0, 2000.0, 80000.0, 5000.0, 6000.0],
            "window_end_s": [9000.0, 4000.0, 2500.0, 90000.0, 4000.0, 6000.0],
            "in_window": [True, True, True, True, True, True],
            "event_type": ["fire_present", "fire_present", "residual", "fire_present", "fire_present", "fire_present"],
            "r60_status": ["direct"] * 5 + ["proxy"],
            "r60_row": [5, 5, 5, 5, 5, 9], "r60_col": [7, 7, 7, 7, 7, 2]})

    def test_smallest_deadline_not_earliest_start(self):
        s = selection(self.cfg, self.reports()).set_index("obs_id")
        # OBS0001 starts first (3000 s) but its window ends at 9000 s. OBS0002 starts later
        # and ends at 4000 s: the tighter bound wins. OBS0005 ties it but comes later in the file.
        self.assertEqual(s.selection["OBS0002"], "primary")
        self.assertEqual(s.selection_rule["OBS0001"], "later_deadline")
        self.assertEqual(s.selection_rule["OBS0005"], "same_deadline")
        self.assertEqual(s.selection_rule["OBS0003"], "residual")            # earliest of all, but residual
        self.assertEqual(s.selection_rule["OBS0004"], "window_after_run")    # ends after 86,400 s
        self.assertEqual(s.selection["OBS0006"], "primary")                  # alone in its cell
        self.assertEqual(s.deadline_s["OBS0002"], 4000.0)
        self.assertTrue((s.selection_reason != "").all())

    def test_allowance_and_overrides(self):
        self.cfg["target"]["deadline_allowance_s"] = 600
        s = selection(self.cfg, self.reports()).set_index("obs_id")
        self.assertEqual(s.deadline_s["OBS0002"], 4600.0)
        self.cfg["target"]["deadline_overrides"] = {"OBS0001": {"u_model_s": 3500.0, "reason": "test"}}
        s = selection(self.cfg, self.reports()).set_index("obs_id")
        self.assertEqual(s.selection["OBS0001"], "primary")                  # now the tightest bound
        self.assertTrue(s.deadline_basis["OBS0001"].startswith("override"))
        self.cfg["target"]["deadline_overrides"] = {"OBS0001": {"u_model_s": 9999.0, "reason": "looser"}}
        with self.assertRaises(StudyError):
            selection(self.cfg, self.reports())

    def test_deadline_residuals(self):
        u = np.array([3600.0, 3600.0, 3600.0, 3600.0])
        T = np.array([[1800.0, 3600.0, 7200.0, np.nan]])
        st = np.array([[REACHED, REACHED, REACHED, NOT_REACHED]], np.int8)
        r, censored = deadline_residuals(T, st, u, np.full(4, 3600.0))
        np.testing.assert_array_equal(r[0, :3], [0.0, 0.0, 1.0])              # early, on time: no penalty; 1 h late: 1
        self.assertTrue(np.isnan(r[0, 3]))                                     # non-arrival: no number,
        self.assertEqual(list(censored[0]), [False, False, False, True])       # flagged as censored
        with self.assertRaises(StudyError):                                    # and never into an update
            eki_update(np.random.default_rng(0).normal(size=(12, 9)), np.tile(r, (12, 1)), np.zeros(4), np.ones(4))
        with self.assertRaises(StudyError):
            deadline_residuals(T, np.array([[REACHED, IGNITION, REACHED, REACHED]], np.int8), u, np.full(4, 3600.0))

    @unittest.skipUnless(HAVE_OBS and (PREPARED / "holdout_v1.csv").exists(), "holdout not run")
    def test_real_constraints(self):
        obs, part = load_table(REAL), load_split(REAL).partition
        c = candidates(REAL, obs, part)
        self.assertTrue((part[c.index] == "fit").all())
        self.assertTrue(c.in_window.all() and (c.deadline_s <= REAL["simulation"]["stop_s"]).all())
        self.assertTrue((c.event_type == "fire_present").all())
        self.assertTrue(c.r60_status.isin(["direct", "proxy"]).all())
        self.assertFalse(c.duplicated(["r60_row", "r60_col"]).any())
        self.assertTrue(c.feature_index.is_monotonic_increasing)
        np.testing.assert_array_equal(c.deadline_s, c.window_end_s)           # no allowance, no overrides
        s = selection(REAL, obs)
        comp = obs[(s.selection != "excluded").to_numpy()].assign(u=s.deadline_s[(s.selection != "excluded").to_numpy()])
        best = comp.groupby(["r60_row", "r60_col"]).u.min()
        got = c.set_index(["r60_row", "r60_col"]).deadline_s
        self.assertTrue((best[got.index].to_numpy() == got.to_numpy()).all())  # smallest bound in every cell

    def fake(self, status):
        J, N = status.shape
        return {"run_id": np.array(["BASELINE"] + [f"I00_M{j:02d}" for j in range(1, J)]),
                "obs_id": np.array([f"OBS{i:04d}" for i in range(1, N + 1)]),
                "theta": np.arange(J * 9, dtype=float).reshape(J, 9),
                "status_primary": status, "arrival_s_primary": np.where(status == REACHED, 3600.0, np.nan)}

    def test_fit_arrays_gate(self):
        ids = ["OBS0001", "OBS0002"]
        st = np.zeros((4, 3), np.int8)
        st[0, :] = NOT_REACHED                                           # BASELINE is not a member
        members, theta, G = fit_arrays(self.fake(st), ids)
        self.assertEqual(list(members), ["I00_M01", "I00_M02", "I00_M03"])
        self.assertEqual(G.shape, (3, 2))
        self.assertTrue(np.isfinite(G).all())
        np.testing.assert_array_equal(theta[0], np.arange(9, 18))
        st2 = st.copy()
        st2[2, 1] = NOT_REACHED                                          # one member misses one candidate
        with self.assertRaisesRegex(StudyError, "not reached by every member"):
            fit_arrays(self.fake(st2), ids)
        st3 = st.copy()
        st3[3, :] = RUN_FAILED
        with self.assertRaisesRegex(StudyError, "failed"):
            fit_arrays(self.fake(st3), ids)
        members, _, _ = fit_arrays(self.fake(st2), ["OBS0001", "OBS0003"])   # its miss is not a candidate
        self.assertEqual(len(members), 3)


class TestPilot(unittest.TestCase):
    def test_distance_to_burned_by_deadline(self):
        toa = np.full((50, 50), -9999.0)
        toa[25, 25] = -1.0                       # ignition cell
        toa[25, 30] = 1000.0
        toa[25, 40] = 5000.0
        rows, cols = np.array([25, 25, 25]), np.array([40, 40, 30])
        d = distance_to_burned(toa, rows, cols, np.array([6000.0, 2000.0, 500.0]), 30.0)
        self.assertEqual(d[0], 0.0)              # its own cell burned by the deadline
        self.assertEqual(d[1], 10 * 30.0)        # nearest burned by 2000 s: (25, 30)
        self.assertEqual(d[2], 5 * 30.0)         # by 500 s only the ignition cell
        far = np.full((400, 400), -9999.0)
        far[0, 0] = -1.0
        self.assertAlmostEqual(distance_to_burned(far, np.array([399]), np.array([399]), np.array([1e5]), 30.0)[0],
                               np.hypot(399, 399) * 30.0)       # found beyond the first search squares


class TestUpdate(unittest.TestCase):
    def problem(self, J=32, seed=0):
        rng = np.random.default_rng(seed)
        A = rng.normal(size=(8, 9))
        R = np.full(8, 0.3 ** 2)
        y = A @ rng.normal(size=9) + rng.normal(0, 0.3, 8)
        B = np.diag(np.full(9, 0.5))
        return A, R, y, B, rng.normal(0, 2, (J, 9))

    def test_one_step_is_the_kalman_analysis_with_the_ensemble_as_prior(self):
        """Exactness of ONE step, not convergence. For linear G = A theta the new mean
        minimizes 1/2 |z - H th|^2/(alpha Gamma) + 1/2 |th - m|^2/C. The second term is a
        penalty from the CURRENT ensemble covariance C, so this reference is the Kalman
        analysis step, not the minimizer of the intended objective (see the next test)."""
        A, R, y, B, th = self.problem()
        for alpha in (1.0, 4.0):
            new = eki_update(th, th @ A.T, y, R, B, alpha)
            H = np.vstack([A, B])
            z = np.concatenate([y, np.zeros(9)])
            W = np.diag(1 / (alpha * np.concatenate([R, np.ones(9)])))
            m, C = th.mean(axis=0), np.cov(th.T)
            Ci = np.linalg.inv(C)
            expect = np.linalg.solve(H.T @ W @ H + Ci, H.T @ W @ z + Ci @ m)
            np.testing.assert_allclose(new.mean(axis=0), expect, atol=1e-9)

    def test_iterating_approaches_the_tikhonov_solution(self):
        """Convergence toward the intended objective 1/2 |y - A th|^2/R + 1/2 |B th|^2,
        whose minimizer has no ensemble-covariance term."""
        A, R, y, B, th = self.problem()
        star = np.linalg.solve(A.T @ (A / R[:, None]) + B.T @ B, A.T @ (y / R))
        errors = []
        for n in range(300):
            th = eki_update(th, th @ A.T, y, R, B)
            errors.append(np.linalg.norm(th.mean(axis=0) - star))
        self.assertLess(errors[-1], errors[0])
        self.assertLess(errors[-1], 0.01 * np.linalg.norm(star))

    def test_weighted_distance_iterates_approach_the_penalized_minimum(self):
        """Distance-like outputs g_i = max(0, c_i + A_i th) (km), unequal block weights,
        outputs sqrt(w) g / L with zero target and Gamma = I: the EKI mean approaches the
        grid minimum of Phi = 1/2 sum w (g/L)^2 + 1/2 |B th|^2 (an independent reference)."""
        rng = np.random.default_rng(1)
        A = rng.normal(size=(6, 2))
        c = np.array([1.0, 2.0, 0.7, -1.0, -1.0, -1.0])
        w = block_weights(["b1", "b1", "b1", "b2", "b3", "b3"])
        B, L = np.diag([1 / 1.5, 1 / 1.5]), 1.0
        g = lambda th: np.maximum(0.0, c + th @ A.T)
        phi = lambda th: 0.5 * (w * (g(th) / L) ** 2).sum(-1) + 0.5 * ((th @ B.T) ** 2).sum(-1)
        grid = np.linspace(-3, 3, 1201)
        P = np.stack(np.meshgrid(grid, grid), -1).reshape(-1, 2)
        star = P[np.argmin(phi(P))]
        th = rng.normal(0, 1.5, (32, 2))
        start = phi(th.mean(axis=0))
        for _ in range(600):
            G = np.sqrt(w) * g(th) / L
            th = eki_update(th, G, np.zeros(6), np.ones(6), B)
        m = th.mean(axis=0)
        self.assertLess(np.linalg.norm(m - star), 0.05)
        self.assertLess(phi(m) - phi(star), 0.01 * (start - phi(star)))

    def test_basic_update_without_B(self):
        A, R, y, _, th = self.problem()
        new = eki_update(th, th @ A.T, y, R)
        g = th @ A.T
        dt, dg = th - th.mean(0), g - g.mean(0)
        K = (dt.T @ dg / 31) @ np.linalg.inv(dg.T @ dg / 31 + np.diag(R))
        np.testing.assert_allclose(new, th + (K @ (y - g).T).T, atol=1e-9)

    def test_refuses_nan_and_too_few_members(self):
        A, R, y, B, th = self.problem()
        G = th @ A.T
        G[3, 2] = np.nan
        with self.assertRaises(StudyError):
            eki_update(th, G, y, R, B)
        with self.assertRaises(StudyError):
            eki_update(th[:9], th[:9] @ A.T, y, R, B)

    def test_controlled_step_limits(self):
        cfg = Config()
        A, R, y, B, th = self.problem()
        G = (th @ A.T - y) / np.sqrt(R)                                 # zero target, Gamma = I
        new, rec = controlled_step(cfg, th, G, B, lambda row: True)
        sd = th.std(axis=0, ddof=1)
        self.assertLessEqual(np.max(np.abs(new.mean(0) - th.mean(0)) / sd), 0.5 + 1e-12)
        self.assertLessEqual(np.max(np.abs(new - th) / sd), 1.0 + 1e-12)
        self.assertEqual(rec["boundary_projections"], 0)
        if rec["alpha"] > 1:                                            # the alpha before it broke a limit
            prev = rec["tried"][-2]
            self.assertTrue(prev["max_mean_move_sd"] > 0.5 or prev["max_member_move_sd"] > 1.0)

    def test_controlled_step_collapsed_and_admissibility(self):
        cfg = Config()
        A, R, y, B, th = self.problem()
        th[:, 4] = 0.25                                                 # d_k2 has no spread
        G = (th @ A.T - y) / np.sqrt(R)
        new, rec = controlled_step(cfg, th, G, B, lambda row: True)
        self.assertEqual(rec["collapsed_coefficients"], ["d_k2"])
        np.testing.assert_allclose(new[:, 4], 0.25, atol=1e-12)          # outside the ensemble span
        # Admissibility is enforced by damping: if a member is inadmissible at every alpha
        # up to the one the limits alone would choose, the next alpha (twice it) is taken.
        n_tries, calls = len(rec["tried"]), []

        def admissible(row):
            calls.append(1)
            return not (len(calls) <= n_tries * len(th) and len(calls) % len(th) == 3)
        _, rec2 = controlled_step(cfg, th, G, B, admissible)
        self.assertEqual(rec2["alpha"], 2 * rec["alpha"])
        self.assertEqual(rec2["tried"][n_tries - 1]["inadmissible_members"], [2])
        self.assertEqual(rec2["tried"][-1]["inadmissible_members"], [])

    @unittest.skipUnless(HAVE_BASELINE, "prepare not run")
    def test_write_ensemble(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = Config()
            path = Path(tmp) / "I01.csv"
            th = np.zeros((32, 9))
            write_ensemble(cfg, path, 1, th, {"from": "test"})
            t = pd.read_csv(path)
            self.assertEqual(list(t.run_id[:2]), ["I01_M01", "I01_M02"])
            np.testing.assert_allclose(t.speed_mult_k1, 1.0)
            with self.assertRaises(StudyError):
                write_ensemble(cfg, path, 1, th, {})                     # never edited
            bad = th.copy()
            bad[5, 7] = -20.0                                            # moisture far below 1 %
            with self.assertRaises(StudyError):
                write_ensemble(cfg, Path(tmp) / "I02.csv", 2, bad, {})


class TestObjective(Base):
    def test_block_weights(self):
        blocks = np.array(["a"] * 5 + ["b"] * 1 + ["c"] * 2)
        w = block_weights(blocks)
        self.assertAlmostEqual(w.sum(), 1.0)
        for b in "abc":
            self.assertAlmostEqual(w[blocks == b].sum(), 1 / 3)          # equal total per block
        np.testing.assert_allclose(w[blocks == "a"], 1 / 15)

    def test_regularization_reproduces_its_definition(self):
        """1/2 |B th|^2 = R_amp + R_smooth, with R_amp the time average over the run of
        q(t)^2 by fine quadrature and R_smooth from the knot differences: independent."""
        B = regularization_matrix(self.cfg)
        self.assertEqual(B.shape, (9 + 6, 9))
        th = np.random.default_rng(4).normal(0, 1, (20, 9)) * np.repeat([0.7, 30.0, 1.0], 3)
        amp, smooth = regularization_terms(self.cfg, th)
        n_amp = 9
        np.testing.assert_allclose(0.5 * ((th @ B[:n_amp].T) ** 2).sum(1), amp, rtol=1e-6)
        np.testing.assert_allclose(0.5 * ((th @ B[n_amp:].T) ** 2).sum(1), smooth, rtol=1e-12)

    def test_reference_values(self):
        s = self.cfg["fit"]["penalty_scales"]
        one = np.zeros(9)
        one[3:6] = s["d"]                                               # direction one scale at every knot
        amp, smooth = regularization_terms(self.cfg, one)
        self.assertAlmostEqual(amp[0], 0.5)                             # time average of 1 is 1
        self.assertAlmostEqual(smooth[0], 0.0)
        step = np.zeros(9)
        step[4:6] = s["d"]                                              # d rises by one scale from 06:00 to 12:00 PST
        _, smooth = regularization_terms(self.cfg, step)
        self.assertAlmostEqual(smooth[0], 0.5)                          # one scale over 6 h -> 1/2
        M = mass_matrix([0, 6, 12], 24)
        np.testing.assert_allclose(M * 24, [[2, 1, 0], [1, 4, 1], [0, 1, 14]])
        self.assertAlmostEqual(np.ones(3) @ M @ np.ones(3), 1.0)

    def test_objective_parts_and_eki_outputs(self):
        rng = np.random.default_rng(5)
        th = rng.normal(0, 1, (4, 9)) * np.repeat([0.7, 30.0, 1.0], 3)
        d = rng.uniform(0, 5000, (4, 7))
        d[:, 0] = 0.0                                                   # a met deadline contributes nothing
        w = block_weights(["a", "a", "b", "b", "b", "c", "c"])
        o = objective(self.cfg, th, d, w)
        np.testing.assert_allclose(o["total"], o["distance"] + o["amplitude"] + o["smoothness"])
        F = np.hstack([distance_outputs(d, w, 1000.0), th @ regularization_matrix(self.cfg).T])
        np.testing.assert_allclose(0.5 * (F ** 2).sum(1), o["total"])   # EKI outputs reproduce Phi exactly
        np.testing.assert_allclose(o["distance"], 0.5 * (w * (d / 1000.0) ** 2).sum(1))
        with self.assertRaises(StudyError):
            distance_outputs(np.where(d > 4000, np.nan, d), w, 1000.0)  # distances are always finite


class TestCausality(Base):
    def test_a_knot_never_changes_weather_before_its_basis_turns_on(self):
        base = self.flat_baseline()
        hours = band_hours(self.cfg)
        B = basis(hours, self.cfg["corrections"]["knots_h"])
        names = theta_names(self.cfg)
        s0, _ = corrected_weather(self.cfg, base, self.theta())
        for i, n in enumerate(names):
            r = i % 3
            s1, _ = corrected_weather(self.cfg, base, self.theta(**{n: {"a": 0.3, "d": 20.0, "m": 0.5}[n[0]]}))
            off = hours < hours[np.argmax(B[:, r] > 0)]                 # before the basis first turns on
            for v in ("ws", "wd", "m1"):
                np.testing.assert_array_equal(s1[v][off], s0[v][off])
            changed = np.any([s1[v] != s0[v] for v in ("ws", "wd", "m1")], axis=0)
            self.assertTrue((B[changed, r] > 0).all())                   # only where its basis is nonzero
        self.assertEqual([float(hours[np.argmax(B[:, r] > 0)]) for r in range(3)], [0.0, 1.0, 7.0])


class TestCheckRuns(Base):
    def test_check_ids_and_coefficients(self):
        self.assertEqual(run_group("CHK01_I00MEAN"), "CHK01")
        self.assertEqual(run_group("CHK01_I00MEAN_REP"), "CHK01")
        self.assertEqual(run_group("I01_M05"), "I01")
        self.assertEqual(run_group(BASELINE), BASELINE)
        self.assertEqual(run_dir(self.cfg, "CHK01_I01MEAN"), self.cfg.scratch_root / "runs/CHK01/CHK01_I01MEAN")
        with tempfile.TemporaryDirectory() as tmp:
            self.cfg["paths"]["ensembles"] = tmp
            names = theta_names(self.cfg)
            pd.DataFrame({"run_id": ["CHK01_I00MEAN", "CHK01_I00MEAN_REP"], **{n: [0.1, 0.1] for n in names}}) \
                .to_csv(Path(tmp) / "CHK01.csv", index=False)
            th, src = coefficients(self.cfg, "CHK01_I00MEAN_REP")
            self.assertEqual(th["d_k2"], 0.1)                           # never the baseline's zeros
            self.assertEqual(selection_runs(self.cfg, "CHK01")[0], ["CHK01_I00MEAN", "CHK01_I00MEAN_REP"])
            with self.assertRaises(StudyError):
                coefficients(self.cfg, "CHK02_I00MEAN")
        for bad in ("CHK1_I00MEAN", "CHK01_MEAN", "CHK01_I00MEAN_REP2"):
            with self.assertRaises(StudyError):
                run_group(bad)


class TestExtract(unittest.TestCase):
    def test_classify(self):
        s, a = classify([3600.0, -9999.0, -1.0])
        self.assertEqual(list(s), [REACHED, NOT_REACHED, IGNITION])
        self.assertEqual(a[0], 3600.0)
        self.assertTrue(np.isnan(a[1:]).all())
        with self.assertRaises(StudyError):
            classify([0.0])

    def test_frozen_cell_stays_unreached(self):
        """The assigned cell is unburned but its neighbour burned: the primary rule keeps
        non-arrival; only the separate 3x3 sensitivity sees the neighbour."""
        toa = np.full((3, 3), -9999.0)
        toa[1, 2] = 5000.0
        s, a = classify(toa[[1], [1]])
        self.assertEqual(s[0], NOT_REACHED)
        self.assertTrue(np.isnan(a[0]))
        blocks = [(np.array([0, 1, 1]), np.array([0, 1, 2])), (np.array([], int), np.array([], int)),
                  (np.array([0]), np.array([0]))]
        s, a = sample_nearby(toa, blocks)
        self.assertEqual(list(s), [REACHED, NOT_MAPPED, NOT_REACHED])
        self.assertEqual(a[0], 5000.0)


@unittest.skipUnless(REAL_SHARED.exists() and HAVE_BASELINE and HAVE_I00, "prepare/ensemble not run")
class TestBuild(Base):
    def setUp(self):
        super().setUp()
        self.cfg.shared_dir.mkdir(parents=True)
        shutil.copy(REAL_SHARED, self.cfg.shared_dir / "manifest.json")

    def test_baseline_build_writes_the_baseline(self):
        out = Path(self.tmp.name) / "b"
        rec = build_case(self.cfg, BASELINE, out, out / "scratch")
        self.assertEqual(rec["theta"], [0.0] * 9)
        self.assertEqual((rec["P_EPS"], rec["SIMULATION_TSTART"], rec["SIMULATION_TSTOP"]), ("0.01", "900.0", "86400.0"))
        base = load_baseline(self.cfg)
        with rasterio.open(out / "inputs/ws.tif") as f:
            ws = f.read()
        self.assertEqual(ws.shape[0], 25)
        np.testing.assert_array_equal(ws[:, 0, 0], base.ws_mph.to_numpy(np.float32))
        text = (out / "elmfire.data.in").read_text()
        self.assertNotIn("{", text)
        self.assertRegex(text, r"(?m)^SIMULATION_TSTART = 900\.0$")
        self.assertRegex(text, r"(?m)^P_EPS = 0\.01$")
        self.assertRegex(text, r"(?m)^WIND_FLUCTUATIONS = \.FALSE\.$")

    def test_member_build_writes_its_correction(self):
        out = Path(self.tmp.name) / "m"
        rec = build_case(self.cfg, "I00_M05", out, out / "scratch")
        theta, _ = coefficients(self.cfg, "I00_M05")
        series, _ = corrected_weather(self.cfg, load_baseline(self.cfg), theta)
        for v in ("ws", "wd", "m1", "m100"):
            with rasterio.open(out / f"inputs/{v}.tif") as f:
                np.testing.assert_array_equal(f.read()[:, 3, 7], series[v].astype(np.float32))
        w = read_json(out / "inputs/weather.json")
        self.assertAlmostEqual(w["bands"][6]["correction"]["speed_mult"], rec["speed_mult_knots"][1])

    def test_run_one_reaches_validation(self):
        """run_one through build, log line and the ELMFIRE call, with /bin/true as
        ELMFIRE. It must get as far as checking the outputs (and fail there: none)."""
        self.cfg["case"]["elmfire"] = "/bin/true"
        saved = os.environ.get("L_SCRATCH_JOB")
        os.environ["L_SCRATCH_JOB"] = str(Path(self.tmp.name) / "node")
        try:
            chk = REAL.path(REAL["paths"]["ensembles"]) / "CHK01.csv"
            for rid in (BASELINE, "I00_M01") + (("CHK01_I01MEAN",) if chk.exists() else ()):
                self.assertEqual(run_one(self.cfg, rid), "failed")
                attempt = read_json(run_dir(self.cfg, rid) / "attempt.json")
                self.assertEqual(attempt["exit_code"], 0)
                self.assertIn("time_of_arrival", attempt["error"])
                self.assertTrue((run_dir(self.cfg, rid) / "run.json").exists())
        finally:
            if saved is None:
                os.environ.pop("L_SCRATCH_JOB")
            else:
                os.environ["L_SCRATCH_JOB"] = saved

    def test_auto_stops_at_failed_runs(self):
        """auto with /bin/true as ELMFIRE: the newest ensemble's runs all fail
        validation, so it must stop there (exit 1), log why, and promote nothing."""
        self.cfg["case"]["elmfire"] = "/bin/true"
        self.cfg["paths"]["batches"] = str(Path(self.tmp.name) / "batches")
        ens = sorted(p.name for p in REAL.path(REAL["paths"]["ensembles"]).iterdir())
        saved = os.environ.get("L_SCRATCH_JOB")
        os.environ["L_SCRATCH_JOB"] = str(Path(self.tmp.name) / "node")
        try:
            self.assertEqual(auto(self.cfg, workers=4), 1)
        finally:
            if saved is None:
                os.environ.pop("L_SCRATCH_JOB")
            else:
                os.environ["L_SCRATCH_JOB"] = saved
        log = read_json(self.cfg.dir("results") / "auto_log.json")
        n = latest_ensemble(self.cfg)
        self.assertIn(f"runs of I{n:02d} failed", log["jobs"][-1]["stopped"])
        self.assertEqual(log["jobs"][-1]["steps"], [])
        self.assertEqual(sorted(p.name for p in REAL.path(REAL["paths"]["ensembles"]).iterdir()), ens)


class TestAuto(Base):
    def test_decide(self):
        self.assertIn("reached I10", decide(10, 10, 0.5, 0.02))
        self.assertIn("1.0%", decide(4, 10, 0.01, 0.02))
        self.assertIsNone(decide(4, 10, 0.05, 0.02))
        self.assertIsNone(decide(0, 10, None, 0.02))          # I00 has no update behind it

    def test_update_gain_matches_verify(self):
        path = REAL.dir("results") / "CHK01_verify.json"
        if not path.exists():
            self.skipTest("no CHK01_verify.json")
        v = read_json(path)
        cur = v["objective"]["CHK01_I00MEAN"]["total"]
        prop = v["objective"]["CHK01_I01MEAN"]["total"]
        self.assertAlmostEqual(update_gain(REAL, 1), (cur - prop) / cur)
        self.assertIsNone(update_gain(REAL, 0))



class TestV2(Base):
    def subset_file(self, points):
        path = Path(self.tmp.name) / "subset.geojson"
        feats = [{"type": "Feature", "geometry": {"type": "Point", "coordinates": [lon, lat]},
                  "properties": {"Date": d, "Time": t, "elapsed": m / 60, "elapsed_hours": m / 60,
                                 "elapsed_minutes": m}} for lon, lat, d, t, m in points]
        path.write_text(json.dumps({"type": "FeatureCollection", "features": feats,
                                    "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"}}}))
        return path

    def test_select_subset(self):
        day = "2018-11-08T00:00:00"
        df = pd.DataFrame({"obs_id": ["OBS0001", "OBS0002", "OBS0003", "OBS0004"],
                           "lon": [-121.5, -121.5, -121.5, -121.6], "lat": [39.7] * 4,
                           "date_original": [day] * 4, "time_original": ["06:30:00"] * 3 + ["07:00:00"],
                           "elapsed_minutes_original": [15.0] * 3 + [45.0]})
        # keeps 2 of the 3 identical points (listed in another order) and the fourth
        path = self.subset_file([(-121.6, 39.7, day, "07:00:00", 45.0), (-121.5, 39.7, day, "06:30:00", 15.0),
                                 (-121.5, 39.7, day, "06:30:00", 15.0)])
        keep, info = select_subset(df, path)
        self.assertEqual(keep, [0, 1, 3])                          # the first 2 in workbook order
        self.assertEqual(info["partly_kept_groups"][0]["dropped_obs_ids"], ["OBS0003"])
        self.assertEqual(info["dropped_obs_ids"], ["OBS0003"])
        for bad in ([(-121.5001, 39.7, day, "06:30:00", 15.0)],     # moved
                    [(-121.5, 39.7, day, "06:31:00", 16.0)],        # retimed
                    [(-121.5, 39.7, day, "06:30:00", 20.0)],        # elapsed disagrees
                    [(-121.6, 39.7, day, "07:00:00", 45.0)] * 2):   # more copies than the original
            with self.assertRaises(StudyError):
                select_subset(df, self.subset_file(bad))

    def test_v2_holdout_keeps_v1_partitions(self):
        v2 = REAL.root / "v2" / "config.yaml"
        if not v2.exists() or not (REAL.root / "v2/prepared/holdout_v2.csv").exists():
            self.skipTest("v2 not prepared")
        cfg2 = Config(v2)
        h1, h2 = load_split(REAL), load_split(cfg2)
        m = h2.merge(h1[["obs_id", "block", "partition"]], on="obs_id", suffixes=("", "_v1"))
        self.assertEqual(len(m), len(h2))
        self.assertTrue((m.partition == m.partition_v1).all() and (m.block == m.block_v1).all())
        obs2 = load_table(cfg2)
        self.assertEqual(len(obs2), read_json(cfg2.path(cfg2["paths"]["prepared"]) / "observations.json")["subset"]["kept"])
        self.assertTrue(set(obs2.obs_id) <= set(load_table(REAL).obs_id))


class TestSensitivity(Base):
    def test_sens_ids_and_overrides(self):
        self.assertEqual(run_group("SEN01_I08MEAN_PEPS0.001"), "SEN01")
        self.assertEqual(parse_run_id("SEN01_I08MEAN_PEPS0.001"), (None, None))
        text = "A = 1\nP_EPS = 0.01 ! cut\nB = 2\n"
        self.assertEqual(apply_overrides(text, {"P_EPS": 0.001}), "A = 1\nP_EPS = 0.001 ! cut\nB = 2\n")
        self.assertEqual(apply_overrides(text, {}), text)
        with self.assertRaises(StudyError):
            apply_overrides(text + "P_EPS = 0.02\n", {"P_EPS": 0.001})      # set twice

    def test_only_sens_files_may_override(self):
        ens = Path(self.tmp.name) / "ensembles"
        ens.mkdir()
        self.cfg["paths"]["ensembles"] = str(ens)
        (ens / "CHK99.csv").write_text("run_id,P_EPS\nCHK99_I00MEAN,0.001\n")
        with self.assertRaises(StudyError):
            namelist_overrides(self.cfg, "CHK99_I00MEAN")
        (ens / "SEN99.csv").write_text("run_id,P_EPS\nSEN99_X,0.001\n")
        self.assertEqual(namelist_overrides(self.cfg, "SEN99_X"), {"P_EPS": 0.001})
        self.assertEqual(namelist_overrides(self.cfg, BASELINE), {})

    def test_score_fit_matches_verify(self):
        path = REAL.dir("results") / "CHK08_verify.json"
        if not path.exists():
            self.skipTest("no CHK08_verify.json")
        from eki.score import one, partition_constraints
        obs, split = load_table(REAL), load_split(REAL)
        c, _, w = partition_constraints(REAL, obs, split, "fit")
        s = one(REAL, "CHK08_I08MEAN", c, w)[0]
        self.assertAlmostEqual(s["distance_term"], read_json(path)["objective"]["CHK08_I08MEAN"]["distance"], places=10)

class TestV3(Base):
    """P_EPS as a coefficient and the over-spread term: present only when configured."""
    def setUp(self):
        super().setUp()
        path = REAL.root / "v3" / "config.yaml"
        if not path.exists():
            self.skipTest("no v3/config.yaml")
        self.cfg3 = Config(path)
        self.cfg3["paths"]["scratch_root"] = str(Path(self.tmp.name) / "scratch3")

    def test_v1_layout_unchanged(self):
        self.assertEqual(len(theta_names(self.cfg)), 9)
        self.assertEqual(regularization_matrix(self.cfg).shape, (15, 9))
        o = objective(self.cfg, np.zeros((1, 9)), np.zeros((1, 3)), np.full(3, 1 / 3))
        self.assertNotIn("overspread", o)
        with self.assertRaises(StudyError):
            objective(self.cfg, np.zeros((1, 9)), np.zeros((1, 3)), np.full(3, 1 / 3), [1.0])
        from eki.case import coefficient_overrides
        self.assertEqual(coefficient_overrides(self.cfg, self.theta()), {})

    def test_v3_regularization(self):
        names = theta_names(self.cfg3)
        self.assertEqual((len(names), names[-1]), (10, "p"))
        B = regularization_matrix(self.cfg3)
        self.assertEqual(B.shape, (9 + 1 + 6, 10))
        th = np.random.default_rng(7).normal(0, 1, (20, 10)) * np.r_[np.repeat([0.7, 30.0, 1.0], 3), 2.0]
        amp, smooth = regularization_terms(self.cfg3, th)
        np.testing.assert_allclose(0.5 * ((th @ B[:10].T) ** 2).sum(1), amp, rtol=1e-6)
        np.testing.assert_allclose(0.5 * ((th @ B[10:].T) ** 2).sum(1), smooth, rtol=1e-12)
        one = np.zeros(10)
        one[9] = np.log(10.0)                                           # a factor of 10 in P_EPS
        amp, smooth = regularization_terms(self.cfg3, one)
        self.assertAlmostEqual(amp[0], 0.5)
        self.assertAlmostEqual(smooth[0], 0.0)
        # the 9 weather rows are v1's
        np.testing.assert_array_equal(B[:9, :9], regularization_matrix(self.cfg)[:9])
        np.testing.assert_array_equal(B[10:, :9], regularization_matrix(self.cfg)[9:])

    def test_v3_overspread_objective(self):
        rng = np.random.default_rng(8)
        th = rng.normal(0, 0.5, (3, 10))
        d = rng.uniform(0, 3000, (3, 5))
        w = block_weights(["a", "a", "b", "c", "c"])
        acres = np.array([0.0, 10000.0, 25000.0])
        o = objective(self.cfg3, th, d, w, acres)
        np.testing.assert_allclose(o["overspread"], [0.0, 0.5, 3.125])
        np.testing.assert_allclose(o["total"], o["distance"] + o["amplitude"] + o["smoothness"] + o["overspread"])
        from eki.objective import overspread_outputs
        F = np.hstack([distance_outputs(d, w, 1000.0), th @ regularization_matrix(self.cfg3).T,
                       overspread_outputs(self.cfg3, acres)])
        np.testing.assert_allclose(0.5 * (F ** 2).sum(1), o["total"])   # the EKI outputs reproduce Phi
        with self.assertRaises(StudyError):
            objective(self.cfg3, th, d, w)                              # areas are required in v3

    def test_burned_outside_acres(self):
        from eki.overspread import burned_outside_acres
        toa = np.array([[-1.0, 100.0, 0.0], [50.0, 90000.0, 86400.0]])
        inside = np.array([[True, False, False], [False, True, True]])
        # outside: (0,1) burned at 100 s, (0,2) never, (1,0) at 50 s -> 2 cells
        self.assertEqual(burned_outside_acres(toa, inside, 0.5, 86400.0), 1.0)
        self.assertEqual(burned_outside_acres(toa, np.zeros_like(inside), 1.0, 86400.0), 4.0)  # 90000 s is after at_s
        self.assertEqual(burned_outside_acres(toa, np.zeros_like(inside), 1.0, 60.0), 2.0)     # -1 and 50 s

    def test_p_eps_floor_and_namelist(self):
        from eki.case import coefficient_overrides
        from eki.corrections import p_eps
        base = self.flat_baseline()
        t = {n: 0.0 for n in theta_names(self.cfg3)}
        self.assertEqual(p_eps(self.cfg3, t), 0.01)
        self.assertEqual(coefficient_overrides(self.cfg3, t), {"P_EPS": "0.01"})   # BASELINE: v1's literal
        corrected_weather(self.cfg3, base, {**t, "p": np.log(0.0006 / 0.01)})
        with self.assertRaises(StudyError):
            corrected_weather(self.cfg3, base, {**t, "p": np.log(0.0004 / 0.01)})
        lit = coefficient_overrides(self.cfg3, {**t, "p": 0.3})["P_EPS"]
        self.assertEqual(float(lit), 0.01 * np.exp(0.3))                           # exact
        text = "A = 1\nP_EPS = 0.01 ! cut\n"
        self.assertEqual(apply_overrides(text, {"P_EPS": lit}), f"A = 1\nP_EPS = {lit} ! cut\n")

    def test_v3_I00_copies_v1_weather(self):
        import csv
        import io
        p3 = self.cfg3.path(self.cfg3["paths"]["ensembles"]) / "I00.csv"
        if not p3.exists():
            self.skipTest("no v3 I00")
        a = list(csv.DictReader(io.StringIO((REAL.path(REAL["paths"]["ensembles"]) / "I00.csv").read_text())))
        b = list(csv.DictReader(io.StringIO(p3.read_text())))
        self.assertEqual(len(a), len(b))
        for ra, rb in zip(a, b):
            self.assertEqual(ra, {k: rb[k] for k in ra})                # every v1 field, as text
            self.assertAlmostEqual(float(rb["p_eps_value"]) / (0.01 * np.exp(float(rb["p"]))), 1.0, places=14)
        pe = np.array([float(r["p_eps_value"]) for r in b])
        self.assertTrue((pe >= 0.001).all() and (pe <= 0.05).all())


if __name__ == "__main__":
    unittest.main()
