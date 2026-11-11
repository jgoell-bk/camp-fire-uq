"""python -m unittest tests.test_baselines  (seconds; no ELMFIRE; needs `prepare` done)"""
import copy
import os
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from baselines.case import build_case, template_text
from baselines.common import Config, read_json
from baselines.forcing import from_uv, to_uv, vector_mean
from baselines.qoi import summarize
from baselines.run import run_one

REAL = Config()


class TestBaselines(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.cfg = copy.deepcopy(REAL)
        self.cfg.file, self.cfg.root, self.cfg.sha256 = REAL.file, REAL.root, REAL.sha256
        self.cfg["paths"] = dict(REAL["paths"], scratch_root=str(Path(self.tmp.name) / "root"))
        # read the real shared inputs, write runs into the temporary root
        Path(self.tmp.name, "root").mkdir()
        os.symlink(REAL.shared_dir, Path(self.tmp.name, "root", "shared"))

    def tearDown(self):
        self.tmp.cleanup()

    def test_uv_round_trip_and_mean(self):
        ws, wd = np.array([10.0, 5.0, 0.0]), np.array([350.0, 90.0, 0.0])
        s, d = from_uv(*to_uv(ws, wd))
        self.assertTrue(np.allclose(s, ws) and np.allclose(d[:2], wd[:2]))
        # north wind at 10 and east wind at 10 average to a 7.07 mph wind from 45 deg
        a = pd.DataFrame({"ws_mph": [10.0], "wd_deg": [0.0]})
        b = pd.DataFrame({"ws_mph": [10.0], "wd_deg": [90.0]})
        s, d = vector_mean([a, b])
        self.assertAlmostEqual(s[0], 10 / np.sqrt(2))
        self.assertAlmostEqual(d[0], 45.0)

    def test_prepared_weather(self):
        w = pd.read_csv(REAL.path(REAL["paths"]["prepared"]) / "weather.csv", index_col="time_utc")
        self.assertEqual(list(w.run_id.unique()), ["JG", "CO", "MEAN"])
        jg, co, mean = (w[w.run_id == r] for r in ("JG", "CO", "MEAN"))
        self.assertEqual(len(mean), REAL.nbands)
        s, d = vector_mean([jg, co])       # pandas' CSV parser can differ in the last bit
        self.assertTrue(np.allclose(s, mean.ws_mph, rtol=0, atol=1e-12))
        self.assertTrue(np.allclose(d, mean.wd_deg, rtol=0, atol=1e-12))
        for col in ("m1_pct", "m10_pct", "m100_pct"):      # one moisture series for every run
            self.assertTrue(np.array_equal(jg[col], co[col]) and np.array_equal(jg[col], mean[col]))

    def test_namelist_is_the_template(self):
        out = Path(self.tmp.name) / "build"
        build_case(self.cfg, "CO", out, out / "scratch")
        text = (out / "elmfire.data.in").read_text()
        lines, template = text.splitlines(), template_text(self.cfg).splitlines()
        self.assertEqual(len(lines), len(template))
        # every line without a {placeholder} is copied unchanged
        self.assertEqual([l for l, t in zip(lines, template) if "{" not in t], [t for t in template if "{" not in t])
        for line in ("SEED = 2023", "SIMULATION_TSTART = 900.0", "SIMULATION_TSTOP = 86400.0",
                     "WIND_FLUCTUATIONS = .FALSE.", "RANDOMIZE_RANDOM_SEED = .FALSE.", "P_EPS = 0.01"):
            self.assertIn(line, text)

    def test_run_one_reaches_validation(self):
        """run_one through build and the ELMFIRE call, with /bin/true as ELMFIRE. It
        must get as far as checking the outputs (and fail there: none)."""
        self.cfg["case"]["elmfire"] = "/bin/true"
        saved = os.environ.get("L_SCRATCH_JOB")
        os.environ["L_SCRATCH_JOB"] = str(Path(self.tmp.name) / "node")
        try:
            for rid in self.cfg.run_ids:
                self.assertEqual(run_one(self.cfg, rid), "failed")
                attempt = read_json(self.cfg.run_dir(rid) / "attempt.json")
                self.assertEqual(attempt["exit_code"], 0)
                self.assertIn("time_of_arrival", attempt["error"])
                self.assertTrue((self.cfg.run_dir(rid) / "run.json").exists())
        finally:
            if saved is None:
                os.environ.pop("L_SCRATCH_JOB")
            else:
                os.environ["L_SCRATCH_JOB"] = saved

    def test_summary(self):
        t = np.array([3600.0, 7200.0, np.nan, np.nan])
        s = summarize(t, 86400.0)
        self.assertEqual((s["cells"], s["burned_cells"]), (4, 2))
        self.assertAlmostEqual(s["burned_share"], 0.5)
        self.assertAlmostEqual(s["earliest_h"], 1.0)
        self.assertAlmostEqual(s["mean_burned_h"], 1.5)
        self.assertAlmostEqual(s["mean_capped_h"], (1 + 2 + 24 + 24) / 4)

    def test_cell_order_is_matlab_find(self):
        """np.nonzero(mask.T) lists cells column by column, as MATLAB's find(mask)."""
        mask = np.array([[0, 1], [1, 1]], bool)
        cols, rows = np.nonzero(mask.T)
        self.assertEqual(list(zip(rows + 1, cols + 1)), [(2, 1), (1, 2), (2, 2)])


if __name__ == "__main__":
    unittest.main()
