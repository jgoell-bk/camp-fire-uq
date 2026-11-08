"""Write a versioned Sobol run matrix (Saltelli layout with second-order blocks).

Construction, identical to the one behind Camp_Fire_Sobol_N128_run_matrix_4608.csv
(--check-original regenerates that file's values exactly):
  u = scipy.stats.qmc.Sobol(d=2D, scramble=True, seed=SEED).random_base2(m=M)
  A = columns 0..D-1, B = columns D..2D-1, each mapped linearly to its bounds.
  Each base row expands to A, AB_k (A with factor k from B), BA_k (B with factor k
  from A), B; stations get identical vectors, in station-major order.

SciPy is not in the elmfire-vnv venv; use Sherlock's module:
  ml math python/3.12.1 py-scipy/1.16.0_py312
  PY=/share/software/user/open/python/3.12.1/bin/python3
  $PY -P tools/make_matrix.py --check-original manifest/Camp_Fire_Sobol_N128_run_matrix_4608.csv
  $PY -P tools/make_matrix.py manifest/Camp_Fire_Sobol_N128_D6_run_matrix_5376.csv
"""
import argparse
import csv

import numpy as np
from scipy.stats import qmc

SEED, M = 20261019, 7                    # N = 2**M = 128 base rows
FACTORS = [("a_U", 0.8, 1.2), ("b_theta_deg", -20.0, 20.0), ("b_1_pp", -1.0, 1.0),
           ("b_10_pp", -1.0, 1.0), ("b_100_pp", -2.0, 2.0),
           ("q_spot", 0.0, 1.0)]         # P_EPS = 0.001 * 50**q_spot (log-uniform 0.001-0.05)
STATIONS = [("JG", "JARBO GAP", "041214"), ("OP", "OPENSHAW", "041215"), ("CO", "COHASSET", "041211")]
MOISTURE = ("JARBO GAP", "041214")
BASE_COLUMNS = ["run_id", "wind_station", "wind_station_id", "moisture_station",
                "moisture_station_id", "design_row", "base_row", "sobol_block", "sampling_seed"]


def p_eps(q):
    return 0.001 * 50 ** q


def saltelli_rows(factors):
    """[(base_row, block, values)] for one station, in design_row order."""
    lo = np.array([f[1] for f in factors])
    hi = np.array([f[2] for f in factors])
    d = len(factors)
    u = qmc.Sobol(d=2 * d, scramble=True, seed=SEED).random_base2(m=M)
    a, b = lo + u[:, :d] * (hi - lo), lo + u[:, d:] * (hi - lo)
    names = [f[0] for f in factors]
    rows = []
    for i in range(len(u)):
        rows.append((i + 1, "A", a[i]))
        for k in range(d):
            r = a[i].copy(); r[k] = b[i, k]
            rows.append((i + 1, f"AB_{names[k]}", r))
        for k in range(d):
            r = b[i].copy(); r[k] = a[i, k]
            rows.append((i + 1, f"BA_{names[k]}", r))
        rows.append((i + 1, "B", b[i]))
    return rows


def write(path):
    names = [f[0] for f in FACTORS]
    rows = saltelli_rows(FACTORS)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(BASE_COLUMNS + names + ["P_EPS"])
        for prefix, station, sid in STATIONS:
            for n, (base, block, v) in enumerate(rows, start=1):
                w.writerow([f"{prefix}_{n:04d}", station, sid, *MOISTURE, n, base, block, SEED,
                            *(f"{x:.17g}" for x in v), f"{p_eps(v[names.index('q_spot')]):.17g}"])
    print(f"wrote {path}: {len(rows) * len(STATIONS)} rows, {len(rows)} per station")


def check_original(path):
    rows = saltelli_rows(FACTORS[:5])
    with open(path) as f:
        old = [r for r in csv.DictReader(f) if r["run_id"].startswith("JG")]
    same = len(old) == len(rows) and all(
        r["sobol_block"] == blk and [float(r[f[0]]) for f in FACTORS[:5]] == list(v)
        for r, (_, blk, v) in zip(old, rows))
    print(f"{path}: regenerated exactly: {same}")
    return same


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("out", nargs="?", help="CSV to write (must not exist)")
    p.add_argument("--check-original", metavar="CSV", help="regenerate the 5-factor matrix and compare")
    args = p.parse_args()
    if args.check_original:
        raise SystemExit(0 if check_original(args.check_original) else 1)
    if args.out:
        import os
        if os.path.exists(args.out):
            raise SystemExit(f"{args.out} exists; matrices are never overwritten")
        write(args.out)
