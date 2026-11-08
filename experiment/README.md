# Camp Fire ELMFIRE Sobol experiment

Runs a 6-factor Sobol design through ELMFIRE and collects the outputs in Sobol
order: `manifest/Camp_Fire_Sobol_N128_D6_run_matrix_5376.csv`, 3 wind stations
x 1,792 rows. Spec: `../Camp_Fire_ELMFIRE_coding_agent_handoff.md`. Every
setting is in `config.yaml`; the ELMFIRE namelist is `elmfire.data.template`,
whose header lists every difference from the baseline case.

The 6-factor matrix adds `q_spot` to the handoff's five factors (decided
2026-10-19). `tools/make_matrix.py` wrote it with the same construction and
seed as the handoff's 5-factor matrix. That script regenerates the original
`Camp_Fire_Sobol_N128_run_matrix_4608.csv` exactly, and the original stays in
`manifest/` unchanged.

## What a run is

- Weather: the station's observed hourly wind (Jarbo Gap, Openshaw or Cohasset)
  and Jarbo Gap NFDRS dead fuel moisture, perturbed by the run's five
  coefficients:
  `ws = a_U * ws_obs`, `wd = (wd_obs + b_theta_deg) mod 360`,
  `m1 = m1_obs + b_1_pp` (percent), and likewise for m10 and m100.
  The values are the same everywhere in the domain, on a 19 x 19 grid of 3 km
  cells. There are 25 hourly bands, from 2018-11-08 14:00Z (06:00 PST) through
  24 h later.
- Fire: one ignition at 39.810283, -121.437222 (fuel cell row 447, col 1186,
  model 182) at t = 1800 s (06:30 PST). The run stops at t = 86,400 s.
- Spotting reach: the sixth factor, `q_spot` ~ Uniform(0, 1), sets
  `P_EPS = 0.001 * 50**q_spot`. That is log-uniform on [0.001, 0.05]; smaller
  `P_EPS` lets embers land farther away. It is constant within a run and
  written into the run's spotting settings. The nominal controls use
  `P_EPS = 0.01` (`q_spot = ln 10 / ln 50`). These are exploratory sensitivity
  bounds, not a calibrated uncertainty distribution. Each run's `run.json`
  records `q_spot`, `P_EPS`, and the exact text written to the namelist.
- Model settings: those of the baseline case except no gusts, a fixed seed
  (2023), and no internal weather offsets. The ADJ (+0.2839) and WAF (-0.0197)
  offsets are kept. Spotting is on, with the baseline's other spotting
  settings. With these settings ELMFIRE's spotting has no randomness, so the
  same row always gives the same fire.
- Each model-setting variant gets its own `paths.scratch_root`, so runs made
  with different settings never mix. The 6-factor design uses
  `$SCRATCH/elmfire-class/camp_sobol_d6`. Earlier 5-factor runs are in
  `camp_sobol` (spotting off) and `camp_sobol_spot` (spotting on, plus the
  spotting study).

## Steps

Run everything from this folder in an `sh_dev` session (or a job). The only
setup is the venv:

```bash
cd /oak/stanford/groups/gorle/jgoell/Wildfire/experiment
ml python/3.12.1     # already loaded in the VS Code session; the venv's Python needs it
source /oak/stanford/groups/gorle/jgoell/venvs/elmfire-vnv/bin/activate
```

ELMFIRE itself also needs `ml gcc/12.4.0 openmpi/5.0.5 physics gdal/3.10.2`. The
Slurm script (`slurm/sobol_batch.slurm`) loads these. You only need them yourself
for `run-one`, which runs a fire in your own session. Like the interiorDataset
scripts, the Slurm script writes its logs to `logs/` here.

1. **Tests** (seconds, no ELMFIRE): `python -m unittest -v tests.test_pipeline`
2. **Validate** (reads everything, writes nothing): `python -m sobolfire validate`
3. **Prepare**: `python -m sobolfire prepare` writes `prepared/forcing_<JG|OP|CO>.csv`
   and `prepared/provenance.json`. It also copies the case's fuel inputs plus
   `ignition.csv` to `$SCRATCH/elmfire-class/camp_sobol/shared/`.
4. **Inspect one case** (optional): `python -m sobolfire build JG_0001 --out $SCRATCH/preview_JG_0001`
   shows the namelist, rasters and `inputs/weather.json` exactly as a run gets them.
5. **Pilot**: `python -m sobolfire make-batch pilot` writes `batches/pilot.txt`
   (46 runs: 3 nominal controls, a repeat of `JG_BASELINE`, and the 14 rows of
   base row 1 for each station) and prints the `sbatch` line. Afterwards:
   - `python -m sobolfire status batches/pilot.txt`
   - `python -m sobolfire compare JG_BASELINE JG_BASELINE__rep2` should report
     `identical: True`.
   - Check `max_rss_mb` and `wall_s` in the runs' `DONE` files, and the
     `seff` of the array tasks. Use them to set CPUs, time and memory in the
     `#SBATCH` lines of `slurm/sobol_batch.slurm`, and `slurm.runs_per_task`
     in `config.yaml`.
6. **Full design**: `python -m sobolfire make-batch all`, then the printed `sbatch` line.
7. **Retries**: after the array has finished,
   `python -m sobolfire make-batch retry --of batches/all.txt`. This reruns the
   same IDs with the same coefficients. The earlier attempt's files move to
   `attempts/<k>/` in the run folder. Validated runs are never rerun unless you
   use `run-one --force`.
8. **Collect**: `python -m sobolfire assemble JG` (and OP, CO). This refuses to
   run until all 1,792 runs of the station are validated.

## Where things are

| Path | Contents |
|---|---|
| `$SCRATCH/elmfire-class/camp_sobol_d6/runs/<JG\|OP\|CO>/b<base_row>/<run_id>/` | One folder per run (the 14 rows of a base row together) |
| `$SCRATCH/elmfire-class/camp_sobol_d6/controls/<run_id>/` | Nominal controls (`JG_BASELINE`, ...) and repeats (`..__rep2`) |
| `results/<station>_runs.csv` | From `assemble`: one row per run in design_row order |
| `prepared/` | Forcing tables and provenance |
| `logs/` | Slurm logs, one per array task |

A run folder holds:

| File | Contents |
|---|---|
| `time_of_arrival_0000001_0086400.tif` | Seconds after t0 (2018-11-08 14:00Z) when each 30 m cell burned. Cells not burned by the end are -9999. The ignition cell itself is -1 (ELMFIRE's marker); its neighbours burn from about 2150 s. |
| `fire_size_stats.csv` | ELMFIRE's summary: burned and crown-fire area (acres), run time |
| `run.json` | Coefficients (as numbers and as the CSV text), station, design row, checksums of every input, seed policy |
| `weather.json` | Observed and written value of every band, with UTC/PST times, and each weather raster's sha256 |
| `inputs/ws.tif`, `wd.tif`, `m1.tif`, `m10.tif`, `m100.tif` | The run's weather rasters (25 hourly bands, uniform in space). Kept by every run since 2026-10-20; for earlier runs, `sbatch slurm/restore_inputs.slurm` rebuilds them and checks them against `weather.json` |
| `elmfire.data.in`, `elmfire.stdout`, `elmfire.stderr` | The exact namelist and ELMFIRE's output |
| `attempt.json`, `DONE` | Status, timing, memory; `DONE` exists only if the outputs passed validation |

Scratch is purged after 90 days without changes. Copy anything you need to
keep (the `results/` files are already on Oak).

## Getting the data out

The rows of every station file are in design_row order, 1 to 1,792. That is
the order SALib's `sobol.analyze(..., calc_second_order=True)` expects, with
the factors in the CSV column order. Analyze each station separately. The
sixth factor is `q_spot` on [0, 1], the variable that was sampled uniformly.
`P_EPS` is a derived column, recorded so you can report actual values.

```python
import numpy as np, pandas as pd, rasterio
from pyproj import Transformer

runs = pd.read_csv("results/JG_runs.csv")          # inputs + outputs, Sobol order
X = runs[["a_U", "b_theta_deg", "b_1_pp", "b_10_pp", "b_100_pp", "q_spot"]].to_numpy()
p_eps = runs["P_EPS"].to_numpy()                     # = 0.001 * 50**q_spot, for plots/tables
area = runs["burned_area_ac"].to_numpy()             # one scalar response per run

# Arrival time at one point, every run (example: one point in central Paradise;
# the analysis itself uses the border cells below)
x, y = Transformer.from_crs("EPSG:4326", "EPSG:26910", always_xy=True).transform(-121.6219, 39.7596)
with rasterio.open(runs.toa_path[0]) as f:
    row, col = f.index(x, y)
toa = np.array([rasterio.open(p).read(1)[row, col] for p in runs.toa_path])
burned = toa > 0                                     # not burned by 24 h otherwise
```

Unburned cells have no arrival time, so turn them into a finite quantity before
any variance-based analysis. One option is `np.where(burned, toa, 86400)`,
which is censored at 24 h; another is `burned.astype(float)`. These answer
different questions.

**The Paradise cells used for the analysis**: `qgis_setup/paradise_border_mask.tif`
(chosen 2026-10-21). It marks the wildland (burnable, non-urban) cells whose centres
lie within 1 km outside Paradise's 2018 town limits: 33,758 cells, on exactly the fuel
grid. `tools/qgis_setup_export.py` makes it; `BORDER_M` sets the width. Any other region
works the same way, as a raster on the fuel grid (same size and transform as
`shared/fuels/fbfm40b.tif`) that is nonzero inside the region.

In MATLAB, `matlab/extract_border_toa.m` collects the arrival times at those cells for
one station into a `.mat` file. It saves `T`, `Tcap` (capped at 24 h), `burned24`,
`X`, the blocks and the cell positions; its header explains them. In Python, run
`python -m sobolfire assemble JG --mask qgis_setup/paradise_border_mask.tif`, which
writes `results/JG_arrival.npz`:

- `arrival_s`: shape (1792, cells), NaN where not burned;
- `row`, `col`, `x`, `y`: the cells;
- `fbfm`, `unburnable`: fuel model, and whether the cell can never burn (urban
  cells, model 91, can burn through the building model).

```python
d = np.load("results/JG_arrival.npz")
A = d["arrival_s"]
```
