# Camp Fire local sensitivity study: arrival time vs wind direction

A one-at-a-time (local) study, separate from the Sobol experiment in `../experiment/`.
Set up 2026-10-26 at the suggestion of Jacob's advisor. The quantity of interest is
still the fire arrival time; the Paradise wildland border band
(`../experiment/qgis_setup/paradise_border_mask.tif`) works here unchanged.

## What a run is

For each wind station (Jarbo Gap JG, Openshaw OP, Cohasset CO):

1. **Aim.** Take the station's observed hourly wind directions over the 24 h window
   (25 bands) and their mean. Find the one rotation, the *aim*, that turns the mean
   into the direction blowing from the ignition straight at Paradise.
2. **Sweep.** Run `<station>_WD<offset>` turns every hour by `aim + offset`, for
   offset = -50, -49, ..., +50 degrees:
   `wd = (wd_obs + aim + offset) mod 360`.
   `WD+00` is aimed at Paradise. Positive offsets turn the wind clockwise: it blows
   toward 246.4 + offset degrees, so the fire passes north of Paradise. Negative
   offsets send it south. For JG, `WD-31` is close to the unturned observed wind.
3. **Nothing else changes.** Wind speed and the 1/10/100-h dead fuel moisture are
   the observed hourly values. `P_EPS` is the baseline case's 0.01. The model
   settings match the Sobol runs: same case, ignition, 24 h window, no gusts,
   fixed seed, spotting on.

That is 101 runs per station and 303 in total. ELMFIRE is deterministic with these
settings, so each run is one fire.

**Numbers** (from `prepared/aim.json`, written by `prepare`):

| | value |
|---|---|
| Target | centroid of Paradise's 2018 town limits (Census TIGER GEOID 0655520), 15.70 km from the ignition |
| Wind must blow from | 66.40 deg (toward 246.40 deg), on the model grid |

| Station | Mean wind from | Aim | Rotations swept |
|---|---|---|---|
| JG | 35.54 deg | +30.86 deg | -19.14 .. +80.86 |
| OP | 351.24 deg | +75.17 deg | +25.17 .. +125.17 |
| CO | 45.43 deg | +20.97 deg | -29.03 .. +70.97 |

**Choices behind these numbers** (all in `config.yaml`):
- **Mean direction**: `study.mean: vector`, the direction of the mean wind vector,
  where each hour counts by its speed. The arithmetic mean of the angles is wrong
  for directions that straddle north: for Openshaw it gives 188 deg instead of
  351 deg. The alternative, `unit`, counts every hour the same. It agrees within
  1.2 deg for JG and OP, but gives 59.58 deg for CO, 14 deg away, because CO's
  wind is light (3.4 mph mean) and variable.
- **Target point**: `study.target_lonlat`, the centroid of the town limits.
- **Grid north**: ELMFIRE reads directions relative to the raster grid (the case sets
  no `GRID_DECLINATION`), so the bearing to Paradise is computed on the grid. It is
  1.0 deg from true north there. Station azimuths are used unchanged, as in the
  Sobol runs.

Changing any of these means rerunning `prepare` and using a new
`paths.scratch_root`, because validated runs are skipped on resubmission.

## Fixed-direction variant (`fixed_direction/`, 2026-10-28)

The same 303 runs, except that the direction doesn't change from hour to hour.
Every hour of run `<station>_WD<offset>` blows from 66.40° + offset, which points
straight at Paradise when the offset is 0. Wind speed and the three moistures are
the recorded hourly values, as before. Jacob asked for it after seeing the hourly
runs, to isolate wind direction from the stations' hour-to-hour swings. The
hourly runs stay the main, more physical set.

- Config: `fixed_direction/config.yaml`. It differs from `config.yaml` only in
  `study.direction: fixed`, its relative paths, its scratch root
  (`$SCRATCH/elmfire-class/local_sensitivity_fixed`) and its Slurm script
  (`slurm/lsa_batch_fixed.slurm`). A test checks that this stays true.
- Every command takes the config, from this folder:
  `python -m lsa --config fixed_direction/config.yaml status` (likewise `collect JG`, ...).
- Its `prepared/`, `batches/`, `results/` and `logs/` are under `fixed_direction/`.
- Submit:
  ```bash
  sbatch --array=0-2%3 slurm/lsa_batch_fixed.slurm fixed_direction/batches/all.txt
  ```
- Analysis: set `RUN_SET = "fixed"` in `analysis/lsa_gradient.py`. Its figures
  and GeoTIFFs go to `fixed_direction/figures/` and `fixed_direction/qgis/`.

## Steps

Run from this folder in a compute session:

```bash
cd /oak/stanford/groups/gorle/jgoell/Wildfire/local_sensitivity
ml python/3.12.1     # already loaded in the VS Code session
source /oak/stanford/groups/gorle/jgoell/venvs/elmfire-vnv/bin/activate
```

1. **Tests** (seconds, no ELMFIRE): `python -m unittest -v tests.test_lsa`
2. **Prepare** (done 2026-10-26): `python -m lsa prepare`. This reads the station CSVs,
   writes `prepared/forcing_<S>.csv` and `prepared/aim.json`, and copies the case's
   fuel inputs and `ignition.csv` to `$SCRATCH/elmfire-class/local_sensitivity/shared/`.
3. **Inspect one case** (optional):
   `python -m lsa build JG_WD+00 --out $SCRATCH/preview_lsa_JG_WD+00`
4. **Batch**: `python -m lsa make-batch` writes `batches/all.txt` and prints the
   sbatch line. There are 3 array tasks, one per station; each runs 101 fires, 24 at
   a time.
   ```bash
   sbatch --array=0-2%3 slurm/lsa_batch.slurm batches/all.txt
   ```
5. **Check**: `python -m lsa status`. Failed or unfinished runs:
   `python -m lsa make-batch retry --of batches/all.txt`, then the printed line.
6. **Collect**: `python -m lsa collect JG` (and OP, CO). This writes `results/<S>_runs.csv`.
   With `--mask ../experiment/qgis_setup/paradise_border_mask.tif`, it also writes
   `results/<S>_arrival.npz`, the arrival times at those cells.

## Where things are

| Path | Contents |
|---|---|
| `config.yaml` | Every setting |
| `elmfire.data.template` | ELMFIRE namelist; its header lists the differences from the baseline case |
| `lsa/` | The code: `forcing.py` (station data, mean, aim), `case.py` (one run's inputs), `weather.py` (rasters), `run.py` (running, validation, batches), `collect.py` |
| `prepared/` | Forcing tables (identical to `../experiment/prepared/`) and `aim.json` |
| `$SCRATCH/elmfire-class/local_sensitivity/runs/<S>/<run_id>/` | One folder per run |
| `results/` | From `collect` |
| `logs/` | Slurm logs |

A run folder has the same files as in the Sobol experiment:
- `time_of_arrival_0000001_0086400.tif`: seconds after t0 (06:00 PST 8 Nov 2018). Unburned cells are -9999; the ignition cell is -1.
- `fire_size_stats.csv`
- `run.json`: offset, aim, rotation, mean and target directions, checksums
- `weather.json`: observed and written values of every band
- `inputs/{ws,wd,m1,m10,m100}.tif`
- `elmfire.data.in` and the ELMFIRE logs
- `attempt.json`, and `DONE` (only once validation has passed)

## Getting the data out

`results/<S>_runs.csv` has one row per run, in offset order from -50 to +50. Its
columns are:
- `run_id`, `offset_deg`, `aim_deg`, `rotation_deg`, `mean_from_deg`, `target_from_deg`;
- the outputs `burned_area_ac`, `crown_area_ac`, `burned_cells`, `first_arrival_s`,
  `tstop_h` and `stalled`;
- `toa_path`, the run's arrival-time raster.

```python
import numpy as np, pandas as pd, rasterio
mask = rasterio.open("../experiment/qgis_setup/paradise_border_mask.tif").read(1).astype(bool)
runs = pd.read_csv("results/JG_runs.csv")                      # offset order
toa = np.stack([rasterio.open(p).read(1)[mask] for p in runs.toa_path])
toa[toa <= 0] = np.nan                                          # not burned by 24 h
```

**Local sensitivity map**: `python analysis/lsa_gradient.py`, after `collect` for all
three stations. It estimates dT/dθ for every cell by central difference at θ0
(`THETA0`, `H` at the top of the script). It writes:
- `figures/lsa_gradient_theta+00_h1.png`, with one contour map per station;
- for QGIS, `qgis/<S>_dTdtheta_theta+00_h1.tif` (minutes per degree, with a matching
  `.qml` style) and `qgis/<S>_toa_theta+00_hours.tif` (arrival time of the θ0 run,
  in hours). These are on exactly the fuel grid, like `../experiment/qgis_setup/`,
  with NoData -9999.

It also prints summary numbers for the border band.

The MATLAB route is the same as for the Sobol tables: read `toa_path` row by row and
index with the mask. The rows are in offset order, not design_row order, and the
input column is `offset_deg`.

Scratch is purged after 90 days without changes. Copy anything that must be kept
to Oak. Groupmates can't read Oak, so share copies through
`$SCRATCH/elmfire-class` as for the Sobol results.
