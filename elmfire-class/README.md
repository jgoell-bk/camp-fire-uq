# Running the 2018 Camp Fire case with ELMFIRE

## What's here

| Path | Contents |
|---|---|
| `~/elmfire/` | ELMFIRE 1.1 (commit `cbf924a`), prebuilt with gcc 12.4.0 + OpenMPI 5.0.5. The executable is `build/linux/bin/elmfire_1.1`. |
| `~/elmfire-class/camp_fire/` | The Camp Fire case (from `ELMFIRE_VnV_Suite`, Validation/landscape_scale) |
| `~/elmfire-class/camp_fire.slurm` | Job script that runs the case |
| `$SCRATCH/elmfire-class/` | Where every run and its Slurm log goes |

Inside the case folder:

| Path | Contents |
|---|---|
| `elmfire.data.in` | ELMFIRE settings (namelist). All paths in it are relative to the case folder. |
| `case.yaml` | Validation settings: comparison time (129,600 s = 36 h), probability threshold (0.5), expected member count (50) |
| `data/archive_payload/fuels_weather/` | All input rasters: terrain, fuels, canopy, weather (37 hourly bands), ignition mask, buildings. 1457 × 1476 cells at 30 m, EPSG:26910. |
| `data/archive_payload/misc/` | `fuel_models.csv`, `building_fuel_models.csv` |
| `data/viirs_observation/` | VIIRS satellite fire detections used for scoring |
| `scripts/postprocess.py` | Scores the run against VIIRS and makes figures |

## One-time setup: Python environment

ELMFIRE itself only needs modules. The scoring step (`postprocess.py`) needs a
venv with numpy, matplotlib, rasterio, pyyaml, pandas, geopandas, pyproj and
shapely ≥ 2.0. The job script uses `$OAK/jgoell/venvs/elmfire-vnv` unless
`VENV` is set.

To make your own (from `sh_dev`). `--only-binary` is required: Sherlock's
glibc 2.17 is too old for current numpy source builds.

```bash
sh_dev
module purge
module load python/3.12.1
python3 -m venv $GROUP_HOME/venvs/elmfire-vnv
source $GROUP_HOME/venvs/elmfire-vnv/bin/activate
export PIP_CACHE_DIR=$SCRATCH/.cache/pip
pip install --upgrade pip
pip install --only-binary=:all: numpy matplotlib rasterio pyyaml pandas geopandas pyproj "shapely>=2.0"
exit
```

## Running the case

Submit from the scratch project folder so the Slurm log lands there:

```bash
mkdir -p $SCRATCH/elmfire-class
cd $SCRATCH/elmfire-class
sbatch ~/elmfire-class/camp_fire.slurm
```

The job:

1. Copies the case to `$SCRATCH/elmfire-class/camp-<jobid>/`.
2. Runs ELMFIRE there on 50 MPI ranks, one fire per rank.
3. Deletes ELMFIRE's temporary converted inputs (`logs/scratch`, ~1.6 GB).
4. Runs `postprocess.py` and prints the main scores at the end of the log.

Resources: one `serc` node, 50 cores, 80 GB, 3 h limit. The first full run
took 45 minutes and used about 56 GB.

Optional overrides, set before `sbatch`:

| Variable | Default | Use |
|---|---|---|
| `ELMFIRE` | `~/elmfire` | A different ELMFIRE copy |
| `CASE` | `~/elmfire-class/camp_fire` | A modified copy of the case |
| `VENV` | `$OAK/jgoell/venvs/elmfire-vnv` | Your own Python venv |

## What you get

In `$SCRATCH/elmfire-class/camp-<jobid>/`:

| File | Contents |
|---|---|
| `outputs/time_of_arrival_<member>_0129600.tif` | One per member. Seconds since t = 0 when each cell burned. **−9999 (NoData) = never burned**; **−1 marks the member's ignition cell** (corrected 2026-10-26; this table used to say −1 = never burned). |
| `outputs/flame_length_<member>_0129600.tif` | Flame length, **feet** |
| `outputs/ember_flux_<member>_0129600.tif` | Ember flux. Has unexplained negative values; don't rely on it. |
| `outputs/fire_size_stats.csv` | One row per member: ignition x/y, burned area and crown-fire area (acres), wall-clock time |
| `outputs/coeffs.csv` | One row per member: the input tweaks it used |
| `outputs/metrics.json` | Scores vs VIIRS: Jaccard, precision, recall, area ratio, Brier score |
| `figures/*.pdf` | `validation_comparison`, `arrival_time_summary`, `burned_area_growth`, `optional_output_fields` |
| `logs/elmfire.stdout`, `logs/elmfire.stderr` | ELMFIRE's own output; look here first if the run fails |
| `elmfire-camp-<jobid>.out` (in `$SCRATCH/elmfire-class/`) | Slurm log, ending with the main scores |

The burn probability of a cell is the fraction of members whose
`time_of_arrival` is between 0 and 129,600 s. `postprocess.py` scores the cells
with probability ≥ 0.5 against the VIIRS footprint.

## Changing the case

Copy the case first and point `CASE` at the copy, so the original stays
intact:

```bash
cp -r ~/elmfire-class/camp_fire ~/elmfire-class/camp_fire_test
# edit ~/elmfire-class/camp_fire_test/elmfire.data.in
cd $SCRATCH/elmfire-class
CASE=~/elmfire-class/camp_fire_test sbatch ~/elmfire-class/camp_fire.slurm
```

Settings you're most likely to change in `elmfire.data.in`:

| Setting | Group | Camp value | Notes |
|---|---|---|---|
| `NUM_ENSEMBLE_MEMBERS` | `&MONTE_CARLO` | 50 | Fewer members than cores wastes cores; more members than cores runs in rounds. If you change it, change `expected_members` in `case.yaml` too (`MONTE_CARLO.NUM_ENSEMBLE_MEMBERS`), or scoring reports NOT EVALUABLE. |
| `SIMULATION_TSTOP` | `&TIME_CONTROL` | 129600 (36 h) | Scoring needs every member to reach 129,600 s |
| `PDF_LOWER_LIMIT(i)` / `PDF_UPPER_LIMIT(i)` | `&MONTE_CARLO` | Equal (fixed tweaks) | Set them apart and each member draws its own value between them. Tweaks are added to the whole raster. Moisture tweaks are in percentage points (−1.0 = −1 point). |
| `WIND_FLUCTUATIONS` | `&SIMULATOR` | `.TRUE.` | Random gusts every 30 s. `.FALSE.` removes them. |
| `RANDOMIZE_RANDOM_SEED` | `&SIMULATOR` | `.TRUE.` | Seeds from the clock, so runs aren't repeatable. `.FALSE.` uses `SEED` on **every** core, so all members draw the same random numbers. |

Repeatable single-fire runs (for sensitivity studies): set
`WIND_FLUCTUATIONS = .FALSE.`, `RANDOMIZE_RANDOM_SEED = .FALSE.` and
`NUM_ENSEMBLE_MEMBERS = 1` (and the member count in `case.yaml` to 1). Every
run then starts at the same ignition cell with no random gusts. Check this by
running the same settings twice and comparing the `time_of_arrival` rasters.

The core count is set by `--ntasks` in the job script (one fire per core). For
a single fire use `sbatch --ntasks=1 --mem=8G ...` (the weather rasters alone
take about 1.3 GB in memory); command-line options override the `#SBATCH`
lines.

## Known issues

- **Report PDF.** The case's LaTeX report isn't built; the job runs the
  simulation and `postprocess.py` only.
- **Time zero.** `case.yaml` labels t = 0 as 2018-11-08 14:00 Pacific, but the
  fire was reported around 06:30 Pacific and the first VIIRS detection
  (19:50 UTC) is before the configured start. This affects the VIIRS scores.
- **Weather provenance.** Not recorded. The wind is blocky at 1–3 km
  (consistent with a ~3 km hourly model such as HRRR), and dead fuel moisture
  is a constant 6 / 7 / 8 % everywhere.
- **Scratch purge.** Files on `$SCRATCH` unmodified for 90 days are deleted.
  Copy results you need to keep to `$HOME` or Oak.
