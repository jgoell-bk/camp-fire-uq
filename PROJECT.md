# Camp Fire ELMFIRE sensitivity study (UQ class project): project state

Start here. This is the single entry point for agents working on this project. It
records what exists, the decisions behind it and their reasons, the current state, and
how to work with Jacob (the user). It was last updated on 2026-10-26. Things not
written here are recorded in the files it links to.

**Before doing anything on Sherlock, read `/etc/agents/AGENTS.md`**, at the start of
every session and again after any context summary or compaction (its rule 12). Read its
topic files (`storage.md`, `slurm.md`, `software.md`, `policy.md`) before the work they
cover. Keep this requirement in any handoff you write.

Companion: `ELMFIRE_CODEBASE.md` (next to this file) covers how ELMFIRE works
internally, with the source-level facts this project depends on.

## 1. What the project is

This is a class project in uncertainty quantification (Jacob's group). It measures how
uncertain wind, dead fuel moisture and spotting reach change the ELMFIRE-simulated
2018 Camp Fire. The quantity of interest is **fire arrival time in the wildland cells
bordering Paradise, CA**. The analysis is a variance-based (Sobol) global sensitivity
analysis, run separately for three wind-station scenarios.

The specification is `Camp_Fire_ELMFIRE_coding_agent_handoff.md` (2026-10-19). Later
decisions extend it (section 4). Older documents:
- `uq_project_plan.md` (2026-10-08) is **superseded**. It describes Morris screening and
  a GP surrogate on gridded weather, and was never carried out.
- `class_project_setup.md` and `elmfire-class/README.md` describe the personal ELMFIRE
  copy and the baseline case setup; they are still accurate. The README's note on
  unburned cells was corrected on 2026-10-26.

**Division of labour: deliver data, not analysis.** Agents prepare runs and extract
clean data, and show how to load it. The group does the Sobol analysis, chooses the
quantity of interest treatment, and interprets results. Don't compute Sobol indices or
write analysis notebooks unless asked.

## 2. Current state (2026-10-26; EKI stage added 2026-11-08)

- **Current stage (from 2026-11-08): EKI parameter estimation** in `eki_calibration/`.
  - It calibrates time-varying corrections to the JG+CO mean weather against observed
    arrival times. The plan is `ParameterEstimationFiles/camp_fire_eki_calibration_plan.md`.
  - The observation table (linked to the NIST TN 2135 records), baseline, frozen
    holdout split and 32-member initial ensemble are prepared, and an EKI update is
    drafted. Batch I00 (33 runs, job 51576874) is done, extracted and reviewed
    (2026-11-08).
  - The target is one-sided presence deadlines, fitted through the distance to the
    burned area at each deadline (block-weighted, L = 1 km).
  - The first controlled update is proposed (`proposals/I01.csv`, α = 32). It waits on
    the CHK01 check runs at the ensemble means, then `verify` and `promote`. I01 is not
    submitted.
  - Ignition is at 06:15 PST (t = 900 s) there, so it differs from the earlier studies.
  - Start from `eki_calibration/README.md`, which has the decisions, status and open
    items.

- **The 6-factor design is complete.** All 5,376 runs passed validation: job 48588168,
  839 core-hours, none stalled. The weather rasters of every run were restored into
  each run folder by job 48701382 (5,374 restored, 2 already present).
- **Results tables are written**: `experiment/results/{JG,OP,CO}_runs.csv`.
- **The QoI region is chosen**: the wildland border band around Paradise (section 6).
- **Sobol analysis under way (group's own work)**: as of 2026-10-26,
  `experiment/matlab/` holds `{JG,OP,CO}_border_toa.mat` (dated 10-22) and Jacob's
  `sobol.m`, `main_sobol.m` and `JG_sobol_results.mat`. The group works in MATLAB with
  a homework `sobol.m` (pick-freeze estimator); data extraction is described in
  section 6.
- **Local sensitivity study (wind direction)**, at the suggestion of Jacob's advisor.
  All 303 runs validated (job 49725936, 2026-10-26). Analysis script:
  `local_sensitivity/analysis/lsa_gradient.py`, a contour map of dT/dθ by central
  difference. Written for Jacob to run; the agent has not run it.
- **Fixed-direction variant of the local study (set up 2026-10-28)**. The first
  submission, job 50023648, failed for all 303 runs before ELMFIRE started: a log
  line in `lsa/run.py` formatted the fixed runs' empty `rotation_deg`. That was an
  agent bug, fixed 2026-10-29, and a test now runs `run_one` end to end for both
  kinds of run. Resubmit the same batch:
  `local_sensitivity/fixed_direction/`, 303 runs. Every hour blows from 66.40° + θ,
  with speed and moisture hourly. It has been prepared and tested, and its inputs
  match the hourly set's. Jacob submits it, from `local_sensitivity/`:
  `sbatch --array=0-2%3 slurm/lsa_batch_fixed.slurm fixed_direction/batches/all.txt`. It is in its own folder,
  `local_sensitivity/`; its `README.md` describes it and is the place to start. Each
  station's hourly directions are turned so that their mean (vector mean) blows from
  the ignition at Paradise. The sweep then turns them -50..+50° in 1° steps, with
  nothing else perturbed: 303 runs. Its inputs match the Sobol pipeline's. The
  vector mean and the town-centroid target are the agent's proposals and are still
  awaiting Jacob's explicit confirmation.
- **Shared copies made**: `results/` and `qgis_setup/` now also exist in
  `$SCRATCH/elmfire-class/camp_sobol_d6/` for groupmates.
- **Optional, not done**:
  - running the three nominal controls (`make-batch controls`);
  - extending N from 128 to 256 (section 5);
  - copying run data to Oak before the scratch purge.

## 3. Where things are

| Path | What |
|---|---|
| `/oak/stanford/groups/gorle/jgoell/Wildfire/` | Project root on Oak (this file). Oak is gorle's space; groupmates **cannot** read it. |
| `Wildfire/experiment/` | The pipeline (`sobolfire` package), config, design matrices, results, QGIS bundles. Its `README.md` is the usage guide. |
| `Wildfire/eki_calibration/` | The EKI calibration stage (2026-11-08): package `eki`, config, ensembles, observation table, README. Run data in `$SCRATCH/elmfire-class/eki_calibration/`. Inputs from `Wildfire/ParameterEstimationFiles/`. |
| `Wildfire/station_baselines/` | Station baselines (2026-11-09): one deterministic run each for JG, CO and the JG+CO mean, at the EKI settings (06:15 ignition), scored on the border band. Package `baselines`; run data in `$SCRATCH/elmfire-class/station_baselines/`. Spotting replicates were dropped because the configured model is deterministic (see its README). |
| `Wildfire/local_sensitivity/` | The local wind-direction study: its own config, namelist template, `lsa` package, Slurm script, tests and README. It is independent of `experiment/`. |
| `$SCRATCH/elmfire-class/local_sensitivity/` | Its run data: `shared/`, `runs/<JG\|OP\|CO>/<run_id>/` |
| `Wildfire/HistoricalDataClass/` | Source data: the hourly RAWS wind CSV for the three stations, the Jarbo Gap NFDRS moisture CSV, and the original 5-factor matrix CSV. The monthly summaries in there are **not** used. |
| `~/elmfire` | ELMFIRE 1.1, commit `cbf924a`, prebuilt (gcc 12.4.0 + OpenMPI 5.0.5). Executable: `~/elmfire/build/linux/bin/elmfire_1.1` |
| `~/elmfire-class/camp_fire` | Baseline Camp Fire case (from `$OAK/jgoell/ELMFIRE_VnV_Suite/.../camp_fire`), read only |
| `$SCRATCH/elmfire-class/camp_sobol_d6/` | **Current run data**: `shared/` (inputs common to all runs), `runs/<JG\|OP\|CO>/b<base_row>/<run_id>/`, `controls/`. 7.9 GB+. |
| `$SCRATCH/elmfire-class/camp_sobol/` | 5-factor spotting-off pilot (superseded) |
| `$SCRATCH/elmfire-class/camp_sobol_spot/` | 5-factor spotting-on pilot plus the spotting study (superseded, but informative) |
| `Wildfire/results/camp_fire/` | The earlier 50-member validation-suite run of the baseline case |

`$SCRATCH` is `/scratch/users/jgoell`. Scratch files are purged after 90 days without
modification. Never refresh files to dodge the purge (Sherlock rule); copy anything
needed to Oak or `$GROUP_HOME` instead.

## 4. The experiment as run

### Model setup (`experiment/config.yaml`, `experiment/elmfire.data.template`)
- **Case**: the baseline Camp Fire case. 30 m grid of 1,457 × 1,476 cells (2,150,532),
  EPSG:26910, about 44 × 44 km. Of these, 84.1% are burnable wildland, 3.5% urban
  (fuel 91), 9.4% non-burnable and 2.9% no-data border.
- **Time**: t0 is 2018-11-08 14:00Z (06:00 PST). Ignition at t = 1800 s (06:30 PST).
  The run stops at t = 86,400 s (24 h). There are 25 hourly weather bands.
- **Ignition**: one fixed cell at 39.810283, −121.437222 (row 447, col 1186, fuel 182
  TL2), via `CSV_FIXED_IGNITION_LOCATIONS`. Jacob supplied the coordinate. It lies about
  0.8 km east of the baseline case's random-ignition mask.
- **Weather**: hourly station observations, made spatially uniform on a 19 × 19 grid of
  3 km cells in the case CRS (margin of 2 cells, see the ELMFIRE doc). Wind comes from
  the scenario station. Dead fuel moisture (1/10/100-h) comes from Jarbo Gap NFDRS
  `fuelModelType=V` for all scenarios.
- **Kept from the baseline**: fuels, terrain, canopy, buildings, live moisture, crown
  fire, the WUI (Hamada) building model, the ADJ +0.2839 and WAF −0.0197 offsets, and the
  other spotting settings.
- **Changed from the baseline**:
  - `WIND_FLUCTUATIONS` off;
  - `RANDOMIZE_RANDOM_SEED` off, `SEED` 2023;
  - one member per run;
  - the baseline's WS +3 mph / WD +45° / M1-M100 offsets removed;
  - spotting on with `P_EPS` from the design;
  - no ember-flux or flame-length outputs.

  With these settings ELMFIRE is fully deterministic: a repeated run is bit-identical,
  as checked in both pilots.

### Uncertain inputs (Sobol factors), all independent
| Factor | Applied as | Distribution | Nominal |
|---|---|---|---|
| `a_U` | WS(t) = a_U · WS_obs(t) | U(0.8, 1.2) | 1 |
| `b_theta_deg` | WD(t) = (WD_obs(t) + b) mod 360 | U(−20°, 20°) | 0 |
| `b_1_pp`, `b_10_pp` | M(t) = M_obs(t) + b (percentage points) | U(−1, 1) | 0 |
| `b_100_pp` | same | U(−2, 2) | 0 |
| `q_spot` | `P_EPS = 0.001·50^q_spot` (log-uniform 0.001-0.05; smaller = longer spotting) | U(0, 1) | 0.588592 (P_EPS 0.01) |

The `P_EPS` bounds came from a literature check Jacob did with another agent. They are
exploratory sensitivity bounds, not a calibrated distribution.

**Scenarios (not factors)**: wind station Jarbo Gap (JG, 041214), Openshaw (OP, 041215)
or Cohasset (CO, 041211). Each scenario uses the same design and is analysed separately.

### Design
- **Construction**: Saltelli sampling with second-order blocks. N = 128 base rows
  (2^7, required for Sobol' balance). Each base row has 2D+2 = 14 runs: A, AB_i (D=6),
  BA_i (D=6), B.
- **Size**: 1,792 runs per station, 5,376 in total. Scrambled Sobol', d = 12, seed
  20261019. A = first 6 columns, B = last 6.
- **Matrix**: `experiment/manifest/Camp_Fire_Sobol_N128_D6_run_matrix_5376.csv`,
  read-only, sha256 `5a378386…`. It has a derived `P_EPS` column, validated exactly
  against the formula.
- **Generator**: `experiment/tools/make_matrix.py`. It reproduces the handoff's original
  5-factor matrix (`…_run_matrix_4608.csv`, sha `6472d320…`, also kept unchanged)
  bit-for-bit.
- **Extensible**: the first 128 points of the N = 256 design with the same seed equal
  the current design (verified). Extending N reuses every existing run.

### Decision log (why things are the way they are)
| Date | Decision | Reason / evidence |
|---|---|---|
| 10-19 | t0 = 14:00Z, not the validation suite's 22:00Z | The case's ignition at t0+1800 s would otherwise be 14:30 PST, after the first VIIRS detection; the fire was reported ~06:30 PST. Later lag correlation of the case's gridded wind with the Jarbo Gap record also peaks at 14:00Z (r 0.61 over 37 h; 0.30 at 22:00Z). |
| 10-19 | 24 h window | Ends exactly at the Brewer & Clements (2020) Jarbo Gap reliability cutoff (06:00 PST 9 Nov). The pipeline refuses later windows unless `allow_after_jarbo_gap_cutoff` is set. |
| 10-19 | Uniform weather on a coarse grid | Station forcing is uniform. ELMFIRE maps fuel cells to weather cells by coordinates, so a small grid works. The Berkeley 2.5 km rasters in the handoff are UTM 11N and unusable. |
| 10-19 | Gusts off, fixed seed, ADJ/WAF kept | Gusts add run-to-run jitter that masquerades as interaction effects. ADJ/WAF are baseline calibration settings, not wind or moisture inputs. |
| 10-19 | Spotting: first off, then on | Off gave fires of 5-13k acres that never came within 8 km of Paradise. On, the fires reach the Paradise wildland border in most runs. |
| 10-19 | Spotting study (5-factor, nominal runs) | Ember rate had no effect; crosswind spread had essentially none. `P_EPS` and `CELL_IGNITION_DELAY` dominate and act alike, with effects larger than the weather/moisture spread. So only one spotting knob became a factor. |
| 10-19 | `q_spot` → `P_EPS` as the 6th factor | Jacob's choice after a literature review. A new versioned matrix; the old one stays immutable. |
| 10-20 | Keep per-run weather rasters | Jacob: "don't delete the rasters"; scratch has room. |
| 10-21 | QoI cells = the wildland border band | Jacob: the project is about the wildland burn, not urban cells (fuel 91 burns via the Hamada building model). |
| 10-28 | Local study keeps the recorded hourly wind variation (each run rotates all 25 hours by one angle) as the main set | Jacob had pictured a fixed direction per run, but chose the hourly variation as "more physical". After an explanation, Jacob said he is "fine with these results": the vector mean and the town-centroid target are accepted. |
| 10-28 | Add a fixed-direction variant (`local_sensitivity/fixed_direction/`) | Jacob: "an additional run with the fixed direction variant (all other inputs in the raster varied hourly)". Built as a config switch (`study.direction`) in the same `lsa` code, with its own config, outputs and scratch root. |
| 10-26 | Local wind-direction study in a separate folder | Advisor's suggestion. Jacob: "a whole new folder", nothing tied to the Sobol pipeline. Aim the mean direction at Paradise, then sweep ±50° in 1° steps, everything else observed/nominal. Defaults chosen by the agent and recorded in `local_sensitivity/README.md`: vector mean, town-limits centroid as the target, bearing on the grid. |

### Known caveats for interpreting results
- **The direction range doesn't reach the observed fire path.** Station winds blow
  toward ~214° (SSW), while Paradise lies at 249° from the ignition. The ±20° range
  cannot reach the baseline's undocumented +45° offset, which aimed the gridded wind at
  Paradise. The fire still reaches the border wildland in most runs (section 6), but
  the modelled spread direction differs from VIIRS (WSW).
- **Burned area is truncated at the domain edge** in OP 1,791/1,792, JG 1,201 and CO
  431 runs. Don't use total burned area without caveats. Paradise metrics are
  unaffected.
- **Censoring at 24 h**: unburned cells have no arrival time and need a finite
  treatment (capped arrival or a burned indicator) before variance-based analysis.
- **Urban cells still affect wildland arrival.** They are excluded from the QoI but
  still spread fire in the model (Hamada), so wildland west of town can be reached
  partly through the urban burn.
- The **nominal controls were never run** in the d6 root.

## 5. The pipeline (`experiment/`)

Package `sobolfire` (Python, run from `experiment/` with the venv). Key commands, all
`python -m sobolfire <cmd>`:

| Command | Does |
|---|---|
| `validate` / `prepare` | Check matrix, sources, case, ignition and grids / also write `prepared/` forcing tables and `<scratch_root>/shared/` |
| `build RUN --out DIR` | Build one case without running it (namelist, rasters, records) |
| `run-one RUN [--force]` | Build, run ELMFIRE and validate one case in the current session (about 4-15 min) |
| `make-batch {pilot,all,controls,spotting,JG,OP,CO,retry --of B}` | Write `batches/<name>.txt` and print the exact `sbatch` line |
| `run-batch` | What the Slurm array calls; `runs_per_task` comes from config |
| `status [BATCH]` | Counts validated, failed, running and not started |
| `assemble STATION [--mask TIF]` | `results/<S>_runs.csv`; with a mask, also `results/<S>_arrival.npz` |
| `restore-inputs` | Rebuild weather rasters into run folders, checked against recorded sha256 |
| `compare A B` | Arrival-time difference between two runs |

**Run IDs**:
- `JG_0001`…`JG_1792`: matrix rows.
- `JG_BASELINE`: nominal control.
- `…__rep2`: a repeat in its own folder.
- `JG_BASELINE__P_EPS-0.001`: a spotting-study override. Refused in d6 because `P_EPS`
  is a factor there.

**Per-run folder contents**:
- `time_of_arrival_0000001_0086400.tif`
- `fire_size_stats.csv`
- `run.json`: coefficients, `P_EPS`, the namelist literal, and checksums of every input
- `weather.json`: observed and written value of every band, plus raster checksums
- `inputs/{ws,wd,m1,m10,m100}.tif`
- `elmfire.data.in` and the ELMFIRE logs
- `attempt.json`; `DONE` exists only after validation passed

Validation checks:
- the grid;
- the value set (>0, −1, −9999);
- coverage of `SIMULATION_TSTOP`;
- the ELMFIRE ignition marker (−1) at the configured cell;
- nothing burning before ignition.

**Slurm**: `experiment/slurm/sobol_batch.slurm` sets serc, `--nodes=1`, sh04
CBASE|CPERF, 24 CPUs, 2G/CPU, 4 h. Each array task runs `runs_per_task` (192) fires, 24
at a time, each as a single-core ELMFIRE singleton (no `mpirun`) in `$L_SCRATCH_JOB`.
Full-run timings: 1.7-25 min per fire (means: JG 11.2, OP 8.5, CO 8.4 min), under
1 GB each. `slurm/restore_inputs.slurm` is the
raster-restore job.

**Tests**: `python -m unittest tests.test_pipeline` (17 tests, seconds, no ELMFIRE).
Run them after any change.

**Other tools**:
- `tools/make_matrix.py`: design generator. Needs SciPy, see section 8.
- `tools/qgis_export.py RUN…`: arrival-time rasters in hours for QGIS.
- `tools/qgis_setup_export.py <census shp>`: the setup-slide layers.

## 6. Quantity of interest and the data for it

**QoI cells: the wildland border band.** These are burnable, non-urban cells whose
centres lie within 1 km outside Paradise's 2018 town limits (US Census TIGER/Line,
GEOID 0655520). That's 33,758 cells.
- **Files**:
  - `experiment/qgis_setup/paradise_border_mask.tif` (1 = QoI cell, exactly on the fuel
    grid);
  - layer `paradise_border_cells` in `qgis_setup/setup.gpkg`.
- **Generator**: `tools/qgis_setup_export.py`, setting `BORDER_M = 1000`.
- **How often fire reaches it** (any burned):

  | Station | Band reached | Earliest arrival |
  |---|---|---|
  | JG | 96% of runs | 3.2 h |
  | OP | 73% | 13.0 h |
  | CO | 90% | 13.2 h |

  Wildland inside the town limits is reached in 78% / 52% / 83% of runs.
- **Earlier metrics were wrong for this.** The centre point and 2 km circle I used
  before were urban-dominated and understated reach; ignore them.

**Extraction** (no pipeline needed):
```python
mask = rasterio.open(".../paradise_border_mask.tif").read(1).astype(bool)
runs = pd.read_csv(".../results/JG_runs.csv")                     # design_row order = SALib order
toa = np.stack([rasterio.open(p).read(1)[mask] for p in runs.toa_path])
```
Arrival times are seconds after t0; values ≤ 0 mean unburned (−9999), or the ignition
marker (−1).

**MATLAB route** (the group computes Sobol indices in MATLAB with its own `sobol.m`):
- **Function**: `experiment/matlab/extract_border_toa.m`. It writes
  `<station>_border_toa.mat` with:
  - `T` (single, 1792 × 33758, NaN = unburned),
  - `Tcap` (unburned set to the horizon; arrival capped at 24 h), `burned24` (1/0),
    `horizon_s` (86,400). Jacob added these three himself;
  - `X` (the 6 inputs), `P_EPS`, `block`, `run_id`,
  - `row`/`col` (1-based, MATLAB column-major cell order, which differs from numpy's
    row-major order).
- **Location**:
  - `experiment/matlab/` on Oak has all three stations (10-22, Jacob's runs);
  - `$SCRATCH/elmfire-class/camp_sobol_d6/matlab/` has only an older
    `JG_border_toa.mat` (10-21, verified against Python). It was written before
    `Tcap`/`burned24` were added, so it lacks them.
- **Learning version**: Jacob prefers to write the extraction himself. The sectioned,
  explained version of the same steps was given in chat on 2026-10-22. Don't run the
  extraction for him unless he asks.
- **Speed**: about 2 min per station with `ml matlab/R2023b`; `imread` reads these
  GeoTIFFs with no toolbox.
- **Blocks**: `reshape(T,14,128,[])` gives (block, base row, cell). In the homework's
  pick-freeze estimator, `C_i` = B with column i from A = the **BA** blocks (layers
  8-13).
- **Centre outputs before computing**: arrival times are large numbers.

**Other QGIS material**:
- `experiment/qgis/`: example runs, VIIRS, context layers.
- `experiment/qgis_setup/`: domain, ignition, stations, town limits, building cells,
  hillshade. Its `README.txt` lists sources.
- Station coordinates (NIFC RAWS): JG 39.73591, −121.48898; OP 39.58983, −121.63516;
  CO 39.87184, −121.76898.

**Sharing with groupmates**: they can't read Oak. Grant ACLs on
`$SCRATCH/elmfire-class` (commands in `class_project_setup.md`, section 3), then copy
`results/` and `qgis_setup/` into `$SCRATCH/elmfire-class/camp_sobol_d6/`. As of
2026-10-22 those copies exist there (Jacob made them); whether the ACLs were granted
hasn't been checked. `controls/` doesn't exist in d6: the nominal runs were never run.

## 7. Pilots and studies (5-factor matrix, kept for reference)
| Job | What | Outcome |
|---|---|---|
| 48557650 | Spotting-off pilot, 40 runs | All ran (21-57 s each). They first failed validation because of a bug in my ignition check (ELMFIRE writes −1 there); now fixed. |
| 48561952 | Spotting-on pilot, 40 runs | All validated; repeat identical; 4-10 min per fire; 96-153k acres. |
| 48577168 | Spotting study, 21 runs (one setting at a time at the nominal runs) | See the decision log. Details are in `camp_sobol_spot/controls/`. |
| 48588168 | Full 6-factor design, 5,376 runs | All validated |
| 48701382 | Restore weather rasters | 5,374 restored, 2 present |

## 8. Working with Jacob on Sherlock (rules for agents)
- **Sherlock policy** is in `/etc/agents/AGENTS.md` and its topic files. Read it at the
  start of every session and after any context compaction (rule 12). Check `hostname`
  and `$SLURM_JOB_ID` first. Jacob normally works inside a compute allocation (an
  OnDemand code-server session on an sh04 node; the node and job ID change between
  sessions).
- **Learning mode**: for class work Jacob wants to understand, such as MATLAB extraction
  and Sobol code, give sectioned, commented code for him to run. Don't run it or
  produce the outputs unless asked.
- **Agents may** run checks, unit tests, `validate`/`prepare`/`build`, and small
  scripts (seconds to a few minutes) in his session. **Jacob runs** every `sbatch` and
  any large Python job. Give him literal, paste-ready commands: no `<placeholders>`;
  use real paths and job IDs.
- **Slurm style**:
  - always `--nodes=1`;
  - `--constraint=CLASS:SH4_CBASE|CLASS:SH4_CPERF`;
  - generous time limits;
  - no `sleep`, no `--exclude`;
  - keep total cores at or below about 600 (his other project runs CharLES arrays at
    the same time);
  - resources go in `#SBATCH` lines in the script, not long `sbatch` command lines;
  - logs go to `experiment/logs/`.
- **Python**: `source /oak/stanford/groups/gorle/jgoell/venvs/elmfire-vnv/bin/activate`.
  It needs `ml python/3.12.1`, which is already loaded in his VS Code session. The venv
  has numpy, pandas, rasterio, pyproj, pyyaml, geopandas, shapely and pyogrio. It has
  **no scipy and no fiona**.
- **SciPy** (only for `make_matrix.py`): `ml math python/3.12.1; ml py-scipy/1.16.0_py312`,
  then call `/share/software/user/open/python/3.12.1/bin/python3 -P` explicitly. A stray
  `/oak/stanford/groups/gorle/jgoell/.venv` (SciPy 1.11.4) shadows `python3` on PATH in
  his VS Code shell.
- **Style preferences**:
  - simple, statable choices (uniform distributions stated in one sentence; no
    hand-tuned heuristics);
  - short scripts;
  - say plainly what was verified and what wasn't.
- **Any model-setting change gets a new `paths.scratch_root`**, so runs never mix,
  because validated runs are skipped on resubmission. **Any design change gets a new,
  versioned matrix file**; never edit an existing one.
- **Git**: this folder is a git repository (branch `main`; private GitHub repo
  `jgoell-bk/camp-fire-uq`, set up 2026-11-08). `.gitignore` keeps data and run output
  out, and `DATA.md` says where it lives. Commit or push only when Jacob asks. Never
  add rasters, `.mat` or `.npz` files, logs or files over ~50 MB. `gh` loads with
  `ml devel gh/2.88.1`.
- **QGIS isn't on Sherlock.** Jacob downloads bundles (VS Code Explorer → Download, or
  `scp` via `dtn.sherlock.stanford.edu`) and views them on his laptop.

## 9. References
- Brewer, M. J., Clements, C. B. (2020). The 2018 Camp Fire: meteorological analysis
  using in situ observations and numerical simulations. *Atmosphere* 11(1), 47 (sec. 3.3:
  Jarbo Gap data quality).
- Saltelli (2002), *Comput. Phys. Commun.* 145, 280-297; Saltelli et al. (2010),
  *Comput. Phys. Commun.* 181, 259-270; Sobol' (2001), *Math. Comput. Simul.* 55,
  271-280.
- Station metadata: NIFC Public View Interagency RAWS. Town limits: Census TIGER/Line
  2018 places (CA).
