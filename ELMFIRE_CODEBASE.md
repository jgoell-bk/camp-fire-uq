# ELMFIRE codebase notes (for agents)

These notes describe how ELMFIRE works internally, as needed for the Camp Fire UQ
project (see `PROJECT.md`). They cover ELMFIRE 1.1 at commit `cbf924a` (2026-09-16),
installed at `~/elmfire`, a copy of `$GROUP_HOME/codes/elmfire`.

Everything marked *verified* was read in the source, or confirmed by a run, during
this project (Oct 2026). Line numbers refer to that commit. Online documentation
(https://elmfire.io/input_reference.html, `docs/archive/`) can differ from this
version: when they disagree, `build/source/elmfire_namelists.f90` (the defaults and
input checks) and the routines below are authoritative.

## 1. What ELMFIRE is

ELMFIRE (Eulerian Level set Model of FIRE spread) is a Fortran/MPI landscape fire
model. It does two things:
- it advances a level-set field over a raster with a narrow-band scheme;
- it computes local spread rates from the Rothermel surface model, plus crown-fire,
  spotting and building (WUI) submodels.

One "case" is one fire. An ensemble is many cases distributed over MPI ranks.

## 2. Repository layout (`~/elmfire`)
| Path | Contents |
|---|---|
| `build/source/*.f90` | All model source (about 18.6k lines, table below) |
| `build/linux/` | `make_gnu.sh`, `make_intel.sh`, `Makefile_elmfire`, `Makefile_elmfire_post` |
| `build/linux/bin/` | `elmfire_1.1`, `elmfire_post_1.1`, plus `_debug`, `_gprof`, `_perf`, `_block` variants. The unversioned names are symlinks. |
| `docs/` | `ELMFIRE_Guide.pdf`, `validation_report.pdf`, `archive/user_guide/*.rst` (older docs; its spotting page describes the **superseded** model) |
| `tutorials/`, `examples/`, `validation/`, `verification/` | Example cases (also copied into `Wildfire/tutorials/`) |
| `cloudfire/` | Python tools the tutorials use to fetch fuels and weather |

**Build**: GNU compiler with OpenMPI. Precompiler flags (`make_gnu.sh`, *verified*):
`-D_SMOKE -D_WUI -D_UMDSPOTTING -D_SUPPRESSION`. The ember-generation code is inside
`#ifdef _UMDSPOTTING`, so a build without that flag silently has no modular spotting.

**Runtime**:
- modules `system gcc/12.4.0 openmpi/5.0.5`, plus `physics gdal/3.10.2`, because
  ELMFIRE shells out to `gdal_translate` (`PATH_TO_GDAL = 'auto'` finds it on PATH);
- invocation: `elmfire_1.1 elmfire.data.in`, where the argument is the namelist file.

Running the binary directly works as an MPI singleton (one process, one fire). This
project runs 24 such singletons per node concurrently (*verified*). With
`mpirun -np N`, the ranks of a host share rasters through MPI shared-memory windows
and split the cases between them.

## 3. Source map (`build/source/`)
| File | Lines | Role |
|---|---|---|
| `elmfire.f90` | 983 | `PROGRAM ELMFIRE`: driver (section 4) |
| `elmfire_namelists.f90` | 936 | Namelist groups INPUTS, OUTPUTS, TIME_CONTROL, SIMULATOR, MONTE_CARLO, SPOTTING, WUI, MISCELLANEOUS, CALIBRATION, SMOKE, SUPPRESSION. Defaults are set immediately before each `READ`; derived setup follows (e.g. `/` appended to directories, `IWX_BAND_*`). |
| `elmfire_vars.f90` | 734 | Global variables and types (`RASTER_TYPE`, `NODE`, `DLL` linked lists) |
| `elmfire_init.f90` | 1603 | Input checks (`CHECK_INPUTS`, `CHECK_RASTER_DIMS`), lookup tables, `INIT_RASTERS` (WAF, 1−cos slope, `ISNONBURNABLE`), shared memory, fuel and building fuel tables, `ROTATE_ASP_AND_WD` |
| `elmfire_io.f90` | 3347 | Raster I/O: GeoTIFF → ENVI `.bsq` via `gdal_translate`, weather slices (`UPDATE_WEATHER_SLICE`, `READ_WEATHER_SLICE[_TILED]`, which apply the unit conversions), output dumps (`MAIN_DUMP_ROUTINE`), fire-size stats |
| `elmfire_level_set.f90` | 3017 | `LEVEL_SET_PROPAGATION` (the per-case time loop), `TAG_BAND`, `UNTAG_CELLS`, `RK2_INTEGRATE`, `CFL_AND_FLUX_LIMITER` (Superbee), `UX_AND_UY_ELLIPTICAL`, `COMPUTE_SPREAD_VELOCITIES`, WUI tagging, `EULERIAN_/LAGRANGIAN_SPOTTING_MAIN` |
| `elmfire_spread_rate.f90` | 963 | `ROTHERMEL_SURFACE_SPREAD_RATE` (line 13), `CROWN_CRITICAL_FLIN`, `CFFDRS_SPREAD_RATE`, building models `HAMADA` and `UMD_UCB_BLDG_SPREAD`, `ELLIPSE_UCB`, `CALC_WUI_HEATFLUX`, `HRR_TRANSIENT` |
| `elmfire_spotting.f90` | 1133 | Modular (UMD) spotting: generation, distance distributions, Eulerian/Lagrangian transport, `EMBER_IGNITION`, `CALC_SPOTTING_DURATION` |
| `elmfire_spotting_superseded.f90` | 285 | Old spotting model (`USE_SUPERSEDED_SPOTTING`) |
| `elmfire_subs.f90` | 2097 | Utilities: MPI broadcasts, `PERTURB_RASTERS`, `MAP_FINE_TO_COARSE`, wind/moisture interpolation (`INTERP_WIND_*_BILINEAR`, `INTERP_RASTER_*`), `APPLY_WIND_FLUCTUATIONS`, linked lists, `SHUTDOWN` |
| `elmfire_ignition.f90` | 456 | Case counts and ignition locations (mask, CSV, line ignitions) |
| `elmfire_suppression.f90` | 1779 | Containment and suppression models (not used here) |
| `elmfire_calibration.f90` | 119 | Per-pyrome calibration tables (not used here) |
| `elmfire_post.f90` | 1116 | Separate `elmfire_post` program for binary outputs (not used here) |

## 4. Program flow (`elmfire.f90`, *verified* call sequence)

1. **Setup**:
   - read the namelist path from argument 1, then all namelist groups;
   - read the fuel model table;
   - seed the random-number generator (section 9);
   - lookup tables; parallel-I/O rank assignment;
   - read raster headers (ASP for the fuel grid, WS for the weather grid);
   - `SETUP_SHARED_MEMORY_1`.
2. **Read inputs**: `READ_FUEL_TOPOGRAPHY`, `UPDATE_WEATHER_SLICE(1, …)`, broadcasts,
   `CHECK_INPUTS`.
3. **Cases**: `LINE_IGN_TO_POINT_IGN`, then `DETERMINE_NUM_CASES_TOTAL[_CSV]`.
4. **Grids and shared state**: `MAP_FINE_TO_COARSE(WS, ASP, …)` (fuel→weather index
   tables), `INIT_RASTERS`, `SETUP_SHARED_MEMORY_2`, `DETERMINE_IGNITION_LOCATIONS`.
5. **Case loop**: for each case, draw the Monte Carlo coefficients
   (`RANDOM_NUMBER(COEFFS)`, lines 527/749), then `CALL LEVEL_SET_PROPAGATION`
   (line 767).
6. **Finish**: statistics and output, then `SHUTDOWN`.

## 5. The solver (`LEVEL_SET_PROPAGATION`)
- **Field**: φ (`PHIP`); φ < 0 means burned. Only a narrow band of cells
  (`LIST_TAGGED`, `BANDTHICKNESS` cells wide) around the front is updated.
- **Each step**:
  - interpolate the weather to time T;
  - compute spread rates for the tagged cells;
  - optionally apply wind fluctuations, spotting, WUI and suppression;
  - compute velocities from the elliptical template;
  - pick DT from `TARGET_CFL`, capped by `SIMULATION_DTMAX`;
  - advance with second-order Runge-Kutta (RK2) plus the Superbee limiter;
  - untag cells that are done;
  - write dumps every `DTDUMP`.
- **Ending a case** (*verified*):
  - the normal end is at T ≥ `TSTOP`, where `TSTOP = SIMULATION_TSTOP + (band−1)·DT_METEOROLOGY`
    or comes from the ignition CSV;
  - **early ends** are a stalled front, fewer than 2 tagged nodes, the acres stop, or
    initial-attack containment. On a stall (all |Δφ| < 0.001 in a step), the stall
    test runs only when spotting is off or no embers are in flight;
  - after any finish the loop advances in DT_METEOROLOGY steps, so the **final dump
    still happens at T ≥ TSTOP**.
- **Ignition timing**: a case ignites when `T ≥ SIMULATION_TSTART + (IWX_BAND−1)·DT_METEOROLOGY`.

## 6. Spread physics
- **Surface**: Rothermel (`ROTHERMEL_SURFACE_SPREAD_RATE`) with the table from
  `FUEL_MODEL_FILE` (CSV: ID, name, dynamic flag, loads, SAV, depth, extinction
  moisture, heat).
  - Midflame wind = 20-ft wind × WAF, plus `PERTURB_WAF`. WAF is computed in
    `INIT_RASTERS` (`CALC_WIND_ADJUSTMENT_FACTOR_EVERYWHERE`); its exact inputs weren't
    read.
  - The ADJ raster multiplies spread rate; values below 0 are reset to 1.
  - `SURFACE_SPREAD_MODEL` defaults to `'ROTHERMEL'`; `'CFFDRS'` is the alternative.
- **Fire shape**: an elliptical template (`UX_AND_UY_ELLIPTICAL`): slope and wind φ
  factors, length-to-width ratio, head and back speeds.
- **Crown fire**: enabled when `CROWN_FIRE_MODEL > 0` and surface intensity reaches
  `CROWN_CRITICAL_FLIN` (from CBD, CBH, foliar moisture). Passive and active crown
  rates follow Cruz, per the code comments.
  - *Verified*: the source only ever tests `> 0`, so values 1 and 2 behave the same in
    this version. The Camp case uses 2.
- **Non-burnable cells** (Rothermel; `INIT_RASTERS`, *verified*): fuel 90-100, 256,
  ≤ 0, and > 303 (which is remapped to 256).
  - Urban fuel 91 is therefore non-burnable for *surface* spread, but burns through the
    building model.
- **WUI / buildings** (`&WUI`): `USE_BLDG_SPREAD_MODEL`; `BLDG_SPREAD_MODEL_TYPE`
  selects the model (level_set line 2387):
  - **1 = `HAMADA`** (Hamada model, contributed by Y. Qin, UMD); the Camp case uses
    this;
  - **2 = `UMD_UCB_BLDG_SPREAD`**, with `ELLIPSE_UCB`, `CALC_WUI_HEATFLUX` and
    `TAG_WUI`.

  The building rasters are `BLDG_AREA`, `BLDG_FOOTPRINT_FRAC`, `BLDG_FUEL_MODEL`,
  `BLDG_NONBURNABLE_FRAC` and `BLDG_SEPARATION_DIST`. Building fuel table columns:
  ID, name, T_1MW, T_EARLY, T_FULLDEV, T_DECAY, FUEL_LOAD, HRRPUA_PEAK, FTP_CRIT,
  Q_CRIT, ABSORPTIVITY, HEIGHT, NONBURNABLE_FRAC, P_IGNITION, HARDENING_FACTOR, TAU_IGN.
  In the Camp case, cells with buildings are exactly the fuel-91 cells (76,307).

## 7. Weather handling (*verified*)
- **Files**: `WS`, `WD`, `M1`, `M10`, `M100` multiband rasters in `WEATHER_DIRECTORY`,
  one band per `DT_METEOROLOGY`.
  - Band k corresponds to t = (k−1)·DT. After the last band, the last band is held.
  - Bands are loaded in slices (`WX_BANDS_KEPT_IN_MEM`, default 30).
  - `NUM_METEOROLOGY_TIMES` must cover the run; a case whose band window exceeds the
    file is skipped ("CYCLING").
- **Units**:
  - WS is mph at 20 ft. `WS_AT_10M` multiplies by 0.87; `WS_IN_KPH` multiplies by
    0.621.
  - WD is the direction the wind comes **from**, in degrees.
  - `DEAD_MC_IN_PERCENT` (default `.TRUE.`) multiplies the moisture rasters by 0.01 at
    read time.
  - `ROTATE_WD` with `GRID_DECLINATION` (only applied if |declination| > 0.1°) rotates
    the directions.
- **Grids**: the weather rasters must match **each other** (`CHECK_RASTER_DIMS(WS, …)`),
  **not** the fuel grid.
  - `MAP_FINE_TO_COARSE` maps each fuel cell centre to a weather cell by coordinates,
    clamped to the grid. There is **no CRS check**.
  - Wind is interpolated bilinearly in space (`GET_BILINEAR_INTERPOLATE_COEFFS`, with
    I1 = 1 + NINT(d)) and linearly in time, **on U/V components**, so direction
    wrap-around is handled.
  - **Gotcha**: if I1 and J1 are both clamped at the high edge (a 1×1 grid, or a fuel
    cell beyond the last weather-cell centre in both directions), `BILINEAR_INTERPOLATE`
    computes 0/0. Pad coarse weather grids by at least 2 cells; `sobolfire.weather.check_coverage`
    repeats this arithmetic.
- **Precision**: ELMFIRE stores grid corners as single-precision REALs, so output
  headers differ from the inputs by up to about 0.5 m (e.g. 4421293.0 vs 4421292.763).
  Compare grids with a tolerance.

## 8. Ignition and Monte Carlo
- **Ignitions** (`&MONTE_CARLO`, `&INPUTS`):
  - `RANDOM_IGNITIONS` with `USE_IGNITION_MASK` samples cells from `IGNITION_MASK_FILENAME`;
  - `CSV_FIXED_IGNITION_LOCATIONS = .TRUE.` with `USE_IGNITION_MASK = .FALSE.` reads
    `FUELS_AND_TOPOGRAPHY_DIRECTORY/IGNITIONS_CSV_FILENAME`. Its columns are
    case, starting band, x, y, acres stop, hours stop (TSTOP if > 0), with a header
    line;
  - `X_IGN`/`Y_IGN`/`T_IGN` in `&SIMULATOR` are a further option, not used here;
  - non-burnable ignition cells give an empty case unless
    `ALLOW_NONBURNABLE_PIXEL_IGNITION`.
- **Starting bands**: `METEOROLOGY_BAND_START/STOP/SKIP_INTERVAL`. The default of 1/1/1
  means one starting band.
- **Perturbations**: `NUM_RASTERS_TO_PERTURB`, `RASTER_TO_PERTURB(i)`.
  - Accepted names: ADJ, CBD, CBH, CC, CH, FBFM, FMC, M1, M10, M100, MLH, MLW, WAF, WD,
    WS.
  - `PDF_TYPE` is uniform, Gaussian or lognormal (`PERTURB_RASTERS`). Values are drawn
    per case from the shared random stream; `PDF_LOWER_LIMIT = PDF_UPPER_LIMIT` gives a
    fixed offset. Offsets are **additive**, which is why a speed *multiplier* can't be
    expressed this way.
  - Only `'GLOBAL'`/`'STATIC'` uniform offsets were used and checked in this project.
    Other spatial/temporal modes are accepted by the input check, but their behaviour
    wasn't verified.
- **Ensemble size**: `NUM_ENSEMBLE_MEMBERS` (mask mode) or the number of CSV rows.

## 9. Randomness (*verified*)
- **One stream per process**: Fortran `RANDOM_NUMBER` with a single stream.
  - With `RANDOMIZE_RANDOM_SEED = .FALSE.`, every rank gets the same seed:
    `K(:) = SEED`, with the rank offset commented out (elmfire.f90 around line 197).
    MPI ensembles then repeat random numbers across ranks.
  - With `.TRUE.`, the seed comes from the clock and runs aren't repeatable.
- **Places that draw random numbers**:
  - ignition sampling;
  - Monte Carlo coefficients;
  - `APPLY_WIND_FLUCTUATIONS` (only if `WIND_FLUCTUATIONS`);
  - initial φ noise (only when `NOISE`);
  - spotting: the spotting-percentage test, ember-count rounding, and the ignition
    draw;
  - initial-attack containment;
  - `RANDOMIZE_SIMULATION_TSTOP`.
- **The configuration used here is deterministic**: fixed CSV ignition, no
  fluctuations, a fixed seed, and spotting settings under which the spotting draws
  can't change the outcome. A repeated run is bit-identical.

## 10. Spotting, the modular model (`&SPOTTING`, *verified*)

The Camp case uses `GENERATION_MODEL='PER-MW'`, `SPOTTING_DISTANCE_MODEL='EMPIRICAL'`,
`ACCUMULATION_MODEL='EULERIAN'`, `IGNITION_MODEL='PHYSICAL'` and
`DIFF_WILDLAND_IGNITION=.TRUE.`. With those settings:

1. **Who emits embers**:
   - every burning cell with `FLIN > CRITICAL_SPOTTING_FIRELINE_INTENSITY(fuel)`, if a
     random draw is below `SURFACE_FIRE_SPOTTING_PERCENT(fuel)`/100;
   - that percentage stays at its **default of 100 unless `STOCHASTIC_SPOTTING`**,
     because the `GLOBAL_*` and `_MULT` keys are only used by `SET_SPOTTING_PARAMETERS`;
   - emission lasts `TAU_EMBERGEN` (default 6 s) after the cell ignites, unless
     `USE_PHYSICAL_SPOTTING_DURATION`.
2. **How many**: N/s = FLIN · cellsize · `EMBER_GR_PER_MW_VEGE` (or `_BLDG` for fuel
   91) / 1000. The intensity used is the **surface** fireline intensity.
3. **How far**: a lognormal landing distance from Sardoy (vegetation) or Himoto
   (buildings) as a function of intensity and 10-m wind (20-ft wind / 0.87).
   - μ is capped at 5.
   - The distribution is truncated at the (1 − `P_EPS`) quantile (default 0.01; a code
     comment wrongly says 0.001).
4. **Transport (Eulerian)**: embers are spread along the straight downwind path at the
   wind speed, with no sampling. `USE_CROSSWIND_DISTRIBUTION` only matters when the
   crosswind width exceeds one cell, which is rare at 30 m.
5. **Ignition**: for vegetation, **any ember flux > 0 ignites**. The local ignition
   time is `LOCAL_IGNITION_TIME`, so P = 1 when it is 0. The cell becomes a full fire
   after `CELL_IGNITION_DELAY` (case: 1000 s).
   - Buildings ignite on a deterministic mass-load criterion (de Beer), using
     `HARDENING_FACTOR` and `TAU_IGN` from the building table.

Consequences, confirmed by the 2026-10-19 spotting study in `PROJECT.md`:
- for vegetation, the ember rate has **no effect**;
- reach (`P_EPS`) and delay (`CELL_IGNITION_DELAY`) dominate, and act almost
  interchangeably;
- the whole model is deterministic.

The keys are listed in the `NAMELIST /SPOTTING/` block of `elmfire_namelists.f90`.

## 11. Outputs (`&OUTPUTS`, *verified*)
- **Rasters**: written as `.bil` in `SCRATCH` (or the outputs folder), then converted
  to GeoTIFF in `OUTPUTS_DIRECTORY` if `CONVERT_TO_GEOTIFF`.
  - Name pattern: `<quantity>_<case:7 digits>_<NINT(T):7 digits>.tif`, e.g.
    `time_of_arrival_0000001_0086400.tif`.
  - `DUMP_TIME_OF_ARRIVAL` writes only at the final dump.
- **Time of arrival**: seconds since t = 0. Unburned cells are **−9999**, the NoData
  value. **The ignition cell itself is −1**: in the baseline, the only −1 cell was
  member 1's ignition cell (*verified*), and in this project −1 appears at the fixed
  ignition cell (*verified*). The README in `elmfire-class/` said −1 meant unburned
  until it was corrected on 2026-10-26; older copies may still say so.
- **`fire_size_stats.csv`**: one row per case. Columns include `tstop (h)`,
  `Wall clock time (s)`, `Total fire area (ac)`, `Crown fire area (ac)` and the
  ignition x/y. Values are padded with spaces.
- **Other files**: `dump_times_<case>.csv`; with Monte Carlo variables, the drawn
  coefficients (as in the baseline's `coeffs.csv`).
- **Ember flux**: the baseline's `ember_flux` rasters contain unexplained negative
  values (`USE_EMBER_CONSUMPTION` is off, so ember consumption isn't the cause); don't
  use them as a quantity of interest.

## 12. Practical gotchas (all hit or checked in this project)
1. **Namelist comments**: gfortran finds a group by `&NAME`, so never write `&` inside
   comments in the namelist file. Comments before the first group are safe.
2. **Directory keys** get a `/` appended; give them without a trailing slash.
3. **Scratch per run**: the `SCRATCH` folder receives converted `.bsq` inputs and
   `.bil` outputs under fixed names. Concurrent runs need separate `SCRATCH` folders.
   Node-local storage is best.
4. **Shared fuel folder**: `FUELS_AND_TOPOGRAPHY_DIRECTORY` can be a read-only shared
   folder. The ignition CSV must be in it too.
5. **Stalled-front exit**: with spotting off, it can stop slow, light-wind fires early.
   They still produce a final output; check stdout for `STOPPED:`.
6. **The validation suite's time zero is wrong**: its `case.yaml` labels t = 0 for the
   Camp case as 14:00 Pacific. The evidence points to 14:00 UTC (see `PROJECT.md`).
7. **Undocumented baseline offsets**: the Camp case's Monte Carlo offsets
   (WS +3 mph, WD +45°, ADJ +0.2839, WAF −0.0197, M offsets) are retained from an
   archive. There is no record of their derivation.
