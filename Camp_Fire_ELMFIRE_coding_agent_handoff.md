# Camp Fire ELMFIRE sensitivity study: coding-agent handoff

Prepared 19 October 2026. This document specifies the workflow to implement; it does not mean the raster-generation or cluster scripts already exist or that simulations have been submitted.

## 1. Objective and scope

Build a reproducible pipeline that prepares hourly wind and dead-fuel-moisture rasters, creates isolated ELMFIRE cases, submits and monitors cluster jobs, retries failures, and collects outputs in the ordering required for Sobol analysis.

Use the existing runnable Camp Fire case as the physical baseline. Start by reading the repository's `AGENTS.md`, case documentation, ELMFIRE version, executable/build information, and any existing submission scripts. Prefer configuration and preprocessing changes. Do not modify ELMFIRE source unless a demonstrated requirement cannot be handled through the installed version's supported inputs; explain any such requirement before changing the experiment.

The primary deliverable is a working preparation/submission pipeline. Preserve outputs suitable for later sensitivity analysis at cells surrounding Paradise's wildland–urban interface (WUI). A WUI geometry and the complete baseline case have not been supplied with this handoff.

## 2. Fixed experimental design

Run three **separate wind-station scenarios**, using the same sampled factor vectors in each scenario:

| Wind scenario | Station ID | Run prefix | Moisture source |
| --- | --- | --- | --- |
| JARBO GAP | `041214` | `JG` | JARBO GAP |
| OPENSHAW | `041215` | `OP` | JARBO GAP |
| COHASSET | `041211` | `CO` | JARBO GAP |

Each scenario uses its station's observed hourly wind speed and direction. Each hourly wind raster is spatially homogeneous. The observed baseline changes with time; the perturbation coefficients remain constant throughout a simulation. Do not introduce WindNinja, blend stations, resample perturbations independently each hour, or replace the event time series with historical seasonal distributions.

All three scenarios use the same hourly Jarbo Gap series for 1-, 10-, and 100-hour dead fuel moisture. Hold ignition mask/location/time, terrain, fuel models, live moisture, suppression, and other physical settings fixed. Preserve the baseline's spotting/crown-fire settings. Distinguish these fixed physical settings from required changes to file paths, weather timing, output names, and reproducibility settings.

### Independent uncertain factors

| CSV column | Distribution | Meaning | Nominal value |
| --- | --- | --- | --- |
| `a_U` | Uniform(0.8, 1.2) | Multiplier on hourly wind speed | 1 |
| `b_theta_deg` | Uniform(-20, 20) | Offset to wind-from azimuth, degrees | 0 |
| `b_1_pp` | Uniform(-1, 1) | Offset to 1-hour moisture, percentage points | 0 |
| `b_10_pp` | Uniform(-1, 1) | Offset to 10-hour moisture, percentage points | 0 |
| `b_100_pp` | Uniform(-2, 2) | Offset to 100-hour moisture, percentage points | 0 |

These are **prescribed pilot sensitivity ranges**, not empirically calibrated measurement-error distributions or confidence intervals. Assume the five factors are mutually independent, as requested. This assumption concerns the uncertain coefficients; it does not erase the temporal structure in the observed weather/moisture series.

For run `r`, station `s`, time `t`, and every valid weather-grid cell:

```text
WS[r,t]   = a_U[r] * WS_observed[s,t]
WD[r,t]   = (WD_observed[s,t] + b_theta_deg[r]) modulo 360
M1[r,t]   = M1_observed[Jarbo Gap,t]   + b_1_pp[r]
M10[r,t]  = M10_observed[Jarbo Gap,t]  + b_10_pp[r]
M100[r,t] = M100_observed[Jarbo Gap,t] + b_100_pp[r]
```

The moisture equations above use **percent**. If ELMFIRE is configured to consume fractions, divide the entire resulting moisture value by 100 before writing. For example, 7% + 1 percentage point becomes 8%, or 0.08 in fraction units; it does not become 7.07%.

## 3. Canonical Sobol run matrix: reuse it unchanged

Companion file: [Camp_Fire_Sobol_N128_run_matrix_4608.csv](Camp_Fire_Sobol_N128_run_matrix_4608.csv).

SHA-256 of the supplied file:

```text
6472d3200e8e03cdf727f74f85d3c3ad95254d881543cf47db1d3683aeadddd6
```

The design has five factors, base sample size `N=128`, and second-order estimation enabled:

```text
D = 5
calc_second_order = True
runs per station = N * (2*D + 2) = 1,536
three stations = 4,608 simulations
```

Three optional nominal controls, one per wind station, are separate from this CSV. Including them gives **4,611 runs**, before any retries, stochastic replicates, or local finite-difference experiments. Use IDs such as `JG_BASELINE`, `OP_BASELINE`, and `CO_BASELINE`; do not insert controls into the Sobol response vectors.

The matrix was generated from `scipy.stats.qmc.Sobol(d=10, scramble=True, seed=20261019).random_base2(m=7)` with SciPy 1.17.0. The first five coordinates form `A`, the last five form `B`; coordinates were mapped to the factor bounds above. **This seed is the sampling seed, not an ELMFIRE random seed.**

Each of the 128 base rows expands in this exact interleaved order:

```text
A,
AB_a_U, AB_b_theta_deg, AB_b_1_pp, AB_b_10_pp, AB_b_100_pp,
BA_a_U, BA_b_theta_deg, BA_b_1_pp, BA_b_10_pp, BA_b_100_pp,
B
```

`AB_k` is A with only factor k taken from B. `BA_k` is B with only factor k taken from A. This matches SALib's second-order sample/analysis layout. Do not replace these rows with arbitrary independent random samples or round coefficients to displayed precision.

The CSV has these columns:

```text
run_id,wind_station,wind_station_id,moisture_station,moisture_station_id,
design_row,base_row,sobol_block,sampling_seed,
a_U,b_theta_deg,b_1_pp,b_10_pp,b_100_pp
```

`design_row` is 1–1536 within each station; `base_row` is 1–128. File ordering is station-major: JG, OP, CO. IDs are `JG_0001` through `JG_1536`, then the corresponding OP and CO IDs. Each station has the identical factor vectors in the identical design-row order. Preserve station IDs as strings so leading zeros survive.

Treat the CSV as immutable. Record its checksum in experiment metadata and retain original coefficients and identifiers in every run manifest. If a future design is regenerated or extended, version it explicitly rather than overwriting this experiment.

## 4. Source files and preprocessing

Configure file locations rather than hardcoding this chat workspace. At handoff creation, the matrix and this document were in `outputs/`; source files were in `upload/`.

### Hourly wind

Use this exact supplied basename, including its literal ellipsis:

```text
COHASSET,JARBO GAP,...OPENSHAW__WeatherHourly__11-06-2018__11-11-2018.csv
```

Read `StationName`, `StationId`, `DateTime`, `ObservationType`, `WindSpeed(mph)`, `WindAzimuth(degrees)`, `WSflag`, and `WAflag`. Baseline forcing uses wind speed and wind azimuth, **not gust speed/direction**. Retain original flags and source values in the cleaned table.

Inspection found 360 records: 120 hourly records per station, from `2018-11-07T00:00:00Z` through `2018-11-11T23:00:00Z`. The November 8–10 PST wind records are complete, unique per station/time, and have zero WS/WA flags. Repeat these checks during preparation; do not assume future exports share these properties. A once-per-hour station wind report need not represent an average over the entire preceding hour.

### Hourly moisture

Use:

```text
JarboGap__NfdrsHourly__09-30-2003__09-30-2026.csv
```

Its `(1).csv` counterpart is a duplicate upload; do not ingest it twice. Read timestamps from `observationTime`; select Jarbo Gap `stationId=041214`, `NFDRType=O`, and **`fuelModelType=V`**. Use:

```text
oneHR_TL_FuelMoisture
tenHR_TL_FuelMoisture
hundredHR_TL_FuelMoisture
```

The export has five NFDRS fuel-model rows per timestamp, V/W/X/Y/Z. These are not independent observations. Verify the three moisture fields agree across those rows in the selected 2018 window, then retain V deterministically. Do not select the minimum across models or map this NFDRS label onto ELMFIRE landscape fuel-model codes. If agreement fails, report the discrepancy before silently changing the series.

These are NFDRS weather-driven moisture estimates, not three direct moisture sensors. Treat values as percent for the perturbation calculation and verify their input-unit conversion against ELMFIRE configuration. The filename's claimed 2003 start is not authoritative; the inspected records begin in 2005. Use actual timestamps.

For this large file, reading selected columns in chunks and retaining the requested time window is sufficient. Produce one sorted, unique moisture record per UTC hour and join it to each station's wind series by timestamp.

Do not use the supplied monthly weather summaries, daily-minimum moisture products, monthly moisture medians, graph images, or digitized graph text as the hourly forcing source. The raw hourly exports supersede those approximations.

### Time origin, coverage, and quality

Keep parsed timestamps timezone-aware in UTC. On the dates at issue, local PST is UTC−8. Examples:

| Local PST | UTC |
| --- | --- |
| November 8, 2018 00:00 | November 8, 2018 08:00Z |
| November 8, 2018 06:00 | November 8, 2018 14:00Z |
| November 11, 2018 00:00 | November 11, 2018 08:00Z |

The simulation start/end and existing case ignition time remain settings to verify. Do not assume simulation start equals ignition or silently move the ignition when changing weather timing. A first-day window is a sensible pilot, but an exact window has not been finalized in this handoff.

Brewer and Clements (2020), section 3.3, identify questionable Jarbo Gap observations after approximately 06:00 PST on November 9, including possible station fire effects. A zero export QA flag does not establish that these observations are physically reliable. Preserve this warning in provenance and make any extension beyond that time an explicit configured choice. Do not automatically repair the record or switch stations mid-simulation.

Validate coverage of every timestamp needed by the model's temporal interpolation, including the end boundary. For a 24-hour interval, hourly snapshots at both boundaries can require 25 bands. The number of intervals and number of snapshots are different quantities. Record the UTC timestamp represented by every band; do not infer band timing from file modification time or band count.

## 5. Raster and ELMFIRE integration

### Establish the correct grid and configuration first

The supplied Berkeley `ws.tif` and `wd.tif` each contain 13 Float32 bands on a 19×19 grid with 2,500 m cells. Their CRS is NAD83 / UTM zone 11N (`EPSG:26911`), top-left corner `(82500, 4432500)`, and bounds `(82500, 4385000, 130000, 4432500)`. Their NoData metadata is `-9999`. They are spatially varying baseline rasters; they are not the desired homogeneous forcing.

Keep these originals unchanged. They may provide a weather-grid template **after** checking the complete case, installed ELMFIRE version, domain coverage, and interpolation settings. Do not assume the weather grid must equal the landscape/fuel grid, or that the coarse template is automatically valid. Its 13 bands do not establish a 13-hour interval, units, or data provenance; timestamps and speed units were not encoded in these TIFFs.

Inspect the actual case to establish:

- Which raster grids ELMFIRE accepts, how they cover the simulation domain, and any interpolation-edge requirements.
- Wind-speed units and measurement height, whether height conversion occurs internally, and whether the case already applies wind adjustment factors. Convert units once; avoid applying height/canopy corrections twice.
- Whether wind azimuths use true north or grid north, and how the case handles projection rotation. Source direction is meteorological **from** azimuth; do not substitute flow-to azimuth.
- Dead-moisture percent versus fraction configuration.
- Band selection, weather time step, interpolation cadence, simulation clock, and ignition clock.

Relevant namelist keys may include `WS_FILENAME`, `WD_FILENAME`, `M1_FILENAME`, `M10_FILENAME`, `M100_FILENAME`, `WEATHER_DIRECTORY`, `DT_METEOROLOGY`, `DEAD_MC_IN_PERCENT`, `WS_AT_10M`, `WS_IN_KPH`, `ROTATE_WD`, `GRID_DECLINATION`, `BAND_ONE_HOUR_OF_YEAR`, `SIMULATION_TSTART`, `SIMULATION_TSTOP`, and `T_IGN`. **Verify names, sections, defaults, and semantics against the installed version.** Documentation versions differ; do not paste an assumed namelist into production. Keep any internal ELMFIRE ensemble/perturbation controls from unintentionally adding a second sampling layer to this externally specified design.

### Generate inputs

Use GDAL/rasterio or the project's established geospatial tooling. A reasonable implementation writes one multiband GeoTIFF per variable (`ws.tif`, `wd.tif`, `m1.tif`, `m10.tif`, `m100.tif`) per run, subject to the installed model's actual format requirements.

1. Read the run's coefficients from the immutable matrix and its station's aligned hourly baseline.
2. Apply the five equations once to the time-series arrays, then convert to the units expected by the case.
3. Fill every valid cell in each band's weather-domain mask with that hour's scalar value. Preserve the template mask; do not perturb or fill NoData cells as valid weather.
4. Preserve verified CRS, affine transform, dimensions, and NoData convention. Use Float32 where supported; verify compression/interleave choices are readable by the executable.
5. Write a sidecar mapping band index to UTC time, source station/record, coefficients, units, grid identity, and source hashes. Use temporary files and atomic rename after a complete successful write.
6. Verify the resulting rasters and namelist before launching ELMFIRE. If fixed live-moisture or other required inputs need matching band coverage, extend their constant values without introducing new uncertain parameters.

Reject missing coverage, duplicate source times, nonfinite valid values, negative speed, or nonphysical moisture rather than silently clipping, redrawing coefficients, filling gaps, or dropping design rows. The agreed fire-day moisture ranges should stay positive under these offsets. If a broader time window creates a problem, expose it explicitly and establish a documented policy.

Direction wrap must be correct both when applying offsets and when ELMFIRE interpolates between bands. For example, 359° to 1° should pass through north rather than 180°. Check the installed implementation and a small boundary-crossing example. Preserve original valid direction values during calms unless a specific handling policy is established; low-speed direction changes are not automatically data errors.

### Numerical checkpoints

These observed values at `2018-11-08T14:00:00Z` (06:00 PST) provide useful extraction/unit checks:

| Series | WS (mph) | Wind-from azimuth (degrees) |
| --- | --- | --- |
| JARBO GAP | 18 | 33 |
| OPENSHAW | 1 | 48 |
| COHASSET | 1 | 76 |

At that time, selected Jarbo Gap moisture is M1 = 7.53%, M10 = 7.18%, and M100 = 10.00%. A nominal control should reproduce these values before any required wind unit/height/direction conversion. For a sampled run, independently calculate a few expected values from its CSV row and verify representative valid pixels, band timestamps, uniformity, and unchanged masks after writing.

## 6. Pipeline, configuration, and case layout

Provide a small Python package or scripts using the repository's conventions. The following CLI stages are suggested interfaces, not commands that already exist:

```text
validate-inputs     Check case, matrix, source records, units, and coverage
prepare-timeseries Build clean station/moisture tables and provenance
build-case         Build one isolated case from --run-id
run-one            Stage, validate, execute, and collect one case
submit --dry-run   Render the exact scheduler script and task mapping
submit             Submit the explicitly selected run set
collect-status     Reconcile scheduler state with validated output state
retry-failed       Select original failed/incomplete runs without resampling
assemble-qoi       Restore ordered per-station output arrays
```

Use one configuration file for input paths, executable/environment, grid template, simulation window, output/scratch roots, scheduler/account/queue, resources, concurrency, output selection, and model random-seed policy. Keep unresolved scientific settings explicit; fail clearly if a required value is absent. Do not put credentials in configuration or logs.

Suggested directory layout:

```text
experiment/
  config.yaml
  manifest/Camp_Fire_Sobol_N128_run_matrix_4608.csv
  prepared/                 # cleaned tables, band times, provenance
  shared/                   # immutable baseline inputs
  runs/JG_0001/             # independent case, inputs or references, metadata
  runs/OP_0001/
  runs/CO_0001/
  controls/                 # nominal cases, outside Sobol ordering
  logs/
  results/
```

Each run needs an immutable coefficient record, effective namelist, source/template hashes, executable version or hash, model seed policy, timing metadata, scheduler job/task ID, resource request, start/end times, exit code, and validation outcome. Link or stage immutable shared fuels/terrain as appropriate; never let concurrent runs overwrite a shared namelist, output directory, or scratch filename.

Generate small weather rasters on node-local scratch if that suits the cluster, or precompute them if storage is inexpensive. Baselines and moisture may be cached using content-addressed immutable keys. A cached moisture raster can be shared across the three matched station runs because their moisture inputs and coefficients match. Publish caches atomically and avoid concurrent partial-file reads.

## 7. Cluster submission, retries, and completion

Discover the scheduler from the project/environment; the scheduler, partition, account, and resource limits are not specified here. Support the site's actual SLURM/PBS/other interface rather than inventing settings. Prefer the repository's existing scheduler wrapper where suitable.

Use job arrays or an equivalent task manifest with configurable concurrency. For example, a combined zero-based array has tasks `0..4607`, each mapped to CSV data row `task_id + 1` after the header. If using per-station or one-based arrays, define and verify that mapping explicitly. Parse CSV with a CSV library; filenames contain spaces/commas, and station IDs require string handling. Every task should log its resolved `run_id`, `wind_station`, `design_row`, and coefficients before execution.

Make CPU count, MPI ranks, OpenMP threads, memory, walltime, scratch, and concurrency configurable. Benchmark a small representative pilot to choose requests and estimate total runtime/storage. Do not oversubscribe cores or assume 4,608 simultaneous jobs are permitted. If runs are short, consider packing independent runs per allocation with isolation; otherwise keep one run per array task. Never start another array unintentionally from a worker task.

Separate preparation from expensive execution. While implementing this handoff, provide dry-run previews and a concrete submission command; preparing scripts alone is not an instruction to launch the entire ensemble. Run or submit simulations within the user's actual authorization and site rules. No production submission has occurred as part of this document.

A successful scheduler exit is insufficient evidence of a completed simulation. Validate required outputs, dimensions/georeferencing, elapsed simulation coverage, and sentinel conventions. Write a completion marker atomically only after validation succeeds. Re-running a completed run should skip it unless explicitly forced; retries should retain the original failed attempt's logs and exact factor vector. A forced rerun should archive the prior attempt rather than overwrite its evidence.

Provide status counts for pending/running/validated/failed/incomplete runs. Reconcile failures and rerun the original IDs. **Never delete a failed row from the Sobol sequence, shift later rows, or draw a replacement sample.** Bulk analysis should fail clearly if any required response is missing.

### Model stochasticity

Check whether the existing case has stochastic spotting, seed randomization, wind fluctuations, or other random processes. Fixed parameter values alone do not guarantee deterministic outputs. The CSV sampling seed does not configure these processes.

For this pilot, document a reproducible ELMFIRE seed policy and check a repeated identical run. Do not silently turn off baseline spotting physics to obtain repeatability. With fixed model randomness, sensitivity results are conditional on that choice. If independent stochastic replicates or expectation-based indices are later desired, specify replication and aggregation explicitly; replication multiplies the simulation budget. Do not add untracked random noise to every Sobol response.

## 8. Output collection and later sensitivity analysis

Retain georeferenced arrival-time outputs and a reliable burned/unburned distinction. Confirm the model's arrival-time units, time origin, NoData values, and unburned sentinels. Configure output volume deliberately after a pilot; saving every intermediate diagnostic for thousands of runs may dominate storage. If the WUI geometry is unavailable, retain full-domain final arrival-time fields so a common mask can be applied later.

For each station separately, assemble `Y` in ascending `design_row` order, with 1,536 responses per scalar quantity of interest. Completion order and lexical directory order are not analysis order. Do not concatenate all 4,608 responses into a single five-factor Sobol analysis. Station choice is a scenario comparison, not a sixth sampled factor in this experiment.

For a fixed WUI mask with P cells, the stored response matrix can have shape `(1536, P)` per station. Each column can yield first-order, total-order, and all ten pairwise second-order Sobol indices from the same simulations; no extra ensemble is needed for each cell. Use the five parameter names/bounds in the CSV order and `calc_second_order=True` in the analysis. Include uncertainty intervals and assess convergence before interpreting weak second-order effects; `N=128` is a starting design, not a guarantee of stable estimates.

Do not analyze raw infinite arrival times or discard runs where a cell does not burn. Define a finite quantity of interest before analysis. Options include capped arrival time `min(T, H)` with unburned cells assigned horizon H, or a burn-by-H indicator; these answer different questions. Keep outside-domain/unburnable NoData distinct from valid cells that remain unburned. Record the chosen horizon/policy consistently across runs and stations. A cell with zero response variance has undefined Sobol indices and should be masked, not reported as zero sensitivity. A capped arrival-time metric is censored arrival, and is not an uncensored arrival-time estimate.

Optional local finite differences can use the same factor definitions and spatial outputs, evaluated around each nominal station control. A central-difference design uses two perturbed runs per factor plus the control: 11 runs per station, with controls reusable. Finite-difference step sizes have not been selected; make them configurable and keep them distinct from the global uncertainty bounds. Derivatives can be computed per WUI cell for a finite, consistently defined quantity of interest; thresholding/censoring and burn/no-burn transitions need careful interpretation.

## 9. Verification and acceptance criteria

Complete inexpensive input/raster checks before large submissions, then use a small actual-model pilot:

1. Validate the matrix checksum, 4,608 unique IDs, 1,536 ordered rows per station, bounds, and identical paired coefficients across stations. Confirm the A/AB/BA/B structure is preserved.
2. Validate source uniqueness, QA records, NFDRS model selection, timezone conversion, required end-band coverage, and the numerical checkpoints above.
3. Verify raster metadata, homogeneous valid cells, NoData masks, percent/fraction handling, one sampled run's expected perturbations, and north-crossing direction behavior.
4. Run the three nominal station controls and a few deliberately selected design rows spanning the factors. Confirm ELMFIRE reads the intended bands and preserves ignition alignment. Repeat one identical case to check the seed policy.
5. Exercise task mapping, isolated work directories, successful resume/skip behavior, and retry of an intentionally incomplete pilot case. Confirm collection restores design-row order regardless of task completion order.
6. Report pilot runtime, memory, output size, and estimated ensemble demands. Provide an exact dry-run scheduler script and the command for the selected full submission.

Deliver the scripts/package, environment specification, example configuration, unchanged run matrix, concise usage instructions, input-quality/provenance report, and a pilot/status report stating what was actually run. Include focused tests for the consequential logic above; do not require thousands of model runs to test the submission wrapper.

## 10. Settings still to resolve from the actual project

| Setting | Required action |
| --- | --- |
| Baseline case and executable | Locate authorized project files; identify ELMFIRE version and a runnable reference configuration. |
| Cluster scheduler/resources | Read existing site/project scripts; configure account, partition/queue, limits, and environment. |
| Simulation interval | Confirm start/end UTC, model clock, ignition timing, and treatment of later Jarbo Gap data. |
| Weather grid | Verify whether the supplied coarse grid is compatible with the full case and domain. |
| Units and directions | Verify wind units/height/rotation and moisture percent/fraction flags. |
| Interpolation and band mapping | Verify installed-version semantics, last-band coverage, and angular interpolation. |
| Model randomness | Preserve baseline physics and establish a reproducible, documented seed policy. |
| WUI mask and output metric | Obtain a fixed Paradise WUI geometry and choose finite unburned-cell treatment. |
| Local finite differences | Choose step sizes only if this optional analysis is implemented. |

Proceed with the implementation and verifiable preparation that do not depend on these answers. Resolve project-specific facts from existing case files and documentation wherever possible; do not silently invent scientific or cluster settings.

## References

- [ELMFIRE input reference](https://elmfire.io/input_reference.html): check against the installed version and source, especially units, weather timing, and reproducibility settings.
- [SALib Sobol sampling implementation](https://salib.readthedocs.io/en/stable/_modules/SALib/sample/sobol.html): second-order design structure.
- [SALib Sobol analysis implementation](https://salib.readthedocs.io/en/stable/_modules/SALib/analyze/sobol.html): response ordering and index estimation.
- [SciPy Sobol sequence documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.qmc.Sobol.html): scrambled digital net and power-of-two base sample sizes.
- [Brewer and Clements (2020), The 2018 Camp Fire: Meteorological Analysis Using In Situ Observations and Numerical Simulations](https://doi.org/10.3390/atmos11010047): event context and Jarbo Gap data-quality discussion, especially section 3.3.
