# Station baselines

One deterministic ELMFIRE run per weather set, with nothing corrected or perturbed. It
gives a baseline for each station to compare against the mean of their data. Set up
2026-11-09 at Jacob's request. It's independent of the EKI calibration and the earlier
studies: it has its own code (`baselines` package), config and scratch root.

| Run | Wind | Dead fuel moisture |
|---|---|---|
| `JG` | Jarbo Gap (041214), hourly as recorded | Jarbo Gap NFDRS (V) |
| `CO` | Cohasset (041211), hourly as recorded | Jarbo Gap NFDRS (V) |
| `MEAN` | equal-weight vector mean of JG and CO (the EKI baseline) | Jarbo Gap NFDRS (V) |

Settings are the EKI calibration's:
- t0 06:00 PST 8 Nov 2018; ignition at 06:15 PST (t = 900 s) in the original cell
  (row 447, col 1186);
- 24 h run;
- P_EPS 0.01, spotting on;
- gusts off, SEED 2023, `RANDOMIZE_RANDOM_SEED` off.

`elmfire.data.template` is a byte-identical copy of the EKI one.

## Decisions (Jacob, 2026-11-09)
- **One run per weather set, deterministic.** The first idea was 5 runs per set, to
  measure the spread caused by spotting. With these settings, though, ELMFIRE's spotting
  has no randomness, so 5 runs would be bit-identical:
  - embers are spread deterministically along the downwind line (Eulerian);
  - every burning cell spots (spotting percentage 100);
  - any ember flux ignites vegetation (`LOCAL_IGNITION_TIME` 0, so P = 1);
  - `ELMFIRE_CODEBASE.md` section 9 has the details.

  Getting a spread would mean adding a random process, such as Lagrangian embers, a
  random ignition wait, a spotting percentage below 100, or gusts. Its size would then
  depend on the setting chosen. Jacob dropped the replicates.
- **Ignition at 06:15 PST**, as in the EKI, so `MEAN` has the same inputs as the EKI
  `BASELINE` run.
- **QoI: arrival time in the Paradise wildland border band.** That's the group's
  QoI cells (`../experiment/qgis_setup/paradise_border_mask.tif`, 33,758 cells).

## Checks done when it was set up
- `prepare`: MEAN matches the groupmate's `weather_MEAN.csv` to 4.9e-5 (its rounding).
  The JG and CO winds equal the station columns in `../eki_calibration/prepared/baseline.csv`.
- The built MEAN case has weather rasters byte-identical to the EKI BASELINE run's. Its
  namelist differs only in the four path lines, and the shared inputs (fuels, ELMFIRE
  binary, weather grid) have the same checksums.
- `python -m unittest tests.test_baselines`: 6 tests, about a second. They include
  `run_one` with `/bin/true` standing in for ELMFIRE.

## Results (job 51693355, 2026-11-09)
All 3 validated (6–11 min per fire). MEAN's arrival raster is byte-identical to the EKI
BASELINE run's.

| Run | Band burned by 24 h | Earliest | Median (burned) | Mean, capped at 24 h | Total burned |
|---|---|---|---|---|---|
| JG | 40.5 % | 10.4 h | 18.8 h | 21.9 h | 103,101 ac |
| CO | 54.1 % | 16.6 h | 19.8 h | 22.0 h | 97,958 ac |
| MEAN | 42.3 % | 14.5 h | 18.7 h | 22.0 h | 97,717 ac |

Hours are after 06:00 PST 8 Nov. Per cell, against MEAN:
- **JG** burns nearly the same cells: 13,609 shared, 64 burned only by JG, 654 only by
  MEAN. On shared cells the median difference is −0.05 h (p10/p90 −1.9/+0.6 h).
- **CO** burns 4,053 more cells. On shared cells its median difference is +1.1 h, with a
  wide spread (p10/p90 −4.1/+5.2 h).

## Use
From this folder, with `ml python/3.12.1` and
`source /oak/stanford/groups/gorle/jgoell/venvs/elmfire-vnv/bin/activate`:

| Command | Does |
|---|---|
| `python -m baselines prepare` | Write `prepared/weather.csv` (each run's hourly inputs) and `$SCRATCH/elmfire-class/station_baselines/shared/` |
| `sbatch slurm/baselines.slurm` | Run all three at once (3 CPUs, 6 h). Validated runs are skipped |
| `python -m baselines status` | Validated / failed / not started, per run |
| `python -m baselines qoi` | Border-band summary in `results/border_qoi.csv`; per-cell arrival times in `$SCRATCH/elmfire-class/station_baselines/border_toa.csv` |
| `python -m baselines build MEAN --out DIR` | Build one case without running it |

`qoi` columns:
- `burned_share`: the share of band cells burned by 24 h;
- `earliest_h`, `median_burned_h`, `mean_burned_h`: over burned cells, in hours after
  06:00 PST;
- `mean_capped_h`: unburned cells count as 24 h, like `Tcap` in
  `../experiment/matlab/extract_border_toa.m`;
- `same_inputs_run`: whether MEAN's arrival raster is identical to the EKI BASELINE run's.

`border_toa.csv` gives `row`/`col` as 1-based indices in MATLAB `find(mask)` order,
then one column per run in seconds after t0 (NaN = unburned).

Run folders: `$SCRATCH/elmfire-class/station_baselines/runs/<JG|CO|MEAN>/`. Each holds
the arrival raster, `run.json`, `weather.json`, the weather rasters and the logs, as in
the EKI.
