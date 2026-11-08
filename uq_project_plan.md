# UQ study of ELMFIRE: 2018 Camp Fire

> **Superseded (kept as history).** This plan (2026-10-08) was replaced by
> `Camp_Fire_ELMFIRE_coding_agent_handoff.md` and the decisions recorded in
> `PROJECT.md`, which describes what was actually run. Read `PROJECT.md` first.

Project plan for the UQ class. Companion to `class_project_setup.md`, which
covers getting ELMFIRE and the job scripts set up.

## 1. Goal

ELMFIRE's Camp Fire case already runs a 50-member ensemble, but every member
uses the same inputs. The spread between members only shows how much the fire
varies by chance. This project asks the questions that ensemble can't answer:

1. **How much does the predicted fire change when the inputs are uncertain**
   (wind, fuel moisture, spread-rate multiplier, ember settings), compared with
   how much it changes by chance alone?
2. **Which inputs matter most?**
3. *(Stretch)* **Which input values are consistent with the satellite
   observations?** In other words, a calibration with uncertainty, instead of
   the single hand-tuned values the case uses now.

## 2. What ELMFIRE's built-in ensemble does, and why it isn't enough

How it works:

- Each member is one fire. Before it starts, it picks a random start cell in
  the ignition mask (a ~650 m patch near Pulga). It can also add a random
  "tweak" to inputs such as wind speed, wind direction, fuel moisture and the
  spread-rate multiplier (ADJ), each drawn from a range set in `&MONTE_CARLO`.
- While it burns, the wind gets random gusts every 30 s (up to ±5 % speed,
  ±27° direction). Ember ignition also makes random draws, but with this
  case's settings the outcome is fixed: vegetation always ignites when any
  embers land, and every building fuel model has an ignition delay of 0, so a
  building either always or never ignites depending on the ember amount.
- So the only real randomness in the Camp case is the start cell, the gusts,
  and ±0.0005 of noise on the initial level-set field.
- Afterwards each cell's burn probability is the fraction of members that
  burned it.

In the Camp case every tweak range has the same number at both ends
(`PDF_LOWER_LIMIT = PDF_UPPER_LIMIT`), so every member gets exactly
ADJ +0.2839, wind speed +3 mph, wind direction +45°, WAF −0.0197. These are
calibration values, not sampled uncertainty. Result: burned area at 36 h is
173,000–193,000 acres across the 50 members (standard deviation 2.4 %).

Limits of the built-in machinery for UQ:

| Limit | Consequence |
|---|---|
| A tweak is one number added to the whole raster for the whole run. Per-cell (`PIXEL`) and time-varying (`DYNAMIC`) tweaks pass the input check but aren't implemented. | No spatial or time-varying input uncertainty unless the rasters are edited outside ELMFIRE |
| Only ADJ, CBD, CBH, FMC, M1, M10, M100, MLH, MLW, WAF, WD, WS can be tweaked. CC, CH, FBFM are accepted but ignored. | Fuel-model parameters and ember settings can't be sampled internally |
| Tweaks are drawn by plain random sampling. You can't hand ELMFIRE your own sample points. | No Latin hypercube, Sobol sequences or Morris designs internally |
| One random-number stream per core is shared by tweaks, gusts and embers. | A member can't be rerun on its own; chance and input effects can't be separated by seed control |
| With `RANDOMIZE_RANDOM_SEED = .FALSE.` (the default) every core gets the same random numbers. | In the first round, every member would draw the same tweak values |

The Monte Carlo error of the case's own validation scores is visible: two
50-member runs gave Jaccard 0.49 vs 0.46 and area ratio 1.82 vs 1.92.

## 3. Approach: ELMFIRE as a black box

Write the model as **Y = f(θ, ω)**:

- **θ**: the uncertain inputs studied (section 5).
- **ω**: everything that varies by chance (start cell and gusts in the Camp
  setup; ember draws too if the ember settings are changed).
- **Y**: a quantity of interest (section 6).

A small driver script, not ELMFIRE's own sampling, chooses the θ values:

1. A design file lists one row per sample point θⱼ.
2. For each row the driver writes a copy of `elmfire.data.in`:
   - raster tweaks set with `PDF_LOWER_LIMIT(i) = PDF_UPPER_LIMIT(i) = θⱼ,i`,
     which shifts the whole raster by exactly that amount;
   - other settings (for example `EMBER_GR_PER_MW_VEGE`) written directly
     into their namelist group;
   - `NUM_ENSEMBLE_MEMBERS = R`, giving R chance-only repeats at that point.
3. A Slurm job array runs one ELMFIRE job per row, with R ranks.
4. A collection script turns every run's outputs into one table:
   design point, repeat number, quantities of interest.

This gets around every limit in section 2: any sampling design works, any
namelist setting can be varied, and chance (the R repeats) is kept separate
from inputs (the design rows).

Keep `RANDOMIZE_RANDOM_SEED = .TRUE.` so the R repeats differ from each other.
To remove the start-cell randomness, fix the ignition with
`CSV_FIXED_IGNITION_LOCATIONS` (columns: case, band, x, y, area stop, time
stop).

## 4. Base case

- The Camp Fire validation case in `ELMFIRE_VnV_Suite`, run through the
  driver rather than `run_case.sh`, so the suite's preprocessing and LaTeX
  report don't run thousands of times.
- Build and debug the pipeline on tutorial 1 or 3 first (seconds per fire).
- Simulation length: 36 h matches the observation comparison
  (`SIMULATION_TSTOP = 129600`). Measure the cost at 12 h and 24 h in phase 0;
  a shorter run may be good enough for screening.

## 5. Uncertain inputs (θ)

Starting list. Ranges are **proposals to be justified** (part of the
project); priors are uniform unless there's a reason otherwise.

| Input | How to set it | Camp value | Proposed range | Why it's uncertain |
|---|---|---|---|---|
| Wind speed offset | `RASTER_TO_PERTURB = 'WS'` | +3 mph | 0 to +6 mph | Source is a ~3 km weather model (likely HRRR) with no terrain downscaling; gap winds are underpredicted |
| Wind direction offset | `'WD'` | +45° | +15° to +75° | Same coarse source; the +45° was hand-tuned |
| Dead fuel moisture offset | `'M1'` (and M10, M100 together) | ≈ 0 | −3 to +2 percentage points | Moisture is a constant 6/7/8 % everywhere, not measured. Enter in percent units (−3.0 = −3 points). |
| Spread-rate multiplier offset | `'ADJ'` | +0.28 | 0 to +0.6 (×1.0 to ×1.6) | Calibration knob; trades off against moisture |
| Ember generation rate, vegetation | `EMBER_GR_PER_MW_VEGE` (&SPOTTING) | 33.3 per MW | e.g. 15–50 | Empirical; spotting drives the Camp Fire |
| Spot-fire development delay | `CELL_IGNITION_DELAY` (&SPOTTING) | 1000 s | e.g. 300–2000 s | Empirical; controls how fast spot fires grow |

Phase 2 screens these six and keeps the 3–4 that matter for phases 3–4.
Other candidates: WAF offset, live woody moisture (MLW), canopy base height
(CBH), `P_EPS` (ember landing cutoff).

## 6. Quantities of interest

| Quantity | Where it comes from | Notes |
|---|---|---|
| Burned area at 36 h (acres) | `fire_size_stats.csv`, one row per member | Simplest; start here |
| Burned area at earlier times | Count cells with 0 < arrival time ≤ t in each `time_of_arrival_*.tif` | Gives growth curves |
| Arrival time at Paradise | Arrival-time raster at the town-center cell | **Censored** if the fire never gets there (cell = −1); handle explicitly |
| Agreement with VIIRS | Jaccard of each member's footprint vs the VIIRS concave hull | Needs the time-zero question resolved first (section 9) |
| Burn-probability map | Fraction of repeats burning each cell | Field quantity; stretch |

## 7. Phases

### Phase 0: pipeline and cost test
- Write the driver, Slurm array script and collection script.
- Check: running with the Camp values (+3 mph, +45°, ADJ +0.2839) reproduces
  the baseline, mean area about 186,000 acres, SD about 4,400.
- Time single fires at 12, 24 and 36 h.

### Phase 1: chance-only variability and number of repeats
- At the Camp values, run 100–200 repeats.
- Measure the spread of each quantity of interest by chance alone, and how
  many repeats R are needed for a stable mean (standard error vs R).
- This sets R for later phases and gives the "chance" part of the variance.

### Phase 2: screening (Morris method)
- 6 inputs, about 10 Morris trajectories → 70 design points × 3 repeats
  ≈ 210 fires.
- Rank inputs by mean absolute effect (μ*) and spread (σ). Keep the
  top 3–4.

### Phase 3: surrogate model, uncertainty propagation and sensitivity
- Latin hypercube design over the kept inputs: about 10 points per input
  (40–80 points) × R repeats.
- Fit a Gaussian process to each quantity of interest, with a noise term
  (the simulator is random). Check it on held-out points.
- Using the GP:
  - propagate the input ranges and get the distribution of each quantity;
  - compute Sobol indices (first-order and total);
  - split the total variance into an input part, Var(E[Y | θ]), and a chance
    part, E[Var(Y | θ)].

### Phase 4 (stretch): calibration against VIIRS
- Find which θ values fit the observed footprint, e.g. approximate Bayesian
  computation on Jaccard and area, or a GP-based Bayesian calibration.
- Compare the result with the hand-tuned Camp values.

### Phase 5 (stretch): wind-field uncertainty
- Replace the single blocky wind stack with perturbed stacks (spatial random
  fields) or a WindNinja/CFD-downscaled stack, one ELMFIRE run per
  realization. This is a model-form comparison, not a parameter study.

## 8. Compute budget

Measured from the full Camp run (36 h, 50 members): about 38 core-minutes and
about 1.1 GB per fire, one core per fire. One 50-core node runs 50 fires in
about 45 minutes.

| Phase | Fires | Core-hours at 36 h |
|---|---|---|
| 1 | 200 | ≈ 130 |
| 2 | 210 | ≈ 135 |
| 3 | 400 | ≈ 255 |
| **Total (phases 1–3)** | **810** | **≈ 520** |

A shorter simulation lowers these (to be measured in phase 0). If one array
task finishes in under about 10 minutes, put several design points in each
task; Sherlock asks that jobs run at least that long.

## 9. Issues to settle first

- **Time zero.** `case.yaml` says t = 0 is 2018-11-08 14:00 Pacific. The fire
  was reported around 06:30 Pacific, and the first VIIRS detection in the
  file (19:50 UTC) is before the configured start. Matching the ignition-point
  wind peak (bands 5–6) to the Jarbo Gap station record should settle it.
  Any VIIRS-based quantity depends on this.
- **Wind height.** ELMFIRE treats the speeds as 20-ft winds; if the source was
  10-m wind, they are about 15 % high.
- **`ember_flux` outputs have negative values.** Don't use them as a quantity
  of interest until that's explained.
- **Seeding.** Keep `RANDOMIZE_RANDOM_SEED = .TRUE.`; with `.FALSE.` every
  core repeats the same random numbers.

## 10. Deliverables

- Table of all runs: design point, repeat, inputs, quantities of interest.
- Chance-only variability and number-of-repeats plot (phase 1).
- Morris ranking (phase 2).
- GP surrogate with validation, output distributions, Sobol indices and the
  input-vs-chance variance split (phase 3).
- *(Stretch)* calibrated input ranges compared with the hand-tuned values.

## 11. Open decisions

- Which 1–2 quantities of interest to commit to.
- Final input list and justified ranges.
- Simulation length (after the phase 0 timing test).
- Random or fixed ignition point.
