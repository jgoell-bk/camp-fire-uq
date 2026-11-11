# Audit of the v1/v2/v3 EKI calibration (2026-11-10)

Scope: the implementation and the saved outputs, before any further calibration runs.
No simulation was run. The holdout was not scored or summarized; its block labels are
read only where the calibration itself reads them (to select the fit partition). v1,
v2 and v3 files are unchanged. New material lives here (`audit/`) and in `../v3_retry_a4/`.

| Script (run from `eki_calibration/`) | Writes `audit/out/` |
|---|---|
| `python audit/early_spread.py` | `early_windows.csv`, `early_by_region.csv`, `early_by_block.csv`, `early_constraints.csv`, `early_members.csv`, `early_summary.json`, `early_map_v3.png`, `early_map_v2.png`, `early_top_contributors.png`, `early_by_region.png`, `early_concentration.png`, `early_tradeoff.png` |
| `python audit/acreage.py` | `acreage.csv`, `acreage.json`, `acreage_map.png` |
| `python audit/history.py` | `history_checks.csv`, `history_spread.csv`, `history_I08_steps.json`, `history_v3_parts.png` |

"Final run" always means the single run at the I08 ensemble-mean parameters
(CHK08_I08MEAN). "Uncorrected" is BASELINE (all corrections zero). Member statistics are
labelled as such. The v2 and v3 fit constraints are identical (1,062 deadlines, 146
blocks).

## Summary

- **No implementation or data error that invalidates the results.** The observation
  operator, timing, time zones, ignition time, CRS and units, weather interpolation,
  perimeter mask and acreage arithmetic all check out against the code, ELMFIRE's
  source and the rasters.
- **Corrections to my earlier reading of v3** (now fixed in `../README.md`):
  - 64 % of v3's added deadline distance comes from the 0–6 h deadlines, not the
    late ones. The late deadlines are where the *met share* fell.
  - v2's 152,264 ac by 24 h are 97,829 ac inside the final perimeter plus 54,435 outside.
    Matching the final total acreage does not mean matching the footprint.
  - P_EPS did not move. The ensemble median was 0.0070 at I00 and 0.0069 at I08, and
    only the spread shrank. That is ensemble concentration, not an inferred interval.
- **The 0–6 h misses are widespread and systematic.** Of the 466 early reports, 97 %
  are missed by more than 1 km and 65 % by more than 5 km. It takes 33 reports (13
  blocks) for half of v3's early term. 99 % of that term comes from reports lying
  straight ahead of the modelled fire along its own spread line. The modelled fire is
  too slow and too narrow, not mis-aimed.
- **Confirmed setup facts that bear on the misses** (none is a coding error):
  - one ignition only, although NIST reports a second one at 06:45;
  - spotting depends on the seed through building ignitions;
  - v1 and v2 reach the western domain edge.
- **The leading explanations are model or parameterization limits, still hypotheses.**
  A spatially uniform wind with three knots, and one global P_EPS, cannot reach Paradise
  early without burning past the southern perimeter edge. The trade-off plot shows this
  clearly.
- **The damped retry is set up.** `../v3_retry_a4/` takes one α = 4 step from the
  unchanged v3 I08 ensemble. It needs **1 new simulation**.

## 1. Observation operator (confirmed correct, unchanged)

A primary fit report i ("fire present here by the end of its window") becomes:

- **Deadline:** u_i = window end in model seconds = (NIST report start − 06:00 PST)
  + 60 × Obs Window (min) + `deadline_allowance_s`, which is 0 here
  (`eki/target.py:31` `deadlines`).
  - One primary report per sampling cell, the one with the smallest deadline
    (`target.py:47`).
  - Residual (after-the-front) reports and windows ending after 24 h are excluded.
- **Sampling cell:** frozen before any run. The report's own cell if it is wildland or
  urban 91; otherwise the nearest such cell within 60 m (a proxy).
- **Residual:** d_i = distance (m) from the sampling-cell centre to the nearest cell of
  the cumulative burned set {0 < T ≤ u_i} ∪ {ignition cell (T = −1)}
  (`eki/pilot.py:37`).
  - d_i = 0 exactly when the cell itself burned by u_i (`pilot.py:42`; T = u_i counts
    as met).
  - Distances are cell-centre to cell-centre on the 30 m grid (EPSG:26910, metres),
    found by expanding squares up to 45 km; it is exact, because a square stops only
    once nothing outside it can be closer.
  - A cell never reached by 24 h (T = −9999) is scored the same way, as its distance to
    the burned set. It is neither dropped nor censored.
- **Tolerance:** none in time beyond the report window, and none in space beyond the
  60 m proxy rule and the cell size. There is no fireline tolerance.
- **Units and weighting:** the EKI output is √w_i · d_i / L, with L = 1000 m.
  - w_i = 1/(K n_b) for the n_b reports in its 1 km block, with K = 146 blocks; the
    weights sum to 1 (`eki/objective.py:48`).
  - Its contribution to Φ is ½ w_i (d_i/L)², dimensionless.
  - With Γ = I, report i's implied error scale is L/√w_i = 12 km (block of 1) to
    93 km (block of 59). This is a block-balanced weighting, not an independent-error
    model.

It is the **burned set, not the instantaneous fireline**: a point the fire passed long
ago has d = 0 at any later deadline. New test
`tests/test_eki.py::TestPilot::test_distance_is_to_the_burned_set_not_the_fireline`
checks five cases on a radially spreading fire:
- a point passed 900 s before its deadline, 90 cells behind the front, scores 0;
- a point reached exactly at its deadline scores 0;
- a point reached 10 s late scores 30 m;
- a point outside the burned set scores its distance to the set's edge;
- a cell that never burned scores its distance to the set.

The suite now has 59 tests, all passing.
`audit/early_spread.py` also recomputes every distance with a nearest-cell variant and
requires bit-equality with the operator. It reproduces the per-window terms in
`CHK08_verify.json`.

## 2. The persistent 0–6 h misses

### Where the deadline objective sits (final runs, fit constraints)

| Window | Reports | Weight | v2 term | v3 term | v3 − v2 | Share of increase | Met v2 → v3 | Mean km v2 → v3 |
|---|---|---|---|---|---|---|---|---|
| 0–6 h | 466 | 0.331 | 4.880 | 6.533 | +1.653 | 64.5 % | 0.011 → 0.009 | 5.08 → 6.15 |
| 6–12 h | 417 | 0.336 | 0.753 | 1.303 | +0.550 | 21.4 % | 0.185 → 0.165 | 1.17 → 1.34 |
| 12–24 h | 179 | 0.333 | 0.220 | 0.581 | +0.361 | 14.1 % | 0.670 → 0.475 | 0.35 → 0.64 |
| All | 1,062 | 1 | 5.853 | 8.418 | +2.565 | | 0.19 → 0.15 | 2.75 → 3.33 |

- The 0–6 h deadlines are 83 % of v2's distance term and 78 % of v3's.
- From v2 to v3, the met share hardly changes in that window (1.1 % → 0.9 %) while its
  distance grows by 1.65.
- In the 12–24 h window the met share drops by 20 points but the distance grows by
  only 0.36. Met share and the weighted distance measure different things.

### By region (NIST sheet) and by block

`early_by_region.png` and `early_by_region.csv` give v3's share of the 0–6 h term by
NIST sheet:

| Sheet | Reports | v3 share | Met share v3 | Mean miss v2 → v3 |
|---|---|---|---|---|
| Pentz Rd | 151 | 26 % | 0 | 5.5 → 6.5 km |
| Concow | 67 | 19 % | 6 % | 2.3 → 2.9 km |
| Skyway | 51 | 13 % | 0 | 8.1 → 9.5 km |
| Clark Rd | 53 | 10 % | 0 | 4.4 → 5.6 km |
| Pearson and Buschmann Rd | 48 | 9 % | 0 | 4.4 → 5.7 km |
| Butte Creek Canyon | 1 | 9 % | 0 | 10.6 → 12.8 km |
| Elliott and Nunneley Rd | 32 | 6 % | 0 | 4.6 → 5.9 km |
| Bille Rd, Wagstaff Rd, Neal Rd, Magalia | 63 together | 9 % | 0 | 5–11 km |

Every region outside Concow has a met share of 0. v3 is worse than v2 in every
region.

`early_by_block.csv` and `early_concentration.png` show how concentrated the term is:
- 62 blocks hold early reports.
- The top 5 blocks carry 30 % of v3's early term; it takes 13 blocks for 50 % and 28
  for 80 %.
- By report: the top 1 carries 8.6 %, the top 10 30 %, and the top 20 41 %; 33 reports
  give 50 % and 124 give 80 %. v2 is similar (30 for 50 %, 112 for 80 %).
- Blocks holding one or two reports carry 30 % of the early term, because the block
  weighting gives every block equal total weight.

**The early objective is not dominated by a few observations: the misses are
widespread.** One report stands out, a weighting effect rather than a data problem:
OBS2307 (Butte Creek Canyon, 10:55, 24 km from the ignition, alone in its block) is
8.6 % of the early term.

Largest contributors (`early_top_contributors.png`, each at its own deadline):

| Report | Due (PST) | Sheet | Block size | km from ignition | Miss v2 / v3 (km) | v3 Φ share | v3 arrival |
|---|---|---|---|---|---|---|---|
| OBS2307 | 10:55 | Butte Creek Canyon | 1 | 23.8 | 10.6 / 12.8 | 0.563 | 20:18 |
| OBS0268 | 08:54 | Pentz Rd (proxy) | 1 | 14.7 | 8.7 / 9.5 | 0.311 | 15:33 |
| OBS0061 | 07:58 | Concow | 1 | 11.5 | 7.2 / 7.9 | 0.212 | 11:42 |
| OBS1093 | 08:30 | Pearson Rd | 3 | 16.6 | 11.2 / 12.0 | 0.163 | 15:38 |
| OBS0200 | 09:00 | Pentz Rd | 1 | 11.2 | 5.5 / 6.5 | 0.143 | 18:58 |
| OBS0058 | 07:54 | Concow | 1 | 9.2 | 5.9 / 6.4 | 0.141 | 20:53 |
| OBS1105 | 09:20 | Pearson Rd | 3 | 16.7 | 10.0 / 10.7 | 0.130 | 16:17 |
| OBS1102 | 09:30 | Pearson Rd | 3 | 16.9 | 9.5 / 10.4 | 0.124 | 16:25 |

**Maps** (`early_map_v3.png`, `early_map_v2.png`): each panel is one deadline hour,
with the area burned by the end of that hour in grey and the reports coloured by their
own-deadline distance.
- By 09:00 PST the reports already cover the Paradise and Pentz Rd area, 10–16 km WSW
  of Pulga. The modelled fire is a narrow plume about 5 km long on the same line.
- For every unmet v3 report, the gap runs from the nearest burned cell to the report.
  99.4 % of the early term comes from gaps within 45° of the ignition → report
  bearing, with a median angle of 7°.
- The modelled front is therefore behind along its own path. This is a lack of speed
  and reach, not of direction.

### Can the calibrated parameters move early spread?

- **Basis:** only k1 (06:00) and k2 (12:00) act before noon (`eki/corrections.py:63`).
  The correction is linear between them, so the morning wind is a fixed blend.
  - v3's 20-ft wind is 9.4 mph uncorrected at 06:00, and corrected to 21.7 mph at
    06:00, 13.8 at 09:00 and 7.8 at 12:00.
- **Bounds:** the only hard limits are dead moisture ≥ 1 % at every band and
  P_EPS ≥ 0.0005 (`corrections.py:53`, `:121`). Members failing them are refused
  (by damping), never clipped.
  - No member was ever inadmissible at any α tried, in any v2 or v3 step, and there
    were no projections. The step limits alone set α.
  - v2's final 06:00 multiplier (×3.51) lies beyond I00's sampling range (0.5–3). That
    range is not a bound.
- **Reach of the ensembles** (`early_members.csv`, `early_tradeoff.png`; 578 v2/v3 runs
  from I00–I08 plus 3 SEN01 runs):
  - no calibration member gets the 0–6 h term below 3.37, and that member burned
    78,880 ac outside the perimeter;
  - 310 of the 466 early reports are never met by any of the 581 runs;
  - only SEN01 (v1's I08 weather at P_EPS 0.003 and 0.001) gets the term down, to 1.73
    and 0.49, but at 70,000–72,000 ac outside.
  - The lower edge of the cloud is a trade-off front: under this parameterization,
    fewer acres outside costs early-deadline fit.
- **Where the over-spread is** (`acreage_map.png`): almost all of v1/v2's outside
  acreage lies SW of the perimeter's southern edge. That is downwind on the same
  Pulga → Paradise line. One spatially uniform wind cannot speed the run into Paradise
  without also pushing fire over that edge.

### Timing, placement, mapping, weather: verified

- **Time origin:**
  - model t = 0 is 14:00Z, 06:00 PST on 8 Nov (PST = UTC−8; DST ended on 4 Nov);
  - arrival rasters hold seconds since t = 0, with an earliest burned cell at 1,190 s
    and the ignition cell set to −1;
  - ELMFIRE starts the fire at T ≥ SIMULATION_TSTART = 900 s, i.e. 06:15
    (`elmfire_level_set.f90:153`).
- **Report times:** they are converted in America/Los_Angeles; the code refuses any
  report not in PST, and requires the observation origin (06:15) to equal t0 +
  ignition_s and the GeoJSON elapsed times to agree (`eki/observations.py:83`).
- **Weather timing:**
  - band k sits at t = (k − 1) × 3600 s; ELMFIRE interpolates linearly between bands
    (`elmfire_level_set.f90:842-846`);
  - there are 25 bands from 06:00 PST on 8 Nov to 06:00 on 9 Nov, with corrections
    exactly at the knots (checked in a run's `weather.json`);
  - **the 06:00 band, before ignition, is only an interpolation anchor.** It weighs
    0.75 at 06:15 and 0 at 07:00. It does not change when the fire starts, which is
    fixed by SIMULATION_TSTART.
- **Ignition:**
  - the configured point (−121.437222, 39.810283) lies 5 m from the centre of its cell
    (row 447, col 1186, fuel 182), and the CSV places the fire there;
  - the 06:25 report under the tower (OBS0001) is 553 m NNE (28°) of that cell, which is
    upwind under the NE wind. No run burns it until the evening. Its Φ share is 0.001;
  - the fire starts 553 m downwind of where the tower report puts it, a head start
    rather than a delay.
- **Second ignition:** OBS0005 at 06:45 reads "second power line ignition determined
  near intersection of Rim Rd and Concow Rd". The model has a single ignition, and it
  misses that site by 2.9 km (v2) and 3.7 km (v3). Its own Φ share is 0.046. Whether a
  second source would speed up Concow or Paradise is untested.
- **Mapping:** the early reports are 301 direct and 165 proxy, all within 60 m on
  eligible fuel. The top contributors are direct except OBS0268. No early misses come
  from mapping: the gaps are kilometres, against at most 60 m of proxy displacement.

### Confirmed vs. hypotheses

| Item | Status |
|---|---|
| Operator, deadlines, time zones, ignition time, CRS and units, weather timing and interpolation | Confirmed correct |
| Single ignition; NIST second ignition at 06:45 not represented | Confirmed setup fact; effect on fit untested |
| Ignition cell 553 m from the 06:25 tower report | Confirmed; negligible in Φ |
| Early misses widespread, radial (the modelled front lags) | Confirmed from the rasters |
| Spatially uniform wind + three knots + one P_EPS can't reach Paradise early without over-burning past the southern edge | Hypothesis, supported by the trade-off front and the location of the outside area |
| Long-range spotting under-represented at P_EPS ≈ 0.007 (P_EPS truncates the spotting-distance distribution at its 1 − P_EPS quantile, `elmfire_spotting.f90:573`) | Hypothesis, supported by SEN01 |
| The baseline wind (JG + CO vector mean) understates the event's strength and terrain channelling | Hypothesis |

## 3. Acreage (verified; `acreage.csv`, `acreage_map.png`)

| Run (single run) | Total by 24 h | Inside final perimeter | Outside | Inside / perimeter | Cells ≤ 900 m from domain edge |
|---|---|---|---|---|---|
| Uncorrected | 97,717 | 57,470 | 40,247 | 37 % | 4,421 (south) |
| v1 final | 149,871 | 97,511 | 52,359 | 64 % | 6,079 (west) |
| v2 final | 152,264 | 97,829 | 54,435 | 64 % | 5,977 (west) |
| v3 final | 109,531 | 90,164 | 19,368 | 59 % | 0 |
| v3 rejected I09 mean | 118,132 | 92,746 | 25,387 | 61 % | 332 (west) |

Your inside figures check out: 57,470, 97,829, and 90,163.5 (90,164) ac for v3.

Checks:
- **Mask:** the polygon's area in the grid CRS is 153,270.7 ac; the mask counts 689,178
  cells × 0.22239 ac = 153,269.6 ac. The attribute GIS_ACRES (153,336) differs by 0.04 %.
  The mask is on the fuel grid and fully inside the domain.
- **Arrival rasters:** they hold only −9999 (NoData: not reached or not burnable) and
  one −1 (the ignition cell), with no NaN.
- **ELMFIRE's own totals** agree within 0.3 ac.

Corrections:
- **Matching total acreage is not matching the footprint.** v2's total is 99 % of the
  final perimeter's area, but 36 % of it is outside the perimeter. Its inside area is
  64 % of the perimeter.
- **No progression evidence.** There is no time-resolved perimeter here, so nothing can
  be said about whether 90,000–98,000 ac inside by 24 h is too much. The NIST reports are
  presence-only. Only the outside area is known to be wrong, as a lower bound.
- **Domain truncation.** v1 and v2 (and the uncorrected run, to the south) reach the
  domain edge, where ELMFIRE stops spread. Their totals and outside areas are lower
  bounds, and the over-spread term cannot see burning beyond the domain. v3 final does
  not touch the edge.
- **Statistics in the logs.** `auto_log.json`'s "outside" and "met / km" columns are
  **member medians**, and its Φ "members mean" is the **mean over members**. The
  verify and compare tables are **single runs at the ensemble-mean parameters**. P_EPS
  at a mean run is 0.01 · exp(mean p), the geometric mean of the members' P_EPS.

## 4. The exact v3 objective

Φ(θ) = D(θ) + O(θ) + R_amp(θ) + R_smooth(θ) (`eki/objective.py:145`):

- **Deadline term:** D = ½ Σ_i w_i (d_i/L)², over the 1,062 primary fit reports.
  - L = 1000 m and w_i = 1/(K n_b), with K = 146 and Σ w_i = 1.
  - The normalization is by block, not by the number of observations: D is half the
    block-balanced mean of (d/L)².
  - The weighted RMS distance is L·√(2D).
- **Perimeter term:** O = ½ (A_out/A_ref)², with A_ref = 10,000 ac.
  - A_out is the area burned by 24 h outside the final perimeter
    (`eki/overspread.py:95`).
  - It is one EKI output with unit noise.
- **Amplitude:** R_amp = ½ Σ_g q_gᵀ M q_g + ½ (p/ln 10)².
  - q_g = θ_g / s_g, with s = ln 2 (log-speed), 30° (direction), 1 pp (moisture) and
    knots at 0, 6, 12 h.
  - M = (1/24)[[2,1,0],[1,4,1],[0,1,14]] is the exact time average of the
    piecewise-linear curve over 0–24 h. The 18:00 knot weighs most, because the curve
    is held constant after it.
- **Smoothness:** R_smooth = ½ Σ_g Σ_r (6 h/Δt_r)(q_g,r+1 − q_g,r)², with Δt = 6 h.
- **EKI form:** F = [√w d/L; A_out/A_ref; Bθ], z = 0, Γ = I (1,079 outputs: 1,062 + 1 +
  16 B rows), with ½|F|² = Φ.
  - Update: θ_j ← θ_j + C_θF (C_FF + αΓ)⁻¹ (0 − F_j), deterministic, with no perturbed
    observations (`eki/update.py:40`).
  - α doubles from `alpha0` until every coefficient's mean move ≤ 0.5 sd, every
    member's move ≤ 1 sd, and every member is admissible (`update.py:70`).

**What A_ref does:**
- Multiplying A_ref by c multiplies O by 1/c². As A_ref → ∞, v3 becomes v2's objective
  plus the p prior.
- At an optimum the marginal costs balance. ∂D/∂RMS = RMS/L² per km and
  ∂O/∂A_out = A_out/A_ref² per acre, so the rate at which acres trade for km is
  (RMS/L²)/(A_out/A_ref²). It depends on where the run is.
- **"10,000 ac equals 1 km of deadline error" is true only for the penalty levels.**
  10,000 ac outside costs 0.5, the same as a weighted-RMS distance of 1 km (every
  deadline missed by exactly 1 km). It is not a fixed exchange rate:
  - at v3's final state (RMS 4.10 km, 19,368 ac), 1 km of RMS is worth about 21,200 ac;
  - at v2's final state (3.42 km, 54,435 ac), it would be worth about 6,300 ac.

**Each component, at the runs at each ensemble mean** (`history_checks.csv`,
`history_v3_parts.png`):

| Mean run | α of the step that made it | Distance | Over-spread | Amplitude | Smoothness | Total | Outside (ac) | Verdict |
|---|---|---|---|---|---|---|---|---|
| I00 | – | 29.991 | 13.512 | 0.308 | 0.016 | 43.827 | 51,985 | – |
| I01 | 64 | 23.901 | 12.226 | 0.310 | 0.015 | 36.452 | 49,448 | accepted |
| I02 | 32 | 18.573 | 6.005 | 0.366 | 0.039 | 24.982 | 34,654 | accepted |
| I03 | 16 | 16.242 | 5.247 | 0.419 | 0.112 | 22.020 | 32,394 | accepted |
| I04 | 8 | 16.531 | 2.996 | 0.411 | 0.220 | 20.157 | 24,477 | accepted |
| I05 | 4 | 14.584 | 2.229 | 0.419 | 0.386 | 17.618 | 21,116 | accepted |
| I06 | 2 | 8.906 | 4.745 | 0.529 | 0.598 | 14.778 | 30,807 | accepted |
| I07 | 4 | 8.717 | 3.007 | 0.562 | 0.746 | 13.032 | 24,523 | accepted |
| **I08 (final)** | 2 | **8.418** | **1.876** | **0.612** | **0.971** | **11.877** | **19,368** | accepted |
| I09 proposal | 2 | 7.937 | 3.222 | 0.671 | 1.079 | 12.910 | 25,387 | **rejected** (+1.033) |

**The accepted final state** is `v3/ensembles/I08.csv` (sha256 `eb9e8c54…6ad750`). It
was promoted from `v3/proposals/I08.csv` after `v3/results/CHK08_verify.json` accepted
it. Its mean-parameter runs CHK08_I08MEAN, CHK09_I08MEAN and CHK09_I08MEAN_REP have the
same arrival raster (sha256 `41b9ba3e…`).

## 5. Stopping, spread and seeds

- **The stopping rule was met; convergence was not demonstrated.** v3 stopped because
  the CHK09 proposal raised Φ at the mean (`auto.py`, `controlled.py:219`).
  - The last accepted step had gained 8.9 %.
  - The ensemble spread is still 0.19–0.58 of the penalty scales.
  - The linearized ensemble prediction failed for that step, predicting −0.85 against
    an actual +1.03.
  - Earlier ratios of actual to predicted change ran 0.64–1.68 (`history_checks.csv`).
  - The over-spread part is not monotone: it rose from 2.23 to 4.75 at I06, in a step
    that was still accepted because the distance part fell by 5.7.
  - v2 stopped on the 2 % rule: CHK08 gained 0.6 % against a predicted 6.5 %, a small
    step, not a test of convergence either.
- **Spread** (sd / penalty scale, I00 → I08):
  - speed 0.66, 0.62, 0.69 → 0.40, 0.32, 0.31;
  - direction 1.19, 1.04, 1.02 → 0.36, 0.19, 0.27;
  - moisture 1.19, 1.07, 1.19 → 0.58, 0.48, 0.47;
  - p 0.52 → 0.11.

  Deterministic EKI shrinks the spread at every step whatever the data say, so this is
  ensemble concentration, not posterior uncertainty.
- **P_EPS:**
  - the member median was 0.0070 (I00), 0.0073–0.0080 (I01–I04) and 0.0069 (I08);
  - at the mean p it was 0.0070 at I00 and 0.0069 at I08;
  - the range narrowed from 0.0011–0.049 to 0.0038–0.012.

  The ensemble concentrated around its starting median. It never shows that the data
  identify P_EPS.
- **Seeds:** every run uses SEED 2023 (`RANDOMIZE_RANDOM_SEED = .FALSE.`, wind
  fluctuations off), so a repeat is bit-identical and the measured variability is 0.
  But ELMFIRE's spotting draws random numbers:
  - one draw per spotting call (`elmfire_spotting.f90:89`);
  - one ember-ignition draw `R0 < P_N` (`:958`). For wildland cells P_N ≈ 1, because
    LOCAL_IGNITION_TIME = 0, so ignition is effectively certain. For building cells
    (fuel 91) it is a real Bernoulli trial.

  The stream is consumed in event order, so any change in θ re-draws building
  ignitions. The bit-identical repeat therefore measures only numerical determinism,
  not spotting-realization variability. The acceptance threshold (0) ignores the
  latter, so a ΔΦ of order 1 could partly be realization noise.
  - Untested hypothesis: the size of that noise. Measuring it needs a few runs of the
    I08 mean at other seeds. That is a separate experiment; the retry keeps the seed
    policy, so it tests damping alone.

## 6. Next experiment: one damped retry (set up, not run)

There's no material implementation or data problem, so no existing result needs
recomputing. The retry is in `../v3_retry_a4/` (its `README.md` has the details):

- **Restart state:** v3's I08 ensemble and its responses (`results/I08_*`), copied with
  checksums checked. A recomputation at α = 2 reproduces v3's rejected I09 to 7e−15.
- **Step:** `fit.step.alpha0` goes from 1 to 4, so α = 4 through the same EKI solve
  (not a halved displacement).
  - Max mean move 0.234 sd and max member move 0.423 sd (α = 2: 0.432 and 0.778); all
    32 admissible.
  - The step points the same way as α = 2 (cosine 0.9998) and is 0.54 times as long;
    per coefficient the ratio is 0.45–0.60.
  - The linearized prediction is −0.47.
- **Unchanged:** the objective (including A_ref), observations, constraints,
  parameterization, penalty scales, seed 2023 and the acceptance rule.
- **Controls:** CHK09_I08MEAN and its repeat are v3's validated runs, reused after a
  fresh build under the retry config matched every input (`REUSE.json`).
- **New simulations: 1** (CHK09_I09MEAN). Submit:

  ```
  cd /oak/stanford/groups/gorle/jgoell/Wildfire/eki_calibration
  sbatch --array=0-0%2 slurm/eki_check.slurm v3_retry_a4/batches/CHK09.txt v3_retry_a4/config.yaml
  ```

  Then run:

  ```
  python -m eki --config v3_retry_a4/config.yaml verify CHK09
  python v3_retry_a4/compare.py
  ```
- **Acceptance (existing rule):** accept iff Φ(α = 4 mean) − 11.8769 ≤ variability, which
  is 0 because the repeat is bit-identical. The comparison reports distance, over-spread,
  amplitude, smoothness and total; the three windows' term, met share and km; and the
  outside acres, for the control, the rejected α = 2 proposal (12.9102) and the retry.
- **Acceptance is not convergence.** It would show that one smaller step lowers Φ at
  the mean for this seed. With the realization noise above unmeasured, a small ΔΦ
  should be read as inconclusive.

## 7. Changes made

- `tests/test_eki.py`: the burned-set vs fireline test above (59 tests pass).
- `slurm/eki_check.slurm`: an optional second argument, the config (default
  `config.yaml`, so existing usage is unchanged). Before this it could only run v1's
  config.
- `eki/run.py` `make_batch`: the printed sbatch line resolves the Slurm script relative
  to the config folder (v2/v3 configs give `../slurm/...`, which was wrong to paste from
  `eki_calibration/`) and adds the config argument for check batches.
- `../README.md`: the v3 results text is corrected (see "Corrections" there).

No calibration code path, objective, operator or output was changed.
