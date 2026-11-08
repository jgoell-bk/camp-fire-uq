# Camp Fire EKI calibration: time-varying weather corrections

This study estimates a few coefficients that correct the hourly JG+CO mean weather over
the day. It fits them by ensemble Kalman inversion (EKI) against the NIST fire
observations. The plan is `../ParameterEstimationFiles/camp_fire_eki_calibration_plan.md`.

This folder holds:
- the observation table, linked to the NIST records;
- the baseline weather;
- the initial diagnostic ensemble, and the code to run it and sample arrival times at
  the observations;
- the frozen holdout split;
- the one-sided presence constraints (deadlines);
- the pilot (I00) diagnostics;
- the regularized EKI update with step control, and a verification step before each
  promotion.

Agents: read `/etc/agents/AGENTS.md` at the start of every Sherlock session (and after
any context summary), and `../PROJECT.md` for how to work with Jacob.

## Status (2026-11-08)

- **I00 done**: job 51576874, 33/33 runs validated, 30 min wall time. It is extracted
  and piloted; see `results/`.
- **First controlled update proposed**, not yet promoted.
  - `proposals/I01.csv`: α = 32, mean direction correction +1/0/+4° → +15/+15/+18°.
  - Check batch CHK01 (3 runs at the ensemble means), ready for Jacob:
    ```
    cd /oak/stanford/groups/gorle/jgoell/Wildfire/eki_calibration
    sbatch --array=0-0%2 slurm/eki_check.slurm batches/CHK01.txt
    ```
  - Then `python -m eki verify CHK01`, and only if accepted, `python -m eki promote I01`.
- The arrival-time (lateness) update stays unused: of the 1,199 fit constraints, only
  21 are reached by every member. The fit uses the distance response instead.
- 44 unit tests pass.

## Decisions

All made by Jacob on 2026-11-08 unless noted.

| Item | Setting | Reason |
|---|---|---|
| Ignition | the original cell (39.810283, −121.437222; row 447, col 1186, fuel 182), at **06:15 PST** (`SIMULATION_TSTART` 900 s) | 06:15 PST is the observations' time origin, so model time = elapsed + 900 s exactly. t0 stays 06:00 PST so band 1 is the 06:00 station hour. |
| Baseline wind | equal-weight vector mean of Jarbo Gap and Cohasset (MEAN; IDW not used) | Plan; Openshaw excluded over data concerns. Reproduces the groupmate's `weather_MEAN.csv` within its rounding (5e-5). |
| Moisture | Jarbo Gap NFDRS 1/10/100-h, as in every earlier study | |
| P_EPS | fixed at **0.01** | Nominal of all earlier controls and the local study. Results are conditional on it. |
| Model | otherwise identical to `../local_sensitivity` (spotting on, gusts off, fixed seed: deterministic) | Namelist body is byte-identical to the local study's template. |
| Knots | 06:00, 12:00, 18:00 PST (bands 1, 7, 13); corrections linear between, **constant after 18:00** | 87 % of reports are before 18:00. |
| Moisture mapping | one m(t), the same percentage-point offset on 1-, 10- and 100-h | Simple to state. |
| Moisture admissibility | every written value ≥ 1 %; below is **refused, never clipped** | ELMFIRE silently raises dead moisture below 1 % (`elmfire.f90` line 550). |
| Initial ensemble | 32 members; speed multiplier U(0.5, 3) → a = ln(mult); direction U(−60°, 60°); moisture U(−2, 2) pp; knots correlated C_rs = 0.7^\|r−s\| (normal scores); families independent | Exploratory initialization ranges, separate from regularization and update constraints. |
| Cell rule | direct / 60 m proxy / excluded; 90 m and 3×3 as sensitivities | See "Observations". |
| Observation source | NIST TN 2135 Appendix F workbook, matched to the GeoJSON | See "NIST records". |
| Review points | OBS0073 proxy (30 m), OBS2440 proxy (44 m, kept with OBS2434 for splitting), OBS2482 excluded | Static decisions; they don't depend on any simulation. |
| Holdout | provisional 1 km spatial blocks, ~20 % held out, seed 20261108, frozen as v1; **kept as is** (not reshuffled for its composition); scored only after settings are fixed | See "Holdout". |
| Non-arrival | recorded per run and observation, never dropped per member, zeroed, moved, set to the end time or fed as NaN; target restricted to the run window | Common in I00; the fit now uses the distance response, which is defined whether or not a cell is reached. |
| Constraints | one-sided presence deadlines: T_i ≤ u_i, with u_i = window end | A report bounds arrival from above; it is not an exact arrival. |
| Fitted response | distance from the fixed cell to the **cumulative burned area at the deadline**, 0 if burned by then; L = 1,000 m (provisional geometric scale); lateness with s = 1 h only as a diagnostic, never a second fitting term | Defined for every member whether or not the cell is reached. |
| Weights | w_i = 1/(K n_b), equal total weight per occupied fit block (K = 161) | Area balance, so dense Paradise blocks don't dominate by count; not an independence claim. |
| Regularization | baseline-centred; scales ln 2 / 30° / 1 pp; R_amp = ½ Σ_g time average of q_g² over the run (exact quadrature); R_smooth = ½ Σ (6 h/Δt)(Δq)² | Penalty scales, not initialization ranges or hard limits. |
| Step control | α doubles from 1 until the mean move ≤ 0.5 sd and every member's move ≤ 1 sd in every coefficient, and every member is admissible (by damping, no projection) | Numerical limits; they do not promise a decrease. |
| Verification | actual runs at the current and proposed ensemble means, plus a repeat; promote only if Φ is no worse beyond the repeat difference | The ensemble-average response is not the response at the mean. |
| One constraint per cell | among comparable reports in a fixed cell, the **smallest deadline** (ties: first in the file); a reason is recorded for every report | Metadata only, never model output. |

## Parameterization (`eki/corrections.py`)

theta = (a_k1, a_k2, a_k3, d_k1, d_k2, d_k3, m_k1, m_k2, m_k3), the corrections at the
three knots. At each hourly band:

    ws = ws0 · exp(a(t)),    wd = (wd0 + d(t)) mod 360,    m_k = m_k0 + m(t)

Zero coefficients reproduce the baseline exactly. Admissibility is checked at every band.

**Initial ensemble** (`ensembles/I00.csv`, details in `I00.json`):
- rank 9 of 9; condition number 4.8 (centred, columns scaled by their standard
  deviations);
- sample lag-1 correlations of the normal scores are 0.57–0.74 against the 0.70 target,
  ordinary sampling noise for 32 members;
- the speed range has mean 1.75, so the ensemble is centred near 1.8×;
- BASELINE (all zeros) is run alongside as a diagnostic.

## NIST records (`eki/nist.py`)

**Source**: `../ParameterEstimationFiles/NIST_TN2135_AppendixF/NIST TN 2135 Appendix F Fire Spread Data.xlsx`
(NIST dataset mds2-2954, downloaded 2026-11-08). The field definitions are in TN 2135
Table 6 and sec. 7.4.1. The report itself is in `../ParameterEstimationFiles/NIST_TN2135_report/`.
- 15 region sheets. A record is identified by **sheet + row**. `Source #` is a source
  reference (TD#, 911 call, photo), not unique: there are 284 values for 2,839 rows.
- `Time` is the observation's start, and `Obs Window (min)` its duration, or the window
  in which it occurred (TN 2135: "end time − start time").
- `Type of Fire`: S spot, V vegetation, R residential structure, C commercial structure,
  O other/undetermined.
- `Residual Fire?` x: fire lingering after the fuels burned, i.e. after the front.
- `SSI` x: documented structure-to-structure ignition.
- The workbook's dates are in 2021 (serial 44508 = 2021-11-08). Only the day (8 or 9
  November) is used.
- Regions overlap, so one record can appear on several sheets. NIST also notes
  near-duplicates from several sources or photos.
- The workbook reads with the standard library: neither the venv nor Sherlock has
  openpyxl.

**Matching** (`observations.match_nist`):
- The GeoJSON is the workbook's located rows in sheet order. Feature i is matched to
  located row i, and the match is accepted only because all 2,836 pairs agree exactly
  on longitude, latitude (6 decimals), day and minute. Any disagreement stops the build:
  there is no nearest match.
- Three workbook rows are not in the GeoJSON. They are TD-082 (Magalia, 9 Nov) with no
  time, and are listed in `observations.json`.
- `nist_match_basis`:

  | Basis | Features |
  |---|---|
  | unique key | 1,611 |
  | key and order, candidates identical | 1,090 |
  | key and order, candidates differ | 135 |

  In the last case the GeoJSON features are identical (same point and minute), so the
  assignment within the group changes no observation.
- `nist_record_repeat_of` names the first copy of a repeated record (556 repeats).

**Event types** come from NIST fields only:
- `residual` (225): the NIST residual flag;
- `fire_present` (2,611): any other report.

NIST does not label arrivals. A report says fire was present at the place during
[start, start + window], so it **bounds arrival from above**; it is not an exact arrival
time. One report has no window.

## Observations (`eki/observations.py` → `prepared/observations.csv`)

- One row per GeoJSON feature, `obs_id` OBS0001 ... OBS2836, with every NIST field
  (`nist_*`).
- Times: `model_t_s` = window start − t0, and `window_end_s`.
- Flags:
  - `in_window`: start inside the run, 2,730 rows;
  - `window_end_after_stop`: 146 rows.
- **Eligible cell**, from the fixed landscape only: a surface-burnable fuel, or fuel 91
  with building area > 0.
- **Frozen mapping**, primary 60 m (`r60_*`) and sensitivity 90 m (`r90_*`):

  | Status | Rule | 60 m | 90 m |
  |---|---|---|---|
  | direct | containing cell eligible | 1,790 | 1,790 |
  | proxy | nonburnable containing cell: nearest eligible cell centre to the reported point, within the radius | 1,023 | 1,044 |
  | excluded | none within the radius, or excluded on review | 23 | 2 |

  For each point the table keeps the original location, both cells and fuels, the
  displacement, the status, the reason, and `reviewed` / `review_decision`.
- **Review decisions** (`observations.review_decisions`):
  - They are checked against the coordinates.
  - OBS0073 (TD-009, "fire burning across Rim Rd", V) gets a 30 m proxy.
  - OBS2440 (TD-092, "fire coming across valley from the east", V) gets a 44 m proxy.
    It shares its proxy cell with OBS2434, which is 8 m away (TD-091, "both sides of
    Neal Rd are already burned north of pull out", residual). They are different sources
    and different statements, so not one event; they stay in one partition. TD-092
    describes fire seen across the valley, not at the point.
  - OBS2482 (TD-008, "fire established in canyon to point where river turns into lake")
    is a landmark on water and is excluded for now.

## Holdout (`eki/holdout.py` → `prepared/holdout_v1.csv`, frozen)

**Rule**:
- Blocks are 1 km squares in EPSG:26910, counted from the fuel grid's top-left corner.
  Each report is located by its frozen sampling cell (else its containing cell).
- These stay in one partition: reports sharing a cell, repeats of a NIST record, and
  `same_event_groups` (OBS2440 + OBS2434).
- The units are shuffled with seed 20261108, and round(0.2 × units) are held out.
- No model output is used.

**Result**: 208 occupied blocks, no link spanning blocks, 42 blocks held out (20 %).

| | fit | holdout |
|---|---|---|
| reports | 2,286 | 550 |
| in run window | 2,210 | 520 |
| after window | 76 | 30 |
| primary direct / proxy / excluded | 1,467 / 803 / 16 | 323 / 220 / 7 |
| residual | 175 | 50 |
| **primary constraints** (proxy) | **1,199** (359) | **231** (88) |
| in window 0–3 / 3–6 / 6–9 / 9–12 / 12–18 / 18–24 h | 309 / 602 / 634 / 356 / 169 / 140 | 103 / 179 / 104 / 74 / 40 / 20 |
| fire type V / R / C / S / O | 574 / 535 / 283 / 370 / 524 | 107 / 123 / 122 / 73 / 125 |
| km from ignition p10 / p50 / p90 | 12.3 / 14.5 / 19.3 | 12.4 / 13.0 / 16.0 |

The holdout is somewhat earlier, closer to the ignition and richer in commercial
structures. This is reported, not corrected. The split tests spatial generalization
within this fire; it is not validation on an independent fire.

## Presence constraints (`eki/target.py` → `prepared/constraints.csv`)

A fire_present report says fire was at its place by the end of its window. The frozen
sampling cell's arrival should therefore satisfy **T_i ≤ u_i**, where u_i = NIST start +
Obs Window, with no allowance (`target.deadline_allowance_s` = 0).

**Tighter bounds** go only in `target.deadline_overrides`, with a reason. There are
none. On 2026-11-08 the 1,013 reports with a window > 0 were searched for clock times
and durations. Five mention one, but none anchors a time for its own location inside
its window. For example, "1 min to 2 min after wind increase" doesn't say when the wind
increased.

**Selection**, from metadata only (never model output), with a rule and a reason for
every report:

| Rule | fit | holdout |
|---|---|---|
| primary: smallest deadline among the comparable reports of its cell | 1,199 | 231 |
| later deadline than the cell's primary | 480 | 142 |
| same deadline, later in the file | 316 | 87 |
| residual (after the front) | 165 | 46 |
| starts outside the run | 76 | 30 |
| window ends after the run | 34 | 6 |
| no sampling cell | 16 | 7 |
| no window | 0 | 1 |

Comparable means eligible, mapped, not residual, and with the window inside the run.

**Residual and objective**:

    r_i = max(0, (T_i − u_i) / s_i)          finite T_i only
    Φ(θ) = Σ r_i² / (2N) + ‖Bθ‖² / 2

This is penalized constraint fitting, not an exact-arrival Gaussian likelihood.
- Meeting a deadline costs nothing.
- Satisfying every deadline does not pin down the weather; B chooses among the
  histories that do.
- s_i is a declared lateness scale. The window duration is not used as a standard
  deviation.

**Non-arrival**: an unreached cell has no residual (NaN, `censored`). It is late by more
than (stop − u_i)/s_i, but that number is never used.

## Observation operator (`eki/extract.py`)

`python -m eki extract batches/I00.txt`:
- refuses while runs are running or not started;
- records failed runs (status 4);
- writes `results/I00_G.npz`: `status_*` and `arrival_s_*` per rule as separate arrays,
  `deadline_s`, the selection and its reasons, the partition and the candidate flags;
- writes `results/I00_reach.csv`, the reach per run for fit direct, fit proxy and fit
  after-window reports.

The holdout is stored but not summarized until the settings are fixed.

| Status | Meaning |
|---|---|
| 0 | reached |
| 1 | eligible cell not reached within the 24 h run, kept as non-arrival |
| 2 | ELMFIRE ignition marker |
| 3 | not mapped (statically nonburnable, excluded or outside) |
| 4 | run failed |

## Pilot I00 (`eki/pilot.py`; fit constraints only, holdout untouched)

`python -m eki pilot batches/I00.txt` → `results/I00_pilot.json`, `I00_pilot_cells.csv`
and `I00_distance.npz`. "Lateness" is T − u: a deadline violation, not an arrival-time
error.

**Reach**: the 1,199 fit constraints are reached by all 32 members in 21 cells, by some
in 1,163, and by none in 15.

| Run | Reached | Met | Late | Median lateness (h, reached) |
|---|---|---|---|---|
| BASELINE | 426 | 27 | 399 | 13.8 |
| members, range | 25–1,161 | 1–214 | 20–1,040 | 3.5–15.1 |
| best (I00_M08) | 1,161 | 214 | 947 | 4.6 |

**Breakdown** (members pooled; reached share / violated share of reached / median
lateness of violated, h):

| Deadline | Cells | Reached | Violated | Lateness (h) |
|---|---|---|---|---|
| 0–3 h | 184 | 0.66 | 0.99 | 11.1 |
| 3–6 h | 322 | 0.56 | 1.00 | 10.3 |
| 6–9 h | 290 | 0.48 | 0.97 | 8.0 |
| 9–12 h | 182 | 0.50 | 0.88 | 5.6 |
| 12–18 h | 127 | 0.48 | 0.62 | 3.7 |
| 18–24 h | 94 | 0.38 | 0.10 | 1.9 |

- Direct vs proxy: 0.53 / 0.87 / 8.4 against 0.50 / 0.94 / 9.3.
- Fire type V is violated least (0.79); S, R, O and C are at 0.92–0.96.
- By NIST region, Magalia (0.31) and Nelson Bar Rd (0.47) are violated least; Pentz,
  Pearson, Bille and Wagstaff are at ≥ 0.98.
- By distance from the ignition there is little pattern (0.84–0.94).

The modelled fire is hours too slow in the first ~12 h. NIST describes spotting up to
6.3 km into the community.

**Parameter-response** (32 members, correlation with standardized coefficients;
|r| > 0.35 ≈ p < 0.05):
- **Direction dominates.** Reach correlates with d_k1/d_k2/d_k3 at +0.48 / +0.65 / +0.92
  (R² 0.91 on all 9 coefficients). Mean distance to the burned area correlates at
  −0.72 / −0.81 / −0.73 (R² 0.94). Per cell, |corr| with each d knot has a median of
  0.6–0.7, and exceeds 0.35 in 97–99 % of cells.
- **Speed barely moves reach** (−0.17 / −0.18 / −0.05). It lowers the violated share
  among reached cells (a_k2 −0.46, a_k3 −0.54).
- **Moisture is weak** everywhere (|r| ≤ 0.3).

**Conditioning**:
- θ ensemble: rank 9, condition number 4.8.
- Arrival on the 21 common cells: effective rank 15.
- Distance on all cells: effective rank 19.
- C_FF + I for the deadline residuals on the common cells (s = 1 h): condition number
  9.3.

**Distance response** (evaluated, not adopted). d_i is the distance from the sampling
cell to the nearest cell burned by u_i, and is 0 if the cell itself burned by then.
- It is defined for every member and cell, with no NaN. It is 0 exactly when the
  deadline is met.
- Where a deadline is violated, it correlates with lateness at 0.68.
- Positive distances for members: p10 / p50 / p90 = 2.9 / 8.9 / 13.7 km.
- A violated deadline finds the fire a median 0.81 km away per hour of lateness (IQR
  0.57–1.06).
- **Spatial scale**: 30 m cells, Euclidean between cell centres, searched to the whole
  domain (45 km).
- **Limitations**:
  - it ignores wind, slope, fuel and barriers between the cell and the fire;
  - it jumps when a spot fire lands nearby;
  - its relation to time depends on the candidate's own spread rate;
  - it needs a declared length scale L to sit in Φ;
  - the cells are spatially correlated (dense in Paradise), but are weighted equally.

## First controlled update (adopted 2026-11-08; config `fit`)

**Response**: d_i is the distance from constraint i's fixed sampling-cell centre to the
nearest cell inside the **cumulative burned area at its deadline** (0 < T ≤ u_i, or the
ignition cell). It is not the final perimeter or the active front, and it is 0 when the
cell itself burned by u_i. The constraint list, deadlines and proxy mapping are fixed
for all members. Standardized lateness (T − u)/1 h is a diagnostic on reached cells
only; it is never fitted.

**Objective** (`eki/objective.py`):

    Phi = 1/2 Σ_i w_i (d_i / L)²  +  R_amp  +  R_smooth,      L = 1,000 m
    w_i = 1 / (K n_b)
    R_amp    = 1/2 Σ_g (1/T) ∫_0^T q_g(t)² dt,                T = 24 h, q_g = c_g / scale_g
    R_smooth = 1/2 Σ_g Σ_adjacent (6 h / Δt_h) (q_next − q_prev)²

- **L** is a provisional geometric weighting scale. It is not an equivalent amount of
  lateness.
- **Weights**: K = 161 occupied fit blocks, n_b = constraints in block b (1 to 60,
  median 2). Each block carries the same total weight and the weights sum to 1. This
  balances area; it does not establish independence. Holdout blocks are not used.
- **Penalty scales**: ln 2, 30° and 1 pp, centred on the baseline (zero correction).
  They are not the initialization ranges or hard limits.
- **Amplitude quadrature**: exact for the piecewise-linear curve, over the whole run,
  including the constant part after the last knot. With knots at 0/6/12 h and T = 24 h,
  the time-average matrix is M = (1/24)[[2,1,0],[1,4,1],[0,1,14]]. So the 18:00 knot,
  which acts for 18 of the 24 h, weighs most.
- **B**: rows L_Mᵀ q_g (with L_M L_Mᵀ = M), then √(6/Δt)·Δq. ½‖Bθ‖² equals
  R_amp + R_smooth exactly (tested against direct quadrature).
- **EKI outputs**: [√w·d/L, Bθ], zero target, Γ = I, with no extra 1/√N.

**Step control** (`update.controlled_step`; numerical limits, not a promise that Φ
decreases): α doubles from 1 until all three hold:
1. every coefficient's mean move is ≤ 0.5 of the current ensemble sd;
2. every member's move in every coefficient is ≤ 1 sd;
3. every member is admissible.

A coefficient with sd < 10⁻⁶ × its penalty scale is listed as collapsed, and the floor
is used in place of its sd. Admissibility is reached by damping, never by projection.

**Proposal from I00** (`proposals/I01.csv` and `.json`; not an ensemble file):

| | |
|---|---|
| α tried | 1 → 2 → 4 → 8 → 16 → 32 |
| largest mean move at α = 16 | 0.72 sd (fails limit 1) |
| largest member move at α = 16 | 1.39 sd (fails limit 2) |
| **chosen α** | **32**: mean ≤ 0.46 sd, member ≤ 0.90 sd |
| inadmissible members / projections / collapsed | 0 / 0 / none |
| mean direction correction (k1 / k2 / k3) | +0.7 / +0.4 / +4.5° → **+14.5 / +14.7 / +17.8°** |
| mean speed multiplier | ×1.66 / 1.66 / 1.71 → ×1.64 / 1.64 / 1.73 |
| mean moisture correction | +0.09 / +0.01 / +0.08 → −0.10 / −0.16 / +0.03 pp |
| I00 members' Φ | mean 48.9 (distance 45.1, amplitude 1.6, smoothness 2.2); min 14.0 |
| BASELINE Φ | 44.7 |

α was recomputed under the new limits and weights. It happens to equal the earlier
dry run's value.

**I00 by deadline window** (BASELINE | members' median [range]). Met share is the share
of deadlines met (d = 0); km is the unweighted mean distance; the weighted term is that
window's part of ½Σw(d/L)².

| Window | Cells | Met | km | Weighted term |
|---|---|---|---|---|
| 0–6 h | 506 | 0.00 \| 0.00 [0.00–0.02] | 10.1 \| 9.8 [6.8–12.3] | 15.1 \| 14.6 |
| 6–12 h | 472 | 0.01 \| 0.01 [0.00–0.11] | 8.7 \| 7.6 [2.8–13.8] | 16.9 \| 13.6 |
| 12–24 h | 221 | 0.10 \| 0.15 [0.01–0.72] | 5.9 \| 4.0 [0.4–12.3] | 12.8 \| 7.6 |

- By fire type, members' met share and mean km: V 0.10 / 6.4, R 0.02 / 7.7, O 0.01 /
  8.2, C 0.00 / 7.4, S 0.01 / 9.6.
- Direct 0.04 / 7.4 against proxy 0.02 / 8.5.
- Median standardized lateness among reached cells: 0–6 h 10.8, 6–12 h 7.0,
  12–24 h −1.2.

**Temporal causality**:
- The basis supports are k1 0–6 h, k2 0–12 h and k3 6–24 h (model hours after 06:00 PST).
- Tested: no coefficient changes weather before its basis turns on, or where its
  basis is zero.
- Standardized regression of each window's mean distance on the nine coefficients:

  | Window | d_k1 | d_k2 | d_k3 | a_k2 |
  |---|---|---|---|---|
  | 0–6 h | −0.62 | −0.34 | −0.15 | −0.37 |
  | 6–12 h | −0.30 | −0.65 | −0.18 | −0.38 |
  | 12–24 h | −0.16 | −0.38 | −0.62 | −0.04 |

  The 0–6 h weight on d_k3 cannot be causal. It reflects correlated initialization
  (C = 0.7^|r−s|) or chance. The earlier reach correlation of +0.92 with d_k3 was a
  late-day effect and does not show control over early progression.

**Early deadlines**:
- No I00 member meets more than 2 % of the 506 deadlines before 12:00 PST.
- The closest (I00_M08: direction about +50°, speed ×1.9–2.0) is still 6.8 km away on
  average.
- No member combines speed > ×2.5 with direction > +30° at k1.
- If early violations persist while the corrections build up toward their scales,
  investigate the model configuration (fixed spotting, P_EPS 0.01; uniform wind) and
  representativeness before widening any range.

The result is a **regularized effective-forcing calibration**, not a unique recovery of
the true weather.

**Verification before I01** (`ensembles/CHK01.csv`, `batches/CHK01.txt`):

| Run | Inputs |
|---|---|
| CHK01_I00MEAN | the current ensemble mean θ̄₀ |
| CHK01_I00MEAN_REP | identical inputs and seed, in its own folder: measures run-to-run variability (earlier repeats were bit-identical) |
| CHK01_I01MEAN | the proposed mean θ̄₁ |

No existing run matches θ̄₀ (BASELINE is θ = 0), so nothing is reused.
- `verify CHK01` compares the actual Φ of the runs at the means. It never uses the
  ensemble average of the member responses.
- The proposal is accepted if Φ(θ̄₁) − Φ(θ̄₀) ≤ the repeat difference.
- `promote I01` copies the proposal to `ensembles/I01.csv` only then, and writes its
  batch.
- If the proposal worsens Φ, increase damping and reconsider. No iteration runs
  unattended.

## EKI update (`eki/update.py`, `eki/objective.py`, `eki/controlled.py`)

    F = [G; B theta],  z = [y; 0],  Gamma = diag(R, I)
    theta_j <- theta_j + C_thetaF (C_FF + alpha Gamma)^-1 (z - F_j)

The update is deterministic (one z for all members) and uses linear solves. Here
G = √w·d/L, y = 0, R = 1 and B = `regularization_matrix`.
- `controlled_step` implements the step control above.
- `write_ensemble` refuses inadmissible members and existing files.
- `controlled.propose`, `verify` and `promote` are the CLI steps.

**Tests**:
- *One-step exactness*: equal to the Kalman analysis with the current ensemble as prior.
  That reference includes the current-ensemble covariance term, so it is not the
  intended objective.
- *Convergence to the intended objective*: the Tikhonov minimizer, and the grid minimum
  of a weighted distance-like Φ with unequal block weights.
- *Objective*: block weights sum to 1 and are equal per block. ½‖Bθ‖² matches direct
  quadrature. Reference values: one scale everywhere gives R_amp = 0.5, and one scale
  over 6 h gives R_smooth = 0.5. The EKI outputs reproduce Φ exactly. Non-finite
  distances are refused.
- *Step control*: both limits hold. A collapsed coefficient is listed and not moved.
  Inadmissible members raise α.
- *Causality*: as described above.
- *Check runs*: their IDs, groups and coefficients never fall back to the baseline's
  zeros.
- *Kept from before*: zero lateness before the deadline and positive after;
  smallest-deadline selection; non-arrival kept as NaN plus a censored flag, and refused
  by the update.

## Open

1. **Jacob**: submit CHK01, then `python -m eki verify CHK01`.
2. If accepted: `python -m eki promote I01` and its sbatch. Then extract I01, recompute
   the distances, and reassess before any further update.
3. If not: increase damping (a larger alpha0 or tighter limits), then `propose` again
   into a new proposal; files are never edited.
4. If the early-window violations persist as the corrections grow: review spotting
   (P_EPS), uniform wind and representativeness before widening ranges.
5. The holdout stays frozen and unscored.

## Commands (from this folder)

```
ml python/3.12.1
source /oak/stanford/groups/gorle/jgoell/venvs/elmfire-vnv/bin/activate
python -m eki prepare          # shared inputs on scratch; prepared/baseline.csv (checked vs MEAN csv)
python -m eki observations     # prepared/observations.csv + .json (GeoJSON + NIST records, mapping)
python -m eki holdout          # prepared/holdout_v1.csv + .json (frozen; refuses to change)
python -m eki constraints      # prepared/constraints.csv: deadlines, primary constraint per cell, reasons
python -m eki ensemble         # ensembles/I00.csv + .json
python -m eki build I00_M05 --out /tmp/eki_preview
python -m eki make-batch I00
python -m eki status batches/I00.txt
python -m eki extract batches/I00.txt   # results/I00_G.npz, I00_reach.csv
python -m eki pilot batches/I00.txt     # results/I00_pilot.json, _pilot_cells.csv, _distance.npz (~3 min)
python -m eki propose I00               # proposals/I01.csv + .json; ensembles/CHK01.csv; batches/CHK01.txt
python -m eki verify CHK01              # results/CHK01_verify.json (after the CHK01 runs)
python -m eki promote I01               # ensembles/I01.csv + batches/I01.txt, only if verify accepted
python -m unittest tests.test_eki
```

## Files

| Path | What |
|---|---|
| `config.yaml` | every setting above |
| `elmfire.data.template` | namelist; body identical to `../local_sensitivity/elmfire.data.template` |
| `eki/` | `forcing`, `corrections`, `ensemble`, `nist`, `observations`, `landscape`, `holdout`, `target`, `case`, `weather`, `run`, `extract`, `pilot`, `objective`, `update`, `controlled`, `cli` |
| `prepared/` | `baseline.csv`, `observations.csv`, `holdout_v1.csv`, `constraints.csv` (+ .json each) |
| `ensembles/` | `I00.csv`, `CHK01.csv` (+ .json) |
| `proposals/` | `I01.csv` (+ .json): the proposed ensemble, before verification |
| `results/` | `I00_G.npz`, `I00_reach.csv`, `I00_pilot.json`, `I00_pilot_cells.csv`, `I00_distance.npz` |
| `slurm/eki_batch.slurm` | serc, 1 node, 24 CPUs, 12 h, 2 G/CPU; logs in `logs/` |
| `slurm/eki_check.slurm` | the same with 3 CPUs, 6 h, for CHK batches |
| `$SCRATCH/elmfire-class/eki_calibration/` | `shared/` and `runs/<BASELINE\|I00\|CHK01\|...>/<run_id>/` |
