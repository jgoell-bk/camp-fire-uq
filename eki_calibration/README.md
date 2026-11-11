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
- **First controlled update verified and promoted** (2026-11-09).
  - CHK01 (job 51581900, 3/3 validated, 12–14 min each): Φ is 34.76 at the I00 mean and
    **23.10 at the proposed mean**. The repeat was bit-identical, so the variability is 0
    and the change of −11.66 counts as an improvement (`results/CHK01_verify.json`).
  - The gain is in the 6–12 h deadlines (met 1 → 9 %, mean distance 7.5 → 5.3 km) and
    the 12–24 h ones (14 → 19 %, 5.1 → 2.5 km). The **0–6 h deadlines are unchanged**:
    0 % met, 9.2 → 8.7 km.
  - `ensembles/I01.csv` = `proposals/I01.csv`. Batch I01 (32 runs) is ready for Jacob:
    ```
    cd /oak/stanford/groups/gorle/jgoell/Wildfire/eki_calibration
    sbatch --array=0-0%2 slurm/eki_batch.slurm batches/I01.txt
    ```
- **I01 done** (job 51676274, 32/32 validated, 29 min). It has been extracted and
  piloted.
  - Members' Φ: mean 48.9 → **34.1**, median 41.2 → 30.9, worst 111.9 → 65.9, best
    14.0 → 13.1. 26 of 32 members improved (median −8.8). The run at the I01 mean scores
    23.1 (CHK01).
  - By deadline window (members' median): 0–6 h is unchanged (0 % met, 9.8 → 9.3 km);
    6–12 h 1 → 4 % met, 7.6 → 5.9 km; 12–24 h 15 → 20 % met, 4.0 → 2.2 km.
  - Constraints reached by all members: 21 → 51; by none: 15 → 14.
  - Direction spread narrowed from 31–36° to 24–31°. Speed (mean ×1.64–1.73) and
    moisture barely moved.
- **Second update proposed** (2026-11-09, `propose I01`): `proposals/I02.csv`.
  - α = 16. At α = 8 a member moved 1.15 sd; at 16 the largest moves are 0.32 sd (mean)
    and 0.73 sd (member).
  - Mean direction correction +14.5/+14.7/+17.8° → **+22.1/+22.2/+25.6°**. Speed
    ×1.64/1.64/1.73 → ×1.73/1.70/1.81. Moisture −0.10/−0.16/+0.03 → −0.33/−0.34/−0.02 pp.
  - Check batch CHK02, ready for Jacob:
    `sbatch --array=0-0%2 slurm/eki_check.slurm batches/CHK02.txt`.
  - CHK02_I01MEAN matches CHK01_I01MEAN up to 5e-15 (rounding in the mean), so it is
    re-run, not reused. Its checksums will be compared with CHK01's as a cross-job
    reproducibility check.
- **Second update verified and promoted** (2026-11-09).
  - CHK02: 3/3 validated, 13–15 min. Φ is 23.10 at the I01 mean and **17.57 at the
    proposed mean** (−5.53). The repeat was bit-identical.
  - Across jobs, CHK02_I01MEAN is byte-identical to CHK01_I01MEAN in both its weather
    and arrival rasters.
  - By deadline window (met share; mean km; weighted term):
    - 0–6 h: 0.00 → 0.01; 8.7 → 7.7 km; 11.8 → 9.6. The first movement in this window.
    - 6–12 h: 0.09 → 0.06; 5.3 → 4.0 km; 7.3 → 5.0. Fewer met, but closer.
    - 12–24 h: 0.19 → 0.37; 2.5 → 1.4 km; 3.5 → 2.2.
  - `ensembles/I02.csv` = `proposals/I02.csv`.
- **I02 ensemble run and reassessed** (2026-11-09). Job 51686250: 32/32 validated, up to
  27 min per fire. `extract` and `pilot` were run on `batches/I02.txt`.
  - Member Φ: mean 34.1 → **27.0**, median 30.9 → 25.4, worst 65.9 → 54.7. The best
    member hasn't improved (I00 14.0, I01 13.1, I02 13.7). 28 of 32 members improved
    over their I01 parent (median −5.0). The run at the I02 mean scores 17.6 (CHK02).
  - By deadline window (members' median): 0–6 h still 0 % met, 9.3 → 8.3 km;
    6–12 h 4 → 5 % met, 5.9 → 4.7 km; 12–24 h 20 → 35 % met, 2.2 → 1.7 km.
  - Constraints reached by all members: 51 → 120; by none: 14 → 14.
  - Coefficient spread: direction sd is now 28/22/20° (from 31/26/24). Speed (sd about
    0.44, mean ×1.73/1.70/1.81) and moisture (sd about 1.0 pp) are still close to their
    initial spread.
  - Causality: in 0–6 h, the largest non-causal (k3) weight is m_k3 −0.22. The causal
    ones are d_k1 −0.46, a_k1 −0.42, a_k2 −0.40 and m_k1 +0.35. This k3 weight is the
    size that correlated initialization gives, as with d_k3 −0.15 in I00.
- **Third update proposed** (2026-11-09, Jacob: "yes next update"; `propose I02`):
  `proposals/I03.csv`.
  - α = 8. At α = 4 a member moved 1.13 sd; at 8 the largest moves are 0.32 sd (mean)
    and 0.69 sd (member). No inadmissible members, no projections.
  - Unlike the first two steps, this one moves speed and moisture more than direction.
    - Speed ×1.73/1.70/1.81 → **×1.99/1.87/1.97**.
    - Moisture −0.33/−0.34/−0.02 → **−0.60/−0.55/−0.08 pp**.
    - Direction +22.1/+22.2/+25.6° → +24.0/+23.4/+29.3°.
    The largest moves are at the first knot, which acts on the 0–6 h window.
  - Check batch CHK03, ready for Jacob:
    `sbatch --array=0-0%2 slurm/eki_check.slurm batches/CHK03.txt`, then `verify CHK03`.
  - CHK03_I02MEAN equals CHK02_I02MEAN to 1e-16 (rounding in the mean). It is re-run,
    and its checksums will be compared with CHK02's.
- **Third update verified and promoted** (2026-11-09).
  - CHK03 (job 51693833): 3/3 validated, 16–20 min. Φ is 17.57 at the I02 mean (equal
    to CHK02) and **14.70 at the proposed mean** (−2.87). The repeat was bit-identical.
    CHK03_I02MEAN is byte-identical to CHK02_I02MEAN in both its weather and arrival
    rasters.
  - By deadline window (met share; mean km; weighted term):
    - 0–6 h: 0.008 → 0.008; 7.7 → 7.3 km; 9.6 → 8.8;
    - 6–12 h: 0.06 → 0.06; 4.0 → 3.4 km; 5.0 → 3.8;
    - 12–24 h: 0.37 → 0.56; 1.4 → 0.9 km; 2.2 → 1.0.
  - Every fire type and both direct and proxy cells improved. Direct: 12.8 → 10.3;
    proxy: 4.0 → 3.3.
  - The penalty terms grew, as expected for larger corrections: amplitude 0.69 → 0.93,
    smoothness 0.06 → 0.14.
  - `ensembles/I03.csv` = `proposals/I03.csv`. I03 ran as job 51698097 (32/32
    validated); it has not been extracted yet. The `auto` job does that first.
- **The loop now runs in one job** (2026-11-09, Jacob: "come up with a way that all of
  this can be run with one command"). See "One job for the whole loop" below. The
  stopping rules are Jacob's: at most I10, and stop once an update gains less than 2 %.
  A replay of the I02 → I03 step through `auto`, in a sandbox, reproduced
  `ensembles/I03.csv` byte for byte. The CHK03 objective and the I02 extract were also
  identical, apart from timestamps.
- **Automated run finished** (job 51732941, 2026-11-10, 6.4 h; exit 0). It scored I03,
  then ran CHK04–CHK08 and I04–I08. It stopped on the 2 % rule: CHK08 gained 1.4 %.
  **I08 is the final ensemble with the original observations.**
  - Φ at the mean: 14.70 → 13.88 → 11.23 → 9.98 → 8.21 → **8.10**. All checks improved.
  - Members' Φ: mean 21.7 (I03) → 10.0 (I08); best 12.1 → 8.0.
  - By deadline window (members' median, I03 → I08):
    - 0–6 h: 0 → 2 % met, 7.6 → 5.3 km;
    - 6–12 h: 6 → 18 % met, 3.8 → 1.4 km;
    - 12–24 h: 52 → 60 % met, 1.3 → 0.5 km.
  - The full table is in `results/auto_log.json` and the end of
    `logs/camp_eki_auto_51732941.out`.
- **Holdout scored once, after the v1 calibration** (2026-11-10, `python -m eki score
  v1_final BASELINE CHK01_I00MEAN CHK08_I08MEAN --members I08` →
  `results/score_v1_final.json`). It isn't used to choose any setting.
  - The scorer's fit numbers for CHK08_I08MEAN equal `verify`'s exactly.

    | Distance term (met share, mean km) | Uncorrected MEAN | I08 mean | I08 members (median) |
    |---|---|---|---|
    | Fit (1,199 constraints, 161 blocks) | 44.73 (2 %, 8.8) | 6.19 (19 %, 2.7) | 6.99 |
    | Holdout (231 constraints, 37 blocks) | 44.35 (4 %, 8.8) | 10.65 (17 %, 3.2) | 10.98 |

  - Holdout by window at the I08 mean: 0–6 h 3 % met, 5.1 km; 6–12 h 10 %, 1.7 km;
    12–24 h 75 %, 0.3 km.
  - The holdout term's excess over the fit is concentrated: two early blocks with three
    reports contribute 5.96 of 10.65; the other 35 blocks contribute 4.7.
    - B22_14 is OBS2306: 09:11 PST, a mile beyond the Honey Run covered bridge,
      16.6 km from the modelled fire.
    - B22_19 is OBS2188/2189: 08:30–09:00 PST, west Paradise, 12.9 km away.
    - All three are in the cleaned data. Each holdout block weighs 1/37, against 1/161 on
      the fit side.
- **Early-window test SEN01** (2026-11-10). Sensitivity runs, not calibration: the I08
  mean (CHK08_I08MEAN's coefficients, copied verbatim) at P_EPS 0.01 (control), 0.003
  and 0.001.
  - Run group `SEN<nn>`: its CSV may override namelist keys (`case.OVERRIDABLE`, P_EPS);
    calibration files may not.
  - Checked: the control builds identical to CHK08_I08MEAN (weather, and namelist apart
    from paths); the others differ only in the P_EPS line.
  - Jacob submits `sbatch --array=0-0%2 slurm/eki_check.slurm batches/SEN01.txt`, then
    `python -m eki score SEN01 CHK08_I08MEAN SEN01_I08MEAN_PEPS0.01 SEN01_I08MEAN_PEPS0.003 SEN01_I08MEAN_PEPS0.001`.
  - **Results** (scored 2026-11-10, `results/score_SEN01.json`). 3/3 validated. The
    control's arrival raster is byte-identical to CHK08_I08MEAN's. Wall time: 1,853 s
    (control), 2,548 s (0.003), 2,638 s (0.001), all in the same job.

    | P_EPS | fit distance term | fit met | fit 0-6 h km | holdout distance term | holdout met |
    |---|---|---|---|---|---|
    | 0.01 | 6.195 | 0.19 | 5.1 | 10.651 | 0.17 |
    | 0.003 | 2.197 | 0.40 | 2.5 | 4.960 | 0.26 |
    | 0.001 | 0.768 | 0.54 | 1.2 | 1.537 | 0.36 |

  - Overshoot check against the final CAL FIRE perimeter (153,336 ac, contained 26 Nov),
    rasterized on the 30 m grid. Burned area by 12 h / 24 h, with the share outside the
    final perimeter: uncorrected 18,602 (7 %) / 97,717 ac (41 %); I08 mean 46,331 (2 %)
    / 149,871 ac (35 %); P_EPS 0.003 101,733 (37 %) / 168,721 ac (42 %); P_EPS 0.001
    113,829 (39 %) / 170,922 ac (42 %).
  - Among reached fit cells, the share reached before its observation window opened:
    5 % uncorrected, 16 % I08 mean, 36 % (0.003), 51 % (0.001).
  - Reading: the deadline objective is one-sided, so a faster fire always scores
    better. Lower P_EPS meets the early deadlines by over-burning: 38-45k ac outside
    the final perimeter by 12 h, against 1k for the I08 mean. P_EPS is not added to the
    EKI until the objective penalizes over-spread (proposal to Jacob, 2026-11-10).
  - Cosmetic: the namelist comment line still says "P_EPS fixed at 0.01" in overridden
    runs; the P_EPS line itself is correct.
- Jacob's labmate is correcting the arrival-time observations (2026-11-09). When the
  corrected file arrives:
  - the observation table, holdout and constraints must be rebuilt;
  - the calibration restarts from I00, in its own result folders, so the current results
    stay intact.
  - The I00 ELMFIRE runs don't depend on the observations, but the updates do.
- The arrival-time (lateness) update stays unused: of the 1,199 fit constraints, only
  21 are reached by every member. The fit uses the distance response instead.
- 52 unit tests pass.

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

1. **Jacob**: submit `slurm/eki_auto.slurm`. It scores I03, then iterates to I10 at most
   (see "One job for the whole loop"). Afterwards, reassess from `results/auto_log.json`.
2. If a check is rejected: increase damping (a larger alpha0 or tighter limits), then
   `propose` again into a new proposal; files are never edited.
3. The corrected arrival-time observations, when they arrive: restart from I00 in a new
   results location (see Status).
4. If the early-window violations persist as the corrections grow: review spotting
   (P_EPS), uniform wind and representativeness before widening ranges.
5. The holdout stays frozen and unscored.

## v2: the cleaned observations (`v2/`, set up 2026-11-10)

The same calibration rerun on the cleaned observations from Jacob's labmate. Run it
from this folder with `--config v2/config.yaml`; v1 (the folders above) is untouched.

- **Data**: `../ParameterEstimationFiles/arrival_observations_cleaned/camp_fire_arrival_times_2.geojson`.
  - It has 2,513 points. All of them are exact copies of original points (location,
    day, minute, elapsed). 323 were removed, none moved or retimed; the order differs.
  - `elapsed` is now hours rather than a duration string.
- **Table** (`sources.observations_subset_geojson`, `observations.select_subset`):
  - built from the full original file, so the NIST matches and obs_ids are v1's, then
    restricted to the kept points;
  - 5 groups of identical points lost only some copies. Their kept copies take the
    group's first records in workbook order (Jacob). The only deadline this decides is
    OBS2772 (06:00 9 Nov, 5-min window): its window ends after the run, so it isn't a
    constraint;
  - the removed points: 250 direct and 72 proxy; by fire type O 169, V 118, C 19, S 10,
    R 7.
- **Review decisions**: OBS0073 (proxy) and OBS2482 (excluded) carry over; they're the
  same points. OBS2440 was removed, so its decision and the same-event group
  [OBS2440, OBS2434] are dropped.
- **Holdout v2** (`holdout.blocks_from`): every 1 km block keeps its v1 role (Jacob), so
  v1 → v2 differences come from the cleaned data only. That gives 40 of 191 blocks
  held out (21 %); every report has its v1 partition.
  - The first attempt re-drew with the v1 seed and matched only 4 of 38 blocks. It is
    kept, unused, as `v2/prepared/*_fresh_draw_unused.*`.
- **Constraints**: 1,062 primary fit constraints (v1 1,199). Of these, 1,051 are shared
  with v1 and 11 are new; 148 of v1's were dropped. Holdout: 208.
- **Same as v1**:
  - the weather, the template, the fit settings and the stopping rules (I10, 2 %);
  - the I00 ensemble: `v2/ensembles/I00.csv` is a copy of v1's (Jacob).

  Built v2 I00 cases have weather rasters and a namelist identical to v1's (apart from
  paths), and the shared inputs have the same checksums. v2's I00 arrival rasters
  should therefore equal v1's.
- Scratch root: `$SCRATCH/elmfire-class/eki_calibration_v2`.
- To run: `sbatch slurm/eki_auto.slurm v2/config.yaml`. Its log is
  `v2/results/auto_log.json`.

## v3: over-spread term and P_EPS as a coefficient (`v3/`, set up 2026-11-10)

Why: SEN01 showed that lower P_EPS meets the early deadlines by over-burning. The
deadline term never penalizes early arrival, so on its own it rewards any faster fire.
v3 is v2 (cleaned observations, holdout v2, the same constraints) plus two changes.
Jacob confirmed each setting by name on 2026-11-10.

- **Over-spread term** (`fit.overspread`, `eki/overspread.py`):
  - R_over = ½ (A_out / 10,000 ac)², where A_out is the area the run burned outside the
    final CAL FIRE perimeter by the end of the run (24 h);
  - burning 10,000 ac outside costs 0.5, as much as missing every deadline by 1 km;
  - it is one extra EKI output, A_out / 10,000, with unit noise, so ½|F|² is still Φ;
  - the perimeter (`sources.final_perimeter`, 153,336 ac, contained 26 Nov) is
    rasterized on the fuel grid by `python -m eki --config v3/config.yaml perimeter`
    into `v3/prepared/final_perimeter.tif`: 153,270 ac, cell-centre rule;
  - it is a lower bound: by 24 h the real fire filled only part of the final perimeter.
    The uncorrected MEAN weather already burns 40,247 ac outside, which adds 8.10.
- **P_EPS as a coefficient** (`corrections.p_eps`; theta has 10 entries, p last):
  - p = ln(P_EPS / 0.01), one constant per run, written to the namelist as
    P_EPS = 0.01 exp(p), an exact literal (`case.coefficient_overrides`);
  - penalty scale ln 10, so a factor of 10 costs as much as ×2 wind or 30°. It adds
    ½ (p / ln 10)² to R_amp and nothing to R_smooth;
  - floor 0.0005: a member below it is inadmissible, so the step is damped until no
    member is. Nothing is clipped;
  - lower P_EPS means slower runs (SEN01: +40 % at 0.001).
- **I00** (`initial_ensemble.weather_from`, `ensemble.make_extended`):
  - v1's I00 rows (= v2's), copied as text, plus p with P_EPS log-uniform on
    0.001–0.05, independent of the weather, from seed 20261110;
  - the drawn P_EPS values are 0.0011 / 0.0070 / 0.049 (min / median / max); rank 10.
- **Same as v2**: observations, holdout v2, constraints, baseline weather, shared inputs
  (all byte-identical to v2's) and the stopping rules (I10, 2 %).
- `v3/elmfire.data.template` differs from `../elmfire.data.template` only in the P_EPS
  comment line.
- **Checks before submission**:
  - built v3 cases have v1's weather rasters; their namelists differ from v1's only in
    paths, the P_EPS value and that comment;
  - a sandbox replay of extract → pilot → propose → verify → promote on stand-in rasters
    (v1's) ran end to end with the over-spread output;
  - 58 tests pass.
- Scratch root: `$SCRATCH/elmfire-class/eki_calibration_v3`.
- To run: `sbatch slurm/eki_auto.slurm v3/config.yaml`. Its log is
  `v3/results/auto_log.json`.
- Comparing v3 with v2: their Φ differ by definition, so compare deadline distance,
  share met, holdout score and over-spread acres. `score` reports over-spread only
  under v3's config.

### Results: v1, v2, v3 (2026-11-10)

Jobs: v2 = 51800200 (8 h 38 min, stopped at I08: gain 0.6 % < 2 %); v3 = 51842711
(4 h 43 min, stopped at CHK09: the I09 mean worsened Φ by +1.03, so I08 is final).
v1 = 51732941 (stopped at I08, gain 1.4 %).

`python -m tools.compare_finals` (→ `v3/results/compare_finals.json`) measures every
final run (CHK08_I08MEAN) on the same yardsticks: the cleaned fit deadlines (1,062, v2 =
v3) and the final perimeter. The holdout is not touched.

Every row is a single run: BASELINE, or the run at the I08 ensemble-mean parameters.
None is an average over members.

| run | RMS km | met | mean km | 0–6 h met / km | 6–12 h | 12–24 h | burned by 24 h (ac) | inside perimeter (ac) | outside perimeter (ac) |
|---|---|---|---|---|---|---|---|---|---|
| uncorrected | 9.47 | 0.02 | 8.94 | 0.00 / 10.2 | 0.01 / 8.8 | 0.09 / 6.0 | 97,717 | 57,470 | 40,247 |
| v1 final | 3.44 | 0.18 | 2.79 | 0.01 / 5.1 | 0.16 / 1.2 | 0.66 / 0.4 | 149,871 | 97,511 | 52,359 |
| v2 final | 3.42 | 0.19 | 2.75 | 0.01 / 5.1 | 0.18 / 1.2 | 0.67 / 0.3 | 152,264 | 97,829 | 54,435 |
| v3 final | 4.10 | 0.15 | 3.33 | 0.01 / 6.2 | 0.17 / 1.3 | 0.47 / 0.6 | 109,531 | 90,164 | 19,368 |

The weighted distance term by deadline window (from `CHK08_verify.json`; the parts sum
to the run's distance term):

| window | v2 | v3 | v3 − v2 | share of the increase |
|---|---|---|---|---|
| 0–6 h | 4.880 | 6.533 | +1.653 | 64 % |
| 6–12 h | 0.753 | 1.303 | +0.550 | 21 % |
| 12–24 h | 0.220 | 0.581 | +0.361 | 14 % |
| all | 5.853 | 8.418 | +2.565 | |

Final coefficients (I08 mean; speed multiplier, direction, moisture at the knots 06, 12,
18 h):

| | speed × | direction ° | moisture pp | P_EPS |
|---|---|---|---|---|
| v1 | 3.55 / 2.66 / 2.55 | +31.1 / +24.7 / +27.3 | −1.11 / −0.80 / −0.28 | 0.01 (fixed) |
| v2 | 3.51 / 2.62 / 2.49 | +31.3 / +24.7 / +27.5 | −1.12 / −0.79 / −0.27 | 0.01 (fixed) |
| v3 | 2.31 / 1.29 / 0.77 | +25.2 / +31.9 / +22.2 | −0.94 / −0.78 / −0.06 | 0.0069 |

- v1 → v2: the cleaned observations barely change the answer. v2's I00 arrival rasters
  are pixel-identical to v1's (the GeoTIFF files differ only in metadata), and BASELINE
  is identical in all three versions.
- v2 → v3: the over-spread term cuts the area outside the perimeter from 54,435 to
  19,368 ac. The cost is 0.7 km of RMS deadline distance (3.42 → 4.10).
  - 64 % of the added distance term is in the 0–6 h deadlines; the 12–24 h deadlines
    add 14 %.
  - The share of deadlines *met* falls mostly at 12–24 h (0.67 → 0.47). At 0–6 h it
    hardly changes (0.011 → 0.009), but those misses grow from 5.1 to 6.2 km.
  - Met share and the weighted distance term measure different things.
  - v3 does this mainly by slowing the wind (×3.51 → ×2.31 at 06 h, ×2.62 → ×1.29 at
    12 h, ×2.49 → ×0.77 at 18 h).
- Acreage: v2's 152,264 ac by 24 h are 97,829 ac inside the final perimeter (64 % of
  it) and 54,435 ac outside.
  - Matching the final fire's total acreage does not mean matching its footprint.
  - Without a time-resolved perimeter, nothing here says whether 90,000–98,000 ac inside
    by 24 h is too much. Only the outside area is known to be wrong.
  - v1 and v2 reach the western domain edge, where ELMFIRE stops spread, so their total
    and outside areas are lower bounds (`audit/README.md` §3).
- P_EPS: the ensemble median was 0.0070 at I00 and 0.0069 at I08. The range narrowed
  from 0.0011–0.049 to 0.0038–0.0119, and p's sd from 1.20 to 0.24.
  - This is ensemble concentration around the starting median, which deterministic EKI
    produces at every step.
  - It is not a posterior interval, and it does not show that the data identify P_EPS.
- v3 stopped on a worsening step, not on the 2 % rule; the last accepted gain was
  8.9 %. The I09 mean lowered the distance term (8.42 → 7.94) but burned more outside
  (19,368 → 25,387 ac). Meeting the stopping rule is not convergence. A more damped
  step is set up in `v3_retry_a4/`.
- The early deadlines (0–6 h) are missed by 5–6 km in every version.
  - The misses are widespread (33 reports give half of v3's early term) and lie
    straight ahead of the modelled fire along its spread line.
  - No calibration member brings the 0–6 h term below 3.4 (`audit/README.md` §2).
- A_ref = 10,000 ac sets equal penalty *levels*: 10,000 ac outside costs as much as a
  weighted-RMS distance of 1 km. Both terms are quadratic, so the marginal rate of
  trading acres for km depends on the state. At v3's final state it is about
  21,000 ac per km. v3 is one point on the trade-off, not a unique answer.

Corrections (2026-11-10 audit, `audit/README.md`). The text above replaces three earlier
statements:
- "almost all of the cost is in the 12–24 h deadlines" (that was the met share);
- "v1 and v2 burn nearly the whole final perimeter" (total, not footprint);
- "the data pin P_EPS near the default" (concentration).

## One job for the whole loop (`eki/auto.py`, `slurm/eki_auto.slurm`)

```
cd /oak/stanford/groups/gorle/jgoell/Wildfire/eki_calibration
sbatch slurm/eki_auto.slurm
```

One Slurm job (24 CPUs, 2 days) runs every step on its compute node. It doesn't poll
the scheduler or resubmit itself, and it leaves nothing running on a login node. It
starts from the newest `ensembles/I<n>.csv`, and each pass does what the manual
commands do, with the same files:

1. run I<n>'s members (validated runs are skipped);
2. `extract` and `pilot` it;
3. `propose` I<n+1> and CHK<n+1>, then run the three check runs;
4. `verify`, then `promote` if accepted, and go back to 1.

It stops, after logging why:
- once I<n> has been run and scored, if I<n> is `auto.max_ensemble` (10);
- once I<n> has been run and scored, if the update that made I<n> lowered Φ at the mean
  by less than `auto.min_relative_gain` (2 %, from `results/CHK<n>_verify.json`);
- if a check isn't accepted. Nothing is promoted, and I<n> is the final ensemble;
- at the first failed run or error. It never retries on its own.

`results/auto_log.json` records every ensemble (members' Φ, and per window the members'
median share met and distance) and every check (Φ at both means, gain, α, verdict). The
end of the Slurm log prints the same table. Resubmitting resumes: validated runs are
skipped, and an existing proposal or promotion is reused (a promotion only if it is
byte-identical to its proposal).

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
python -m eki auto --workers 24         # the whole loop (what slurm/eki_auto.slurm runs)
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
| `tools/compare_finals.py` | v1/v2/v3 final runs on the cleaned fit deadlines and the final perimeter (`python -m tools.compare_finals`) |
| `audit/` | 2026-11-10 audit: operator, early misses (maps), acreage, objective, history; `audit/README.md` is the report |
| `v3_retry_a4/` | one damped (α = 4) retry of v3's rejected I08 → I09 step, set up, not run; see its `README.md` |
| `$SCRATCH/elmfire-class/eki_calibration/` | `shared/` and `runs/<BASELINE\|I00\|CHK01\|...>/<run_id>/` |
