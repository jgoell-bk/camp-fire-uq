# v3_retry_a4: one damped retry of v3's rejected I08 → I09 step

Set up 2026-11-10 (audit, `../audit/README.md` §6). **Not yet run.**

v3 stopped because its I09 proposal (α = 2) raised Φ at the ensemble mean from 11.8769
to 12.9102 (`../v3/results/CHK09_verify.json`). This folder retries that one step with
more damping, and nothing else.

## What is the same as v3

The objective: deadline term, over-spread term (A_ref = 10,000 ac at 24 h), amplitude
and smoothness penalties, and their scales. Also the same:
- observations, constraints (1,062 fit deadlines) and holdout labels;
- the parameterization (10 coefficients, knots 0/6/12 h);
- the admissibility checks and step limits (mean move ≤ 0.5 sd, member move ≤ 1 sd);
- SEED 2023 for every run;
- the acceptance rule.

`config.yaml` is v3's with three lines changed (`setup_retry.py` checks that nothing
else differs):

| Setting | v3 | here | why |
|---|---|---|---|
| `fit.step.alpha0` | 1.0 | 4.0 | α doubles from alpha0 until the limits hold, so α ≥ 4 (v3's step used 2) |
| `paths.scratch_root` | `…/eki_calibration_v3` | `…/eki_calibration_v3_retry_a4` | own runs; v3's are never touched |
| `auto.max_ensemble` | 10 | 7 | guard: `python -m eki auto` refuses here (I08 > I07); this is one manual check |

A fresh folder is needed because `auto.py` reuses an existing `proposals/I09.csv` and
`ensembles/CHK09.csv`. Resubmitting v3's config would reuse the rejected proposal. v3's
proposal, check runs and rejection record are preserved as they were.

## Restart state and provenance (`PROVENANCE.json`, `REUSE.json`)

- **Copied byte for byte from v3** (checksums compared):
  - `ensembles/I08.csv` and `.json`, the accepted final state (sha256 `eb9e8c54…`);
  - `results/I08_*`, the I08 members' responses;
  - `prepared/*` and `elmfire.data.template`.
- **Shared inputs:** this folder's `shared/` in scratch was built fresh, and every file's
  checksum equals v3's.
- **The proposal** (`proposals/I09.csv`, from `python -m eki --config
  v3_retry_a4/config.yaml propose I08`):
  - α = 4 was the first α tried and met the limits: max mean move 0.234 sd, max member
    move 0.423 sd, no inadmissible member;
  - an independent EKI step at α = 4 matches it to 7e−15, and the same code at α = 2
    reproduces v3's rejected proposal to 7e−15;
  - the mean step points the same way as α = 2 (cosine 0.9998) and is 0.54 times as
    long. It is the EKI solve at α = 4, not a halved displacement.
- **Proposed mean** vs I08 mean and the rejected α = 2 mean:

| | I08 mean | α = 4 (this) | α = 2 (rejected) |
|---|---|---|---|
| speed × at 06/12/18 h | 2.31 / 1.29 / 0.77 | 2.47 / 1.32 / 0.78 | 2.60 / 1.35 / 0.78 |
| direction ° | +25.2 / +31.9 / +22.2 | +25.6 / +31.7 / +23.1 | +25.8 / +31.5 / +23.7 |
| moisture pp | −0.94 / −0.78 / −0.06 | −0.99 / −0.77 / −0.08 | −1.03 / −0.76 / −0.09 |
| P_EPS | 0.00687 | 0.00707 | 0.00724 |

- **Controls:** CHK09_I08MEAN and CHK09_I08MEAN_REP are v3's validated runs at the I08
  mean. They are bit-identical to each other and to CHK08_I08MEAN, with Φ 11.8769.
  - `reuse_controls.py` built both runs under this config and compared them with v3's:
    θ, every weather-raster checksum, the `weather.json` bands, P_EPS, SEED, start and
    stop times, template, baseline, ELMFIRE binary, shared files, and the namelist with
    its four folder paths masked. All were equal, so it copied the run folders here.
    v3's originals are untouched, and each copy carries a `REUSED.json`.
  - Expected differences, recorded only: the config checksum, the CHK09.csv path and
    checksum (its third row is new), the shared-manifest timestamp, build time and host.

## Run

**1 new simulation** (CHK09_I09MEAN, about 6 min on one core; v3's alpha = 2 run took 351 s):

```
cd /oak/stanford/groups/gorle/jgoell/Wildfire/eki_calibration
sbatch --array=0-0%2 slurm/eki_check.slurm v3_retry_a4/batches/CHK09.txt v3_retry_a4/config.yaml
```

Then (small, a minute):

```
cd /oak/stanford/groups/gorle/jgoell/Wildfire/eki_calibration
python -m eki --config v3_retry_a4/config.yaml verify CHK09
python v3_retry_a4/compare.py
```

## Comparison and acceptance (defined before the run)

All values come from single runs at ensemble-mean parameters, never from averages over
members. `compare.py` prints them for the I08 control, the rejected α = 2 proposal and
this α = 4 proposal:
- distance, over-spread, amplitude, smoothness and total Φ;
- the outside acres;
- for each deadline window (0–6, 6–12, 12–24 h): the weighted distance term, the met
  share and the mean km.

**Acceptance** is the existing rule (`eki/controlled.py` `verify`): accepted iff
Φ(α = 4 mean) − Φ(I08 mean) ≤ variability = |Φ(repeat) − Φ(I08 mean)|. That variability
is 0, so the rule means Φ(α = 4 mean) ≤ 11.8769.

**Acceptance is not convergence.** It would show that one smaller step lowers Φ at the
mean for SEED 2023.
- The bit-identical repeat does not measure spotting-realization variability, because
  building ignitions are random draws that change whenever θ changes. A small ΔΦ is
  inconclusive.
- If the step is accepted, `promote I09` and running I09's 32 members would be the next
  decision. That is not part of this setup, and `auto` is blocked here until
  `auto.max_ensemble` is raised.
- If it is rejected, v3's I08 stays final. A further retry (α = 8, max mean move
  0.122 sd) or a seed-variability check would be the options.

This setup does not cover a full calibration, an acreage-weight sweep or any holdout
evaluation.

## Files

| File | What |
|---|---|
| `config.yaml` | v3's config, three lines changed (above) |
| `setup_retry.py` | builds this folder: config, copies, shared inputs, `PROVENANCE.json` |
| `reuse_controls.py` | checks and reuses v3's two control runs, `REUSE.json` |
| `compare.py` | the side-by-side table after `verify` |
| `proposals/I09.csv`, `ensembles/CHK09.csv`, `batches/CHK09.txt` | the α = 4 proposal and its check batch |
| `results/auto_log.json` | one record: the guard refusing `auto` (tested 2026-11-10) |
