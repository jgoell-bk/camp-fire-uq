# Camp Fire ELMFIRE calibration with time-varying weather corrections

**Project:** CFD Wildfire — uncertainty quantification and parameter estimation  
**Prepared:** November 8, 2026  
**Status:** Proposed research and implementation plan; no calibration has been run for this document.

## 1. Proposed approach

Estimate a small set of coefficients describing time-varying corrections to the average of the Jarbo Gap and Cohasset station records. Calibrate those coefficients against observed fire-arrival times using regularized ensemble Kalman inversion (EKI). Evaluate predictive performance on withheld observations, and quantify uncertainty through a separate, explicitly specified procedure.

The attraction of EKI is computational: an iteration evaluates an ensemble of candidate weather histories in parallel and uses their input–output relationships to propose the next ensemble. ELMFIRE remains a black box; derivatives of its source code are unnecessary. This is a method to investigate, not an established guarantee of successful inversion for this Camp Fire configuration. [1–3]

The scientific target is **effective weather forcing that improves ELMFIRE's reconstruction, conditional on the fixed model configuration**. Estimated weather corrections can also compensate for fuel-map errors, missing physics, ignition errors, or spotting assumptions. They should not automatically be interpreted as uniquely recovered atmospheric conditions.

### Established project choices and provisional recommendations

| Item | Status |
| --- | --- |
| Case study | Camp Fire |
| Forward model | Existing ELMFIRE configuration; approximately 45 minutes per run, as reported by the researcher |
| Baseline weather | Time-varying average of Jarbo Gap and Cohasset records |
| Openshaw | Excluded from the proposed baseline because the researcher has concerns about its data |
| Observational data | An Excel file containing arrival locations and associated timing information; its contents still need inspection |
| Inferred input families | Wind speed, wind direction, and fuel moisture |
| Temporal behavior | Corrections may vary with physical time |
| P_EPS | Fixed during this calibration; its chosen value must be recorded |
| Inversion method | Regularized EKI is the proposed first method to test |
| Initial temporal resolution | Provisionally three correction coefficients per input family, giving nine unknowns |
| Ensemble and budget | Provisionally 32 members and up to eight evaluated batches for an initial fitting budget |
| Moisture-class mapping, knot times, bounds, error scales | Unresolved; determine from the data and pilot before production calibration |

Neither a full hourly inversion nor a full Cartesian grid is the starting proposal. Use a compact temporal representation and increase complexity only when recovery tests and predictive validation justify it.

## 2. Research questions and success criteria

### Primary questions

1. Can low-dimensional, time-varying weather corrections improve arrival predictions beyond the uncorrected station-average baseline?
2. Which correction coefficients are informed by the arrival observations, and which remain weakly constrained or interchangeable?
3. Do improvements persist at withheld locations or times and across plausible observation-error assumptions?
4. How much do results depend on temporal resolution, regularization, initial ensemble, and stochastic simulation variability?

### Evidence of a useful result

- Improved error on withheld observation groups, alongside training improvement.
- Physically admissible corrected weather histories.
- Comparable predictive behavior from different initial ensembles.
- Transparent reporting of parameter tradeoffs and poorly constrained periods.
- A quantified distinction between observation uncertainty, stochastic ELMFIRE variability, and conditional parameter uncertainty.

A small training error or a tightly clustered final ensemble alone is insufficient evidence of successful inference.

## 3. Audit and freeze the existing simulation

Before calibration, record the ELMFIRE revision, build settings, run command, resource allocation, namelist/configuration files, input checksums, ignition specification, simulation horizon, and all fixed physical parameters.

Freeze fuel, terrain, canopy, ignition, suppression, spread-adjustment, and spotting settings for the main experiment. Record any internal Monte Carlo perturbations or randomized wind fluctuations. Avoid unknowingly adding random weather perturbations on top of the explicitly calibrated weather histories.

Fixing P_EPS removes it from the optimization but does not establish its correctness. The primary results are conditional on its fixed value and the remaining spotting configuration. An optional later robustness experiment can test alternative fixed configurations without adding them to the first inversion.

Verify the units and conventions of the installed ELMFIRE version against its documentation and existing scripts: wind measurement height, speed units, meteorological direction, moisture percentage versus fraction, raster-band times, and arrival-output units. The official guide documents time-varying weather rasters and relevant input switches. [4]

## 4. Prepare the observations

### Canonical observation table

Convert the Excel data into a documented analysis table without changing the source workbook. Suggested fields are:

| Field | Meaning |
| --- | --- |
| observation_id | Stable unique identifier |
| x, y | Location in a declared coordinate reference system |
| crs | Source and transformed CRS metadata |
| observed_time_original | Original timestamp or time description |
| observed_time_utc | Unambiguous absolute timestamp where possible |
| arrival_seconds | Elapsed time relative to the simulation's documented time origin |
| lower_seconds, upper_seconds | Arrival interval, if reported |
| timing_sigma_seconds | Justified timing-error scale, if available |
| source_id, source_group | Provenance and shared-report grouping |
| location_uncertainty | Available positional uncertainty |
| quality_flag | Inclusion or exclusion decision and reason |
| validation_group | Spatial/temporal/source grouping used for validation |

### Required checks

- Resolve time zones and daylight-saving conventions explicitly; do not infer them from the spreadsheet's display format.
- Verify the ignition time and simulation time origin against observational timestamps.
- Map every point onto the simulation domain and check CRS transformations.
- Identify duplicates, rounded timestamps, reconstructed values, and observations derived from the same report or perimeter.
- Distinguish an exact arrival estimate from statements such as "burned by this time" or "arrival occurred within this interval."
- Plot observation coverage through physical time and space before choosing correction times.

An interval or one-sided report is not automatically a point observation with independent Gaussian error. Retain its original information. A Gaussian approximation using an interval midpoint requires a stated justification and a sensitivity check; otherwise use an observation formulation that respects the interval or inequality.

### Holdout design

Reserve groups of observations for final predictive assessment before tuning the model. Prefer spatial, temporal, or source-based groups over a purely random point split when nearby observations share information. Use training-group cross-validation for tuning if sample size permits; do not repeatedly tune against the final holdout.

## 5. Construct the baseline weather histories

Let the baseline records be U_0(t), D_0(t), and M_{k,0}(t), where k identifies a moisture class.

1. Quality-control the two station records independently.
2. Align timestamps, sampling intervals, units, direction convention, and measurement-height treatment.
3. Define missing-data rules in advance; avoid silently changing station weights as records disappear.
4. Construct the station-average record using a documented rule.
5. Preserve the observed physical-time variation when preparing the model inputs.

Wind direction requires circular handling. Two defensible but different constructions are scalar-mean speed with circular-mean direction, or a mean of full horizontal wind vectors followed by reconstruction of magnitude and direction. Opposing station winds can make the vector-mean speed small, while their circular mean direction may also be poorly defined. Inspect such periods explicitly. [5]

For the main baseline, retain the researcher's intended averaging convention once verified. Do not silently change between scalar and vector averaging. Station disagreement is a diagnostic of representativeness as well as possible error; it is not automatically an estimate of measurement-noise variance.

A domain-wide corrected average is an effective forcing assumption. It cannot reproduce arbitrary spatial differences in local wind. Preserve any established spatial mapping in the existing workflow, and document whether corrections act on station records before mapping or on the final model fields.

## 6. Parameterize time-varying corrections

### 6.1 Proposed equations

Use the following starting parameterization:

$$
U(t)=U_0(t)\exp[a(t)],
$$

$$
D(t)=\operatorname{wrap}_{360}\{D_0(t)+d(t)\},
$$

$$
M_k(t)=M_{k,0}(t)+w_k m(t).
$$

Here a(t) is a dimensionless log-speed correction, d(t) is a direction correction in degrees, and m(t) is a moisture correction in percentage points if the baseline is represented in percent. The fixed weights w_k determine how the shared moisture correction acts on different fuel classes.

Zero corrections reproduce the baseline exactly. The exponential keeps the speed multiplier positive. It cannot create a nonzero wind from an exactly zero baseline; if meaningful calm periods require that flexibility, test an additive horizontal-vector correction instead.

The moisture mapping is an explicit modeling decision. Applying equal percentage-point offsets to M1, M10, and M100 is one option, not an established project choice. Using only M1 or a different fixed class relationship is another. Do not assign independent class corrections without updating the parameter count and the identifiability assessment. Live moisture remains fixed unless deliberately added to scope.

### 6.2 Temporal basis

Represent each correction with piecewise-linear basis functions at selected physical times:

$$
a(t)=\sum_{r=1}^{K_U}a_r b_r^{U}(t),\quad
d(t)=\sum_{r=1}^{K_D}d_r b_r^{D}(t),\quad
m(t)=\sum_{r=1}^{K_M}m_r b_r^{M}(t).
$$

Provisionally set K_U=K_D=K_M=3. The parameter vector then contains nine coefficients:

$$
\theta=(a_1,a_2,a_3,d_1,d_2,d_3,m_1,m_2,m_3)^T.
$$

Choose knot times from the observation coverage and meteorological transitions, including the calibration-window endpoints. Different input families may eventually need different knots. Do not estimate changes after the latest informative arrival merely because the simulation continues.

Interpolate an unwrapped direction-correction curve and wrap only the resulting physical direction. Interpolation and station averaging must handle the 0/360-degree boundary consistently. Very large directional uncertainty may favor a horizontal-component representation over a local angular correction model.

### 6.3 Bounds and scaling

Set correction bounds from event information, quality-controlled station records, and relevant literature. Record the evidence and mark any provisional choices. Do not copy broad sensitivity-study distributions as if they measured event-specific station errors.

Standardize coefficients before updating them so degrees, percentage points, and log multipliers have comparable numerical scales. Check admissibility at every generated weather timestamp, not only at knots. If using projection onto bounds, identify the algorithm as a projected EKI variant, log every projection, and inspect persistent boundary accumulation. Frequent clipping can conceal a poor parameterization.

## 7. Define the forward observation operator

For deterministic ELMFIRE, define

$$
G_i(\theta)=T_{\mathrm{ELMFIRE}}(x_i;\theta),\qquad
G(\theta)=(G_1,\ldots,G_N)^T.
$$

The complete operation is: coefficients → correction curves → model weather inputs → ELMFIRE → arrival-time extraction at fixed observation locations.

Use a documented raster-sampling rule. Avoid interpolation across unburned/NoData cells, verify raster orientation and cell alignment, and quantify whether positional uncertainty is important relative to the grid. Keep the observation-vector order identical across all runs.

### Unreached locations and failures

Distinguish job failure, out-of-domain points, insufficient simulation horizon, and genuine non-arrival under a candidate history. Never drop observations separately for each member, substitute zero arrival, or pass NaN into the covariance update.

The initial Gaussian arrival-time formulation is suitable only when the pilot demonstrates adequately defined finite predictions on its calibration observation set. If plausible histories regularly fail to reach observed locations, resolve this before production: extend an insufficient horizon where appropriate, or adopt a fixed observation operator that supports reachability/inequality information. Candidate alternatives include front-position or signed-distance information at specified observation times, after checking available outputs. Arbitrary capped arrival times can create flat artificial responses and must not be introduced without testing and disclosure.

## 8. Handle stochastic ELMFIRE behavior

If active spotting or other physics make the simulator stochastic, write it as G(theta, xi), where xi denotes the random realization.

Run a small replication pilot at baseline and selected perturbed histories. Compare seed-to-seed variation with changes induced by plausible weather corrections. Determine whether the target is a representative realization, mean arrival response, or a richer stochastic observation model.

If targeting a mean response, average over a documented fixed seed set where needed and account for residual Monte Carlo error. Common seeds across parameter candidates can improve comparisons in some models, but changes in random-number consumption can reduce that benefit. Verify it empirically. Replicate final predictions with fresh seeds.

A single fixed seed produces a conditional calibration. It is not an estimate of stochastic-model uncertainty. If seed variability dominates the weather signal, increasing ensemble size alone does not solve the observation-model problem.

## 9. Statistical objective and regularization

For the initial finite, point-arrival formulation, propose

$$
\Phi(\theta)=
\frac12[y-G(\theta)]^T R^{-1}[y-G(\theta)]
+\frac12\|B\theta\|^2.
$$

R is a fixed, justified residual-error covariance for the first implementation. It may include an observation component and an explicitly justified allowance for other residual sources. It must not be tuned arbitrarily to make the fitted error look statistically acceptable. If its value depends on parameters, the objective and likelihood require corresponding changes.

Shared reports, nearby locations, and similar arrival histories can produce dependence. A diagonal R is an approximation that needs assessment, not proof of independent data. Use grouped observations or an appropriate structured covariance when supported. A flexible fitted discrepancy process can confound the weather corrections; avoid adding one without an identifiability strategy.

### Construction of B

Build B from two types of rows:

- **Amplitude rows:** a coefficient divided by its plausible correction scale, multiplied by the square root of an amplitude-penalty weight.
- **Temporal rows:** differences between neighboring correction coefficients, divided by an appropriate scale and time-spacing normalization, multiplied by the square root of a smoothness-penalty weight.

Thus B can encode distinct scales and temporal behavior for speed, direction, and moisture. Penalize correction curves rather than smoothing away the measured baseline weather. Document the chosen normalization so changing knot spacing does not silently redefine the penalty.

Literature-informed scales may be used here without claiming a full Bayesian posterior. With suitable Gaussian assumptions, a related penalized estimator also has a maximum-a-posteriori interpretation. The estimation algorithm and its uncertainty interpretation remain separate decisions.

## 10. EKI formulation

### 10.1 Basic ensemble relationships

Let theta_n^(j) be member j at calibration iteration n, with J members. Evaluate g_n^(j)=G(theta_n^(j)). Physical simulation time t is distinct from optimization iteration n.

Define the ensemble means and sample covariances:

$$
\bar\theta_n=\frac1J\sum_j\theta_n^{(j)},\qquad
\bar g_n=\frac1J\sum_jg_n^{(j)},
$$

$$
C_n^{\theta g}=\frac1{J-1}\sum_j
(\theta_n^{(j)}-\bar\theta_n)(g_n^{(j)}-\bar g_n)^T,
$$

$$
C_n^{gg}=\frac1{J-1}\sum_j
(g_n^{(j)}-\bar g_n)(g_n^{(j)}-\bar g_n)^T.
$$

A basic deterministic update is

$$
\theta_{n+1}^{(j)}=\theta_n^{(j)}+
C_n^{\theta g}(C_n^{gg}+R)^{-1}
[y-g_n^{(j)}].
$$

This equation illustrates how ensemble input–output covariance translates arrival residuals into coefficient changes. By itself it does not implement the explicit B penalty above. [1]

### 10.2 Explicit regularization by augmentation

For the proposed penalized problem, use augmented outputs and targets:

$$
F(\theta)=\begin{bmatrix}G(\theta)\\B\theta\end{bmatrix},\quad
z=\begin{bmatrix}y\\0\end{bmatrix},\quad
\Gamma=\begin{bmatrix}R&0\\0&I\end{bmatrix}.
$$

Apply the ensemble update to F and z:

$$
\theta_{n+1}^{(j)}=\theta_n^{(j)}+
C_n^{\theta F}(C_n^{FF}+\alpha_n\Gamma)^{-1}
[z-F(\theta_n^{(j)})].
$$

The positive numerical parameter alpha_n controls update damping; alpha_n=1 gives the corresponding basic augmented update. Choose and document a tested adaptation or step-control rule rather than claiming this schematic formula is a complete production solver. The augmentation follows the Tikhonov-EKI construction. [2]

Use linear solves, not an explicitly formed matrix inverse. For large observation vectors, exploit ensemble-rank structure rather than forming unnecessary dense matrices. Keep consistent covariance normalizations throughout the implementation.

### 10.3 Initialization and convergence diagnostics

- Include the baseline as a diagnostic candidate and initialize a diverse ensemble around zero correction.
- Sample plausible, temporally correlated coefficient vectors; independently jittering every timestamp is not required.
- With nine unknowns, ensure the centered initial ensemble spans all intended parameter directions. J greater than nine is necessary but not sufficient for full rank.
- Basic unmodified EKI updates remain in the initial ensemble's affine span; constraints, rejuvenation, or other modifications can alter that property. [1]
- Track objective, training and tuning-validation errors, coefficient spread, ensemble rank, boundary activity, covariance conditioning, and accepted step sizes.
- Stop using a documented combination of residual adequacy, validation behavior, objective improvement, and budget. A collapsing ensemble alone is not a success criterion.
- Re-evaluate the final proposed coefficient vector with ELMFIRE; nonlinear predictions at the ensemble mean need not equal the mean of ensemble predictions.

EKI is not guaranteed to find a global optimum or accurately recover multimodal solutions. Test alternative initial ensembles and a modest alternative fitting method when the results suggest multiple solutions. Repeated use of the same observations is an optimization device, not the acquisition of fresh independent evidence. [1–3]

## 11. Staged experiment and decision gates

### Stage A — Data and baseline audit

Prepare the observation table, plot station records, verify conventions, and run the uncorrected baseline. If resources permit, compare separate Jarbo Gap and Cohasset forcings as diagnostic scenarios. These runs assess the consequences of the average; they do not turn station disagreement into a formal error distribution.

### Stage B — Replication and local-response pilot

Use roughly 4–8 baseline realizations if stochastic variability is active. Perturb each of the nine coefficients positively and negatively at justified scales: 18 additional evaluations, before any replication. Inspect which arrivals respond, how strongly they respond, and whether finite predictions remain available.

Compare the scaled response patterns of different coefficients. Near-parallel patterns suggest confounding, such as wind-speed and moisture adjustments producing similar delays. This pilot is diagnostic; the production EKI still uses ensemble covariance rather than these finite differences as required derivatives.

**Gate:** revise parameterization, observation treatment, or stochastic averaging if the useful signal is weak, strongly discontinuous, or dominated by simulation noise.

### Stage C — Synthetic recovery

Choose known smooth correction curves, generate synthetic arrivals at the actual observation locations, add a justified noise realization, and attempt recovery using the same algorithm and data geometry. Include a case with changing direction corrections. Keep a reference case hidden from solver tuning where feasible.

Assess recovery of curves and predictions, not just the objective. Repeat selected noise/initialization cases if the budget permits. A failure may indicate weak information or algorithm limitations; it is not automatically evidence against time-varying weather corrections.

This is a separate inversion with its own simulation cost. It does not establish ELMFIRE's physical realism or automatically validate later confidence intervals.

### Stage D — Real-data calibration

Run regularized EKI using training observations. Treat 32 members and eight evaluated batches as an initial budget cap. Compare a few predeclared regularization choices or use grouped training cross-validation, accounting for the additional fits. Inspect the simpler constant-correction model as a benchmark for the value of temporal flexibility.

### Stage E — Verification and prediction

Run the selected estimate through actual ELMFIRE, assess the untouched holdout, and use fresh random seeds if relevant. Compare baseline and calibrated arrival-error maps, mean bias, MAE, RMSE, and timing-group errors. State whether evaluation is interpolation within the hindcast window or genuine prediction into a later period.

### Stage F — Conditional uncertainty analysis

Select the uncertainty procedure in Section 13. Do not report the final EKI ensemble as a calibrated confidence distribution.

## 12. Cluster execution and reproducibility

Use a controller that writes a batch manifest, submits one isolated run per member/seed, waits for completion, validates outputs, and only then performs the ensemble update. A scheduler job array is a suitable pattern where available; adapt to the actual cluster rather than assuming a specific scheduler.

Every run needs isolated weather inputs, configuration, scratch, logs, and output paths. Read-only landscape data may be shared. Avoid overwriting shared weather files while other members are running.

Record at minimum: experiment ID, iteration, member, coefficient vector, transformed physical inputs, seed, ELMFIRE revision, input hashes, resource allocation, scheduler ID, status, wall time, extracted arrivals, valid-output flags, and objective components.

Checkpoint the evaluated ensemble, outputs, random-generator state, fixed settings, observation order, scaling, B, R, and solver state after every batch. Reuse an evaluation only when all relevant inputs and seeds match; use an explicit cache key. Retry infrastructure failures with identical inputs, while retaining the original failure record. A physically poor simulation is not an infrastructure failure.

### Scheduling illustration

Let P be the number of full ELMFIRE runs that can execute concurrently at the resources used to obtain the 45-minute runtime. For J members, I evaluated batches, and S realizations per member, the idealized fitting time is

$$
T_{\mathrm{fit}}\approx I\left\lceil\frac{JS}{P}\right\rceil(45\;\mathrm{minutes}).
$$

| Example | Evaluations | Idealized time |
| --- | ---: | ---: |
| J=32, I=8, S=1, P=32 | 256 | 6 hours |
| J=32, I=8, S=1, P=16 | 256 | 12 hours |
| J=32, I=8, S=4, P=32 | 1,024 | 24 hours |

These figures exclude queues, I/O, retries, final verification, pilots, synthetic recovery, tuning fits, restarts, and uncertainty calculations. Concurrency means full simulations, not CPU cores. Benchmark actual concurrent throughput before choosing the ensemble size.

## 13. Uncertainty quantification after fitting

Separate four questions: uncertainty in arrival observations, random simulation variability, uncertainty in fitted correction coefficients, and structural/model-configuration uncertainty. The final EKI ensemble alone does not distinguish these. [2]

### Proposed frequentist route

Use a bootstrap/refitting procedure conditional on the selected temporal basis, constraints, regularization, and fixed model configuration. Generate observation perturbations under a justified covariance model, or resample suitably defined groups where that sampling interpretation is defensible. Avoid independently resampling correlated neighboring points by default.

Repeat the complete fitting rule for each bootstrap dataset. If regularization is selected from data, either repeat that selection inside each replicate or clearly state that uncertainty is conditional on the fixed selected value. Include simulation variability in a deliberate way when it matters.

Report conditional parameter/curve intervals and arrival-prediction intervals, with their construction stated. Penalization bias and nonlinear weak identification can affect coverage; use synthetic repeated-data experiments to assess the procedure where feasible. The general combination of surrogate-based frequentist calibration and bootstrap inference has methodological precedent, but its validity for this application still needs assessment. [6]

### Computational acceleration

A surrogate may approximate arrival responses near the fitted region using accumulated simulation evaluations plus additional designed runs. EKI evaluations alone may cover a narrow region. Validate against held-out actual simulations near candidate uncertainty limits, and assess surrogate error relative to the residual scales. A mean-only deterministic emulator can understate uncertainty if stochastic or emulator variability is ignored.

Bootstrap refits can then use the validated surrogate, with selected refits verified using ELMFIRE. If the surrogate is unreliable near arrival discontinuities, use more direct runs, simplify the target, or report a limited uncertainty analysis rather than nominal confidence claims.

Bayesian inference remains a possible later extension over the same coefficients, potentially using a validated emulator. EKI optimization does not provide MCMC samples and does not resolve posterior multimodality.

## 14. Proposed software responsibilities

The implementation can extend the existing project language and scheduler tooling; a language rewrite is unnecessary.

| Component | Responsibility |
| --- | --- |
| Observation preparation | Read workbook, retain provenance, normalize times/CRS, assign groups |
| Weather preparation | QC station records, average consistently, build baseline inputs |
| Parameter transform | Coefficients to admissible correction histories and rasters |
| Simulation adapter | Isolated ELMFIRE run, seed/config capture, status handling |
| Observation operator | Extract a fixed-order output vector and explicit validity flags |
| Objective | R, B, observation weighting, separate diagnostic metrics |
| EKI controller | Ensemble generation, covariance update, constraints, checkpoints |
| Experiment runner | Pilots, synthetic recovery, calibration, holdout evaluation |
| UQ driver | Bootstrap datasets/refits and optional surrogate validation |

### Controller pseudocode

```text
load frozen simulation settings and audited observations
load baseline weather and documented correction specification
construct coefficient scaling, bounds, B, and R
initialize a reproducible, diverse, full-rank ensemble

for each allowed evaluated batch:
    transform each member into a full admissible weather history
    submit isolated ELMFIRE evaluations for required seeds
    collect outputs and validate fixed-order observation predictions
    resolve failures using the predeclared policy
    compute data misfit, penalty, and diagnostic metrics
    checkpoint the evaluated state
    if stopping criteria are met:
        break
    compute regularized ensemble update with linear solves
    apply documented step control and admissibility treatment

select a candidate using training/tuning information
verify it with actual ELMFIRE and assess final holdout observations
perform the separately specified uncertainty analysis
```

## 15. Required validation before a large experiment

- Zero correction reproduces the baseline weather inputs and predictions under identical seeds.
- Direction averaging, correction, wrapping, and time interpolation handle the circular boundary correctly.
- Moisture-class mapping and unit conversions are explicit and reversible where applicable.
- Observation extraction agrees with manually checked raster cells and model times.
- Failure and non-arrival cases cannot silently alter the observation vector.
- The covariance/update code works on a small linear inverse problem with an independently calculated regularized solution.
- The synthetic experiment tests the actual observation layout and proposed changing-direction corrections.
- Checkpoint/restart reproduces the next batch without duplicating completed evaluations.
- Concurrent runs cannot change each other's input files.

## 16. First decisions needed from the actual data

1. Observation timing, location, provenance, and uncertainty structure.
2. Simulation/calibration window and exact time origin.
3. Scalar/circular versus full-vector station averaging and handling of missing records.
4. Existing spatial weather mapping and measurement-height treatment.
5. Which dead-fuel-moisture classes share the correction, and with which fixed weights.
6. Correction knot times, bounds, amplitude scales, and temporal scales.
7. Active stochastic features, seed control, and the target simulation response.
8. Available full-run concurrency and replication budget.

The immediate implementation milestone is an audited observation table, reproducible baseline, and a small response/replication pilot. Production inversion should follow the pilot's decision gate.

## 17. References and scope of evidence

The EKI papers establish general algorithms and properties. They do not demonstrate successful calibration of this specific Camp Fire ELMFIRE setup. The observation design, coefficient counts, budgets, validation gates, and workflow above are project-specific proposals.

1. Iglesias, M. A., Law, K. J. H., and Stuart, A. M. (2013). *Ensemble Kalman methods for inverse problems*. Inverse Problems, 29, 045001. [Paper](https://arxiv.org/abs/1209.2736). Foundational ensemble inversion formulation and its relation to derivative-free optimization.
2. Chada, N. K., Stuart, A. M., and Tong, X. T. (2020). *Tikhonov regularization within ensemble Kalman inversion*. SIAM Journal on Numerical Analysis, 58(2), 1263–1294. [Paper](https://arxiv.org/abs/1901.10382). Explicit regularization and discussion of ensemble collapse and uncertainty limitations.
3. Huang, D. Z., Schneider, T., and Stuart, A. M. (2022). *Iterated Kalman methodology for inverse problems*. Journal of Computational Physics, 463, 111262. [Paper](https://arxiv.org/abs/2102.01580). Related iterative inversion formulations and regularized optimization perspectives.
4. ELMFIRE. *Computational Guide*. [Official documentation](https://elmfire.io/user_guide.html). Verify details against the installed revision rather than assuming the latest documentation exactly matches the existing case.
5. NOAA National Data Buoy Center. *What averaging procedures are performed on the wind measurements?* [Official explanation](https://www.ndbc.noaa.gov/faq/wndav.shtml). Scalar/unit-vector and full-vector averaging are distinct operations; the page does not prescribe a station-combination model for this case.
6. Wong, R. K. W., Storlie, C. B., and Lee, T. C. M. (2017). *A frequentist approach to computer model calibration*. Journal of the Royal Statistical Society: Series B, 79(2), 635–648. [Manuscript](https://arxiv.org/abs/1411.4723). Frequentist calibration with emulation, discrepancy, and bootstrap uncertainty methods.

The previously discussed fire-spread sensitivity paper can help motivate exploration of input distributions. Its assumed sensitivity distributions do not by themselves establish event-specific correction-error distributions for Jarbo Gap, Cohasset, or the Camp Fire.
