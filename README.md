# Camp Fire uncertainty quantification with ELMFIRE

A class project in uncertainty quantification. It studies how uncertain wind, dead fuel
moisture and spotting change the ELMFIRE-simulated 2018 Camp Fire, then calibrates
time-varying weather corrections against the NIST fire-progression observations.

**Start with [`PROJECT.md`](PROJECT.md)**. It records the project state, the decisions
behind it and the conventions, and it links to everything else. ELMFIRE internals that
the project relies on are in [`ELMFIRE_CODEBASE.md`](ELMFIRE_CODEBASE.md).

| Folder | Study | Guide |
|---|---|---|
| `experiment/` | Variance-based (Sobol) sensitivity study, 6 factors × 3 wind stations, 5,376 runs (`sobolfire` package, MATLAB analysis) | `experiment/README.md` |
| `local_sensitivity/` | Local wind-direction sensitivity, ±50° around the bearing to Paradise (`lsa` package) | `local_sensitivity/README.md` |
| `eki_calibration/` | Ensemble Kalman inversion of time-varying weather corrections against NIST arrival observations (`eki` package) | `eki_calibration/README.md` |
| `HistoricalDataClass/`, `ParameterEstimationFiles/` | Source data: RAWS weather, Sobol design matrices, NIST TN 2135 observations, fire perimeter, the JG+CO weather blend | `DATA.md` |
| `tutorials/`, `slurm/`, `elmfire-class/` | ELMFIRE tutorials and the baseline-case setup | `class_project_setup.md` |

The code runs on Stanford's Sherlock cluster with ELMFIRE 1.1 and the `elmfire-vnv`
Python environment. Large data and all run output are kept out of git. [`DATA.md`](DATA.md)
lists what is excluded and how to get or rebuild it.
