# Data not in this repository

The repository holds code, configuration, notes, figures and small source data. The
files below are left out by `.gitignore`: they are large, can be rebuilt by the code,
or can be downloaded again. Paths are relative to the project folder on Oak
(`/oak/stanford/groups/gorle/jgoell/Wildfire`) unless they start with `$SCRATCH`.

## Run data (never in git)

| What | Where | Notes |
|---|---|---|
| Sobol study runs (5,376) | `$SCRATCH/elmfire-class/camp_sobol_d6/` | Per-run arrival rasters, weather rasters, records. `experiment/README.md`. |
| Local wind-direction study runs | `$SCRATCH/elmfire-class/local_sensitivity/` and `local_sensitivity_fixed/` | `local_sensitivity/README.md`. |
| EKI calibration runs | `$SCRATCH/elmfire-class/eki_calibration/` | `eki_calibration/README.md`. |
| Earlier pilots | `$SCRATCH/elmfire-class/camp_sobol/`, `camp_sobol_spot/` | Superseded. |

`$SCRATCH` is purged 90 days after files were last modified. Copy anything worth
keeping to Oak.

## Derived files (rebuilt by the code in this repository)

| Pattern | Made by |
|---|---|
| `experiment/matlab/*_border_toa.mat`, `JG_sobol_results.mat` | `experiment/matlab/extract_border_toa.m`, `main_sobol.m` |
| `experiment/qgis/*.tif` | `experiment/tools/qgis_export.py` |
| `experiment/qgis_setup/*.tif`, `setup.gpkg` (incl. the QoI mask `paradise_border_mask.tif`) | `experiment/tools/qgis_setup_export.py` (needs the Census shapefile below) |
| `local_sensitivity/qgis/*.tif` | `local_sensitivity/analysis/lsa_gradient.py` |
| `eki_calibration/results/*.npz` | `python -m eki extract` and `python -m eki pilot` |
| `*.out`, `*/logs/*` | Slurm job logs |
| `results/` (top level, 496 MB) | ELMFIRE tutorial runs and the 50-member validation-suite run of the baseline case (`slurm/`, `tutorials/`) |

## Source files to download again

| File (where the code expects it) | Source |
|---|---|
| `HistoricalDataClass/drive-download-20261019T162718Z-1-001/JarboGap__NfdrsHourly__09-30-2003__09-30-2026.csv` (97 MB) | Jarbo Gap RAWS (041214) NFDRS hourly record, from the class's shared Google Drive download of 2026-10-19 |
| `ParameterEstimationFiles/NIST_TN2135_report/NIST.TN.2135.pdf` (75 MB) | Maranghides et al. (2021), *A Case Study of the Camp Fire – Fire Progression Timeline*, NIST TN 2135, https://doi.org/10.6028/NIST.TN.2135 |
| `experiment/external/tl_2018_06_place/` (24 MB) | US Census TIGER/Line 2018 places, California: https://www2.census.gov/geo/tiger/TIGER2018/PLACE/tl_2018_06_place.zip (Paradise = GEOID 0655520) |

The NIST Appendix F workbook (`ParameterEstimationFiles/NIST_TN2135_AppendixF/`,
NIST dataset mds2-2954) is small and is in the repository.

## Not in this folder at all

| What | Where |
|---|---|
| ELMFIRE 1.1 (commit `cbf924a`), built with gcc 12.4.0 + OpenMPI 5.0.5 | `~/elmfire` |
| Baseline Camp Fire case (fuels, terrain, buildings), read only | `~/elmfire-class/camp_fire` (from `$OAK/jgoell/ELMFIRE_VnV_Suite`) |
| Python environment | `/oak/stanford/groups/gorle/jgoell/venvs/elmfire-vnv` |
