#!/bin/bash
#SBATCH --job-name=camp_restore_inputs
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --constraint=CLASS:SH4_CBASE|CLASS:SH4_CPERF
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH --output=/oak/stanford/groups/gorle/jgoell/Wildfire/experiment/logs/%x_%j.out

# Rebuild the weather rasters (ws, wd, m1, m10, m100) of every validated run into
# <run folder>/inputs/, using the same code that built them for the run, and keep
# each file only if its sha256 matches the one recorded in the run's weather.json.
# Runs that already have matching rasters are skipped, so this can be rerun.
# No ELMFIRE is run; this only writes about 180 KB per run.
#   cd /oak/stanford/groups/gorle/jgoell/Wildfire/experiment
#   sbatch slurm/restore_inputs.slurm
# The log ends with counts (restored / present / not_validated / mismatch); the
# job exits non-zero if any run's rebuilt rasters differ from the recorded ones.

cd /oak/stanford/groups/gorle/jgoell/Wildfire/experiment || exit 1
module purge
module load system python/3.12.1
source /oak/stanford/groups/gorle/jgoell/venvs/elmfire-vnv/bin/activate

python -m sobolfire restore-inputs --workers "$SLURM_CPUS_PER_TASK"
