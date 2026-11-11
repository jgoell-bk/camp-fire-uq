#!/bin/bash
#SBATCH --job-name=camp_eki_chk
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --constraint=CLASS:SH4_CBASE|CLASS:SH4_CPERF
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/oak/stanford/groups/gorle/jgoell/Wildfire/eki_calibration/logs/%x_%A_%a.out

# Check runs (CHK batches): a few runs at ensemble means, before a proposed
# ensemble is promoted. The same command as eki_batch.slurm; only the size
# differs: one task runs every line of the batch (runs_per_task is 33), with up
# to --cpus-per-task fires at once, each a single-core ELMFIRE process in
# node-local scratch. `python -m eki propose` writes the batch file and prints
# this sbatch line. With the I00 settings, fires took 2-25 min and under 1 GB.
# The optional second argument is the config (default config.yaml, i.e. v1), e.g.
#   sbatch --array=0-0%2 slurm/eki_check.slurm v3_retry_a4/batches/CHK09.txt v3_retry_a4/config.yaml

BATCH=$1
CONFIG=${2:-config.yaml}

cd /oak/stanford/groups/gorle/jgoell/Wildfire/eki_calibration || exit 1

# ELMFIRE needs gcc/openmpi and gdal_translate; the pipeline needs the venv.
module purge
module load system gcc/12.4.0 openmpi/5.0.5
module load physics gdal/3.10.2 python/3.12.1
source /oak/stanford/groups/gorle/jgoell/venvs/elmfire-vnv/bin/activate

python -m eki --config "$CONFIG" run-batch "$BATCH" --index "$SLURM_ARRAY_TASK_ID" --workers "$SLURM_CPUS_PER_TASK"
