#!/bin/bash
#SBATCH --job-name=camp_lsa_fixed
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --constraint=CLASS:SH4_CBASE|CLASS:SH4_CPERF
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/oak/stanford/groups/gorle/jgoell/Wildfire/local_sensitivity/fixed_direction/logs/%x_%A_%a.out

# The fixed-direction variant: same as lsa_batch.slurm, but with
# fixed_direction/config.yaml. Task i runs lines 101*i+1 .. 101*(i+1) of the
# batch file (one station per task), 24 fires at a time, each a single-core
# ELMFIRE process in node-local scratch.
#   writes  $SCRATCH/elmfire-class/local_sensitivity_fixed/runs/<station>/<run_id>/
# `python -m lsa --config fixed_direction/config.yaml make-batch` prints the sbatch line.
# Resources: the hourly-direction runs (job 49725936) took 35-55 min per task.

BATCH=$1

cd /oak/stanford/groups/gorle/jgoell/Wildfire/local_sensitivity || exit 1

module purge
module load system gcc/12.4.0 openmpi/5.0.5
module load physics gdal/3.10.2 python/3.12.1
source /oak/stanford/groups/gorle/jgoell/venvs/elmfire-vnv/bin/activate

python -m lsa --config fixed_direction/config.yaml run-batch "$BATCH" \
    --index "$SLURM_ARRAY_TASK_ID" --workers "$SLURM_CPUS_PER_TASK"
