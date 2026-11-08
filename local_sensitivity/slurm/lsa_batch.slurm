#!/bin/bash
#SBATCH --job-name=camp_lsa
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --constraint=CLASS:SH4_CBASE|CLASS:SH4_CPERF
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/oak/stanford/groups/gorle/jgoell/Wildfire/local_sensitivity/logs/%x_%A_%a.out

# One array task runs one chunk of a batch file: task i runs lines
# R*i+1 .. R*(i+1), where R is slurm.runs_per_task in config.yaml (101: one
# station per task), with up to --cpus-per-task fires at once, each a
# single-core ELMFIRE process in node-local scratch ($L_SCRATCH_JOB).
#   reads   local_sensitivity/{config.yaml, batches/<batch>.txt, prepared/}   (Oak)
#           <scratch_root>/shared/ from config.yaml                          (prepare)
#   writes  <scratch_root>/runs/<station>/<run_id>/
# and one status line per run to the job log above.
#
# `python -m lsa make-batch` writes the batch file and prints the sbatch line.
# Validated runs are skipped, so resubmitting a batch only reruns the others.
#
# Resources: the Sobol runs (same model settings, job 48588168) took 1.7-25 min
# and under 1 GB per fire, so 101 fires at 24 at a time take about 1-2 h.

BATCH=$1

cd /oak/stanford/groups/gorle/jgoell/Wildfire/local_sensitivity || exit 1

# ELMFIRE needs gcc/openmpi and gdal_translate; the pipeline needs the venv.
module purge
module load system gcc/12.4.0 openmpi/5.0.5
module load physics gdal/3.10.2 python/3.12.1
source /oak/stanford/groups/gorle/jgoell/venvs/elmfire-vnv/bin/activate

python -m lsa run-batch "$BATCH" --index "$SLURM_ARRAY_TASK_ID" --workers "$SLURM_CPUS_PER_TASK"
