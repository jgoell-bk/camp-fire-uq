#!/bin/bash
#SBATCH --job-name=camp_eki
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --constraint=CLASS:SH4_CBASE|CLASS:SH4_CPERF
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/oak/stanford/groups/gorle/jgoell/Wildfire/eki_calibration/logs/%x_%A_%a.out

# One array task runs one chunk of a batch file: task i runs lines
# R*i+1 .. R*(i+1), where R is slurm.runs_per_task in config.yaml (33: the
# baseline and the 32 iteration-0 members in one task), with up to
# --cpus-per-task fires at once, each a single-core ELMFIRE process in
# node-local scratch ($L_SCRATCH_JOB).
#   reads   eki_calibration/{config.yaml, batches/<batch>.txt, prepared/, ensembles/}   (Oak)
#           <scratch_root>/shared/ from config.yaml                                  (prepare)
#   writes  <scratch_root>/runs/<BASELINE|I00>/<run_id>/
# and one status line per run to the job log above.
#
# `python -m eki make-batch` writes the batch file and prints the sbatch line.
# Validated runs are skipped, so resubmitting a batch only reruns the others.
#
# Resources: with the same model settings the Sobol runs took 1.7-25 min and
# under 1 GB per fire. Speed multipliers up to 3x make larger, slower fires, so
# the time limit allows for runs several times longer than that (33 fires, 24
# at a time, is two rounds).

BATCH=$1

cd /oak/stanford/groups/gorle/jgoell/Wildfire/eki_calibration || exit 1

# ELMFIRE needs gcc/openmpi and gdal_translate; the pipeline needs the venv.
module purge
module load system gcc/12.4.0 openmpi/5.0.5
module load physics gdal/3.10.2 python/3.12.1
source /oak/stanford/groups/gorle/jgoell/venvs/elmfire-vnv/bin/activate

python -m eki run-batch "$BATCH" --index "$SLURM_ARRAY_TASK_ID" --workers "$SLURM_CPUS_PER_TASK"
