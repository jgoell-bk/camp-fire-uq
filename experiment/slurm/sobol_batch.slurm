#!/bin/bash
#SBATCH --job-name=camp_sobol
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --constraint=CLASS:SH4_CBASE|CLASS:SH4_CPERF
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/oak/stanford/groups/gorle/jgoell/Wildfire/experiment/logs/%x_%A_%a.out

# One array task runs one chunk of a batch file: task i runs lines
# R*i+1 .. R*(i+1), where R is slurm.runs_per_task in config.yaml, with up to
# --cpus-per-task fires at once, each a single-core ELMFIRE process in
# node-local scratch ($L_SCRATCH_JOB).
#   reads   experiment/{config.yaml, batches/<batch>.txt, prepared/}         (Oak)
#           <scratch_root>/shared/ from config.yaml                           (prepare)
#   writes  <scratch_root>/{runs,controls}/.../<run_id>/
# and one status line per run to the job log above.
#
# `python -m sobolfire make-batch <selection>` writes the batch file and prints
# the array range, e.g.
#   sbatch --array=0-1%24 slurm/sobol_batch.slurm batches/pilot_spot.txt
# Validated runs are skipped, so resubmitting a batch only reruns the others.
#
# One node, sh04 CPU classes only, as for the CharLES runs (interiorDataset/slurm).
#
# Resources:
# - Spotting-off pilot (2026-10-19, 24 at once): 21-57 s and ~430 MB per fire.
#   The full run used 24 CPUs, 1G per CPU and 768 runs per task (6 tasks,
#   15-30 min each).
# - Spotting-on pilot (job 48561952, 20 at once): 4-10 min and <1 GB per fire
#   (mean JG 508 s, OP 406 s, CO 261 s). Now 24 CPUs, 2G per CPU, 4 h: one wave
#   of 24 for the spotting study. For the full design, raise runs_per_task in
#   config.yaml to about 192 (8 waves, about 1 h per task) and keep 4 h.

BATCH=$1

cd /oak/stanford/groups/gorle/jgoell/Wildfire/experiment || exit 1

# ELMFIRE needs gcc/openmpi and gdal_translate; the pipeline needs the venv.
module purge
module load system gcc/12.4.0 openmpi/5.0.5
module load physics gdal/3.10.2 python/3.12.1
source /oak/stanford/groups/gorle/jgoell/venvs/elmfire-vnv/bin/activate

python -m sobolfire run-batch "$BATCH" --index "$SLURM_ARRAY_TASK_ID" --workers "$SLURM_CPUS_PER_TASK"
