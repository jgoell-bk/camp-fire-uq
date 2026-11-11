#!/bin/bash
#SBATCH --job-name=camp_eki_auto
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --constraint=CLASS:SH4_CBASE|CLASS:SH4_CPERF
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/oak/stanford/groups/gorle/jgoell/Wildfire/eki_calibration/logs/%x_%j.out

# The whole EKI loop in this one job (python -m eki auto; eki/auto.py): run the newest
# ensemble, extract and pilot it, propose the next update, run its three check runs,
# verify, promote, repeat. Stopping rules are in config.yaml (auto). Fires run 24 at a
# time, each a single-core ELMFIRE process in node-local scratch; an ensemble pass
# plus its check took about 1.5 h with the I02/I03 settings. Every step is logged in
# results/auto_log.json. If the job ends early (time limit, failed run), resubmitting
# it resumes where it stopped: validated runs are skipped.
#
# Optional argument: the config file (default config.yaml, the v1 calibration), e.g.
#   sbatch slurm/eki_auto.slurm v2/config.yaml
# runs the v2 calibration (cleaned observations) in its own folders and scratch root.
#   sbatch slurm/eki_auto.slurm v3/config.yaml
# runs v3 (v2 + the over-spread term + P_EPS as a coefficient).

CONFIG=${1:-config.yaml}

cd /oak/stanford/groups/gorle/jgoell/Wildfire/eki_calibration || exit 1

# ELMFIRE needs gcc/openmpi and gdal_translate; the pipeline needs the venv.
module purge
module load system gcc/12.4.0 openmpi/5.0.5
module load physics gdal/3.10.2 python/3.12.1
source /oak/stanford/groups/gorle/jgoell/venvs/elmfire-vnv/bin/activate

python -m eki --config "$CONFIG" auto --workers "$SLURM_CPUS_PER_TASK"
