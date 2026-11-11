#!/bin/bash
#SBATCH --job-name=camp_baselines
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --constraint=CLASS:SH4_CBASE|CLASS:SH4_CPERF
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/oak/stanford/groups/gorle/jgoell/Wildfire/station_baselines/logs/%x_%j.out

# The three station baselines (JG, CO, MEAN) at once, each a single-core ELMFIRE
# process in node-local scratch. With these settings, EKI runs took 2-25 min and
# under 1 GB. Runs that already validated are skipped.

cd /oak/stanford/groups/gorle/jgoell/Wildfire/station_baselines || exit 1

# ELMFIRE needs gcc/openmpi and gdal_translate; the code needs the venv.
module purge
module load system gcc/12.4.0 openmpi/5.0.5
module load physics gdal/3.10.2 python/3.12.1
source /oak/stanford/groups/gorle/jgoell/venvs/elmfire-vnv/bin/activate

python -m baselines run-all --workers "$SLURM_CPUS_PER_TASK"
