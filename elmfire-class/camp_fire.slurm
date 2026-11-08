#!/bin/bash
#SBATCH --job-name=elmfire-camp
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --time=03:00:00
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --output=%x-%j.out

# 2018 Camp Fire case: 50 ensemble members, one MPI rank each, all on one node.
# Copies the case to $SCRATCH/elmfire-class/camp-<jobid>, runs ELMFIRE there,
# then scores it against VIIRS. Everything it writes stays in that folder.
#   cd $SCRATCH/elmfire-class && sbatch ~/elmfire-class/camp_fire.slurm
# First full run: 45 min, ~56 GB.

ELMFIRE=${ELMFIRE:-$HOME/elmfire}
CASE=${CASE:-$HOME/elmfire-class/camp_fire}
VENV=${VENV:-$OAK/jgoell/venvs/elmfire-vnv}
RUN=$SCRATCH/elmfire-class/camp-$SLURM_JOB_ID

module purge
module load system gcc/12.4.0 openmpi/5.0.5
module load physics gdal/3.10.2 python/3.12.1
source $VENV/bin/activate
export PYTHONNOUSERSITE=1 PROJ_NETWORK=OFF

mkdir -p $SCRATCH/elmfire-class
cp -r $CASE $RUN
mkdir -p $RUN/outputs $RUN/figures $RUN/logs/scratch
cd $RUN

mpirun -np $SLURM_NTASKS $ELMFIRE/build/linux/bin/elmfire_1.1 elmfire.data.in \
    > logs/elmfire.stdout 2> logs/elmfire.stderr
rm -rf logs/scratch    # ~1.6 GB of ELMFIRE's converted inputs

python scripts/postprocess.py
grep -E '"(status|jaccard|precision|recall|area_ratio)"' outputs/metrics.json
