#!/bin/bash
#SBATCH --job-name=elmfire-camp
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --time=2-00:00:00
#SBATCH --mem=128G
#SBATCH --mail-type=ALL
#SBATCH --output=elmfire-camp-%j.out

# 2018 Camp Fire validation case from ELMFIRE_VnV_Suite: 50 ensemble members,
# one MPI rank each, all on one node (ranks on a node share the weather rasters).
# Runs a copy on scratch, then copies outputs, figures, report and ELMFIRE logs
# to Wildfire/results/camp_fire.
#
#   cd Wildfire/slurm && sbatch elmfire_camp_fire.slurm
#
# Time and memory are first guesses -- the earlier (Savio) run recorded no
# runtime. Check seff <jobid> afterwards.

VAL=/oak/stanford/groups/gorle/jgoell/ELMFIRE_VnV_Suite/cases/Validation
RESULTS=/oak/stanford/groups/gorle/jgoell/Wildfire/results/camp_fire
RUN=$SCRATCH/elmfire/camp-$SLURM_JOB_ID
CASE=$RUN/landscape_scale/camp_fire

module purge
module load system gcc/12.4.0 openmpi/5.0.5
module load physics gdal/3.10.2 python/3.12.1 texlive/2019
source /oak/stanford/groups/gorle/jgoell/venvs/elmfire-vnv/bin/activate

export ELMFIRE_BIN=/home/groups/gorle/codes/elmfire/build/linux/bin/elmfire_1.1
export ELMFIRE_MPI_RANKS=$SLURM_NTASKS

# preprocess.py checksums ../../3_FIRES/Camp_VAL.rar, so keep that layout.
# Leave out the old run's outputs/figures/report PDF so nothing stale survives.
mkdir -p $RUN/landscape_scale $RUN/3_FIRES
rsync -a --exclude=outputs --exclude=figures --exclude=case_report.pdf \
    $VAL/landscape_scale/camp_fire $RUN/landscape_scale/
cp $VAL/3_FIRES/Camp_VAL.rar $RUN/3_FIRES/

$CASE/run_case.sh

mkdir -p $RESULTS
cd $CASE
rsync -a --exclude='*.bil' --exclude='*.hdr' --exclude=scratch/ outputs figures report logs $RESULTS/
grep -E '"(status|jaccard|precision|recall|area_ratio)"' $RESULTS/outputs/metrics.json
