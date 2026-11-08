#!/bin/bash
#SBATCH --job-name=elmfire-tut01
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --partition=serc
#SBATCH --mail-type=ALL
#SBATCH --output=elmfire-tut01-%j.out

module purge
module load system gcc/12.4.0 openmpi/5.0.5
module load physics gdal/3.10.2 python/3.12.1 py-gdal-utils/3.10.2_py312

ELMFIRE=/home/groups/gorle/codes/elmfire
export PATH=$ELMFIRE/build/linux/bin:$PATH
export ELMFIRE_VER=$(cat $ELMFIRE/VERSION)

TUT=/oak/stanford/groups/gorle/jgoell/Wildfire/tutorials
RESULTS=/oak/stanford/groups/gorle/jgoell/Wildfire/results/tutorial01
RUN=$SCRATCH/elmfire/tutorials

# Run on scratch; the tutorial wipes and rebuilds its own inputs/ and outputs/
mkdir -p $RUN
cp -r $TUT/functions $TUT/01-constant-wind $RUN/
cd $RUN/01-constant-wind
./01-run.sh

mkdir -p $RESULTS
cp -r outputs/* $RESULTS/
ls -l $RESULTS
