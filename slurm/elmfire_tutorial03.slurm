#!/bin/bash
#SBATCH --job-name=elmfire-tut03
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --partition=serc
#SBATCH --mail-type=ALL
#SBATCH --output=elmfire-tut03-%j.out

module purge
module load system gcc/12.4.0 openmpi/5.0.5
module load physics gdal/3.10.2 python/3.12.1 py-gdal-utils/3.10.2_py312

# grpcio/protobuf for the cloudfire fuel download
source /oak/stanford/groups/gorle/jgoell/venvs/cloudfire/bin/activate

export ELMFIRE_BASE_DIR=/home/groups/gorle/codes/elmfire
export PATH=$ELMFIRE_BASE_DIR/build/linux/bin:$PATH
export ELMFIRE_VER=$(cat $ELMFIRE_BASE_DIR/VERSION)

TUT=/oak/stanford/groups/gorle/jgoell/Wildfire/tutorials
CASE=/oak/stanford/groups/gorle/jgoell/Wildfire/results/03-real-fuels
RUN=$SCRATCH/elmfire/tutorials

# Run on scratch; the tutorial wipes and rebuilds its own inputs/ and outputs/
mkdir -p $RUN
cp -r $TUT/functions $CASE $RUN/
cd $RUN/03-real-fuels
./01-run.sh

cp -r outputs $CASE/
ls -l $CASE/outputs
