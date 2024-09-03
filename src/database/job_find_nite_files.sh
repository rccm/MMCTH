#! /usr/bin/env bash
#SBATCH --job-name="DB_Populate"
#SBATCH -n 380
#SBATCH --constraint=sfp
#SBATCH -p sesempi
#SBATCH --time=160:00:00
#SBATCH --mem=100gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=gzhao1@illinois.edu
module load gnu/openmpi-4.1.2-gnu-9.3.0  
cd  /data/keeling/a/gzhao1/f/mmcth/src/database/
mpirun -n 380 python find_night_files.py 