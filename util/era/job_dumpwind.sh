#!/usr/bin/env bash
#SBATCH --job-name="dump_monthly_wind"
#SBATCH -n 12
#SBATCH --constraint=sfp
#SBATCH -p sesempi
#SBATCH --time=56:00:00
#SBATCH --mem=200gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=gzhao1@illinois.edu
module load gnu/openmpi-3.1.6-gnu-4.8.5  
cd /data/keeling/a/gzhao1/f/mmcth/util/era/
mpirun -n 12 python  dump_monthly_wind.py
 