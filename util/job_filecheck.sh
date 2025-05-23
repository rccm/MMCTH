#!/usr/bin/env bash
#SBATCH -n 200
#SBATCH --constraint=sfp
#SBATCH -p sesempi
#SBATCH --time=46:00:00
#SBATCH --mem=20gb
#SBATCH --job-name="filechecking"
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=gzhao1@illinois.edu
module unload openmpi
module load gnu/openmpi-3.1.6-gnu-4.8.5  
cd /data/keeling/a/gzhao1/f/mmcth/util
mpirun -n 200 python  filecheck.py
 