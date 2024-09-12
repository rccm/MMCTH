#!/usr/bin/env bash
#SBATCH -n 200
#SBATCH --constraint=sfp
#SBATCH -p sesempi
#SBATCH --time=14:00:00
#SBATCH --mem=20gb
#SBATCH --job-name="filechecking"
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=gzhao1@illinois.edu
module load gnu/openmpi-4.1.2-gnu-9.3.0  
cd /data/keeling/a/gzhao1/f/mmcth/util
mpirun -n 200 python  filecheck.py
 