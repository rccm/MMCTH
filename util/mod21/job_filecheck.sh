#!/usr/bin/env bash
#SBATCH -n 200
#SBATCH --constraint=sfp
#SBATCH -p sesebig
#SBATCH --time=46:00:00
#SBATCH --mem=20gb
#SBATCH --job-name="filechecking_mod21"
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=gzhao1@illinois.edu
module load gnu/openmpi-4.1.2-gnu-9.3.0  
cd /data/keeling/a/gzhao1/f/mmcth/util/mod21
mpirun -n 200 python  filecheck.py
 