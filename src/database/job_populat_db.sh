#! /usr/bin/env bash
#SBATCH --job-name="DB_Populate"
#SBATCH -n 400
#SBATCH --constraint=sfp
#SBATCH -p sesebig
#SBATCH --time=48:00:00
#SBATCH --mem=80gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=gzhao1@illinois.edu
module load gnu/openmpi-4.1.2-gnu-9.3.0  
cd  /data/keeling/a/gzhao1/f/mmcth/src/database/
mpirun -n 400 python populate_db.py