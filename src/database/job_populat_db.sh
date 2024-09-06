#!/usr/bin/env bash
#SBATCH --array=0-4               # Array with 5 jobs, index from 0 to 4
#SBATCH -n 400
#SBATCH --constraint=sfp
#SBATCH -p sesebig
#SBATCH --time=48:00:00
#SBATCH --mem=80gb
#SBATCH --job-name="DB_Populate_{%A}_{%a}"
#SBATCH --output="populate_{%A}_{%a}.out"
#SBATCH --error="populate_{%A}_{%a}.err"
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=gzhao1@illinois.edu

module load gnu/openmpi-4.1.2-gnu-9.3.0  
cd /data/keeling/a/gzhao1/f/mmcth/src/database/

# Define an array of years corresponding to the job array indices
years=(2000 2001 2002 2003 2004)

# Use SLURM_ARRAY_TASK_ID to get the current year based on the array index
year=${years[$SLURM_ARRAY_TASK_ID]}

mpirun -n 400 python populate_db.py $year
