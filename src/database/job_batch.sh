#!/usr/bin/env bash

# List of years to process
years=("2000" "2001" "2002" "2003" "2004")  # Add all the years you want to process

# Loop through each year and submit a job
for year in "${years[@]}"; do

    sbatch <<EOF
#!/usr/bin/env bash
#SBATCH -n 400
#SBATCH --constraint=sfp
#SBATCH -p sesebig
#SBATCH --time=48:00:00
#SBATCH --mem=80gb
#SBATCH --job-name="DB_Populate_$year"
#SBATCH --output="populate_$year.out"
#SBATCH --error="populate_$year.err"
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=gzhao1@illinois.edu

module load gnu/openmpi-4.1.2-gnu-9.3.0  
cd /data/keeling/a/gzhao1/f/mmcth/src/database/

mpirun -n 400 python populate_db.py $year
EOF

done