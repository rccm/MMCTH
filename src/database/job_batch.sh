#!/usr/bin/env bash

# List of years to process
years=("2017" "2018" "2019" "2020" "2021")  # Add all the years you want to process
years=("2021" "2022") 
# Loop through each year and submit a job
for year in "${years[@]}"; do

    sbatch <<EOF
#!/usr/bin/env bash
#SBATCH -n 350
#SBATCH --constraint=sfp
#SBATCH -p sesebig
#SBATCH --time=48:00:00
#SBATCH --mem=130gb
#SBATCH --job-name="DB_Populate_$year"
#SBATCH --output="populate_$year.out"
#SBATCH --error="populate_$year.err"
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=gzhao1@illinois.edu

module load gnu/openmpi-4.1.2-gnu-9.3.0  
cd /data/keeling/a/gzhao1/f/mmcth/src/database/

mpirun -n 350 python populate_db.py $year
EOF

done