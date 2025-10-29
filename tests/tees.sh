#!/bin/bash
#SBATCH --job-name=stereo_monthly
#SBATCH --array=1-200%10    # Limits to 10 concurrent jobs
#SBATCH --time=60:00:00
#SBATCH -p sesebig
#SBATCH --output=output_%a.out
#SBATCH --error=error_%a.err

input_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" file_list.txt)
ech  "$(basename "${input_file}")"