#!/bin/bash
#SBATCH --job-name=gpu_job         # Job name
#SBATCH --partition=gpu,sesempi            # Partition (queue) to submit to
#SBATCH --gres=gpu:1               # Request 1 GPU
#SBATCH --nodes=1                  # Request 1 node
#SBATCH --ntasks=1                 # Run 1 task (process)
#SBATCH --time=01:00:00            # Maximum runtime of 1 hour
#SBATCH --output=output.log        # Output file

# Load GPU module if necessary
module load GPU

# Run the Python script
python my_cupy_script.py