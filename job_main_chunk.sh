#!/bin/bash -l
#SBATCH -p sesebig
#SBATCH -J mmcth
#SBATCH -t 48:00:00
#SBATCH --ntasks=90
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=1
#SBATCH --mem=0                 # all memory on each node
#SBATCH --exclusive
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=gzhao1@illinois.edu
#SBATCH -o slurm_create_dataset_%A_%a.out
#SBATCH -e slurm_create_dataset_%A_%a.err
#SBATCH --array=0-2             # Adjust based on how many 5-month chunks you need

module purge
module load mpi/openmpi-x86_64
source ~/miniconda3/etc/profile.d/conda.sh
conda activate k9

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=4
export NUMEXPR_NUM_THREADS=4
export HDF5_USE_FILE_LOCKING=FALSE

# Calculate start and end dates for each job array task
BASE_DATE="2015-01-01"
MONTHS_PER_CHUNK=12

# Calculate start date for this array task
START_DATE=$(date -d "$BASE_DATE + $((SLURM_ARRAY_TASK_ID * MONTHS_PER_CHUNK)) months" +%Y-%m-%d)

# Calculate end date for this array task (5 months later, minus 1 day to avoid overlap)
END_DATE=$(date -d "$START_DATE + $MONTHS_PER_CHUNK months - 1 day" +%Y-%m-%d)

echo "Job array task: $SLURM_ARRAY_TASK_ID"
echo "Processing data from: $START_DATE to: $END_DATE"
#srun --mpi=pmi2  python /data/gdi/f/gzhao1/mmcth/main.py -s $START_DATE -e $END_DATE
mpirun -n 90  python /data/gdi/f/gzhao1/mmcth/main.py -s $START_DATE -e $END_DATE