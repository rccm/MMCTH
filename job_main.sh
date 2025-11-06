#!/bin/bash -l
#SBATCH -p sesebig
#SBATCH -J mmcth
#SBATCH -t 2:00:00
#SBATCH --ntasks=20
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=0                 # all memory on each node
#SBATCH --exclusive
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=gzhao1@illinois.edu
#SBATCH -o slurm_create_dataset.out
#SBATCH -e slurm_create_dataset.err

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

srun --mpi=pmix_v3 python /data/gdi/f/gzhao1/mmcth/main.py -y 2003 -d 01-02