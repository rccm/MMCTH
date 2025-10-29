#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --partition=sesempi
#SBATCH --time=55:00:00
#SBATCH --nodes=8               # sesebig MinNodes=9
#SBATCH --ntasks=60              # total MPI ranks
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G             # total MPI ranks
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=gzhao1@illinois.edu

module load mpi/openmpi-x86_64
source ~/miniconda3/etc/profile.d/conda.sh
conda activate k9
export HDF5_USE_FILE_LOCKING=FALSE
cd /data/gdi/f/gzhao1/mmcth/util/era
mpirun -n 60 python merge_era5.py