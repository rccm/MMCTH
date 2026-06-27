#!/bin/bash -l
#SBATCH -p sesebig
#SBATCH -J mmcth-by-orbit
#SBATCH -t 48:00:00
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=0                 # all memory on each node
#SBATCH --exclusive
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=gzhao1@illinois.edu
#SBATCH -o slurm_create_dataset_%A_%a.out
#SBATCH -e slurm_create_dataset_%A_%a.err
 

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
export OMPI_MCA_tmpdir_base=/data/gdi/f/gzhao1/tmp

 # List of orbits to process, sequentially
for ORBIT in 82545; do
    echo "=== Processing orbit ${ORBIT} ==="
    mpirun -n "$SLURM_NTASKS" python /data/gdi/f/gzhao1/mmcth/main.py -o "${ORBIT}"
done