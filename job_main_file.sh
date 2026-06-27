#!/bin/bash -l
#SBATCH -p sesebig
#SBATCH -J mmcth-file
#SBATCH -t 48:00:00
#SBATCH --ntasks=180
#SBATCH --ntasks-per-node=6
#SBATCH --constraint=m256
#SBATCH --cpus-per-task=1
#SBATCH --mem=0                 # all memory on each node
#SBATCH --exclusive
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=gzhao1@illinois.edu
#SBATCH -o slurm_create_dataset_%A_%a.out
#SBATCH -e slurm_create_dataset_%A_%a.err


##SBATCH --array=0-2             # Adjust based on how many 5-month chunks you need

source /etc/profile.d/modules.sh
module purge
module load mpi/openmpi-x86_64
source ~/miniconda3/etc/profile.d/conda.sh
conda activate k9

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
# export NUMEXPR_MAX_THREADS=4
# export NUMEXPR_NUM_THREADS=4
export HDF5_USE_FILE_LOCKING=FALSE
echo "Processing job"
#mpirun -n "$SLURM_NTASKS"  python /data/gdi/f/gzhao1/mmcth/main.py --modisid-file /data/gdi/f/gzhao1/mmcth_validate/calipso/unique_mm_files.txt --output-dir cali_output
mpirun -n "$SLURM_NTASKS"  python /data/gdi/f/gzhao1/mmcth/main.py --modisid-file /data/gdi/f/gzhao1/mmcth_validate/cats/modisid_list.txt --output-dir cats_output