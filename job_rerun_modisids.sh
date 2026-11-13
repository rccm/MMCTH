#!/bin/bash -l
#SBATCH -p sesempi
#SBATCH -J mmcth-rerun
#SBATCH -t 12:00:00
# The MM retrieval can approach the full memory of an m64 node when 12
# granules reach the ERA5/retrieval stage together (job 986104).  Keep the
# existing rank layout, but place it on the project's high-memory nodes.
#SBATCH --constraint=m256
#SBATCH --ntasks=60
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=gzhao1@illinois.edu
#SBATCH -o slurm_rerun_%j.out
#SBATCH -e slurm_rerun_%j.err

source /etc/profile.d/modules.sh
module purge
module load mpi/openmpi-x86_64
source ~/miniconda3/etc/profile.d/conda.sh
conda activate k9

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

PROJECT_DIR="${PROJECT_DIR:-/data/gdi/f/gzhao1/mmcth}"
MAIN_SCRIPT="${MAIN_SCRIPT:-$PROJECT_DIR/main.py}"
: "${MODISID_FILE:?Set MODISID_FILE to the targeted rerun list}"
: "${OUTPUT_DIR:?Set OUTPUT_DIR to the existing yearly product directory}"

if [[ ! -r "$MODISID_FILE" ]]; then
  echo "ERROR: MODISID_FILE is not readable: $MODISID_FILE" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "Targeted rerun list: $MODISID_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "MPI ranks: $SLURM_NTASKS"
echo "Required node class: m256"

mpirun -np "$SLURM_NTASKS" \
  --mca orte_base_help_aggregate 0 \
  python "$MAIN_SCRIPT" \
  --modisid-file "$MODISID_FILE" \
  --output-dir "$OUTPUT_DIR" \
  --co2-debug
