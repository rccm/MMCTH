#!/bin/bash -l
#SBATCH -p sesebig
#SBATCH -J mmcth-time-chunk
#SBATCH -t 20:00:00
#SBATCH --ntasks=224
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --constraint=m256
#SBATCH --exclude=keeling-g17 #keeling-g18,keeling-j22
#SBATCH --mem=0                 # all memory on each node
#SBATCH --exclusive
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=gzhao1@illinois.edu
#SBATCH -o slurm_create_dataset_%A_%a.out
#SBATCH -e slurm_create_dataset_%A_%a.err
#SBATCH --array=0-1           # Adjust based on how many 5-month chunks you need

source /etc/profile.d/modules.sh
module purge
module load mpi/openmpi-x86_64
source ~/miniconda3/etc/profile.d/conda.sh
conda activate k9

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
# export NUMEXPR_MAX_THREADS=1
# export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

PROJECT_DIR="/data/gdi/f/gzhao1/mmcth"
MAIN_SCRIPT="$PROJECT_DIR/main.py"
OUTPUT_ROOT="$PROJECT_DIR/ds_output"

if [[ ! -r "$MAIN_SCRIPT" ]]; then
  echo "ERROR: main.py is not readable at $MAIN_SCRIPT on $(hostname)" >&2
  exit 1
fi

# Calculate start and end dates for each job array task
BASE_DATE="2015-01-01"
MONTHS_PER_CHUNK=6

# Calculate start date for this array task
START_DATE=$(date -d "$BASE_DATE + $((SLURM_ARRAY_TASK_ID * MONTHS_PER_CHUNK)) months" +%Y-%m-%d)
OUTPUT_YEAR="${START_DATE%%-*}"
OUTPUT_DIR="$OUTPUT_ROOT/$OUTPUT_YEAR"
mkdir -p "$OUTPUT_DIR"

# Calculate end date for this array task (5 months later, minus 1 day to avoid overlap)
END_DATE=$(date -d "$START_DATE + $MONTHS_PER_CHUNK months - 1 day" +%Y-%m-%d)

echo "Job array task: $SLURM_ARRAY_TASK_ID"
echo "Processing data from: $START_DATE to: $END_DATE"
echo "Output directory: $OUTPUT_DIR"
#srun --mpi=pmi2  python /data/gdi/f/gzhao1/mmcth/main.py -s $START_DATE -e $END_DATE
# mpirun -n 90  python /data/gdi/f/gzhao1/mmcth/main.py -s $START_DATE -e $END_DATE
echo "$SLURM_NTASKS" 
mpirun -np "$SLURM_NTASKS" \
  --mca orte_base_help_aggregate 0 \
  python "$MAIN_SCRIPT" -s "$START_DATE" -e "$END_DATE" --output-dir "$OUTPUT_DIR"
