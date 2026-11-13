#!/bin/bash -l
#SBATCH --partition=node
#SBATCH --job-name=mmcth-json
#SBATCH --array=0-22%2
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=slurm_metadata_json_%A_%a.out
#SBATCH --error=slurm_metadata_json_%A_%a.err

# Submit all yearly directories: sbatch job_metadata_json.sh
# Tasks 0-22 map to years 2000-2022 (cluster MaxArraySize is 1001).
# One year (2011): sbatch --array=11 job_metadata_json.sh
# Small trial: LIMIT=2 sbatch --array=11 job_metadata_json.sh
# Regenerate existing sidecars: OVERWRITE=1 sbatch job_metadata_json.sh
# Optional environment settings: PROJECT_DIR, OUTPUT_ROOT, MOD03_ROOT,
# METADATA_SCRIPT, WORKERS, PATTERN, COLLECTION_SHORTNAME, COLLECTION_VERSION.

set -eo pipefail

source /etc/profile.d/modules.sh
module purge
source ~/miniconda3/etc/profile.d/conda.sh
conda activate k9

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

PROJECT_DIR="${PROJECT_DIR:-/data/keeling/a/gzhao1/f/mmcth}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$PROJECT_DIR/ds_output}"
METADATA_SCRIPT="${METADATA_SCRIPT:-$PROJECT_DIR/metadata/creat_json.py}"
MOD03_ROOT="${MOD03_ROOT:-/data/gdi/satellite/TerraDataArchive/MODIS/MOD03}"
OUTPUT_YEAR="${OUTPUT_YEAR:-}"
WORKERS="${WORKERS:-${SLURM_CPUS_PER_TASK:-4}}"

if [[ -z "$OUTPUT_YEAR" && "${SLURM_ARRAY_TASK_ID:-}" =~ ^[0-9]+$ ]]; then
  OUTPUT_YEAR=$((2000 + 10#${SLURM_ARRAY_TASK_ID}))
fi
if [[ ! "$OUTPUT_YEAR" =~ ^[0-9]{4}$ ]]; then
  echo "ERROR: submit as array tasks 0-22 (years 2000-2022), or set OUTPUT_YEAR" >&2
  exit 1
fi

LOG_DIR="$PROJECT_DIR/logs/metadata_json/slurm_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}_$OUTPUT_YEAR"
mkdir -p "$LOG_DIR"

args=(
  --input-dir "$OUTPUT_ROOT/$OUTPUT_YEAR"
  --metadata-script "$METADATA_SCRIPT"
  --mod03-root "$MOD03_ROOT"
  --workers "$WORKERS"
  --pattern "${PATTERN:-MODMISR_L2_CP*.nc}"
  --collection-shortname "${COLLECTION_SHORTNAME:-MODMISR_L2_CP}"
  --collection-version "${COLLECTION_VERSION:-V00}"
  --log-dir "$LOG_DIR"
)
if [[ "${OVERWRITE:-0}" == 1 ]]; then
  args+=(--overwrite)
fi
if [[ -n "${LIMIT:-}" ]]; then
  args+=(--limit "$LIMIT")
fi

echo "Creating JSON metadata for $OUTPUT_ROOT/$OUTPUT_YEAR"
echo "Generator: $METADATA_SCRIPT"
echo "Sidecars will be written beside the original NetCDF files."
python -u "$PROJECT_DIR/metadata/create_json_batch.py" "${args[@]}"
