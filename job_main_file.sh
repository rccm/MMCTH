#!/bin/bash -l
#SBATCH -p sesebig
#SBATCH -J mmcth-file
#SBATCH -t 48:00:00
#SBATCH --ntasks=170
#SBATCH --ntasks-per-node=6
#SBATCH --constraint=m256
#SBATCH --cpus-per-task=1
#SBATCH --mem=0                 # all memory on each node
#SBATCH --exclusive
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=gzhao1@illinois.edu
#SBATCH -o slurm_create_dataset_%A_%a.out
#SBATCH -e slurm_create_dataset_%A_%a.err


set -eo pipefail

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

PROJECT_DIR="${PROJECT_DIR:-/data/gdi/f/gzhao1/mmcth}"
MAIN_SCRIPT="${MAIN_SCRIPT:-$PROJECT_DIR/main.py}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$PROJECT_DIR/ds_output}"
MODISID_FILE="${MODISID_FILE:-/data/gdi/f/gzhao1/mmcth_validate/diagnostics/source_gap_causes_20261028/catalog_repair_applied_20130318_20261028T141629Z/repaired_modis_ids_2013.txt}"
: "${SLURM_NTASKS:?Submit this script with sbatch}"

if [[ ! -r "$MAIN_SCRIPT" || ! -r "$MODISID_FILE" ]]; then
  echo "ERROR: MAIN_SCRIPT and MODISID_FILE must be readable: $MAIN_SCRIPT, $MODISID_FILE" >&2
  exit 1
fi

# Keep temporary ID lists on shared storage so every MPI rank can read them.
mkdir -p "$OUTPUT_ROOT"
YEAR_LIST_DIR=$(mktemp -d "$OUTPUT_ROOT/.modisid-rerun.XXXXXX")
trap 'rm -rf -- "$YEAR_LIST_DIR"' EXIT

# A rerun list can contain several years; derive output years from its IDs.
python - "$MODISID_FILE" "$YEAR_LIST_DIR" <<'PY'
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

source, destination = map(Path, sys.argv[1:])
by_year = defaultdict(set)
for number, line in enumerate(source.read_text().splitlines(), start=1):
    modis_id = line.strip()
    if not modis_id or modis_id.startswith("#"):
        continue
    match = re.fullmatch(r"\.?A(\d{4})(\d{3})\.(\d{4})\.?", modis_id)
    if match is None:
        raise SystemExit(f"Invalid MODIS ID at {source}:{number}: {modis_id}")
    year, doy, hhmm = match.groups()
    try:
        date = datetime.strptime(year + doy + hhmm, "%Y%j%H%M")
        if date.year != int(year):
            raise ValueError("day of year is outside this year")
    except ValueError as exc:
        raise SystemExit(f"Invalid MODIS ID at {source}:{number}: {modis_id}: {exc}")
    by_year[year].add(f".A{year}{doy}.{hhmm}.")

if not by_year:
    raise SystemExit(f"No MODIS IDs found in {source}")
for year, ids in sorted(by_year.items()):
    (destination / f"{year}.txt").write_text("\n".join(sorted(ids)) + "\n")
(destination / "years.txt").write_text("\n".join(sorted(by_year)) + "\n")
PY

while IFS= read -r OUTPUT_YEAR; do
  YEAR_ID_FILE="$YEAR_LIST_DIR/$OUTPUT_YEAR.txt"
  OUTPUT_DIR="$OUTPUT_ROOT/$OUTPUT_YEAR"
  mkdir -p "$OUTPUT_DIR"
  echo "Processing $(wc -l < "$YEAR_ID_FILE") MODIS IDs for $OUTPUT_YEAR"
  echo "Output directory: $OUTPUT_DIR"
  # MPI must not consume the remaining years from the loop's standard input.
  mpirun -n "$SLURM_NTASKS" python "$MAIN_SCRIPT" \
    --modisid-file "$YEAR_ID_FILE" --output-dir "$OUTPUT_DIR" < /dev/null
done < "$YEAR_LIST_DIR/years.txt"
