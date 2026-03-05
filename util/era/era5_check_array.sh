#!/bin/bash -l
#SBATCH -p sesempi
#SBATCH -J era5mulchk
#SBATCH -t 12:00:00

#SBATCH --nodes=2              # sesempi MinNodes=2
#SBATCH --ntasks-per-node=10    # start conservative (2 ranks total)
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G

#SBATCH -o era5mulchk_%j.out
#SBATCH -e era5mulchk_%j.err

module purge
source ~/miniconda3/etc/profile.d/conda.sh
conda activate k9

LIST=era5_multi_files.txt
OUTDIR=era5_multi_check_${SLURM_JOB_ID}
mkdir -p "$OUTDIR"

echo "Job $SLURM_JOB_ID: nodes=$SLURM_JOB_NUM_NODES ntasks=$SLURM_NTASKS"
echo "List: $LIST  (N=$(wc -l < "$LIST"))"
echo "Out:  $OUTDIR"

srun bash -lc '
  rank=$SLURM_PROCID
  nranks=$SLURM_NTASKS
  outcsv="'"$OUTDIR"'/part_${rank}.csv"
  echo "file,var,where,error" > "$outcsv"

  n=$(wc -l < "'"$LIST"'")

  for ((i=1; i<=n; i++)); do
    if ((( (i-1) % nranks ) != rank)); then
      continue
    fi

    f=$(sed -n "${i}p" "'"$LIST"'")

    # This script should: open file, enumerate vars, probe reads
    /data/keeling/a/gzhao1/miniconda3/envs/k9/bin/python era5_multi_check.py "$f" --out "$outcsv" --append --probes 400
  done

  echo "Rank $rank done: $(wc -l < "$outcsv") lines"
'

FINAL="$OUTDIR/era5_multi_corruption_all.csv"
head -n 1 "$OUTDIR/part_0.csv" > "$FINAL"
tail -n +2 -q "$OUTDIR"/part_*.csv >> "$FINAL"
echo "Merged report: $FINAL"