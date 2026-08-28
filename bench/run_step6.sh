#!/bin/bash
# Step 6 final benchmark: all four BCF modes on chr20 and chr10, under limitram,
# capped at 32 threads. Meant to be run detached (babysit-run).
set -uo pipefail

cd "$(dirname "$0")/.."

CHR20=/data/proj/2bfpbwt/one_thousands/phase3/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf
CHR10=/data/proj/2bfpbwt/one_thousands/phase3/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf
MEMLIMIT="${MEMLIMIT:-16G}"
OMP_NUM_THREADS=32
export OMP_NUM_THREADS

OUT=bench/step6_results.txt
: > "$OUT"

for pair in "chr20:$CHR20" "chr10:$CHR10"; do
  name="${pair%%:*}"
  path="${pair#*:}"
  echo "===== $name =====" | tee -a "$OUT"
  BCF_INPUT="$path" MODES="linear sampled blockpar stagpar" OMP_NUM_THREADS=32 \
    limitram "$MEMLIMIT" bash bench/run.sh 2>&1 | tee -a "$OUT"
  echo "" | tee -a "$OUT"
done

echo "DONE" | tee -a "$OUT"
