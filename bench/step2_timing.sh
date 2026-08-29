#!/bin/bash
# Step 2 timing: base (pre-plan master) vs new (Step1+Step2), interleaved,
# single-threaded, 3 runs/mode, chr21, under limitram.
set -u
cd "$(dirname "$0")/.."
export LD_LIBRARY_PATH=/data/proj/2bfpbwt/htslib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
CHR21=/data/proj/2bfpbwt/exp/phase3/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf

for i in 1 2 3; do
  for m in linear sampled blockpar stagpar; do
    for b in bench/sp-pbwt-bcf.base ./sp-pbwt-bcf; do
      echo "=== $m $b run $i ==="
      limitram 16G "$b" "$m" "$CHR21" 2>&1 >/dev/null | grep 'time:'
    done
  done
done
echo "=== STEP2 TIMING DONE ==="
