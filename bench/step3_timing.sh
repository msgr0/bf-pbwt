#!/bin/bash
set -u
cd "$(dirname "$0")/.."
export LD_LIBRARY_PATH=/data/proj/2bfpbwt/htslib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
CHR21=/data/proj/2bfpbwt/exp/phase3/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf

for i in 1 2 3; do
  for b in bench/sp-pbwt-bcf.base ./sp-pbwt-bcf; do
    echo "=== sampled $b run $i ==="
    limitram 16G "$b" sampled "$CHR21" 2>&1 >/dev/null | grep 'time:'
  done
done
echo "=== STEP3 TIMING DONE ==="
