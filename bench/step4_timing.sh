#!/bin/bash
# Step 4 (rrsortx experiments) interleaved timing sweep:
# bench/sp-pbwt-bcf.step3 (Step-3 baseline) vs the current ./sp-pbwt-bcf build,
# chr21, `sampled` mode only, single-threaded.
#
# Usage: bench/step4_timing.sh [runs]
set -u
cd "$(dirname "$0")/.."
export LD_LIBRARY_PATH=/data/proj/2bfpbwt/htslib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
CHR21=/data/proj/2bfpbwt/exp/phase3/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf
RUNS="${1:-3}"

for i in $(seq 1 "$RUNS"); do
  for b in bench/sp-pbwt-bcf.step3 ./sp-pbwt-bcf; do
    echo "=== sampled $b run $i ==="
    limitram 16G "$b" sampled "$CHR21" 2>&1 >/dev/null | grep 'time:'
  done
done
echo "=== STEP4 TIMING DONE ==="
