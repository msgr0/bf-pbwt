#!/bin/bash
# Step 0.3: baseline timing for the linear/sampled BCF optimization plan.
# 3 runs each of linear/sampled on chr21 and chr20, single-threaded.
set -u
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"
export LD_LIBRARY_PATH="/data/proj/2bfpbwt/htslib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1

CHR20=/data/proj/2bfpbwt/exp/phase3/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf
CHR21=/data/proj/2bfpbwt/exp/phase3/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf

OUT="bench/baseline_timings.txt"
: > "$OUT"

for chr_name in chr21 chr20; do
  bcf_var="CHR${chr_name#chr}"
  bcf="${!bcf_var}"
  for mode in linear sampled; do
    for i in 1 2 3; do
      echo "=== $mode $chr_name run $i ===" | tee -a "$OUT"
      limitram 16G ./bench/sp-pbwt-bcf.base "$mode" "$bcf" 2>&1 >/dev/null | tee -a "$OUT"
    done
  done
done

echo "DONE" | tee -a "$OUT"
