#!/bin/bash
set -u
cd "$(dirname "$0")/.."
export LD_LIBRARY_PATH=/data/proj/2bfpbwt/htslib:${LD_LIBRARY_PATH:-}
CHR21=/data/proj/2bfpbwt/exp/phase3/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf
CHR20=/data/proj/2bfpbwt/exp/phase3/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf

# linear/sampled: single-threaded, per plan's Step 0 protocol
for mode in linear sampled; do
  for chr in CHR21 CHR20; do
    f=${!chr}
    for i in 1 2 3; do
      for b in bench/sp-pbwt-bcf.base bench/sp-pbwt-bcf.final; do
        echo "=== $mode $chr $b run $i ==="
        OMP_NUM_THREADS=1 limitram 16G "$b" "$mode" "$f" 2>&1 >/dev/null | grep 'time:'
      done
    done
  done
done

# blockpar/stagpar: 8 threads, matching the manuscript/exp convention (bpr/spr)
for mode in blockpar stagpar; do
  for chr in CHR21 CHR20; do
    f=${!chr}
    for i in 1 2 3; do
      for b in bench/sp-pbwt-bcf.base bench/sp-pbwt-bcf.final; do
        echo "=== $mode $chr $b run $i ==="
        OMP_NUM_THREADS=8 limitram 16G "$b" "$mode" "$f" 2>&1 >/dev/null | grep 'time:'
      done
    done
  done
done
echo "=== STEP6 TIMING DONE ==="
