#!/bin/bash
set -u
cd /data/proj/2bfpbwt/sp-pbwt
export LD_LIBRARY_PATH=/data/proj/2bfpbwt/htslib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
CHR21=/data/proj/2bfpbwt/exp/phase3/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf
for i in 1 2 3 4 5 6 7; do
  for b in bench/sp-pbwt-bcf.base ./sp-pbwt-bcf; do
    echo "=== linear $b run $i ==="
    limitram 16G "$b" linear "$CHR21" 2>&1 >/dev/null | grep 'linc('
  done
done
echo DONE_EXTRA
