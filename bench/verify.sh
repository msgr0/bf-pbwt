#!/bin/bash
# Step-gate for the linear/sampled BCF optimization plan.
#
# Rebuilds sp-pbwt-bcf and re-dumps all four BCF modes (linear, sampled,
# blockpar, stagpar) on both bench/panel.small.bcf and bench/panel.mid.bcf,
# then byte-compares every dump against the pre-change reference in
# bench/ref/ (produced from bench/sp-pbwt-bcf.base, the frozen master build).
#
# Usage: bench/verify.sh
# Exit status: 0 iff build succeeded and every dump is byte-identical.

set -u

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

export HTSLIB="${HTSLIB:-/data/proj/2bfpbwt/htslib}"
export LD_LIBRARY_PATH="$HTSLIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

echo "== build =="
HTSLIB="$HTSLIB" make sp-pbwt-bcf || { echo "BUILD FAILED"; exit 1; }

MODES="linear sampled blockpar stagpar"
PANELS="small mid"
OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
export OMP_NUM_THREADS

fail=0
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

for mode in $MODES; do
  for panel in $PANELS; do
    ref="bench/ref/$mode.$panel.dump"
    if [ ! -f "$ref" ]; then
      echo "MISSING REF: $ref"
      fail=1
      continue
    fi
    out="$tmpdir/$mode.$panel.dump"
    ./sp-pbwt-bcf "$mode" "bench/panel.$panel.bcf" DUMP 2>/dev/null > "$out"
    if cmp -s "$out" "$ref"; then
      echo "OK   $mode.$panel"
    else
      echo "DIFF $mode.$panel  <-- output changed vs reference"
      fail=1
    fi
  done
done

if [ "$fail" -eq 0 ]; then
  echo "== ALL BIT-IDENTICAL =="
else
  echo "== VERIFY FAILED =="
fi
exit $fail
