#!/bin/bash
# Baseline benchmark harness for sp-pbwt BCF modes
# Runs each mode on a fixed input and records wall time, CPU time, peak RSS

set -e

# Defaults
MODES="${MODES:-linear sampled blockpar stagpar}"
BCF_INPUT="${BCF_INPUT:-/data/proj/2bfpbwt/one_thousands/phase3/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf}"
OMP_NUM_THREADS="${OMP_NUM_THREADS:-32}"

# Ensure BCF exists
if [ ! -f "$BCF_INPUT" ]; then
    echo "Error: BCF input not found at $BCF_INPUT" >&2
    exit 1
fi

# Go to repo root
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

# Ensure binary exists
if [ ! -f sp-pbwt-bcf ]; then
    echo "Error: sp-pbwt-bcf not found in $REPO_ROOT" >&2
    exit 1
fi

# htslib is a shared library; make sure the loader can find it (not needed
# for sp-pbwt-bm, only sp-pbwt-bcf).
export LD_LIBRARY_PATH="${HTSLIB:-/data/proj/2bfpbwt/htslib}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

echo "Benchmark harness: BCF modes on $(basename $BCF_INPUT)"
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo ""
echo "Mode         Wall(s)      CPU(s,user+sys)  Peak RSS(MB)"
echo "=================================================="

export OMP_NUM_THREADS

for mode in $MODES; do
    # Run under /usr/bin/time -v to capture wall time, CPU time, and peak RSS
    # Filter to extract key metrics
    output=$( { /usr/bin/time -v ./sp-pbwt-bcf "$mode" "$BCF_INPUT" 2>&1 1>/dev/null; } | grep -E "Elapsed|User|System|Maximum resident" || true)

    # Parse the time output
    wall=$(echo "$output" | grep "Elapsed" | awk '{print $NF}' | head -1)
    user=$(echo "$output" | grep "User" | awk '{print $4}')
    sys=$(echo "$output" | grep "System" | awk '{print $4}')
    rss=$(echo "$output" | grep "Maximum resident" | awk '{print $6}')

    cpu=$(echo "scale=2; ${user:-0} + ${sys:-0}" | bc 2>/dev/null || echo "N/A")

    # Convert RSS from KB to MB
    if [ ! -z "$rss" ]; then
        rss_mb=$(echo "scale=2; $rss / 1024" | bc)
    else
        rss_mb="N/A"
    fi

    printf "%-12s %10s %16s %15s\n" "$mode" "$wall" "$cpu" "$rss_mb"
done
