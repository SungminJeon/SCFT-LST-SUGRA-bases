#!/bin/bash
# Run full SUGRA pipeline for a given tensor (T) range.
#
# Usage:   ./run_pipeline.sh T_MIN T_MAX
# Example: ./run_pipeline.sh 0 10        # T = 0..10  → outputs in cat_*_t010
#          ./run_pipeline.sh 11 20       # T = 11..20 → outputs in cat_*_t1120
#          ./run_pipeline.sh 21 30       # T = 21..30 → outputs in cat_*_t2130
#
# Suffix rule: t<T_MIN><T_MAX> with no separator (always 2 digits each).
#
# Steps:
#   1. phase1     — single-curve external attach
#   2. nhc_ext    — NHC chain external attach
#   3. merge      — combine phase1 + nhc_ext for phase2 input
#   4. phase2 chain (run_chain.sh) — stack more externals through 9 rounds

set -e
cd "$(dirname "$0")"

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 T_MIN T_MAX"
    echo "Example: $0 0 10"
    exit 1
fi

T_MIN=$1
T_MAX=$2

# Build a 4-char suffix from the two values, e.g. T_MIN=0  T_MAX=10 → "010"
#                                                T_MIN=11 T_MAX=20 → "1120"
SUFFIX=$(printf "t%d%02d" "$T_MIN" "$T_MAX")

COMMON_FLAGS="--use-lst-T --no-su2"

echo "=========================================="
echo "SUGRA pipeline: T = ${T_MIN}..${T_MAX}  (suffix: ${SUFFIX})"
echo "=========================================="
date

# --- Build (no-op if up to date) ---
echo ""
echo "=== Build ==="
make 2>&1 | tail -5

# --- Step 1: Phase 1 (single-curve external attach) ---
echo ""
echo "=== Step 1: phase1 (T=${T_MIN}..${T_MAX}) ==="
rm -rf cat_1ext_nosu2_${SUFFIX} cat_1ext_nosu2_${SUFFIX}_nonsugra
date
time ./gen_sugra_phase1 unified.cat ${T_MAX} ${T_MIN} ${COMMON_FLAGS} --save-nonsugra \
    --out cat_1ext_nosu2_${SUFFIX} 2>&1 | tail -5

# --- Step 2: NHC ext (NHC chain attach) ---
echo ""
echo "=== Step 2: nhc_ext (T=${T_MIN}..${T_MAX}) ==="
rm -rf cat_nhc_ext_${SUFFIX} cat_nhc_ext_${SUFFIX}_nonsugra
# gen_sugra_nhc_ext writes to "cat_nhc_ext"; rename after
rm -rf cat_nhc_ext cat_nhc_ext_nonsugra
date
time ./gen_sugra_nhc_ext unified.cat ${T_MAX} ${T_MIN} ${COMMON_FLAGS} --save-nonsugra 2>&1 | tail -5
mv cat_nhc_ext cat_nhc_ext_${SUFFIX}
mv cat_nhc_ext_nonsugra cat_nhc_ext_${SUFFIX}_nonsugra

# --- Step 3: Merge phase1 + nhc_ext for phase2 input ---
echo ""
echo "=== Step 3: merge phase1 + nhc_ext ==="
rm -rf cat_phase1_nhc_merged_${SUFFIX} cat_phase1_nhc_merged_${SUFFIX}_nonsugra
mkdir -p cat_phase1_nhc_merged_${SUFFIX} cat_phase1_nhc_merged_${SUFFIX}_nonsugra
for f in cat_1ext_nosu2_${SUFFIX}/*.cat; do cp "$f" "cat_phase1_nhc_merged_${SUFFIX}/$(basename "$f")"; done
for f in cat_nhc_ext_${SUFFIX}/*.cat; do cp "$f" "cat_phase1_nhc_merged_${SUFFIX}/$(basename "$f")"; done
for f in cat_1ext_nosu2_${SUFFIX}_nonsugra/*.cat; do [ -f "$f" ] && cp "$f" "cat_phase1_nhc_merged_${SUFFIX}_nonsugra/$(basename "$f")"; done
for f in cat_nhc_ext_${SUFFIX}_nonsugra/*.cat; do [ -f "$f" ] && cp "$f" "cat_phase1_nhc_merged_${SUFFIX}_nonsugra/$(basename "$f")"; done
SUGRA_COUNT=$(find cat_phase1_nhc_merged_${SUFFIX} -name '*.cat' -exec grep -c '^ENTRY' {} + | awk -F: '{s+=$2}END{print s}')
NONSUGRA_COUNT=$(find cat_phase1_nhc_merged_${SUFFIX}_nonsugra -name '*.cat' -exec grep -c '^ENTRY' {} + | awk -F: '{s+=$2}END{print s}')
echo "  merged SUGRA:     $SUGRA_COUNT"
echo "  merged non-SUGRA: $NONSUGRA_COUNT"

# --- Step 4: NHC ext Phase 2 chain (multi-round external stacking) ---
echo ""
echo "=== Step 4: nhc_ext_phase2 chain (max 9 rounds) ==="
date
time ./run_chain.sh "./gen_sugra_nhc_ext_phase2" \
    cat_phase1_nhc_merged_${SUFFIX} \
    nhc_ext_unified_${SUFFIX} \
    9 \
    "${COMMON_FLAGS}" \
    cat_phase1_nhc_merged_${SUFFIX}_nonsugra 2>&1 | tail -40

echo ""
echo "=========================================="
echo "Pipeline complete (T = ${T_MIN}..${T_MAX})"
echo "=========================================="
date

# --- Summary ---
echo ""
echo "=== Output directories (round-by-round) ==="
ls -d cat_nhc_ext_nhc_ext_unified_${SUFFIX}_r* 2>/dev/null | while read d; do
    n=$(find "$d" -name '*.cat' -exec grep -c '^ENTRY' {} + 2>/dev/null | awk -F: '{s+=$2}END{print s+0}')
    echo "  $d: $n entries"
done
