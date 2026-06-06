#!/bin/bash
# Run SUGRA pipeline WITH su2 externals enabled, with explicit OpenMP thread count.
#
# Usage:   ./run_pipeline_su2_omp.sh T_MIN T_MAX [N_THREADS]
# Example: ./run_pipeline_su2_omp.sh 0 10 8     # T=0..10 with 8 threads
#          ./run_pipeline_su2_omp.sh 21 30 32   # T=21..30 with 32 threads
#
# If N_THREADS is omitted, OpenMP picks the default (typically all hardware
# threads). Set N_THREADS to a smaller value on shared servers to keep the
# load reasonable.
#
# Same pipeline steps as run_pipeline.sh — only difference is OMP_NUM_THREADS
# is exported up front.

set -e
cd "$(dirname "$0")"

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 T_MIN T_MAX [N_THREADS]"
    echo "Example: $0 0 10 8"
    exit 1
fi

T_MIN=$1
T_MAX=$2
N_THREADS=${3:-}

if [ -n "$N_THREADS" ]; then
    export OMP_NUM_THREADS="$N_THREADS"
    echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
else
    echo "OMP_NUM_THREADS not set — using OpenMP default (all hardware threads)."
fi

SUFFIX=$(printf "t%d%02d_su2" "$T_MIN" "$T_MAX")

COMMON_FLAGS="--use-lst-T"

echo "=========================================="
echo "SUGRA pipeline: T = ${T_MIN}..${T_MAX}  (suffix: ${SUFFIX})"
echo "=========================================="
date

echo ""
echo "=== Build ==="
make 2>&1 | tail -5

echo ""
echo "=== Step 1: phase1 (T=${T_MIN}..${T_MAX}) ==="
rm -rf cat_1ext_nosu2_${SUFFIX} cat_1ext_nosu2_${SUFFIX}_nonsugra
date
time ./gen_sugra_phase1 unified.cat ${T_MAX} ${T_MIN} ${COMMON_FLAGS} --save-nonsugra \
    --out cat_1ext_nosu2_${SUFFIX} 2>&1 | tail -5

echo ""
echo "=== Step 2: nhc_ext (T=${T_MIN}..${T_MAX}) ==="
rm -rf cat_nhc_ext_${SUFFIX} cat_nhc_ext_${SUFFIX}_nonsugra
rm -rf cat_nhc_ext cat_nhc_ext_nonsugra
date
time ./gen_sugra_nhc_ext unified.cat ${T_MAX} ${T_MIN} ${COMMON_FLAGS} --save-nonsugra 2>&1 | tail -5
if [ -d cat_nhc_ext ]; then mv cat_nhc_ext cat_nhc_ext_${SUFFIX}; fi
if [ -d cat_nhc_ext_nonsugra ]; then mv cat_nhc_ext_nonsugra cat_nhc_ext_${SUFFIX}_nonsugra; fi

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

OUTDIR="T_${T_MIN}_${T_MAX}_su2"
echo ""
echo "=== Step 5: collect catalogs into ${OUTDIR}/ ==="
mkdir -p "${OUTDIR}"
mv cat_1ext_nosu2_${SUFFIX} cat_1ext_nosu2_${SUFFIX}_nonsugra "${OUTDIR}/" 2>/dev/null
mv cat_nhc_ext_${SUFFIX} cat_nhc_ext_${SUFFIX}_nonsugra "${OUTDIR}/" 2>/dev/null
mv cat_phase1_nhc_merged_${SUFFIX} cat_phase1_nhc_merged_${SUFFIX}_nonsugra "${OUTDIR}/" 2>/dev/null
mv cat_nhc_ext_nhc_ext_unified_${SUFFIX}_* "${OUTDIR}/" 2>/dev/null
echo "  moved to ${OUTDIR}/  ($(ls "${OUTDIR}" | wc -l | xargs) directories)"

echo ""
echo "=== Output directories (round-by-round) ==="
ls -d "${OUTDIR}"/cat_nhc_ext_nhc_ext_unified_${SUFFIX}_r* 2>/dev/null | while read d; do
    n=$(find "$d" -name '*.cat' -exec grep -c '^ENTRY' {} + 2>/dev/null | awk -F: '{s+=$2}END{print s+0}')
    echo "  $(basename "$d"): $n entries"
done
