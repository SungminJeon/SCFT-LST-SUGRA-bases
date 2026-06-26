#!/bin/bash
# Finalize a SUGRA pipeline output directory into the SINGLE canonical,
# deduplicated SUGRA catalog, and archive everything else.
#
#   <OUT>/final_byT/T<n>/<combo>.cat   = THE blessed SUGRA catalog (deduped)
#   <OUT>/_intermediate/               = all seeds / rounds / _nonsugra / merged
#
# Why this exists (read carefully — this is the part everyone gets wrong):
#   The pipeline leaves a ZOO of dirs. Only some are final SUGRA blocks:
#     KEEP  cat_1ext_*            phase1 seeds (SUGRA)
#     KEEP  cat_nhc_ext_*         nhc seeds + ALL chain rounds (rN sext & nsext)
#     DROP  *_nonsugra            intermediates that feed the next round's ns lane
#     DROP  cat_phase1_nhc_merged_*   redundant (= phase1 + nhc seed, phase2 input)
#   TRAP: 'rN nsext' (no _nonsugra suffix) is a SUGRA dir (non-SUGRA->SUGRA
#         recovery lane); KEEP it. Only the literal *_nonsugra dirs are dropped.
#   The s-lane (sext) and ns-recovery-lane (nsext) reach some of the SAME bases,
#   so they MUST be deduped together per (T,combo) -- consolidate_dedup.py does
#   exactly that. Plain concatenation leaves cross-stream duplicates.
#
# T is taken from each entry's PHYSICS anomaly.T (--byphys), so this works
# the same for flat single-OUT layouts and per-T (byT/ + chain/) layouts.
#
# Usage: ./finalize_sugra.sh [--archive] <OUT_DIR> [SRC_GLOB ...]
#   With no SRC_GLOB, the source dir set is auto-detected:
#     - per-T layout  (OUT/byT + OUT/chain) -> OUT/byT/T*  OUT/chain/T*
#     - flat layout                          -> OUT/cat_1ext_*  OUT/cat_nhc_ext_*
#   Pass SRC_GLOB(s) explicitly for non-standard layouts (e.g. chain_su2/), so
#   mutually-exclusive experiment variants in one master dir are NOT merged.
#
#   --archive : after producing final_byT, MOVE every other top-level entry of
#               OUT into OUT/_intermediate/. Use ONLY for throwaway flat-pipeline
#               zoos. Do NOT use on per-T family masters: byT/ and chain/ are the
#               durable, .done-resumable slices and must stay in place.
#   Without --archive (default) the sources are left untouched; final_byT is just
#   produced alongside. finalize is then idempotent and safe to re-run anytime.
set -u
cd "$(dirname "$0")"
SELF="$(pwd)"
ARCHIVE=0
if [ "${1:-}" = "--archive" ]; then ARCHIVE=1; shift; fi
OUT="${1:?usage: finalize_sugra.sh [--archive] OUT_DIR [SRC_GLOB ...]}"
shift
[ -d "$OUT" ] || { echo "ERROR: no such dir: $OUT"; exit 1; }

FINAL="$OUT/final_byT"
INT="$OUT/_intermediate"

# --- Pick the SUGRA-clean input dir set (consolidate_dedup auto-skips *_nonsugra) ---
shopt -s nullglob
GLOBS=()
if [ "$#" -gt 0 ]; then
    # explicit source globs (already shell-expanded by the caller)
    for d in "$@"; do [ -e "$d" ] && GLOBS+=("$d"); done
elif [ -d "$OUT/byT" ] || [ -d "$OUT/chain" ]; then
    # per-T family layout (run_so7hat2_byT.sh etc.): byT/T* seeds + chain/T* rounds
    for d in "$OUT"/byT/T* "$OUT"/chain/T*; do GLOBS+=("$d"); done
else
    # flat pipeline layout (run_pipeline*.sh): cat_1ext_* + cat_nhc_ext_*
    #   (cat_phase1_nhc_merged_* is NOT matched by either glob -> correctly excluded)
    for d in "$OUT"/cat_1ext_* "$OUT"/cat_nhc_ext_*; do
        case "$d" in *_nonsugra) continue;; esac
        GLOBS+=("$d")
    done
fi

if [ "${#GLOBS[@]}" -eq 0 ]; then
    echo "ERROR: no SUGRA source dirs found under $OUT"; exit 1
fi

echo "=== finalize: consolidate_dedup (--byphys) over ${#GLOBS[@]} source dirs ==="
rm -rf "$FINAL"
python3 "$SELF/consolidate_dedup.py" --byphys "$FINAL" "${GLOBS[@]}"
rc=$?
[ "$rc" -eq 0 ] || { echo "ERROR: consolidate_dedup failed (rc=$rc) — NOT archiving."; exit 1; }

# --- Archive everything except final_byT / _intermediate into _intermediate/ ---
if [ "$ARCHIVE" -eq 1 ]; then
    echo ""
    echo "=== finalize: archiving intermediates -> $INT/ ==="
    mkdir -p "$INT"
    for e in "$OUT"/*; do
        b="$(basename "$e")"
        case "$b" in final_byT|_intermediate) continue;; esac
        mv "$e" "$INT/"
    done
else
    echo ""
    echo "=== finalize: sources left in place (no --archive) ==="
fi

n_combo=$(find "$FINAL" -name '*.cat' | wc -l | xargs)
n_entry=$(find "$FINAL" -name '*.cat' -exec grep -h '^ENTRY' {} + 2>/dev/null | wc -l | tr -d ' ')
echo ""
echo "=== DONE ==="
echo "  canonical SUGRA catalog : $FINAL/  ($n_combo combo files, $n_entry entries)"
echo "  archived intermediates  : $INT/"
