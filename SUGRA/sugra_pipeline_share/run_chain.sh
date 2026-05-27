#!/bin/bash
# Chain with non-SUGRA handling
# Usage: ./run_chain.sh <phase2_cmd> <input_dir> <tag> [max_round] [flags] [initial_nonsugra_dir]

set -e
PHASE2="$1"
INPUT="$2"
TAG="$3"
MAX=${4:-12}
FLAGS="$5"
NS_INIT="${6:-}"

# Merge cat files: append entries from src/*.cat into dst/ (same filename → concatenate)
merge_cat_dirs() {
    local src="$1" dst="$2"
    for f in "$src"/*.cat; do
        [ -f "$f" ] || continue
        local base=$(basename "$f")
        if [ -f "$dst/$base" ]; then
            # Append entries (skip header comments from src)
            grep -v '^#' "$f" >> "$dst/$base"
        else
            cp "$f" "$dst/$base"
        fi
    done
}

echo "=== Chain: $TAG (max $MAX rounds) ==="

SUGRA_IN="$INPUT"
NS_IN="$NS_INIT"

for round in $(seq 2 $MAX); do
    echo ""
    echo ">>> Round $round"

    ANY_OUTPUT=0

    # --- SUGRA chain ---
    if [ -d "$SUGRA_IN" ]; then
        S_N=$(find "$SUGRA_IN" -name "*.cat" -exec grep -c '^ENTRY' {} + 2>/dev/null | awk -F: '{s+=$2}END{print s+0}')
        if [ "$S_N" -gt 0 ]; then
            echo "  SUGRA chain ($S_N from $SUGRA_IN)"
            OUTPUT=$($PHASE2 "$SUGRA_IN" ${TAG}_r${round}s $FLAGS --save-nonsugra 2>&1)
            echo "$OUTPUT" | grep -E 'Written|Non-SUGRA'
            # Parse output dir names
            S_OUT=$(echo "$OUTPUT" | grep 'Written.*bases.*to ' | head -1 | sed 's/.*to //' | sed 's|/$||')
            S_NS_OUT=$(echo "$OUTPUT" | grep 'Non-SUGRA.*to ' | sed 's/.*to //' | sed 's|/$||')
            ANY_OUTPUT=1
        else
            S_OUT=""
            S_NS_OUT=""
        fi
    else
        S_OUT=""
        S_NS_OUT=""
    fi

    # --- non-SUGRA chain ---
    # Merge: previous SUGRA chain's nonsugra + previous NS chain's nonsugra
    NS_MERGE="tmp_ns_merge_${TAG}_${round}"
    rm -rf "$NS_MERGE"; mkdir -p "$NS_MERGE"

    # Previous SUGRA's nonsugra (from round-1 SUGRA chain)
    if [ -n "$PREV_S_NS" ] && [ -d "$PREV_S_NS" ]; then
        merge_cat_dirs "$PREV_S_NS" "$NS_MERGE"
    fi
    # Previous NS chain's nonsugra
    if [ -n "$PREV_NS_NS" ] && [ -d "$PREV_NS_NS" ]; then
        merge_cat_dirs "$PREV_NS_NS" "$NS_MERGE"
    fi
    # Also add current NS_IN if first round
    if [ -n "$NS_IN" ] && [ -d "$NS_IN" ]; then
        merge_cat_dirs "$NS_IN" "$NS_MERGE"
    fi

    NS_N=$(find "$NS_MERGE" -name "*.cat" -exec grep -c '^ENTRY' {} + 2>/dev/null | awk -F: '{s+=$2}END{print s+0}')

    NS_S_OUT=""
    NS_NS_OUT=""
    if [ "$NS_N" -gt 0 ]; then
        echo "  non-SUGRA chain ($NS_N)"
        OUTPUT=$($PHASE2 "$NS_MERGE" ${TAG}_r${round}ns $FLAGS --save-nonsugra 2>&1)
        echo "$OUTPUT" | grep -E 'Written|Non-SUGRA'
        NS_S_OUT=$(echo "$OUTPUT" | grep 'Written.*bases.*to ' | head -1 | sed 's/.*to //' | sed 's|/$||')
        NS_NS_OUT=$(echo "$OUTPUT" | grep 'Non-SUGRA.*to ' | sed 's/.*to //' | sed 's|/$||')
        ANY_OUTPUT=1
    fi
    rm -rf "$NS_MERGE"

    if [ "$ANY_OUTPUT" -eq 0 ]; then
        echo "  No output. Done."
        break
    fi

    # Merge SUGRA for next round: S_OUT + NS_S_OUT
    NEXT_S="tmp_next_s_${TAG}_$((round+1))"
    rm -rf "$NEXT_S"; mkdir -p "$NEXT_S"
    [ -n "$S_OUT" ] && [ -d "$S_OUT" ] && merge_cat_dirs "$S_OUT" "$NEXT_S"
    [ -n "$NS_S_OUT" ] && [ -d "$NS_S_OUT" ] && merge_cat_dirs "$NS_S_OUT" "$NEXT_S"

    NEXT_N=$(find "$NEXT_S" -name "*.cat" -exec grep -c '^ENTRY' {} + 2>/dev/null | awk -F: '{s+=$2}END{print s+0}')
    if [ "$NEXT_N" -eq 0 ]; then
        rm -rf "$NEXT_S"
        SUGRA_IN=""
    else
        SUGRA_IN="$NEXT_S"
    fi

    # Track nonsugra for next round
    PREV_S_NS="$S_NS_OUT"
    PREV_NS_NS="$NS_NS_OUT"
    NS_IN=""
done

# Clean tmp dirs
rm -rf tmp_ns_merge_${TAG}_* tmp_next_s_${TAG}_*

echo ""
echo "=== Chain $TAG complete ==="
