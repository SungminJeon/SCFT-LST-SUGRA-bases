#!/usr/bin/env python3
"""Count blocks by external gauge algebra (+ self-intersection), with % of total.

Each external (single-curve OR NHC cluster) supports gauge factors on curves of a
given self-intersection. We map each external TAG to the set of (gauge, self_int)
it carries, then count how many blocks contain AT LEAST ONE curve of each
(gauge, self_int).  su3n2/su3n2mix -> g2 (user request). NHC clusters are counted
by the gauges their curves support. Gauge-less externals 2/1 -> (none).

Rows OVERLAP (a block with e6 and su2 is counted in both rows), so percentages
do not sum to 100. percent = n_blocks / (total blocks).

Output: <out>.csv  with columns gauge, self_int, n_blocks, percent

Usage:  python count_by_gauge.py "<catalog>" out.csv
"""
import os
import sys
import glob
import csv
from collections import defaultdict

# external TAG -> list of (gauge, self_int) it supports
TAG_GAUGE = {
    # single-curve externals
    "su2": [("su2", -2)], "su2n3": [("su2", -2)], "su2n3mix": [("su2", -2)],
    "su3": [("su3", -3)],
    "su3n2": [("g2", -3)], "su3n2mix": [("g2", -3)],       # user: -> g2
    "so7": [("so7", -3)], "so7mix": [("so7", -3)],
    "so8": [("so8", -4)], "f4": [("f4", -5)], "e6": [("e6", -6)],
    "e7p": [("e7", -7)], "e7": [("e7", -8)], "e8": [("e8", -12)],
    "su8": [("su8", -2)], "hat1m1": [("su8", -1)], "hat1m2": [("su16", -1)],
    "so16n2": [("so16", -4)],
    "2": [("(none)", -2)], "1": [("(none)", -1)],
    # NHC clusters: gauges their curves support (lone -2 of -2-2-3 carries none)
    "nhc_2_3":   [("su2", -2), ("g2", -3)],
    "nhc_2_2_3": [("sp1", -2), ("g2", -3)],
    "nhc_2_3_2": [("su2", -2), ("so7", -3)],
}

# display order for the table
ORDER = [("su2", -2), ("sp1", -2), ("su3", -3), ("g2", -3), ("so7", -3),
         ("so8", -4), ("f4", -5), ("e6", -6), ("e7", -7), ("e7", -8),
         ("e8", -12), ("su8", -1), ("su8", -2), ("su16", -1), ("so16", -4),
         ("(none)", -1), ("(none)", -2)]


def run(cat, outcsv, no_nhc=False):
    tag_gauge = {k: v for k, v in TAG_GAUGE.items()
                 if not (no_nhc and k.startswith("nhc_"))}
    count = defaultdict(int)            # (gauge, self_int) -> n_blocks
    total = 0
    unknown = set()
    for f in glob.glob(os.path.join(cat, "T*", "*.cat")):
        combo = os.path.basename(f)[:-4]
        n = sum(1 for l in open(f) if l.startswith("ENTRY"))
        if n == 0:
            continue
        total += n
        pairs = set()
        for tag in combo.split("+"):
            if tag in tag_gauge:
                pairs.update(tag_gauge[tag])
            elif tag not in TAG_GAUGE:
                unknown.add(tag)
        for gp in pairs:
            count[gp] += n

    rows = []
    seen = set()
    for gp in ORDER:
        if gp in count:
            rows.append(gp); seen.add(gp)
    for gp in sorted(count, key=lambda x: -count[x]):     # any extras not in ORDER
        if gp not in seen:
            rows.append(gp)

    with open(outcsv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gauge", "self_int", "n_blocks", "percent_of_total"])
        for g, si in rows:
            n = count[(g, si)]
            w.writerow([g, si, n, f"{100.0*n/total:.3f}"])

    print(f"total blocks: {total:,}")
    if unknown:
        print(f"UNMAPPED tags (ignored): {sorted(unknown)}")
    print(f"{'gauge':8} {'s.i.':>5} {'n_blocks':>12} {'%':>8}")
    for g, si in rows:
        n = count[(g, si)]
        print(f"{g:8} {si:>5} {n:>12,} {100.0*n/total:>7.2f}%")
    print(f"\n-> {outcsv}")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--no-nhc"]
    if len(args) < 2:
        print(__doc__); sys.exit(1)
    run(args[0], args[1], no_nhc="--no-nhc" in sys.argv)
