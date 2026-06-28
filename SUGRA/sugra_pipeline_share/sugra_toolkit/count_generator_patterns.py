#!/usr/bin/env python3
"""Distinct GENERATOR (external combo) patterns, by external gauge content.

Each external tag maps to the gauge factor(s) it supports; a block's generator
pattern is the sorted multiset of those factors over its whole combo:
    su2 / su2n3 / su2n3mix -> su2
    su3                    -> su3
    su3n2 / su3n2mix       -> g2          (su3n2 IS treated as g2)
    so7 / so7mix           -> so7
    so8/f4/e6/e7p/e7/e8    -> so8/f4/e6/e7/e7/e8
    su8 / hat1m1           -> su8 ;  hat1m2 -> su16 ;  so16n2 -> so16
    nhc_2_3   -> su2+g2 ;  nhc_2_2_3 -> sp1+g2 ;  nhc_2_3_2 -> su2+so7+su2
    2 / 1 (gauge-less)     -> (nothing; pattern "(none)" if a block has only these)

CAVEAT (tag-based): a `nhc_2_3` whose -3 acquires a second -2 is geometrically a
`nhc_2_3_2` (su2+so7+su2), but here it is still scored su2+g2. An exact count that
honours that needs the IF-based gauge decomposition (gauge_analysis).

Output CSV: generator_pattern, n_blocks, percent_of_total, n_combos
Usage:  python count_generator_patterns.py "<catalog>" out.csv
"""
import os
import sys
import glob
import csv
from collections import defaultdict

TAG = {
    "su2": ["su2"], "su2n3": ["su2"], "su2n3mix": ["su2"], "su3": ["su3"],
    "su3n2": ["g2"], "su3n2mix": ["g2"], "so7": ["so7"], "so7mix": ["so7"],
    "so8": ["so8"], "f4": ["f4"], "e6": ["e6"], "e7p": ["e7"], "e7": ["e7"],
    "e8": ["e8"], "su8": ["su8"], "hat1m1": ["su8"], "hat1m2": ["su16"],
    "so16n2": ["so16"], "2": [], "1": [],
    "nhc_2_3": ["su2", "g2"], "nhc_2_2_3": ["sp1", "g2"],
    "nhc_2_3_2": ["su2", "so7", "su2"],
}


def run(cat, outcsv):
    n_blocks = defaultdict(int)
    n_combos = defaultdict(int)
    total = 0
    unknown = set()
    for f in glob.glob(os.path.join(cat, "T*", "*.cat")):
        combo = os.path.basename(f)[:-4]
        n = sum(1 for l in open(f) if l.startswith("ENTRY"))
        if n == 0:
            continue
        facs = []
        for t in combo.split("+"):
            if t in TAG:
                facs += TAG[t]
            else:
                unknown.add(t)
        pat = "+".join(sorted(facs)) if facs else "(none)"
        n_blocks[pat] += n
        n_combos[pat] += 1
        total += n

    rows = sorted(n_blocks, key=lambda p: -n_blocks[p])
    with open(outcsv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["generator_pattern", "n_blocks", "percent_of_total", "n_combos"])
        for p in rows:
            w.writerow([p, n_blocks[p], f"{100.0 * n_blocks[p] / total:.4f}", n_combos[p]])

    print(f"distinct generator patterns: {len(rows):,}")
    print(f"total blocks: {total:,}")
    if unknown:
        print(f"UNMAPPED tags: {sorted(unknown)}")
    print(f"-> {outcsv}")
    print("top 10:")
    for p in rows[:10]:
        print(f"  {n_blocks[p]:>9,}  {p}")


def run_exact(cat_dir, outcsv):
    """IF-based: external generator gauge from the actual intersection form, so a
    -2-3 grown into a -2-3-2 is scored su2+so7+su2 (not su2+g2). Slow (full scan)."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import sugra
    import gauge_analysis
    cat = sugra.Catalog(cat_dir)
    n_blocks = defaultdict(int)
    total = 0
    for th in cat.tensors():
        for cb in cat.combos(th):
            for b in cat.load(th, cb):
                g = gauge_analysis.gauge_per_curve(b.IF, b.externals)
                ext = {e["curveIdx"] for e in b.externals
                       if 0 <= e.get("curveIdx", -1) < len(g)}
                facs = sorted(g[i] for i in ext if g[i])
                pat = "+".join(facs) if facs else "(none)"
                n_blocks[pat] += 1
                total += 1
        print(f"  T_H={th} done ({total:,})", flush=True)
    rows = sorted(n_blocks, key=lambda p: -n_blocks[p])
    with open(outcsv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["generator_pattern", "n_blocks", "percent_of_total"])
        for p in rows:
            w.writerow([p, n_blocks[p], f"{100.0 * n_blocks[p] / total:.4f}"])
    print(f"distinct generator patterns (IF-exact): {len(rows):,}")
    print(f"total blocks: {total:,}  -> {outcsv}")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--exact"]
    if len(args) < 2:
        print(__doc__)
        sys.exit(1)
    if "--exact" in sys.argv:
        run_exact(args[0], args[1])
    else:
        run(args[0], args[1])
