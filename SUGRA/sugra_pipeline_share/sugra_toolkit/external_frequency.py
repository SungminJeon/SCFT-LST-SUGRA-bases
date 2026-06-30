#!/usr/bin/env python3
"""Frequency of each individual external (not combo frequency).

Reads a stats `by_combo.csv` (combo, n_bases) and decomposes each combo into
its external tags, accumulating how often each external appears overall.

Under the current label policy an NHC cluster is ONE tag per cluster, so the
counts here are per external object (one nhc_2_2_3 cluster = +1), not per curve.

Outputs (in the same dir as by_combo.csv, or <out>):
  external_frequency.csv   external, total_occurrences, n_bases_containing, occ_per_base
  external_frequency.png   bar chart (log scale)

Usage:  python external_frequency.py <by_combo.csv> [out_dir]
"""
import os
import sys
import csv
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def run(by_combo_csv, out=None):
    out = out or os.path.dirname(os.path.abspath(by_combo_csv))
    os.makedirs(out, exist_ok=True)

    occ, inbase, tot = Counter(), Counter(), 0
    with open(by_combo_csv) as f:
        r = csv.reader(f); next(r)
        for row in r:                       # by_combo may carry extra cols (percent_of_total)
            if len(row) < 2 or not row[1].isdigit():
                continue
            combo, n = row[0], int(row[1]); tot += n
            for tag, c in Counter(combo.split("+")).items():
                occ[tag] += c * n        # total occurrences (per-object)
                inbase[tag] += n         # bases containing >=1
    rows = sorted(occ.items(), key=lambda kv: -kv[1])

    with open(os.path.join(out, "external_frequency.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["external", "total_occurrences", "n_bases_containing", "occ_per_base"])
        for t, o in rows:
            w.writerow([t, o, inbase[t], round(o / tot, 4) if tot else 0])

    labels = [t for t, _ in rows]
    vals = [o for _, o in rows]
    plt.figure(figsize=(max(8, len(labels) * 0.45), 5))
    plt.bar(range(len(labels)), vals, color="#36a")
    plt.yscale("log")
    plt.xticks(range(len(labels)), labels, rotation=60, ha="right", fontsize=8)
    plt.ylabel("total occurrences (log)")
    plt.title("Frequency of each external (over all bases)")
    plt.tight_layout()
    plt.savefig(os.path.join(out, "external_frequency.png"), dpi=150)
    plt.close()

    print(f"{tot:,} bases, {len(rows)} external types -> {out}/external_frequency.csv|.png")
    for t, o in rows:
        print(f"  {t:12} {o:>14,} {inbase[t]:>14,}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    run(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
