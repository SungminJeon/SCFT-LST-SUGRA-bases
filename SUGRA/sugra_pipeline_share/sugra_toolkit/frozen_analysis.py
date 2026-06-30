#!/usr/bin/env python3
"""Frozen-block analysis: T(=neg eigenvalues, sig_neg) distribution + co-gauge.

Reads the collected frozen blocks (collect_frozen.py output:
<FROZEN>/by_external/<tag>/T*/<combo>.cat) and writes, into <FROZEN>/analysis/:

  frozen_T_distribution.csv       pivot  T  vs  n_blocks  per frozen external
  frozen_T_distribution_long.csv  long   (external, T, n_blocks)
  frozen_T_vs_blocks.png/.pdf     one figure, one colour per frozen external
  hat1m1_cogauge.csv              gauge combos co-occurring with hat1m1, by freq
  hat1m2_cogauge.csv              gauge combos co-occurring with hat1m2, by freq

T = sig_neg = number of negative eigenvalues of the block's intersection form
(PHYSICS field 7).

Usage:  python frozen_analysis.py "<Frozen blocks dir>"
"""
import os
import sys
import glob
import csv
from collections import defaultdict, Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TAGS = ["hat1m1", "hat1m2", "su8", "so16n2"]
COLOR = {"hat1m1": "#1414ff", "hat1m2": "#ef7a00", "su8": "#d62728", "so16n2": "#1f9e44"}
NHC_PER = {"nhc_2_3": 2, "nhc_2_3_2": 3, "nhc_2_2_3": 3, "nhc_2_4": 2}


def blocks_of(froot, tag):
    """Yield (sig_neg, combo) for every block under by_external/<tag>/."""
    for f in glob.glob(os.path.join(froot, "by_external", tag, "T*", "*.cat")):
        combo = os.path.basename(f)[:-4]
        for line in open(f):
            if line.startswith("PHYSICS"):
                p = line.split()
                if len(p) >= 8:
                    yield int(p[7]), combo            # p[7] = sigNeg = T


def run(froot):
    out = os.path.join(froot, "analysis")
    os.makedirs(out, exist_ok=True)

    # ---- T distribution per external ------------------------------------
    dist = {t: Counter() for t in TAGS}          # tag -> T -> n
    for t in TAGS:
        for T, _ in blocks_of(froot, t):
            dist[t][T] += 1

    allT = sorted({T for t in TAGS for T in dist[t]})
    # pivot CSV
    with open(os.path.join(out, "frozen_T_distribution.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["T"] + TAGS)
        for T in allT:
            w.writerow([T] + [dist[t].get(T, 0) for t in TAGS])
    # long CSV
    with open(os.path.join(out, "frozen_T_distribution_long.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["external", "T", "n_blocks"])
        for t in TAGS:
            for T in sorted(dist[t]):
                w.writerow([t, T, dist[t][T]])

    # ---- combined plot (one colour per external) ------------------------
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                         "axes.linewidth": 1.1, "savefig.dpi": 200})
    fig, ax = plt.subplots(figsize=(11, 6))
    fig.subplots_adjust(left=0.10, right=0.97, top=0.95, bottom=0.13)
    for t in TAGS:
        if not dist[t]:
            continue
        xs = sorted(dist[t])
        ys = [dist[t][x] for x in xs]
        tot = sum(ys)
        ax.scatter(xs, ys, s=26, c=COLOR[t], edgecolors="none", zorder=4,
                   label=f"{t}  ({tot:,})")
    ax.set_yscale("log")
    ax.grid(True, which="major", color="#cfcfcf", linewidth=0.7, zorder=0)
    ax.set_xlabel(r"$T_\mathrm{block}$", fontsize=17)
    ax.set_ylabel("Number of gravity blocks (log)", fontsize=15)
    ax.set_title("Frozen gravity blocks per $T_\mathrm{block}$", fontsize=16)
    ax.legend(title="frozen external (total)", fontsize=12, title_fontsize=12,
              frameon=True, framealpha=0.95)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"frozen_T_vs_blocks.{ext}"), bbox_inches="tight")
    plt.close(fig)

    # ---- co-occurring gauge combos for hat1m1 / hat1m2 ------------------
    def co_combo(combo, tag):
        """Remaining tags after removing ONE occurrence of `tag`, as a sorted combo."""
        tags = combo.split("+")
        tags.remove(tag)                          # drop one hat1mN
        return "+".join(sorted(tags)) if tags else "(none)"

    for tag in ("hat1m1", "hat1m2"):
        co = Counter()                            # co-combo -> n_blocks
        for T, combo in blocks_of(froot, tag):
            co[co_combo(combo, tag)] += 1
        with open(os.path.join(out, f"{tag}_cogauge.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["co_occurring_gauge_combo", "n_blocks"])
            for combo, n in co.most_common():
                w.writerow([combo, n])

    # ---- console summary ------------------------------------------------
    print(f"DONE -> {out}/")
    for t in TAGS:
        print(f"  {t:8}: {sum(dist[t].values()):>6,} blocks, T range "
              f"{min(dist[t]) if dist[t] else '-'}..{max(dist[t]) if dist[t] else '-'}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    run(sys.argv[1])
