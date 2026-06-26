#!/usr/bin/env python3
"""Bar chart: number of SUGRA blocks containing >=1 of each gauge algebra.

Single-curve external gauges (from blocks_by_gauge_nonhc.csv) and NHC clusters
(from blocks_by_nhc.csv) are shown as SEPARATE bars (NHC considered separately,
not folded into the su2/g2/so7 bars). Gauge-less (none) is omitted.

Reads the two CSVs from <stats_dir> (default: ../stats_final next to this script)
and writes blocks_per_gauge.png/.pdf there.

Usage:  python plot_gauge_census.py [stats_dir]
"""
import os
import sys
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
SD = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "stats_final")

SINGLE_COLOR = "#3b6fb6"     # single-curve external gauge
NHC_COLOR = "#e07b00"        # NHC cluster (separate)


def _rows(name):
    with open(os.path.join(SD, name)) as f:
        r = csv.reader(f)
        next(r)
        return [x for x in r if x]


# single-curve external gauges (drop gauge-less (none))
single = [(f"{g}({si})", int(n)) for g, si, n, _ in _rows("blocks_by_gauge_nonhc.csv")
          if g != "(none)"]
single.sort(key=lambda x: -x[1])

# NHC clusters, shown separately
nhc = [(r[0], int(r[3])) for r in _rows("blocks_by_nhc.csv") if r[0] != "ANY_nhc"]
nhc.sort(key=lambda x: -x[1])

labels = [l for l, _ in single] + [l for l, _ in nhc]
vals = [v for _, v in single] + [v for _, v in nhc]
colors = [SINGLE_COLOR] * len(single) + [NHC_COLOR] * len(nhc)

def _fmt(n):
    if n >= 1_000_000:
        return f"{n / 1e6:.1f}M"
    if n >= 1_000:
        return f"{n / 1e3:.0f}k"
    return str(n)


plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                     "axes.linewidth": 1.1, "savefig.dpi": 200})
fig, ax = plt.subplots(figsize=(max(11, len(labels) * 0.62), 6.8))
fig.subplots_adjust(left=0.085, right=0.985, top=0.90, bottom=0.27)
ax.bar(range(len(labels)), vals, color=colors, edgecolor="none", width=0.72, zorder=3)
ax.set_yscale("log")
ax.set_xlim(-0.8, len(labels) - 0.2)                 # side breathing room
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=50, ha="right", fontsize=10.5)
ax.set_ylabel("Number of SUGRA blocks (log)", fontsize=14)
ax.set_title(r"SUGRA blocks containing each external ($\geq 1$)",
             fontsize=15, pad=12)
ax.grid(True, axis="y", color="#d8d8d8", linewidth=0.7, zorder=0)
ax.set_ylim(0.7, max(vals) * 12)                     # generous top headroom
ax.legend(handles=[Patch(color=SINGLE_COLOR, label="single-curve external gauge"),
                   Patch(color=NHC_COLOR, label="NHC cluster (counted separately)")],
          fontsize=11.5, frameon=True, framealpha=0.95, loc="upper right")
# compact horizontal count just above each bar
for i, v in enumerate(vals):
    ax.text(i, v * 1.45, _fmt(v), ha="center", va="bottom", fontsize=8.5,
            color="#333333")

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(SD, f"blocks_per_gauge.{ext}"), bbox_inches="tight")
plt.close(fig)
print(f"saved blocks_per_gauge.png/.pdf -> {SD}")
print(f"  {len(single)} single-curve gauges + {len(nhc)} NHC clusters")
