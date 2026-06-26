#!/usr/bin/env python3
"""Reproduce the 'Number of Gravity Blocks vs T' figure style for all CSVs.

Reads the gauge-analysis summary CSVs and writes publication-style scatter plots
(.png + .pdf). Runs locally with no pandas (csv + numpy only).

Paths default to this script's own directory (where the CSVs live and the plots
are written). Override with:
    python make_plots.py [SRC_DIR] [OUT_DIR]
"""
import os
import sys
import csv
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, LogLocator, NullFormatter

# ---------------------------------------------------------------- paths
_HERE = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1] if len(sys.argv) > 1 else _HERE       # where the CSVs live
OUT = sys.argv[2] if len(sys.argv) > 2 else _HERE       # where the plots go
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- style
plt.rcParams.update({
    "text.usetex": False,
    "mathtext.fontset": "cm",
    "font.family": "serif",
    "font.serif": ["Latin Modern Roman", "CMU Serif", "DejaVu Serif"],
    "axes.linewidth": 1.1,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 5,
    "ytick.major.size": 5,
    "xtick.minor.size": 2.8,
    "ytick.minor.size": 2.8,
    "savefig.dpi": 220,
})

# per-plot colours so the figures are distinguishable side by side in a paper
BLUE   = "#1414ff"   # gravity blocks vs T^H (linear & log — same data)
ORANGE = "#ef7a00"   # gravity blocks vs T (= sigma_-)
RED    = "#d62728"   # maximum total gauge rank vs T
GREEN  = "#1f9e44"   # gravity blocks vs number of external curves
PURPLE = "#8a3ffc"   # distinct gauge patterns vs T^H+1
GRID = "#cfcfcf"
DASH = "#8a8a8a"


# ── tiny CSV helpers (pandas-free) ────────────────────────────────────
def _rows(name):
    with open(os.path.join(SRC, name)) as f:
        r = csv.reader(f)
        next(r)                       # header
        return [row for row in r if row]


def _col(rows, i):
    return np.array([int(row[i]) for row in rows])


# ── plotting primitives (unchanged style) ─────────────────────────────
def style_axes(ax, x_major=25, x_minor=5):
    ax.grid(True, which="major", color=GRID, linewidth=0.7, zorder=0)
    for s in ax.spines.values():
        s.set_color("black")
        s.set_zorder(5)
    ax.tick_params(which="both", top=False, right=False, labelsize=15,
                   color="black", width=1.0)
    if x_major:
        ax.xaxis.set_major_locator(MultipleLocator(x_major))
    if x_minor:
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor))


def base_fig():
    fig, ax = plt.subplots(figsize=(11.4, 6.0))
    fig.subplots_adjust(left=0.105, right=0.975, top=0.965, bottom=0.135)
    return fig, ax


def dashed_marker(ax, xval, label, logy=False):
    ax.axvline(xval, color=DASH, linestyle=(0, (6, 4)), linewidth=1.3, zorder=2)
    if logy:
        ylo, yhi = ax.get_ylim()
        ypos = ylo * (yhi / ylo) ** 0.045
    else:
        ypos = ax.get_ylim()[1] * 0.045
    ax.annotate(label, xy=(xval, ypos), xytext=(-6, 0),
                textcoords="offset points", ha="right", va="bottom",
                fontsize=13)


def scatter(ax, x, y, color=BLUE, logy=False):
    ax.scatter(x, y, s=15, c=color, edgecolors="none", zorder=4)
    if logy:
        ax.set_yscale("log")


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("saved", name)


def log_yaxis(ax):
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=12))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=range(2, 10), numticks=80))
    ax.yaxis.set_minor_formatter(NullFormatter())


# ================================================================ FILE 1
r1 = _rows("by_T_H.csv")
TH, NB = _col(r1, 0), _col(r1, 1)
xmax1 = int(TH.max())

fig, ax = base_fig()                                    # 1a — linear
scatter(ax, TH, NB)
style_axes(ax)
ax.set_xlim(0, 200); ax.set_ylim(bottom=0)
ax.set_xlabel(r"$T^H$", fontsize=21)
ax.set_ylabel("Number of Gravity Blocks", fontsize=18)
dashed_marker(ax, xmax1, rf"$T^H = {xmax1}$")
save(fig, "Number_GBs_by_TH_linear")

fig, ax = base_fig()                                    # 1b — log y
scatter(ax, TH, NB, logy=True)
style_axes(ax)
ax.set_xlim(0, 200); ax.set_ylim(0.7, NB.max() * 1.6)
log_yaxis(ax)
ax.set_xlabel(r"$T^H$", fontsize=21)
ax.set_ylabel("Number of Gravity Blocks", fontsize=18)
dashed_marker(ax, xmax1, rf"$T^H = {xmax1}$", logy=True)
save(fig, "Number_GBs_by_TH_log")

# ================================================================ FILE 2
r2 = _rows("by_T_summary.csv")
T, NB2, MR = _col(r2, 0), _col(r2, 1), _col(r2, 2)
xmax2 = int(T.max())

fig, ax = base_fig()                                    # 2a — n_bases vs T (log)
scatter(ax, T, NB2, color=ORANGE, logy=True)
style_axes(ax)
ax.set_xlim(0, 200); ax.set_ylim(0.7, NB2.max() * 1.6)
log_yaxis(ax)
ax.set_xlabel(r"$T$", fontsize=21)
ax.set_ylabel("Number of Gravity Blocks", fontsize=18)
dashed_marker(ax, xmax2, rf"$T = {xmax2}$", logy=True)
save(fig, "Number_GBs_by_T_log")

fig, ax = base_fig()                                    # 2b — max rank vs T
scatter(ax, T, MR, color=RED)
style_axes(ax)
ax.set_xlim(0, 200); ax.set_ylim(bottom=0)
ax.set_xlabel(r"$T$", fontsize=21)
ax.set_ylabel("Maximum total gauge rank", fontsize=18)
dashed_marker(ax, xmax2, rf"$T = {xmax2}$")
save(fig, "maxrank_by_T")

# ================================================================ FILE 3
g3 = defaultdict(int)                                    # sum n_bases by n_ext
for row in _rows("by_TH_nExt_gauge.csv"):
    g3[int(row[1])] += int(row[3])
gx = np.array(sorted(g3)); gy = np.array([g3[k] for k in gx])

fig, ax = base_fig()
scatter(ax, gx, gy, color=GREEN, logy=True)
style_axes(ax, x_major=1, x_minor=None)
ax.set_xlim(0.4, gx.max() + 0.6)
ax.set_ylim(gy.min() * 0.6, gy.max() * 1.7)
ax.set_xlabel("Number of external curves", fontsize=18)
ax.set_ylabel("Number of Gravity Blocks", fontsize=18)
save(fig, "Number_GBs_by_nExternal")

# ================================================================ FILE 4
p4 = defaultdict(set)                                    # distinct patterns by T_H+1
for row in _rows("gauge_patterns_by_THplus1.csv"):
    p4[int(row[0])].add(row[1])
px = np.array(sorted(p4)); py = np.array([len(p4[k]) for k in px])
xmax4 = int(px.max())

fig, ax = base_fig()
scatter(ax, px, py, color=PURPLE)
style_axes(ax)
ax.set_xlim(0, 200); ax.set_ylim(bottom=0)
ax.set_xlabel(r"$T^H + 1$", fontsize=21)
ax.set_ylabel("Number of distinct gauge patterns", fontsize=18)
dashed_marker(ax, xmax4, rf"$T^H+1 = {xmax4}$")
save(fig, "Number_distinct_patterns")

print("ALL DONE ->", OUT)
