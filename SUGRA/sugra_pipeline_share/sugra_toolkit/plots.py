#!/usr/bin/env python3
"""Regenerate the gauge-analysis plots from the summary CSVs — NO catalog scan.

Edit the STYLE block (and the per-plot sections) below, then re-run:
    python plots.py
It reads the already-computed CSVs and rewrites the PNGs in seconds, so you can
tweak colours / sizes / labels / scales without the ~20-min full scan.

Inputs:
    gauge_out/by_T_summary.csv   (T, n_bases, max_total_rank)
    stats_final/by_T_H.csv       (T_H, n_bases)
Outputs (in gauge_out/):
    blocks_per_T.png, maxrank_per_T.png, blocks_per_THplus1.png
"""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ─────────────────────── STYLE (edit here) ───────────────────────
GAUGE_DIR = "gauge_out"          # where the plots are written / by_T_summary.csv lives
STATS_DIR = "stats_final"        # where by_T_H.csv lives
FIGSIZE   = (10, 5)              # (width, height) inches
DPI       = 150                  # output resolution
DOT_SIZE  = 16                   # scatter marker size (s=)
COLOR_BLOCKS = "orange"          # block-count plots
COLOR_RANK   = "red"             # rank plot
GRID_ALPHA   = 0.3               # 0 = no grid


def _read(path, key_col, val_col):
    xs, ys = [], []
    with open(path) as f:
        r = csv.reader(f); next(r)
        for row in r:
            xs.append(int(row[key_col])); ys.append(int(row[val_col]))
    return xs, ys


def scatter(x, y, *, color, log, xlabel, ylabel, title, fname):
    plt.figure(figsize=FIGSIZE)
    plt.scatter(x, y, color=color, s=DOT_SIZE)
    if log:
        plt.yscale("log")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, alpha=GRID_ALPHA)
    plt.tight_layout()
    plt.savefig(f"{GAUGE_DIR}/{fname}", dpi=DPI)
    plt.close()
    print("wrote", f"{GAUGE_DIR}/{fname}")


# ─────────────── plot 1: blocks per T (= sig_neg) ───────────────
T, nb = _read(f"{GAUGE_DIR}/by_T_summary.csv", 0, 1)
scatter(T, nb, color=COLOR_BLOCKS, log=True,
        xlabel=r"$T_\mathrm{block}$",
        ylabel="Number of non-Higgsable gravity blocks (log)",
        title=r"Non-Higgsable gravity blocks per $T_\mathrm{block}$",
        fname="blocks_per_T.png")

# ─────────────── plot 2: maximal total gauge rank per T ─────────
T, mr = _read(f"{GAUGE_DIR}/by_T_summary.csv", 0, 2)
scatter(T, mr, color=COLOR_RANK, log=False,
        xlabel=r"$T_\mathrm{block}$",
        ylabel="Maximal rank of gauge algebras",
        title=r"Maximal rank of gauge algebras per $T_\mathrm{block}$",
        fname="maxrank_per_T.png")

# ─────────────── plot 3: blocks per T_H + 1 ─────────────────────
TH, nbh = _read(f"{STATS_DIR}/by_T_H.csv", 0, 1)
scatter([t + 1 for t in TH], nbh, color=COLOR_BLOCKS, log=True,
        xlabel=r"$T_H + 1$",
        ylabel="Number of non-Higgsable gravity blocks (log)",
        title=r"Non-Higgsable gravity blocks per $T_H+1$",
        fname="blocks_per_THplus1.png")
