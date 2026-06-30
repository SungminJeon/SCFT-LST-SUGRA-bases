#!/usr/bin/env python3
"""Separate statistics + graph for the gauge-less unit externals '1' and '2'.

'1' (a gauge-less -1) and '2' (a gauge-less -2) carry no gauge algebra, so they
are invisible in the gauge census; this tool tracks them on their own.

Usage:
    python external_1_2_analysis.py <catalog> [out_dir]   # default out_dir: stats_final

Writes into out_dir:
    external_1_2_T_distribution.csv   T_block -> #blocks containing '1' / '2'
    external_1_2_cogauge.csv          per unit external, co-occurring externals (by #blocks)
    external_1_2_vs_Tblock.png/.pdf   blocks vs T_block, one curve per unit external
"""
import sys, os, glob, csv
from collections import Counter
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

TAGS = ["1", "2"]
COLOR = {"1": "tab:green", "2": "tab:purple"}


def run(cat, out="stats_final"):
    os.makedirs(out, exist_ok=True)
    distT = {t: Counter() for t in TAGS}   # T_block -> #blocks containing tag
    cog   = {t: Counter() for t in TAGS}   # co-external -> #blocks
    total = {t: 0 for t in TAGS}
    single = {t: 0 for t in TAGS}          # blocks whose combo is exactly the tag

    for f in glob.glob(os.path.join(cat, "T*", "*.cat")):
        combo = os.path.basename(f)[:-4]
        tags = combo.split("+")
        present = [t for t in TAGS if t in tags]
        if not present:
            continue
        co_tags = set(tags)
        with open(f) as fh:
            for line in fh:
                if line.startswith("PHYSICS"):
                    T = int(line.split()[1])          # PHYSICS <T> Hc V Hn ...
                    for t in present:
                        distT[t][T] += 1
                        total[t] += 1
                        if combo == t:
                            single[t] += 1
                        for o in co_tags:
                            if o != t:
                                cog[t][o] += 1

    allT = sorted(set().union(*[set(distT[t]) for t in TAGS]) or {0})
    with open(os.path.join(out, "external_1_2_T_distribution.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["T_block", "n_with_1", "n_with_2"])
        for T in allT:
            w.writerow([T, distT["1"].get(T, 0), distT["2"].get(T, 0)])

    with open(os.path.join(out, "external_1_2_cogauge.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["unit_external", "co_external", "n_blocks"])
        for t in TAGS:
            for o, n in cog[t].most_common():
                w.writerow([t, o, n])

    plt.figure(figsize=(10, 5))
    for t in TAGS:
        xs = sorted(distT[t])
        ys = [distT[t][x] for x in xs]
        plt.scatter(xs, ys, s=16, color=COLOR[t],
                    label=f"external '{t}'  ({total[t]:,} blocks)")
    plt.yscale("log")
    plt.xlabel(r"$T_\mathrm{block}$")
    plt.ylabel("Number of gravity blocks (log)")
    plt.title(r"Blocks containing gauge-less unit external '1' / '2' per $T_\mathrm{block}$")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out, "external_1_2_vs_Tblock.png"), dpi=150)
    plt.savefig(os.path.join(out, "external_1_2_vs_Tblock.pdf"))
    plt.close()

    for t in TAGS:
        rng = f"{min(distT[t])}..{max(distT[t])}" if distT[t] else "-"
        print(f"external '{t}': {total[t]:,} blocks ({single[t]:,} as the only external), "
              f"T_block {rng}, top co-external: {cog[t].most_common(1)}")


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "stats_final")
