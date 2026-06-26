#!/usr/bin/env python3
"""Collect 'frozen' gravity blocks (su8 / su16 / so16 gauge) into one folder.

Frozen externals in the catalog (frozen-flux gauges):
    hat1m1  -> su8   (gauge-less (-1) hat1, target -1)
    hat1m2  -> su16  (gauge-less (-1) hat1, target -2)
    su8     -> su8   (frozen (-2))
    so16n2  -> so16  (frozen (-4) adjacent to a (-2))
(The "-2-4 frozen cluster" appears as SEPARATE su8 / so16n2 externals — they never
co-occur in one block.)

Output:
    <OUT>/by_external/<tag>/T<n>_<n>/<combo>.cat   explicit blocks per frozen external
    <OUT>/SUMMARY.txt                              counts + per-T breakdown
    <OUT>/manifest.csv                             combo, T_H, n_blocks, frozen_tags

Usage:  python collect_frozen.py "<catalog>" "<out dir>"
"""
import os
import sys
import glob
import shutil
import csv
from collections import defaultdict

FROZEN = {
    "hat1m1": "su8  (hat1, -1 target)",
    "hat1m2": "su16 (hat1, -2 target)",
    "su8":    "su8  (frozen -2)",
    "so16n2": "so16 (frozen -4)",
}


def n_entries(path):
    return sum(1 for l in open(path) if l.startswith("ENTRY"))


def run(cat, out):
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out, exist_ok=True)

    cnt = defaultdict(int)                       # tag -> blocks
    combos = defaultdict(set)                    # tag -> {combo}
    per_t = defaultdict(lambda: defaultdict(int))  # tag -> T_H -> blocks
    manifest = []                                # (combo, T_H, n, tags)
    frozen_combos = set()

    for f in sorted(glob.glob(os.path.join(cat, "T*", "*.cat"))):
        tdir = os.path.basename(os.path.dirname(f))      # T<n>_<n>
        try:
            th = int(tdir[1:].split("_")[0])
        except ValueError:
            continue
        combo = os.path.basename(f)[:-4]
        tags = combo.split("+")
        hits = [k for k in FROZEN if k in tags]
        if not hits:
            continue
        n = n_entries(f)
        if n == 0:
            continue
        frozen_combos.add((tdir, combo))
        manifest.append((combo, th, n, "+".join(hits)))
        for k in hits:
            cnt[k] += n
            combos[k].add(combo)
            per_t[k][th] += n
            dst = os.path.join(out, "by_external", k, tdir)
            os.makedirs(dst, exist_ok=True)
            shutil.copyfile(f, os.path.join(dst, combo + ".cat"))

    total_blocks = sum(n_entries(os.path.join(cat, td, cb + ".cat"))
                       for td, cb in frozen_combos)

    # manifest.csv
    with open(os.path.join(out, "manifest.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["combo", "T_H", "n_blocks", "frozen_tags"])
        for row in sorted(manifest, key=lambda r: (r[0], r[1])):
            w.writerow(row)

    # SUMMARY.txt
    with open(os.path.join(out, "SUMMARY.txt"), "w") as fh:
        fh.write("FROZEN GRAVITY BLOCKS — summary\n")
        fh.write(f"catalog: {cat}\n\n")
        fh.write("blocks per frozen external (blocks in any combo containing the tag):\n")
        for k, desc in FROZEN.items():
            fh.write(f"  {k:8} [{desc:24}]: {cnt[k]:>7,} blocks  ({len(combos[k])} combos)\n")
        fh.write(f"\n  total distinct frozen blocks (>=1 frozen tag): {total_blocks:,}\n")
        fh.write("\n  note: su8 and so16n2 (the -2-4 frozen cluster pieces) never co-occur.\n")
        for k in FROZEN:
            if not per_t[k]:
                continue
            fh.write(f"\n{k} — blocks per T_H:\n")
            for th in sorted(per_t[k]):
                fh.write(f"  T_H={th:<4} {per_t[k][th]:>6,}\n")
            fh.write(f"  combos: {', '.join(sorted(combos[k]))}\n")

    print(f"DONE -> {out}/")
    for k, desc in FROZEN.items():
        print(f"  {k:8}: {cnt[k]:>7,} blocks ({len(combos[k])} combos)")
    print(f"  total frozen blocks: {total_blocks:,}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    run(sys.argv[1], sys.argv[2])
