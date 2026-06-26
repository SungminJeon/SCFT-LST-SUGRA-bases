#!/usr/bin/env python3
"""Full SUGRA gauge-algebra analysis: tables and plots vs T_H+1.

The full gauge algebra of a base = NHC decomposition of the WHOLE intersection
form (LST + external curves) + su2 on standalone -2 externals + hat1 (su8/su16).
The lone -2 of an nhc_2_2_3 cluster carries NO gauge (correct convention; the
generator's stored V over-counts these by su2 on some bases).

Outputs (in <out>/):
  by_TH_nExt_gauge.csv          (T_H, n_external_curves, gauge_algebra) -> n_bases
  gauge_patterns_by_THplus1.csv (T_H+1, gauge_algebra)                  -> n_bases
  by_T_summary.csv              (T = sig_neg)  -> n_bases, max_total_rank
  blocks_per_T.png              scatter: non-Higgsable gravity blocks vs T (orange)
  maxrank_per_T.png             scatter: maximal total gauge rank vs T (red)
  blocks_per_THplus1.png        scatter: non-Higgsable gravity blocks vs T_H+1 (orange)

T = sig_neg (extended-base tensor count, from PHYSICS), not T_H+1.

Usage:  python gauge_analysis.py "<catalog>" "<out dir>"
"""
import os
import sys
import csv
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stats                                    # noqa: E402  (_bases, _nhc_lookup)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                 # noqa: E402

# rank of each simple factor (4th field of GaugeInfo in sugra_generator.h)
RANK = {"su2": 1, "sp1": 1, "su3": 2, "g2": 2, "so7": 3, "so8": 4, "f4": 4,
        "e6": 6, "e7": 7, "e8": 8, "so16": 8, "su8": 7, "su16": 15}

_SINGLE = {-3: "su3", -4: "so8", -5: "f4", -6: "e6", -7: "e7", -8: "e7", -12: "e8"}
_PERC = {"-2-3": ["su2", "g2"], "-2-3-2": ["su2", "so7", "su2"],
         "-2-2-3": [None, "sp1", "g2"], "-2-4": ["su8", "so16"]}
_NSI = {"-2-3": (-2, -3), "-2-3-2": (-2, -3, -2),
        "-2-2-3": (-2, -2, -3), "-2-4": (-2, -4)}


def _walk(comp, IF):
    if len(comp) <= 1:
        return list(comp)
    adj = {a: [b for b in comp if b != a and IF[a][b] != 0] for a in comp}
    cur = next((c for c in comp if len(adj[c]) <= 1), comp[0])
    o, seen = [], set()
    while True:
        o.append(cur); seen.add(cur)
        nxt = next((b for b in adj[cur] if b not in seen), None)
        if nxt is None:
            break
        cur = nxt
    return o


def gauge_factors(IF, externals):
    """List of gauge-algebra factors for the full base (correct convention)."""
    n = len(IF)
    if n == 0:
        return []
    g = [None] * n
    diag = [int(IF[i][i]) for i in range(n)]
    m1 = {i for i in range(n) if diag[i] == -1}
    pure2, seen = set(), set()
    for s in range(n):
        if s in m1 or s in seen:
            continue
        comp, st = [], [s]; seen.add(s)
        while st:
            c = st.pop(); comp.append(c)
            for j in range(n):
                if j not in m1 and j not in seen and IF[c][j] != 0:
                    seen.add(j); st.append(j)
        if all(diag[c] == -2 for c in comp):           # pure -2: no gauge
            pure2.update(comp); continue
        if len(comp) == 1:
            g[comp[0]] = _SINGLE.get(diag[comp[0]]); continue
        name = stats._nhc_lookup(IF[np.ix_(sorted(comp), sorted(comp))])
        if name in _PERC:
            order = _walk(comp, IF)
            pc = _PERC[name]
            if tuple(diag[c] for c in order) == _NSI[name][::-1]:
                pc = list(reversed(pc))
            for c, a in zip(order, pc):
                g[c] = a
        else:
            for c in comp:
                g[c] = _SINGLE.get(diag[c])
    for e in externals:                                # external enhancements
        c = e.get("curveIdx", -1)
        if not (0 <= c < n):
            continue
        if e.get("isHat1"):                            # hat1 -> su8 / su16
            # NOTE: must guard on isHat1, NOT extSI==-1 — the special_m1 "1"
            # external is also extSI=-1 but GAUGE-LESS (isHat1=False); treating
            # every -1 external as hat1 wrongly su8/su16-tags those bases.
            a = {-1: "su8", -2: "su16"}.get(e.get("targetSI"))
            if a:
                g[c] = a
        elif e.get("extSI") == -2 and e.get("targetSI") == -1 and \
                g[c] is None and c in pure2:            # standalone su2 only
            g[c] = "su2"
    return [x for x in g if x]


def run(root, out):
    os.makedirs(out, exist_ok=True)
    table = Counter()                       # (T_H, n_ext, gauge) -> n
    patterns = Counter()                    # (T_H, gauge) -> n
    count_t = Counter()                     # T (= sig_neg) -> n_bases
    maxrank_t = defaultdict(int)            # T (= sig_neg) -> max total rank
    count_th = Counter()                    # T_H -> n_bases
    total = 0

    th_dirs = sorted(
        (d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d))
         and d.startswith("T")),
        key=lambda s: int(s[1:].split("_")[0]) if s[1:].split("_")[0].isdigit() else 0,
    )
    for d in th_dirs:
        tp = os.path.join(root, d)
        for fn in os.listdir(tp):
            if not fn.endswith(".cat"):
                continue
            for ext, IF, T_H, ct, cid, sig_neg in stats._bases(os.path.join(tp, fn)):
                if IF is None:
                    continue
                total += 1
                facs = gauge_factors(IF, ext)
                gauge = "+".join(sorted(facs)) or "(none)"
                n_ext = len(ext)
                rank = sum(RANK.get(f, 0) for f in facs)
                table[(T_H, n_ext, gauge)] += 1
                patterns[(T_H, gauge)] += 1
                count_t[sig_neg] += 1
                count_th[T_H] += 1
                if rank > maxrank_t[sig_neg]:
                    maxrank_t[sig_neg] = rank
        print(f"  {d} done", flush=True)

    # ── CSVs ────────────────────────────────────────────────────────────
    def _csv(name, header, rows):
        with open(os.path.join(out, name), "w", newline="") as f:
            w = csv.writer(f); w.writerow(header); w.writerows(rows)

    _csv("by_TH_nExt_gauge.csv", ["T_H", "n_external_curves", "gauge_algebra", "n_bases"],
         [(th, ne, g, n) for (th, ne, g), n in
          sorted(table.items(), key=lambda kv: (kv[0][0], kv[0][1], -kv[1]))])
    _csv("gauge_patterns_by_THplus1.csv", ["T_H_plus_1", "gauge_algebra", "n_bases"],
         [(th + 1, g, n) for (th, g), n in
          sorted(patterns.items(), key=lambda kv: (kv[0][0], -kv[1]))])
    ts = sorted(count_t)
    _csv("by_T_summary.csv", ["T", "n_bases", "max_total_rank"],
         [(t, count_t[t], maxrank_t[t]) for t in ts])

    # ── plots (x = T = sig_neg), scatter ────────────────────────────────
    plt.figure(figsize=(10, 5))
    plt.scatter(ts, [count_t[t] for t in ts], color="orange", s=16)
    plt.yscale("log")
    plt.xlabel(r"$T\ (=\sigma_-)$")
    plt.ylabel("Number of non-Higgsable gravity blocks (log)")
    plt.title(r"Non-Higgsable gravity blocks per $T$")
    plt.grid(True, alpha=0.3)
    plt.tight_layout(); plt.savefig(os.path.join(out, "blocks_per_T.png"), dpi=150)
    plt.close()

    plt.figure(figsize=(10, 5))
    plt.scatter(ts, [maxrank_t[t] for t in ts], color="red", s=16)
    plt.xlabel(r"$T\ (=\sigma_-)$")
    plt.ylabel("Maximal total gauge rank")
    plt.title(r"Maximal total gauge-algebra rank per $T$")
    plt.grid(True, alpha=0.3)
    plt.tight_layout(); plt.savefig(os.path.join(out, "maxrank_per_T.png"), dpi=150)
    plt.close()

    ths = sorted(count_th)
    plt.figure(figsize=(10, 5))
    plt.scatter([th + 1 for th in ths], [count_th[th] for th in ths], color="orange", s=16)
    plt.yscale("log")
    plt.xlabel(r"$T_H + 1$")
    plt.ylabel("Number of non-Higgsable gravity blocks (log)")
    plt.title(r"Non-Higgsable gravity blocks per $T_H+1$")
    plt.grid(True, alpha=0.3)
    plt.tight_layout(); plt.savefig(os.path.join(out, "blocks_per_THplus1.png"), dpi=150)
    plt.close()

    print(f"DONE -> {out}/  ({total:,} bases, {len(ts)} T values)", flush=True)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    run(sys.argv[1], sys.argv[2])
