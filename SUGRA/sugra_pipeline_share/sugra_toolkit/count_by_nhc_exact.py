#!/usr/bin/env python3
"""IF-exact NHC-cluster census — the exact counterpart of blocks_by_nhc.csv.

blocks_by_nhc.csv counts blocks by the NHC *tag* attached (nhc_2_3 / nhc_2_2_3 /
nhc_2_3_2). But an attached -2-3 (nhc_2_3) cluster can GROW in a later chain step
into a -2-3-2 (su2+so7+su2) or -2-2-3 (sp1+g2). This tool reads the actual
intersection form for the nhc_2_3 external (its connected non -1 component,
identified via the eigenvalue spectrum, stats._nhc_lookup) and SPLITS the nhc_2_3
row into its three IF-exact outcomes. nhc_2_2_3 and nhc_2_3_2 are separate attach
modes (whole clusters glued directly), so they are left as their tag counts.

Usage:
    python count_by_nhc_exact.py <catalog> [out_dir]     # default out_dir: stats_final

Writes:
    blocks_by_nhc_exact.csv          nhc_2_3 split into stays-2-3 / ->2-3-2 / ->2-2-3
                                     (#blocks); nhc_2_2_3, nhc_2_3_2 kept as tag counts
    nhc_2_3_outcome_by_THplus1.csv   for nhc_2_3 attachments: clusters per T_H+1, by outcome
    nhc_exact_migration.csv          nhc_2_3 -> exact cluster -> #clusters
    nhc_2_3_outcome_by_THplus1.png/.pdf
"""
import sys, os, glob, csv
from collections import Counter, defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sugra
import stats as statsmod

NHC_TAGS = {"nhc_2_3", "nhc_2_2_3", "nhc_2_3_2"}
GAUGE = {"-2-3": "su2+g2", "-2-3-2": "su2+so7+su2", "-2-2-3": "sp1+g2", None: "(other/grown)"}
CURVES = {"-2-3": "(-2,-3)", "-2-3-2": "(-2,-3,-2)", "-2-2-3": "(-2,-2,-3)", None: "-"}


def component_of(IF, start, m1, n):
    seen = {start}; st = [start]; comp = []
    while st:
        c = st.pop(); comp.append(c)
        for j in range(n):
            if j not in m1 and j not in seen and IF[c][j] != 0:
                seen.add(j); st.append(j)
    return comp


def run(cat_dir, out="stats_final"):
    os.makedirs(out, exist_ok=True)

    grand_total = 0   # all blocks in the catalog (for % of total)
    for f in glob.glob(os.path.join(cat_dir, "T*", "*.cat")):
        with open(f) as fh:
            grand_total += sum(1 for ln in fh if ln.startswith("ENTRY"))

    cat = sugra.Catalog(cat_dir)
    nhc23_out = Counter()                    # outcome of nhc_2_3 -> #blocks with >=1 such
    tag_blocks = Counter()                   # nhc_2_2_3 / nhc_2_3_2 -> #blocks (by combo tag)
    tag_exact = Counter()                    # (origin tag, exact name) -> #clusters
    thp1 = defaultdict(Counter)              # T_H+1 -> Counter(outcome) for nhc_2_3 origin
    n_nhc_blocks = 0

    for th in cat.tensors():
        for cb in cat.combos(th):
            cb_tags = set(cb.split("+"))
            if not (cb_tags & NHC_TAGS):
                continue
            other_tags = [t for t in ("nhc_2_2_3", "nhc_2_3_2") if t in cb_tags]
            has_23 = "nhc_2_3" in cb_tags
            for b in cat.load(th, cb):
                n_nhc_blocks += 1
                for t in other_tags:                       # separate attach modes: tag-based, untouched
                    tag_blocks[t] += 1
                if not has_23:
                    continue
                IF = b.IF; n = len(IF)
                diag = [int(IF[i][i]) for i in range(n)]
                m1 = {i for i in range(n) if diag[i] == -1}
                nhc23_curves = [e["curveIdx"] for e in b.externals
                                if e.get("tag") == "nhc_2_3" and 0 <= e.get("curveIdx", -1) < n]
                seen, outcomes = set(), set()
                for cidx in nhc23_curves:
                    if cidx in seen:
                        continue
                    comp = component_of(IF, cidx, m1, n)
                    seen.update(comp)
                    name = statsmod._nhc_lookup(IF[np.ix_(sorted(comp), sorted(comp))])
                    tag_exact[("nhc_2_3", name)] += 1
                    outcomes.add(name)
                    thp1[th + 1][name] += 1
                for nm in outcomes:                        # block counted once per distinct outcome
                    nhc23_out[nm] += 1
        print(f"  T_H={th} done", flush=True)

    # blocks_by_nhc.csv with the nhc_2_3 row SPLIT into its 3 IF-exact outcomes;
    # nhc_2_2_3 / nhc_2_3_2 are separate attach modes -> left as tag-based counts.
    def pct(x):
        return f"{100*x/grand_total:.4f}"
    with open(os.path.join(out, "blocks_by_nhc_exact.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["nhc_cluster", "curves", "gauge", "n_blocks", "percent_of_total", "source"])
        for nm, label in [("-2-3", "nhc_2_3 (stays -2-3)"),
                          ("-2-3-2", "nhc_2_3 -> -2-3-2"),
                          ("-2-2-3", "nhc_2_3 -> -2-2-3")]:
            v = nhc23_out.get(nm, 0)
            w.writerow([label, CURVES[nm], GAUGE[nm], v, pct(v), "nhc_2_3 attach (IF-exact)"])
        if nhc23_out.get(None):
            w.writerow(["nhc_2_3 -> (other)", "-", "-", nhc23_out[None], pct(nhc23_out[None]), "nhc_2_3 attach (IF-exact)"])
        for t, nm in [("nhc_2_2_3", "-2-2-3"), ("nhc_2_3_2", "-2-3-2")]:
            v = tag_blocks.get(t, 0)
            w.writerow([t, CURVES[nm], GAUGE[nm], v, pct(v), "separate attach mode (tag)"])
        w.writerow(["ANY_nhc", "-", "-", n_nhc_blocks, pct(n_nhc_blocks), ""])

    ths = sorted(thp1)
    with open(os.path.join(out, "nhc_2_3_outcome_by_THplus1.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["T_H_plus_1", "as_2_3", "as_2_3_2", "as_2_2_3", "as_other"])
        for t in ths:
            c = thp1[t]
            w.writerow([t, c.get("-2-3", 0), c.get("-2-3-2", 0), c.get("-2-2-3", 0), c.get(None, 0)])

    with open(os.path.join(out, "nhc_exact_migration.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["origin_tag", "exact_cluster", "gauge", "n_clusters"])
        for (tag, nm), k in sorted(tag_exact.items(), key=lambda kv: (kv[0][0], -kv[1])):
            w.writerow([tag, nm if nm else "(other)", GAUGE.get(nm, "?"), k])

    plt.figure(figsize=(10, 5))
    series = [("-2-3", "stays 2-3 (su2+g2)", "tab:blue"),
              ("-2-3-2", "grew to 2-3-2 (su2+so7+su2)", "tab:red"),
              ("-2-2-3", "grew to 2-2-3 (sp1+g2)", "tab:green")]
    for nm, lab, col in series:
        ys = [thp1[t].get(nm, 0) for t in ths]
        plt.scatter([t for t in ths], ys, s=16, color=col, label=lab)
    plt.yscale("log")
    plt.xlabel(r"$T_H + 1$")
    plt.ylabel("Number of nhc_2_3 clusters (log)")
    plt.title(r"Fate of attached $-2{-}3$ clusters per $T_H+1$")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out, "nhc_2_3_outcome_by_THplus1.png"), dpi=150)
    plt.savefig(os.path.join(out, "nhc_2_3_outcome_by_THplus1.pdf"))
    plt.close()

    print(f"\nblocks with >=1 NHC external: {n_nhc_blocks:,} / {grand_total:,}")
    print("nhc_2_3 outcome (blocks, >=1):", dict(nhc23_out))
    print("nhc_2_2_3 / nhc_2_3_2 (tag blocks):", dict(tag_blocks))
    print("nhc_2_3 outcome (clusters):", {k[1]: v for k, v in tag_exact.items()})


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "stats_final")
