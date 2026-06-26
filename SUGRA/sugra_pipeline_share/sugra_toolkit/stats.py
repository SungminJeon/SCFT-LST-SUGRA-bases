"""stats.py — publication summary tables for the SUGRA base catalog.

LSTGauge (the base's gauge algebra) is a faithful Python port of
catalog_analysis.m / SUGRACatalog.m (NHCLookup): BFS the non-(-1),
NON-external curves into connected components and match each (by eigenvalue
spectrum) to the NHC table; pure-(-2) chains carry none. Result is the sorted
"+"-joined list of gauge algebras (-3->su3, -2-3->su2+g2, -2-2-3->g2+sp1,
-2-3-2->su2+so7+su2, -2-4->su8+so16, ...).

Tables written (CSV):
  by_T_H.csv               T_H                              -> n_bases
  by_combo.csv             external combo (file name)       -> n_bases
  by_LST.csv               LST (type, T_H, id)              -> n_bases
  by_LSTGauge.csv          LSTGauge (gauge algebra)         -> n_bases
  classification.csv       (T_H, LSTGauge, combo)           -> n_bases

USAGE
-----
    python stats.py <catalog_root> [out_dir]
"""
from __future__ import annotations

import csv
import os
import re
import sys
from collections import Counter

import numpy as np

# ── NHC reference table (SUGRACatalog.m nhcTable) ─────────────────────────
def _nhc_if(self_ints, bonds):
    n = len(self_ints)
    M = np.zeros((n, n))
    for i in range(n):
        M[i, i] = self_ints[i]
        if i + 1 < n:
            b = bonds[i] if i < len(bonds) else 1
            M[i, i + 1] = M[i + 1, i] = b
    return M

_NHC_DEFS = [
    ("-3", [-3], []), ("-4", [-4], []), ("-5", [-5], []), ("-6", [-6], []),
    ("-7", [-7], []), ("-8", [-8], []), ("-12", [-12], []),
    ("-2-3", [-2, -3], []), ("-2-4", [-2, -4], [2]),
    ("-2-3-2", [-2, -3, -2], []), ("-2-2-3", [-2, -2, -3], []),
]
_NHC_EIGS = [(name, np.sort(np.linalg.eigvalsh(_nhc_if(si, b))))
             for name, si, b in _NHC_DEFS]

# NHC cluster name -> gauge algebra factors (the lone -2 in -2-2-3 carries none)
_NHC_ALG = {
    "-3": ["su3"], "-4": ["so8"], "-5": ["f4"], "-6": ["e6"],
    "-7": ["e7"], "-8": ["e7"], "-12": ["e8"],
    "-2-3": ["su2", "g2"], "-2-4": ["su8", "so16"],
    "-2-3-2": ["su2", "so7", "su2"], "-2-2-3": ["sp1", "g2"],
}


def _nhc_lookup(subIF: np.ndarray) -> str | None:
    eigs = np.sort(np.linalg.eigvalsh(subIF.astype(float)))
    for name, ref in _NHC_EIGS:
        if len(ref) == len(eigs) and np.allclose(eigs, ref, atol=1e-6):
            return name
    return None


def lst_gauge(IF: np.ndarray, externals: list[dict]) -> str:
    n = len(IF)
    if n == 0:
        return ""
    ext_idx = {e["curveIdx"] for e in externals if e.get("curveIdx", -1) >= 0}
    diag = [int(IF[i][i]) for i in range(n)]
    non_m1 = [i for i in range(n) if diag[i] != -1 and i not in ext_idx]
    seen, factors = set(), []
    for s in non_m1:
        if s in seen:
            continue
        comp, queue = [], [s]; seen.add(s)
        while queue:
            cur = queue.pop(0); comp.append(cur)
            for nb in non_m1:
                if nb not in seen and IF[cur][nb] != 0:
                    seen.add(nb); queue.append(nb)
        if all(diag[c] == -2 for c in comp):       # pure -2: no gauge
            continue
        sub = IF[np.ix_(sorted(comp), sorted(comp))]
        name = _nhc_lookup(sub)
        factors.extend(_NHC_ALG[name] if name is not None else ["?"])
    return "+".join(sorted(factors))


# ── display relabeling of external tags in combo names (output only) ──────
# Data is never modified; only how combos are written in the tables.
DISPLAY_RENAME = {"su3n2": "g2", "su3n2mix": "g2mix"}


def _relabel_combo(combo: str) -> str:
    return "+".join(sorted(DISPLAY_RENAME.get(t, t) for t in combo.split("+")))


# ── catalog scan ──────────────────────────────────────────────────────────
ENTRY_RE = re.compile(r"ENTRY\b.*?\nEND", re.S)


def _bases(path: str):
    text = open(path, encoding="utf-8", errors="replace").read()
    for blk in ENTRY_RE.findall(text):
        lines = blk.split("\n")
        ext, IF = [], None
        T_H = cat_id = sig_neg = 0
        cat_type = ""
        i = 0
        while i < len(lines):
            ln = lines[i]
            if ln.startswith("EXTERNALS"):
                k = int(ln.split()[1])
                for j in range(k):
                    p = lines[i + 1 + j].split()
                    if len(p) >= 7:
                        ext.append({"curveIdx": int(p[0]), "extSI": int(p[3]),
                                    "targetSI": int(p[4]), "isHat1": p[6] == "1"})
                i += k
            elif ln.startswith("SPEC"):              # single-external (e.g. hat1)
                p = ln.split()
                if len(p) >= 7:
                    ext.append({"curveIdx": -1, "extSI": int(p[3]),
                                "targetSI": int(p[4]), "isHat1": p[6] == "1"})
            elif ln.startswith("ATTACH"):            # fills the last SPEC's curveIdx
                p = ln.split()
                if len(p) >= 3 and ext and ext[-1]["curveIdx"] == -1:
                    ext[-1]["curveIdx"] = int(p[2])
            elif ln.startswith("BASE"):
                p = ln.split()
                if len(p) >= 5:
                    cat_id, cat_type, T_H = int(p[1]), p[2], int(p[3])
            elif ln.startswith("PHYSICS"):           # T Hc V Hn det sig+ sig- sig0
                p = ln.split()
                if len(p) >= 8:
                    sig_neg = int(p[7])              # T = sig_neg (extended base)
            elif ln.startswith("IF"):
                m = int(ln[3:].strip())
                IF = np.array([[int(x) for x in lines[i + 1 + r].split()]
                               for r in range(m)], dtype=np.int64)
                i += m
            i += 1
        yield ext, IF, T_H, cat_type, cat_id, sig_neg


def run(root: str, out: str):
    os.makedirs(out, exist_ok=True)
    by_th, by_combo, by_lst = Counter(), Counter(), Counter()
    by_lg, classif = Counter(), Counter()
    total = 0
    tdirs = sorted((d for d in os.listdir(root) if re.match(r"T\d+(_\d+)?$", d)),
                   key=lambda s: int(re.match(r"T(\d+)", s).group(1)))
    for td in tdirs:
        tp = os.path.join(root, td)
        for fn in sorted(os.listdir(tp)):
            if not fn.endswith(".cat"):
                continue
            combo = _relabel_combo(fn[:-4])
            for ext, IF, T_H, ct, cid, sig_neg in _bases(os.path.join(tp, fn)):
                total += 1
                lg = lst_gauge(IF, ext) if IF is not None else ""
                by_th[T_H] += 1
                by_combo[combo] += 1
                by_lst[(ct, T_H, cid)] += 1
                by_lg[lg or "(none)"] += 1
                classif[(T_H, lg or "(none)", combo)] += 1
        print(f"  {td} done", flush=True)

    _csv(out, "by_T_H.csv", ["T_H", "n_bases"], sorted(by_th.items()))
    _csv(out, "by_combo.csv", ["combo", "n_bases"],
         sorted(by_combo.items(), key=lambda kv: -kv[1]))
    _csv(out, "by_LST.csv", ["LST_type", "T_H", "LST_id", "n_bases"],
         [(t, th, i, n) for (t, th, i), n in sorted(by_lst.items())])
    _csv(out, "by_LSTGauge.csv", ["LSTGauge", "n_bases"],
         sorted(by_lg.items(), key=lambda kv: -kv[1]))
    _csv(out, "classification.csv", ["T_H", "LSTGauge", "combo", "n_bases"],
         [(th, lg, c, n) for (th, lg, c), n in
          sorted(classif.items(), key=lambda kv: (kv[0][0], -kv[1]))])
    print(f"\nDONE -> {out}/  ({total:,} bases)")


def _csv(out, name, header, rows):
    with open(os.path.join(out, name), "w", newline="") as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)


if __name__ == "__main__":
    if not 2 <= len(sys.argv) <= 3:
        print(__doc__); raise SystemExit("usage: python stats.py <catalog_root> [out_dir]")
    run(sys.argv[1], sys.argv[2] if len(sys.argv) == 3 else "stats_out")
