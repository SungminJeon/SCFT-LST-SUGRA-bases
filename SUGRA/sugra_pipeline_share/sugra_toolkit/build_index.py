"""build_index.py — build CSV summary tables from a clean catalog tree.

Lets people filter/query the catalog with no Python and no parsing of the
intersection-form matrices — just open the CSV in pandas, Excel, or any tool.

Produces two files in the catalog root:

  index.csv     one row per base:
                T_H, T, combo, n_ext, Hc, V, Hn, det, sig_pos, sig_neg,
                sig_zero, n_curves, file            [+ t_min with --tmin]
                (large: one row per base; best for pandas.)

  summary.csv   one row per (T_H, combo):
                T_H, combo, n_bases, T   (small; Excel-friendly overview.)

USAGE
-----
    python build_index.py catalog                 # -> catalog/index.csv, summary.csv
    python build_index.py catalog --tmin          # also compute T_min (slower)

Then, for example:
    import pandas as pd
    df = pd.read_csv("catalog/index.csv")
    df[(df.T_H == 8) & (df.combo == "e6+su3")]
    df[df.combo.str.contains("e8")].groupby("T_H").size()
"""
from __future__ import annotations

import csv
import os
import re
import sys

ENTRY_RE = re.compile(r"ENTRY\b.*?\nEND", re.S)


def _iter_bases(path: str):
    """Yield dicts of the cheap (already-stored) fields for each base."""
    text = open(path, encoding="utf-8", errors="replace").read()
    for blk in ENTRY_RE.findall(text):
        rec = {"combo": "", "n_ext": 0, "T": 0, "Hc": 0, "V": 0, "Hn": 0,
               "det": 0, "sig_pos": 0, "sig_neg": 0, "sig_zero": 0,
               "T_H": 0, "n_curves": 0, "IF": None, "hat1": []}
        lines = blk.split("\n")
        i = 0
        while i < len(lines):
            ln = lines[i]
            if ln.startswith("COMBO"):
                rec["combo"] = ln[6:].strip()
            elif ln.startswith("EXTERNALS"):
                k = int(ln.split()[1]); rec["n_ext"] = k
                for j in range(k):
                    p = lines[i + 1 + j].split()
                    if len(p) >= 7 and p[6] == "1":
                        rec["hat1"].append(int(p[0]))
                i += k
            elif ln.startswith("BASE"):
                p = ln.split()
                if len(p) >= 5:
                    rec["T_H"] = int(p[3]); rec["n_curves"] = int(p[4])
            elif ln.startswith("PHYSICS"):
                p = ln.split()
                if len(p) >= 9:
                    (rec["T"], rec["Hc"], rec["V"], rec["Hn"], rec["det"],
                     rec["sig_pos"], rec["sig_neg"], rec["sig_zero"]) = (int(x) for x in p[1:9])
            elif ln.startswith("IF") and rec["IF"] is None:
                n = int(ln[3:].strip()); mat = []
                for r in range(n):
                    mat.append([int(x) for x in lines[i + 1 + r].split()])
                rec["IF"] = mat; i += n
            i += 1
        if not rec["combo"]:
            rec["combo"] = "?"
        yield rec


def _t_min(IF, hat1):
    import numpy as np
    n = len(IF)
    if n == 0:
        return ""
    M = np.zeros((n + 1, n + 1)); M[:n, :n] = IF
    h1 = set(hat1)
    for i in range(n):
        b = -1.0 if i in h1 else float(IF[i][i] + 2)
        M[i, n] = M[n, i] = b
    M[n, n] = 0.0; B = float(np.linalg.det(M))
    M[n, n] = 1.0; A = float(np.linalg.det(M)) - B
    if abs(A) < 1e-10:
        return ""
    return int(__import__("math").ceil(9.0 - (-B / A) - 1e-9))


def build(root: str, with_tmin: bool = False):
    tdirs = sorted((d for d in os.listdir(root) if re.match(r"T\d+(_\d+)?$", d)),
                   key=lambda s: int(re.match(r"T(\d+)", s).group(1)))
    cols = ["T_H", "T", "combo", "n_ext", "Hc", "V", "Hn", "det",
            "sig_pos", "sig_neg", "sig_zero", "n_curves", "file"]
    if with_tmin:
        cols.insert(2, "t_min")
    index_path = os.path.join(root, "index.csv")
    summ: dict = {}
    nrows = 0
    with open(index_path, "w", newline="") as fout:
        w = csv.writer(fout); w.writerow(cols)
        for td in tdirs:
            tpath = os.path.join(root, td)
            for fn in sorted(os.listdir(tpath)):
                if not fn.endswith(".cat"):
                    continue
                combo_name = fn[:-4]          # filename IS the canonical combo
                for rec in _iter_bases(os.path.join(tpath, fn)):
                    rec["combo"] = combo_name
                    row = [rec["T_H"], rec["T"], rec["combo"], rec["n_ext"],
                           rec["Hc"], rec["V"], rec["Hn"], rec["det"],
                           rec["sig_pos"], rec["sig_neg"], rec["sig_zero"],
                           rec["n_curves"], fn]
                    if with_tmin:
                        row.insert(2, _t_min(rec["IF"], rec["hat1"]))
                    w.writerow(row)
                    nrows += 1
                    key = (rec["T_H"], rec["combo"])
                    s = summ.setdefault(key, [0, rec["T"]])
                    s[0] += 1
            print(f"  {td} done", flush=True)
    with open(os.path.join(root, "summary.csv"), "w", newline="") as fout:
        w = csv.writer(fout); w.writerow(["T_H", "combo", "n_bases", "T"])
        for (t_h, combo), (n, T) in sorted(summ.items()):
            w.writerow([t_h, combo, n, T])
    print(f"\nDONE: {index_path} ({nrows:,} rows), summary.csv ({len(summ):,} rows)")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    if len(args) != 1:
        print(__doc__); raise SystemExit("usage: python build_index.py <catalog_root> [--tmin]")
    build(args[0], with_tmin=("--tmin" in sys.argv))
