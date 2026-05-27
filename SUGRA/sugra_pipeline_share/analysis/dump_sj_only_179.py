#!/usr/bin/env python3
"""Dump intersection forms for 179 SJ-only entries (true c_ext ≤ 16).
Matches the prior intersection_form_sj_only_no223_yuta.txt format."""
import csv, glob, os, re, sys
csv.field_size_limit(sys.maxsize)
import numpy as np
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def gauge_content(s):
    s = (s or "").strip()
    if s in ("(none)", "0", "", "none"): return {}
    s = re.sub(r"su\((\d+)\)", r"su\1", s); s = re.sub(r"so\((\d+)\)", r"so\1", s)
    s = re.sub(r"sp\((\d+)\)", r"sp\1", s)
    s = s.replace("sp1","su2"); s = re.sub(r"\bf5\b","f4",s); s = s.replace("e7'","e7p")
    def expand(s):
        result = {}; tokens, depth, cur = [], 0, ""
        for c in s + "+":
            if c == "(": depth += 1; cur += c
            elif c == ")": depth -= 1; cur += c
            elif c == "+" and depth == 0:
                if cur.strip(): tokens.append(cur.strip())
                cur = ""
            else: cur += c
        for tok in tokens:
            m = re.match(r"^(\d+)\s*\((.+)\)$", tok)
            if m:
                for g, c in expand(m.group(2)).items(): result[g] = result.get(g, 0) + c * int(m.group(1))
                continue
            m = re.match(r"^(\d+)\s+(.+)$", tok)
            if m: result[m.group(2).strip()] = result.get(m.group(2).strip(), 0) + int(m.group(1)); continue
            result[tok] = result.get(tok, 0) + 1
        return result
    return expand(s)

def content_key(s):
    d = gauge_content(s)
    if not d: return "(none)"
    return " + ".join(f"{c} {n}" if c > 1 else n for n, c in sorted(d.items()))

GAUGE_DIM_H = {
    "su2": (3, 2), "su3": (8, 3), "g2": (14, 4), "so7": (21, 5),
    "so8": (28, 6), "f4": (52, 9), "e6": (78, 12), "e7": (133, 18),
    "e7'": (133, 18), "e7p": (133, 18), "e8": (248, 30),
    "su8": (63, 8), "so16": (120, 14),
}

def find_kernel_positive_integer(A):
    A = np.array(A, dtype=float)
    n = A.shape[0]
    U, s, Vt = np.linalg.svd(A)
    tol = max(A.shape) * np.spacing(np.max(s)) if len(s) > 0 else 1e-10
    ker_vecs = Vt[s < tol]
    if len(ker_vecs) == 0:
        if s[-1] < 1e-8: ker_vecs = Vt[-1:]
        else: return None
    for v in ker_vecs:
        v = np.array(v, dtype=float)
        if np.all(v < -1e-9): v = -v
        if not np.all(v > -1e-9): continue
        v = np.abs(v)
        nonzero = v[v > 1e-9]
        if len(nonzero) == 0: continue
        scale = 1.0 / nonzero.min()
        scaled = v * scale
        if not np.all(np.abs(scaled - np.round(scaled)) < 1e-6): continue
        x = np.round(scaled).astype(int)
        if np.gcd.reduce(x[x > 0]) != 1: x = x // np.gcd.reduce(x[x > 0])
        if np.all(x > 0):
            return x
    return None

def per_curve_gauge_from_nhc(M, curves):
    out = {}
    non1 = [c for c in curves if M[c][c] != -1]
    adj = {c: [c2 for c2 in non1 if c != c2 and M[c][c2] != 0] for c in non1}
    visited = set()
    for c in non1:
        if c in visited: continue
        cluster = []
        q = [c]; visited.add(c)
        while q:
            x = q.pop()
            cluster.append(x)
            for y in adj[x]:
                if y not in visited:
                    visited.add(y); q.append(y)
        cluster.sort()
        sis = [M[x][x] for x in cluster]
        if len(cluster) == 1:
            si = sis[0]
            if si == -3: out[cluster[0]] = "su3"
            elif si == -4: out[cluster[0]] = "so8"
            elif si == -5: out[cluster[0]] = "f4"
            elif si == -6: out[cluster[0]] = "e6"
            elif si == -7: out[cluster[0]] = "e7"
            elif si == -8: out[cluster[0]] = "e7"
            elif si == -12: out[cluster[0]] = "e8"
        elif len(cluster) == 2 and sorted(sis) == [-3, -2]:
            for x in cluster:
                if M[x][x] == -2: out[x] = "su2"
                elif M[x][x] == -3: out[x] = "g2"
        elif len(cluster) == 2 and sorted(sis) == [-4, -2]:
            for x in cluster:
                if M[x][x] == -2: out[x] = "su8"
                elif M[x][x] == -4: out[x] = "so16"
        elif len(cluster) == 3 and sorted(sis) == [-3, -2, -2]:
            m3 = [x for x in cluster if M[x][x] == -3][0]
            deg = sum(1 for y in cluster if M[m3][y] != 0 and y != m3)
            if deg == 2:
                out[m3] = "so7"
                for x in cluster:
                    if M[x][x] == -2: out[x] = "su2"
            else:
                out[m3] = "g2"
                mid = [y for y in cluster if M[m3][y] != 0 and y != m3 and M[y][y] == -2]
                if mid: out[mid[0]] = "su2"
    return out

def compute_true_cext(M, ext_idx, lst_idx):
    """Returns (true_c_ext, breakdown_list)."""
    lst_idx_sorted = sorted(lst_idx)
    A = np.array([[M[i][j] for j in lst_idx_sorted] for i in lst_idx_sorted], dtype=int)
    X = find_kernel_positive_integer(A)
    if X is None: return None, []
    x_map = {lst_idx_sorted[i]: int(X[i]) for i in range(len(lst_idx_sorted))}
    all_idx = sorted(set(ext_idx) | set(lst_idx))
    gauges = per_curve_gauge_from_nhc(M, all_idx)
    total = 0.0
    breakdown = []
    for e in sorted(ext_idx):
        g = gauges.get(e)
        if g not in GAUGE_DIM_H: continue
        k_eff = sum(M[e][i] * x_map.get(i, 0) for i in lst_idx_sorted)
        if k_eff <= 0: continue
        dim, h = GAUGE_DIM_H[g]
        c = k_eff * dim / (k_eff + h)
        total += c
        breakdown.append((e, g, k_eff, round(c, 3)))
    return total, breakdown

sys.path.insert(0, BASE)
from gen_classification_csv import parse_cat_file, get_lst_gauge, get_sugra_gauge

T_RANGE = (1,2,3,4,5,6,7)
SCAN = []
_seen = set()
for root in [BASE, os.path.join(BASE, "backup_t010_chain_v4")]:
    if not os.path.isdir(root): continue
    for name in sorted(os.listdir(root)):
        p = os.path.join(root, name)
        if not os.path.isdir(p) or not name.startswith("cat_"): continue
        if "nonsugra" in name or "backup" in name or "merged" in name: continue
        if name in _seen: continue
        if name == "cat_1ext_nosu2_t010" or name == "cat_nhc_ext":
            SCAN.append((root, name)); _seen.add(name)
        elif (name.startswith("cat_nhc_ext_unified_r") or name.startswith("cat_nhc_ext_nhc_ext_unified_r")) and "sext" in name:
            SCAN.append((root, name)); _seen.add(name)

yh_keys = set()
for fn in sorted(glob.glob(os.path.join(BASE, "data/T*/Ext*.tsv"))):
    with open(fn) as f:
        r = csv.DictReader(f, delimiter="\t", quotechar='"')
        for row in r:
            T_H = int(row["TH+1"]) - 1
            if T_H not in T_RANGE: continue
            key = (T_H, content_key(row["LST gauge algebra"]), content_key(row["full gauge algebra"]), int(row["Delta"]))
            yh_keys.add(key)

sj_repr = {}
sj_count = {}
for root, dname in SCAN:
    for cf in os.listdir(os.path.join(root, dname)):
        if not cf.endswith(".cat"): continue
        try:
            for entry in parse_cat_file(os.path.join(root, dname, cf)):
                if any(ext["tag"] in ("hat1m1","hat1m2","su8","su8mix","so16n2","so16n2mix","nhc_2_4") or ext["is_hat1"] for ext in entry["externals"]): continue
                T_H = entry.get("base_T", 0)
                if T_H not in T_RANGE: continue
                T = entry["T"]; Hc = entry["H_charged"]; V = entry["V"]
                delta = Hc - V + 29 * T
                key = (T_H, content_key(get_lst_gauge(entry)), content_key(get_sugra_gauge(entry)), delta)
                if key in yh_keys: continue
                sj_count[key] = sj_count.get(key, 0) + 1
                if key in sj_repr: continue
                sj_repr[key] = (entry, dname, cf)
        except Exception:
            pass

# Header
print(f"# SJ-only IFs (vs Hamada updated data) — T_H ∈ {{1,2,3,4,5,6,7}}")
print(f"# Excluded externals: hat1, su8, so16, nhc_2_4")
print(f"# c_ext computed as fiber-weighted (X ∈ ker(LST IF), k_eff·dim/(k_eff+h))")
print(f"# Filter: true c_ext ≤ 16")
print(f"# Total unique keys: {len(sj_repr)}")
print()

idx = 0
for key in sorted(sj_repr.keys()):
    entry, dname, cf = sj_repr[key]
    M = entry.get("if_matrix")
    if M is None: continue
    ext_idx_set = set()
    for e in entry["externals"]:
        if e.get("curve_idx") is not None:
            ext_idx_set.add(e["curve_idx"])
    if not ext_idx_set and "attach_if_size" in entry:
        ext_idx_set.add(len(M) - 1)
    n = len(M)
    lst_idx_set = [i for i in range(n) if i not in ext_idx_set]
    cval, bd = compute_true_cext(M, list(ext_idx_set), lst_idx_set)
    if cval is None or cval > 16 + 1e-6: continue

    idx += 1
    T_H, lst_ck, sg_ck, delta = key
    T = entry["T"]; Hc = entry["H_charged"]; V = entry["V"]; Hn = entry["H_neutral"]
    bases = sj_count[key]
    ext_sorted = sorted(ext_idx_set)
    lst_sorted = sorted(lst_idx_set)
    bd_str = ", ".join(f"c{e}({g},k={k})→{c}" for e, g, k, c in bd)

    print(f"# T_H={T_H} T={T} Δ={delta} LST={lst_ck} SUGRA={sg_ck} c_ext={cval:.2f} H={Hc} V={V} Hn={Hn} bases={bases}")
    print(f"# ext_idx={ext_sorted} lst_idx={lst_sorted}  c_ext breakdown: {bd_str}")
    print(f"# source: {dname}/{cf}")
    print(n)
    for row in M:
        print(" ".join(str(x) for x in row))
    print()

# (Extended IF with b0 row/col — omitted; depends on hat1 set which we don't track here)
