#!/usr/bin/env python3
"""Delta comparison T_H ∈ {1..5}."""
import csv, os, sys, re, glob
from collections import defaultdict
csv.field_size_limit(sys.maxsize)
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

sys.path.insert(0, BASE)
from gen_classification_csv import parse_cat_file, get_lst_gauge, get_sugra_gauge

T_RANGE = (1,2,3,4,5,6,7)
SCAN = []  # list of (root, name) tuples
# Prefer BASE (LDLT-fixed) over backup for same-named dirs
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
print(f"Scanning {len(SCAN)} dirs..."); sys.stdout.flush()

sj = defaultdict(int); total = 0
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
                sj[key] += 1; total += 1
        except: pass
print(f"  SJ entries: {total}, unique keys: {len(sj)}"); sys.stdout.flush()

yh = defaultdict(int)
for fn in sorted(glob.glob(os.path.join(BASE, "data/T*/Ext*.tsv"))):
    with open(fn) as f:
        r = csv.DictReader(f, delimiter="\t", quotechar='"')
        for row in r:
            T_H = int(row["TH+1"]) - 1
            if T_H not in T_RANGE: continue
            key = (T_H, content_key(row["LST gauge algebra"]), content_key(row["full gauge algebra"]), int(row["Delta"]))
            yh[key] += 1
print(f"  YH rows: {sum(yh.values())}, unique keys: {len(yh)}"); sys.stdout.flush()

common = set(sj) & set(yh); sj_only = set(sj) - set(yh); yh_only = set(yh) - set(sj)
print(f"\n## Per T_H (T_H, LST, SUGRA, Delta) match")
print(f"  T_H | SJ keys | YH keys | common | SJ-only | YH-only")
for t in T_RANGE:
    s = {k for k in sj if k[0]==t}; y = {k for k in yh if k[0]==t}
    print(f"  {t}   | {len(s):>7} | {len(y):>7} | {len(s&y):>6} | {len(s-y):>7} | {len(y-s):>7}")
print(f"\nTotal: common={len(common)}, SJ-only={len(sj_only)}, YH-only={len(yh_only)}")
if yh_only:
    print(f"\nYH-only examples (first 8):")
    for k in sorted(yh_only)[:8]:
        T_H, lst, sg, d = k
        print(f"  T_H={T_H} Delta={d:>5}  LST='{lst[:20]}'  SUGRA='{sg[:50]}'")
