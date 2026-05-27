#!/usr/bin/env python3
"""Compare SJ pipeline catalog with Hamada (YH) TSV reference data.

Match key:  (T_H, LST gauge content, SUGRA gauge content, Delta)
            where Delta = H_charged - V + 29 T.

Usage:
    python3 compare_hamada.py                  # default T_H range = 1..7
    python3 compare_hamada.py T_MIN T_MAX      # restrict to specific range

The SJ side scans the parent directory's `cat_*` catalogs (Phase 1 / NHC ext /
Phase 2 chain output). The YH side reads `data/T*/Ext*.tsv` (Hamada's TSV
files). Provide Hamada's data files in `<pipeline_root>/data/T<NNN>/Ext<N>.tsv`
to enable the YH comparison.
"""
import csv, glob, os, re, sys
from collections import defaultdict
csv.field_size_limit(sys.maxsize)

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ── CLI / parameters ────────────────────────────────────────────
if len(sys.argv) == 1:
    T_RANGE = tuple(range(1, 8))            # default 1..7
elif len(sys.argv) == 3:
    T_RANGE = tuple(range(int(sys.argv[1]), int(sys.argv[2]) + 1))
else:
    print(__doc__); sys.exit(1)

# Externals that we exclude from the comparison (Hamada's tables don't
# enumerate these as separate single-curve externals).
EXCLUDED_EXT_TAGS = {"hat1m1", "hat1m2", "su8", "su8mix",
                      "so16n2", "so16n2mix", "nhc_2_4"}

# ── Gauge content parsing ────────────────────────────────────────
def gauge_content(s):
    s = (s or "").strip()
    if s in ("(none)", "0", "", "none"): return {}
    s = re.sub(r"su\((\d+)\)", r"su\1", s); s = re.sub(r"so\((\d+)\)", r"so\1", s)
    s = re.sub(r"sp\((\d+)\)", r"sp\1", s)
    s = s.replace("sp1", "su2"); s = re.sub(r"\bf5\b", "f4", s); s = s.replace("e7'", "e7p")
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
                for g, c in expand(m.group(2)).items():
                    result[g] = result.get(g, 0) + c * int(m.group(1))
                continue
            m = re.match(r"^(\d+)\s+(.+)$", tok)
            if m:
                result[m.group(2).strip()] = result.get(m.group(2).strip(), 0) + int(m.group(1))
                continue
            result[tok] = result.get(tok, 0) + 1
        return result
    return expand(s)

def content_key(s):
    d = gauge_content(s)
    if not d: return "(none)"
    return " + ".join(f"{c} {n}" if c > 1 else n for n, c in sorted(d.items()))

# ── SJ catalog scan ──────────────────────────────────────────────
sys.path.insert(0, BASE)
from gen_classification_csv import parse_cat_file, get_lst_gauge, get_sugra_gauge

def is_sj_catalog(name):
    """Pick directories that hold pipeline SUGRA output (exclude nonsugra,
       merged-temp, and other internal staging)."""
    if not name.startswith("cat_"): return False
    if "nonsugra" in name or "merged" in name: return False
    if name.startswith("cat_1ext_nosu2_"): return True
    if name.startswith("cat_nhc_ext_") and "_unified_" not in name: return True
    # phase2 chain output dirs end in "sext" or "nsext"
    if "_unified_" in name and name.endswith("sext"): return True
    return False

SCAN = [name for name in sorted(os.listdir(BASE)) if is_sj_catalog(name)]
print(f"Scanning {len(SCAN)} SJ catalog directories (T_H range {T_RANGE[0]}..{T_RANGE[-1]})...")

sj = defaultdict(int); total = 0
for dname in SCAN:
    for cf in os.listdir(os.path.join(BASE, dname)):
        if not cf.endswith(".cat"): continue
        try:
            for entry in parse_cat_file(os.path.join(BASE, dname, cf)):
                if any(ext["tag"] in EXCLUDED_EXT_TAGS or ext["is_hat1"]
                       for ext in entry["externals"]):
                    continue
                T_H = entry.get("base_T", 0)
                if T_H not in T_RANGE: continue
                T = entry["T"]; Hc = entry["H_charged"]; V = entry["V"]
                delta = Hc - V + 29 * T
                key = (T_H,
                       content_key(get_lst_gauge(entry)),
                       content_key(get_sugra_gauge(entry)),
                       delta)
                sj[key] += 1; total += 1
        except Exception:
            pass
print(f"  SJ entries: {total}, unique keys: {len(sj)}")

# ── YH (Hamada) TSV scan ─────────────────────────────────────────
yh = defaultdict(int)
yh_files = sorted(glob.glob(os.path.join(BASE, "data/T*/Ext*.tsv")))
if not yh_files:
    print(f"\n(No data/T*/Ext*.tsv found under {BASE}. "
          "Place Hamada's TSV files there to enable comparison.)")
    sys.exit(0)

for fn in yh_files:
    with open(fn) as f:
        r = csv.DictReader(f, delimiter="\t", quotechar='"')
        for row in r:
            try: T_H = int(row["TH+1"]) - 1
            except Exception: continue
            if T_H not in T_RANGE: continue
            try:
                key = (T_H,
                       content_key(row["LST gauge algebra"]),
                       content_key(row["full gauge algebra"]),
                       int(row["Delta"]))
            except Exception:
                continue
            yh[key] += 1
print(f"  YH rows: {sum(yh.values())}, unique keys: {len(yh)}")

# ── Compare ──────────────────────────────────────────────────────
common  = set(sj) & set(yh)
sj_only = set(sj) - set(yh)
yh_only = set(yh) - set(sj)

print(f"\n## Per T_H match (key = T_H, LST gauge, SUGRA gauge, Delta)")
print(f"  T_H | SJ keys | YH keys | common | SJ-only | YH-only")
for t in T_RANGE:
    s = {k for k in sj if k[0] == t}
    y = {k for k in yh if k[0] == t}
    print(f"  {t:>3} | {len(s):>7} | {len(y):>7} | {len(s & y):>6} | {len(s - y):>7} | {len(y - s):>7}")
print(f"\nTotal: common={len(common)}, SJ-only={len(sj_only)}, YH-only={len(yh_only)}")

if yh_only:
    print("\nYH-only examples (first 8):")
    for k in sorted(yh_only)[:8]:
        T_H, lst, sg, d = k
        print(f"  T_H={T_H} Delta={d:>5}  LST='{lst[:25]}'  SUGRA='{sg[:50]}'")
if sj_only:
    print("\nSJ-only examples (first 8):")
    for k in sorted(sj_only)[:8]:
        T_H, lst, sg, d = k
        print(f"  T_H={T_H} Delta={d:>5}  LST='{lst[:25]}'  SUGRA='{sg[:50]}'")
