#!/usr/bin/env python3
"""Find 8 YH-only entries (after no-223 + LDLT-fix pipeline)."""
import csv, glob, os, re, sys
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

TARGETS = [
    (6, "so8 + 2 su3", "2 g2 + so7 + 2 so8 + 4 su2", 118),
    (6, "so8 + 2 su3", "2 so7 + so8 + 4 su2 + 3 su3", 129),
    (7, "(none)", "f4 + 2 so8", 124),
    (7, "2 g2 + 2 su2", "f4 + g2 + so7 + 3 su2 + 3 su3", 136),
    (7, "2 so8 + su3", "2 g2 + 3 so8 + 2 su2 + 2 su3", 114),
    (7, "f4 + 2 su3", "2 f4 + g2 + so8 + su2 + 3 su3", 67),
    (7, "f4 + 2 su3", "e6 + f4 + 3 g2 + 3 su2 + su3", 67),
    (7, "f4 + so7 + 2 su2", "f4 + g2 + so7 + 3 su2 + 3 su3", 136),
]
target_keys = {(t[0], content_key(t[1]), content_key(t[2]), t[3]): t for t in TARGETS}

found = {k: [] for k in target_keys}
for fn in sorted(glob.glob(os.path.join(BASE, "data/T*/Ext*.tsv"))):
    rel = os.path.relpath(fn, BASE)
    with open(fn) as f:
        r = csv.DictReader(f, delimiter="\t", quotechar='"')
        for ridx, row in enumerate(r, start=2):
            try:
                T_H = int(row["TH+1"]) - 1
                key = (T_H, content_key(row["LST gauge algebra"]),
                            content_key(row["full gauge algebra"]),
                            int(row["Delta"]))
                if key in target_keys:
                    found[key].append((rel, ridx, row))
            except Exception:
                pass

for i, (k, infos) in enumerate(found.items(), 1):
    t_H, _, _, delta = k
    orig = target_keys[k]
    print(f"=== #{i}  T_H={t_H}, Delta={delta} ===")
    print(f"    LST   : '{orig[1]}'")
    print(f"    SUGRA : '{orig[2]}'")
    if not infos:
        print("    (NOT FOUND in Hamada data!)\n")
        continue
    for rel, ridx, row in infos:
        print(f"    {rel}:{ridx}  Ext={row.get('index','?')}  Tmin={row.get('Tmin','?')}")
        print(f"      LST_raw   : '{row['LST gauge algebra']}'")
        print(f"      SUGRA_raw : '{row['full gauge algebra']}'")
    print()
