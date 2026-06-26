#!/usr/bin/env python3
"""
Consolidate chain output into a per-(T,combo) catalog WITH real dedup.

The pipeline dedups only WITHIN each phase2 invocation; the SUGRA-chain
(rNsext) and non-SUGRA-recovery-chain (rNnsext) are separate invocations,
so the same base reached via both paths appears in both. Plain concatenation
(the old consolidate*.sh) leaves those cross-stream duplicates in final_byT.

This script reproduces the pipeline's dedup identity
    (combo, catalog_id, dedup_key_v2)
where dedup_key_v2 = quantized eigenvalues + sorted vertex descriptors
(self_int, degree, sorted neighbour self_ints) + sorted edge multiset.
Vertex descriptors and edges are exact permutation invariants; the eigenvalue
part is coarse-quantized (1e6) so floating noise never splits a true duplicate.

Usage:
    consolidate_dedup.py OUT_DIR  SEED_GLOB  CHAIN_GLOB...
      OUT_DIR     : output root; writes OUT_DIR/T<t>/<combo>.cat
      SEED_GLOB   : glob for 1-ext seed dirs, e.g. 'cat_so7hat2_no2mix/byT/T*'
      CHAIN_GLOB  : glob(s) for SUGRA chain output dirs holding *sext/ subdirs,
                    e.g. 'cat_so7hat2_no2mix/chain/T*'
    Per-T grouping is taken from the 'T<t>' path component of each source.
    Pass '--flat OUT_DIR  GLOB...' to merge a flat (non-per-T) layout instead.
    Pass '--byphys OUT_DIR GLOB...' to group per-T by each entry's PHYSICS
    anomaly.T (layout-agnostic; works for both flat and per-T source trees).
"""
import sys, os, glob, re
import numpy as np
from collections import defaultdict

def dedup_key(rows):
    n = len(rows)
    A = np.array(rows, dtype=float)
    ev = np.linalg.eigvalsh(A)
    evq = tuple(0 if abs(v) < 1e-9 else int(round(v * 1e6)) for v in sorted(ev))
    vd = tuple(sorted(
        (rows[i][i],
         sum(1 for j in range(n) if j != i and rows[i][j] != 0),
         tuple(sorted(rows[j][j] for j in range(n) if j != i and rows[i][j] != 0)))
        for i in range(n)))
    edges = tuple(sorted(
        (min(rows[i][i], rows[j][j]), max(rows[i][i], rows[j][j]), rows[i][j])
        for i in range(n) for j in range(i + 1, n) if rows[i][j] != 0))
    return (evq, vd, edges)

def iter_entries(fn):
    """Yield (cat_id, key, full_block_text, phys_T) per ENTRY...END block, skip '#' headers.
    phys_T = anomaly.T (PHYSICS field 2), used for --byphys per-T grouping."""
    block = []; rows = []; rif = 0; cid = None; inside = False; phys_T = '?'
    for line in open(fn):
        if line.startswith('#'):
            continue
        if line.startswith('ENTRY'):
            block = [line]; rows = []; rif = 0; cid = None; inside = True; phys_T = '?'
            continue
        if not inside:
            continue
        block.append(line)
        s = line.split()
        if line.startswith('BASE'):
            cid = s[1] if len(s) > 1 else '?'
        elif line.startswith('PHYSICS'):
            phys_T = s[1] if len(s) > 1 else '?'
        elif line.startswith('IF'):
            rif = int(s[1]); rows = []
        elif line.startswith('END'):
            if rows:
                yield (cid, dedup_key(rows), ''.join(block), phys_T)
            inside = False
        elif rif > 0 and len(s) == rif and all(x.lstrip('-').isdigit() for x in s):
            rows.append([int(x) for x in s])

def combo_of(fn):
    return os.path.basename(fn)[:-4] if fn.endswith('.cat') else os.path.basename(fn)

def main():
    args = sys.argv[1:]
    flat = False; byphys = False
    while args and args[0] in ('--flat', '--byphys'):
        if args[0] == '--flat':   flat = True
        if args[0] == '--byphys': byphys = True
        args = args[1:]
    out = args[0]
    src_globs = args[1:]
    os.makedirs(out, exist_ok=True)

    # group source .cat files -> (T, combo); T from /T<n>/ path component.
    #   --flat   : ignore path-T, single bucket per combo (T='')
    #   --byphys : ignore path-T at bucket time; assign per-entry T from PHYSICS
    #              anomaly.T at read time (layout-agnostic; matches our T convention)
    # buckets[(T,combo)] = list of file paths  (T='' when flat/byphys)
    buckets = defaultdict(list)
    no_path_T = flat or byphys
    def add_dir(d):
        if d.endswith('_nonsugra'):
            return
        m = re.search(r'/(T\d+)(?:/|$)', d + '/')
        T = m.group(1) if (m and not no_path_T) else ''
        # SUGRA chain dirs: descend into *sext subdirs (skip *_nonsugra)
        sext_dirs = [s for s in glob.glob(os.path.join(d, '*sext')) if not s.endswith('_nonsugra')]
        search = sext_dirs if sext_dirs else [d]
        for sd in search:
            if sd.endswith('_nonsugra'):
                continue
            for cf in glob.glob(os.path.join(sd, '*.cat')):
                buckets[(T, combo_of(cf))].append(cf)
    for g in src_globs:
        for d in glob.glob(g):
            if os.path.isdir(d):
                add_dir(d)
            elif d.endswith('.cat'):
                m = re.search(r'/(T\d+)/', d)
                T = m.group(1) if (m and not no_path_T) else ''
                buckets[(T, combo_of(d))].append(d)

    tot_in = tot_out = 0
    per_T = defaultdict(lambda: [0, 0])
    # In --byphys mode a single (path) bucket can split across several phys-T;
    # dedup is per (final-T, combo), so keep per-T seen-sets / kept-lists.
    for (T0, combo), files in sorted(buckets.items()):
        seen = defaultdict(set); kept = defaultdict(list)
        for fn in files:
            for cid, key, text, phys_T in iter_entries(fn):
                T = ('T' + phys_T) if byphys else T0
                tot_in += 1; per_T[T][0] += 1
                k = (cid, key)
                if k in seen[T]:
                    continue
                seen[T].add(k); kept[T].append(text)
        for T, kept_T in kept.items():
            tot_out += len(kept_T); per_T[T][1] += len(kept_T)
            odir = os.path.join(out, T) if T else out
            os.makedirs(odir, exist_ok=True)
            with open(os.path.join(odir, combo + '.cat'), 'w') as f:
                f.write(''.join(kept_T))

    for T in sorted(per_T):
        i, o = per_T[T]
        print(f"  {T or '(flat)':8} in {i:9} -> unique {o:9}  (-{i-o})")
    print(f"TOTAL: {tot_in} -> {tot_out} unique  (removed {tot_in-tot_out} duplicates, "
          f"{100*(tot_in-tot_out)//max(1,tot_in)}%)")

if __name__ == '__main__':
    main()
