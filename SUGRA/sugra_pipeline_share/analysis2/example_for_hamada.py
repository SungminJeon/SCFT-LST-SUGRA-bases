"""example_for_hamada.py — Python port of example_for_hamada.m.

6D SUGRA base selection by T_H, T_min, and external types.

HOW TO RUN
----------
    python3 example_for_hamada.py

The first load takes a few seconds; subsequent queries are instant since the
entries stay in memory.  Edit the CATALOG / T_H_RANGE constants below to point
at another decade or to enlarge the slice.

PHYSICS CHEAT-SHEET
-------------------
    T_H   = tensor count of the underlying Little String Theory (LST) base.
            Stored on each entry as `base_t`.
    T     = sig_neg of the extended base (= T_H + externals).
    T_min = minimum T that admits a SUGRA b0 lift.
            Returned by entry_t_min(e); None when undefined (unimodular IF).
    external = an extra (-n) curve glued to an LST base; algebra labels are
               e6, e7, e8, su2, su3, so8, f4, ...
    IF    = intersection-form matrix of the extended base.


EXTERNAL TAG REFERENCE
----------------------
The strings below are the legal `tag` values you can pass to

    pick_bases(entries, externals=[...])      # multiset, order-independent
    find_by_combo (entries, "tag1+tag2+...")  # exact "+"-joined string
    find_containing_tag(entries, "tag")       # at least one occurrence

  Simple Lie-algebra externals (extSI = value in entry.externals[i]["extSI"]):
    su2   (-2)    su3 (-3)    so8 (-4)    f4 (-5)
    e6    (-6)    e7p (-7)    e7  (-8)    e8 (-12)

  Mixed / special embeddings:
    su2n3, su2n3mix, su3n2, su3n2mix, so7, so7mix, su8, so16n2
    -- trailing "n<k>" = attachment multiplicity; "mix" = mixed embedding.

  Non-Higgsable clusters (glued in as one external unit):
    nhc_2_3      = (-2,-3) cluster        (gauge: su2 + g2)
    nhc_2_3_2    = (-2,-3,-2) cluster     (gauge: su2 + so7 + su2)
    nhc_2_2_3    = (-2,-2,-3) cluster     (gauge: sp1 + g2)

  Hat-1 externals (isHat1 = 1):
    hat1m1, hat1m2

CATALOG STATISTICS (T1-10, 9,566,360 bases total)
-------------------------------------------------
"Base count" = number of distinct base entries matching the row.
Counts are configuration-based, not "contains at least one tag".

  By T_H:
    T_H= 2:   160,773    T_H= 3: 1,064,708    T_H= 4: 1,793,609
    T_H= 5: 1,808,446    T_H= 6: 1,278,650    T_H= 7: 1,218,321
    T_H= 8: 1,033,171    T_H= 9:   673,698    T_H=10:   534,984

  By number of externals:
    1: 135,070  | 2:   764,994 | 3: 2,640,042 | 4: 3,317,646
    5: 1,940,934| 6:   625,118 | 7:   126,050 | 8:    15,518
    9:     975 | 10:       13

  By external configuration (complete multiset): 4,026 distinct combos.
  Top 5: su2+su2+su2 (436,286), su2+su2+su3 (403,455),
         su2+su2+su2+su3 (349,216), su2+su2+su2+su2 (298,056),
         su2+su2+su2n3mix (291,670).  See README.md for the full top-30.

USAGE EXAMPLES
--------------
    # Exact multiset (any ordering, two NHC clusters + one su2)
    pick_bases(entries, externals=["nhc_2_3", "nhc_2_3", "su2"])

    # Exact combo string (ordering matters; mirrors entry.combo)
    find_by_combo(entries, "su2+su2+su3")
    find_by_combo(entries, "nhc_2_3+nhc_2_3+su2+su2")

    # All bases containing at least one e8 external
    find_containing_tag(entries, "e8")
"""
from __future__ import annotations

import os
import time

from sugra_catalog import (
    entry_t_min,
    find_by_combo,
    find_by_lst,
    find_containing_tag,
    group_by_lst,
    load_catalog,
    lst_if,
    lst_key,
    pick_bases,
    show_entry,
    show_if,
    table_form,
)


# ─── Choose the catalog slice ──────────────────────────────────────────
# Catalogs are bundled per decade in the sibling `catalogs/` folder:
#   catalogs_T1-10.zip   ( T_H =  2 .. 10 )
#   catalogs_T11-20.zip  ( T_H = 11 .. 20 )
#   catalogs_T21-30.zip  ( T_H = 21 .. 30 )
#   ...                  ...
#   catalogs_T141-150.zip( T_H = 145 )
#
# Pass either a .zip (no extraction needed) or an extracted directory.
HERE       = os.path.dirname(os.path.abspath(__file__))
CATALOG    = os.path.join(HERE, "catalogs", "catalogs_T1-10.zip")
T_H_RANGE  = (2, 4)   # start small; widen once your query works.


def main():
    # ─── 1. Load entries ───────────────────────────────────────────────
    t0 = time.time()
    entries = load_catalog(CATALOG, t_h_range=T_H_RANGE)
    print(f"loaded {len(entries)} entries in {time.time()-t0:.1f}s")

    # ─── 2. Entry anatomy ──────────────────────────────────────────────
    # Each element is a `sugra_catalog.Entry` dataclass:
    #   id            - catalog-local identifier
    #   combo         - "+"-joined external tag string (e.g. "e6+su3")
    #   externals     - list of dicts: {curveIdx, specId, tag, extSI, targetSI,
    #                                    intNum, isHat1}
    #   catalog_id    - sub-catalog numeric id
    #   catalog_type  - LST family label ("NN", "DM:A(3)", ...)
    #   base_t        - T_H
    #   T             - sig_neg of the extended base
    #   Hc, V, Hn     - hypers (charged), vectors, hypers (neutral)
    #   det           - det of IF
    #   sig_pos / sig_neg / sig_zero  - IF signature
    #   IF            - extended-base intersection form (numpy ndarray).
    #                   Reading entry.IF always returns a fresh COPY, so
    #                   mutating it cannot poison the catalog.
    #   remaining     - residue curve list
    #   source        - originating .cat path (useful for debugging)
    e = entries[0]
    print(show_entry(e))

    # ─── 3. Convenience selectors are methods on Entry ─────────────────
    print("\nfirst entry ext tags:", e.ext_tags, "  sorted:", e.ext_tag_set)

    # ─── 3b. Intersection form & LST recovery ──────────────────────────
    # `entry.IF` is the *extended* base (LST + glued externals).
    # `lst_if(entry)` drops the external rows/cols, giving the underlying
    # LST's own IF matrix.  Both return safe copies.
    print(f"\nextended IF (shape {e.IF.shape}):")
    print(show_if(e))
    print(f"\nLST IF (shape {lst_if(e).shape}):")
    print(lst_if(e))

    # `lst_key(entry)` is the (catalog_type, base_t, catalog_id) tuple
    # that uniquely identifies the underlying LST.  Entries sharing this
    # key sit on the SAME base; only their externals differ.
    print(f"\nLST key of first entry: {lst_key(e)}")

    # Group all entries by their underlying LST:
    groups = group_by_lst(entries)
    print(f"distinct LSTs in this slice: {len(groups)}")
    print("top 5 LSTs by entry count:")
    for key, lst_entries in sorted(groups.items(), key=lambda kv: -len(kv[1]))[:5]:
        print(f"  {key}  -> {len(lst_entries):>7,} entries")

    # All entries sitting on DM:A(3), T_H=2:
    on_a3 = find_by_lst(entries, "DM:A(3)", base_t=2)
    print(f"\nDM:A(3), T_H=2 -> {len(on_a3):,} entries")

    # ─── 4. Pick by T_H ────────────────────────────────────────────────
    # All bases whose underlying LST has T_H = 3.
    sel = [x for x in entries if x.base_t == 3]
    print(f"\nT_H = 3   -> {len(sel)} entries")

    # ─── 5. Pick by T_min ──────────────────────────────────────────────
    # Bases with T_min == 7. entry_t_min returns None for unimodular IFs.
    sel = [x for x in entries if entry_t_min(x) == 7]
    print(f"T_min = 7 -> {len(sel)} entries")

    # ─── 6. Pick by external types ─────────────────────────────────────

    # 6a) Combo-string match — fastest, order-sensitive.
    sel = find_by_combo(entries, "e6+su3")
    print(f"combo == 'e6+su3'      -> {len(sel)} entries")

    # 6b) Multiset match — order-independent.
    target = sorted(["e6", "su3"])
    sel = [x for x in entries if x.ext_tag_set == target]
    print(f"ext_tag_set == {target} -> {len(sel)} entries")

    # 6c) Contains at least one e8 external (others allowed).
    sel = find_containing_tag(entries, "e8")
    print(f"has e8 external        -> {len(sel)} entries")

    # 6d) Exactly two externals, both e6.
    sel = [x for x in entries if x.ext_tag_set == ["e6", "e6"]]
    print(f"ext_tag_set == [e6,e6] -> {len(sel)} entries")

    # ─── 7. Combined filter (T_H + T_min + externals) ──────────────────
    # Note: T_H = 5 is *outside* T_H_RANGE = (2,4) above, so this will be 0.
    # Widen T_H_RANGE to (2, 8) and rerun to see hits.
    sel = pick_bases(entries, t_h=5, t_min_val=8, externals=["e7", "e7", "su3"])
    print(f"\nT_H=5, T_min=8, ext=[e7,e7,su3] -> {len(sel)} entries")
    for s in sel:
        print(show_entry(s))

    # ─── 8. Reusable filter is `pick_bases`. ───────────────────────────
    # Any kwarg you omit acts as "match all" for that dimension.
    pick_bases(entries, t_h=3)                          # T_H = 3
    pick_bases(entries, t_h=3, t_min_val=6)             # T_H = 3, T_min = 6
    pick_bases(entries, t_h=3, t_min_val=6, externals=["e6", "su3"])
    pick_bases(entries, t_min_val=8, externals=["e8"])  # T_min = 8, single e8

    # ─── 9. Tabulate a result list ─────────────────────────────────────
    # Swap the call below for whichever query you want tabulated.
    print("\n" + table_form(
        pick_bases(entries, t_h=3, externals=["e6", "su3"]),
        columns=("id", "catalog_type", "base_t", "t_min", "combo", "T"),
    ))

    # ─── 10. Inspect a single entry ────────────────────────────────────
    if sel:
        print("\n" + show_entry(sel[0]))
        print("\nIF matrix:")
        print(show_if(sel[0]))


if __name__ == "__main__":
    main()
