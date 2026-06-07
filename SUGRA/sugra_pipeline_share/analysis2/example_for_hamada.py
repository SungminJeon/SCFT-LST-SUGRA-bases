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
    find_by_combo (entries, "tag1+tag2+...")   # exact "+"-joined string
    find_containing_tag(entries, "tag")       # at least one occurrence

Counts below are from the **T1-10** group as loaded by `load_catalog`
(canonical sub-catalogs only). Total bases: 9,566,360.
`#bases` = "number of base entries that contain at least one external
with this tag" (a base whose externals repeat a tag is counted once).
The same tag vocabulary is used in every decade.

  (1) Simple Lie-algebra externals. `extSI` is the value stored in
      `entry.externals[i]["extSI"]`.

        tag      extSI    #bases        %
        ----     -----    ---------   -----
        su2       -2      8,920,079   93.24
        su3       -3      4,012,697   41.95
        so8       -4      1,083,753   11.33
        f4        -5        404,213    4.23
        e6        -6        298,127    3.12
        e7p       -7        138,198    1.44   (e7 "prime" / broken e7)
        e7        -8        146,483    1.53
        e8       -12         41,815    0.44

  (2) Mixed / special embeddings.  The trailing "n<k>" denotes the
      attachment multiplicity of the external curve, and "mix" marks a
      configuration that mixes two simple factors.

        tag         #bases       %    meaning
        ----        ---------  -----  -------
        su2n3mix    2,762,458  28.88  su2 ext, n=3, mixed embedding
        su2n3       1,714,072  17.92  su2 ext with n=3 attachment
        su3n2         202,399   2.12  su3 ext with n=2 attachment
        su3n2mix       38,299   0.40  su3 ext, n=2, mixed embedding
        so7               815   0.01  so7 external
        so7mix             63  <0.01  so7 external, mixed embedding
        su8                 9  <0.01  su8 external
        so16n2              3  <0.01  so16 ext with n=2 attachment

  (3) Non-Higgsable clusters (NHC) glued in as a single external unit.
      The number suffix encodes the self-intersection sequence of the
      cluster's component curves.

        tag           cluster        gauge algebra        #bases
        ----          -------        -------------        ---------
        nhc_2_3       (-2,-3)         su2 + g2            1,699,729
        nhc_2_3_2     (-2,-3,-2)      su2 + so7 + su2       396,287
        nhc_2_2_3     (-2,-2,-3)      sp1 + g2              184,442

  (4) Hat-1 externals (isHat1 = 1).  These are (-1)-curves used in the
      hat-1 / unhiggsing construction.

        tag       #bases
        ----      ------
        hat1m1     425
        hat1m2       6

  See README.md for the top-30 combination ranking and the
  per-#externals distribution.

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
    find_containing_tag,
    load_catalog,
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
    # Each element is a `sugra_catalog.Entry` dataclass with these attributes:
    #   id            - catalog-local identifier
    #   combo         - "+"-joined external tag string (e.g. "e6+su3")
    #   externals     - list of dicts: {curveIdx, specId, tag, extSI, targetSI,
    #                                    intNum, isHat1}
    #   catalog_id    - sub-catalog numeric id
    #   catalog_type  - human-readable LST label
    #   base_t        - T_H
    #   T             - sig_neg of the extended base
    #   Hc, V, Hn     - hypers (charged), vectors, hypers (neutral)
    #   det           - det of IF
    #   sig_pos / sig_neg / sig_zero  - IF signature
    #   IF            - numpy ndarray of the intersection form
    #   remaining     - residue curve list
    #   source        - originating .cat path (useful for debugging)
    e = entries[0]
    print("\nfirst entry keys:", list(vars(e).keys()))
    print(show_entry(e))

    # ─── 3. Convenience selectors are methods on Entry ─────────────────
    print("\nfirst entry ext tags:", e.ext_tags, "  sorted:", e.ext_tag_set)

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
