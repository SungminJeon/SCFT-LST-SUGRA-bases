# Analysis scripts

Python utilities for inspecting catalog output and comparing with reference data.

## Layout

All scripts read the pipeline's catalog directories from the **parent directory**
(the same place that `unified.cat` and the compiled binaries live). They expect:

```
<pipeline_root>/
    cat_1ext_nosu2_<suffix>/           # phase1 output
    cat_nhc_ext_<suffix>/              # NHC ext output
    cat_nhc_ext_nhc_ext_unified_<suffix>_r<N>sext/   # phase2 chain output
    ...
    data/T002/Ext1.tsv, ...            # Hamada's reference TSV files (optional, for comparison)
    analysis/
        <these scripts>
```

The `data/` directory is **not** included in this share — supply your own copy of
the Hamada TSV files to run the comparison scripts, or skip those scripts entirely.

## Requirements

- Python 3.8+
- `numpy` (`pip install numpy`)
- For the kernel-finding step (true c_ext): `numpy.linalg.svd`

## Scripts

### `gen_classification_csv.py`
Shared utility. Provides `parse_cat_file()`, `get_lst_gauge()`, `get_sugra_gauge()`
used by the other scripts. Not run directly.

### `compare_delta_T7.py` — basic Hamada vs SJ comparison
Compares our pipeline (SJ) keys against Hamada (YH) keys, using
`(T_H, LST gauge content, SUGRA gauge content, Δ)` as the matching key,
where `Δ = H_charged − V + 29·T`.

```bash
python3 compare_delta_T7.py
```

Edit the `T_RANGE` tuple near the top to restrict to specific T_H values.

### `compare_filtered.py` — Hamada vs SJ with true c_ext ≤ 16 filter
Same comparison, but additionally applies the fiber-weighted true c_ext ≤ 16
filter to the SJ side. Reports two tables (BEFORE and AFTER the filter) so
you can see which SJ-only entries the filter removes.

```bash
python3 compare_filtered.py
```

### `check_sj_only_cext.py` — c_ext distribution of SJ-only entries
For each SJ-only key (= in our catalog, not in Hamada's), compute the
fiber-weighted true c_ext and report the count of entries with `c_ext ≤ 16` vs
`> 16`.

```bash
python3 check_sj_only_cext.py
```

### `dump_sj_only_179.py` — IF dump for SJ-only (c_ext ≤ 16)
Same SJ-only entries that pass the c_ext filter, but writes a human-readable
file listing each entry's intersection form, ext indices, and c_ext breakdown.
Useful for inspecting "what we have that Hamada doesn't" entry by entry.

```bash
python3 dump_sj_only_179.py > sj_only_IFs.txt
```

### `find_yh_only_8.py` — locate YH-only entries in Hamada's TSV
For a list of `(T_H, LST, SUGRA, Δ)` tuples, find the matching rows in
`data/T*/Ext*.tsv` and report file/row/Ext index. The current list contains
the 8 YH-only entries observed in our last comparison; edit the `TARGETS`
list in the script for other targets.

```bash
python3 find_yh_only_8.py
```

## Methodology notes

- The comparison key is **per-curve gauge content**, not the raw `.cat` IF.
  Cospectral non-isomorphic IFs and same-IF-different-labeling pairs may
  collapse to one entry under this key — that is intentional.
- We always exclude SJ entries whose externals include `hat1`, `su8`, `so16`,
  or `nhc_2_4` (these are not part of the standard NHC ext enumeration we
  compare).
- True c_ext uses the LST IF's positive integer kernel X (from SVD null space)
  to weight ext curve intersections: `k_eff = Σ IF[ext, i] · X[i]`, then
  `c = k_eff · dim(g) / (k_eff + h(g))`. Returns `-1.0` when no valid X exists
  (treat that case as "filter passes by default").
