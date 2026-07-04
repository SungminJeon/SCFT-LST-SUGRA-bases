# Analysis scripts

Python utilities for inspecting the pipeline catalog and comparing against
Hamada's TSV reference data.

## Layout assumption

Each script lives in `analysis/` next to the parent pipeline root. The pipeline
root must contain the catalog directories that `run_pipeline.sh` produced:

```
<pipeline_root>/
    unified.cat
    cat_1ext_nosu2_<suffix>/
    cat_nhc_ext_<suffix>/
    cat_nhc_ext_nhc_ext_unified_<suffix>_r<N>sext/
    ...
    data/T<NNN>/Ext<N>.tsv               # Hamada's reference TSV (optional)
    analysis/
        compare_hamada.py
        gen_classification_csv.py
        README.md
```

Hamada's `data/T*/Ext*.tsv` files are not bundled — supply them at this path to
enable the comparison.

## Requirements

- Python 3.8+
- Standard library only (no extra packages needed for these scripts)

## Scripts

### `gen_classification_csv.py`
Shared utility. Provides `parse_cat_file()`, `get_lst_gauge()`,
`get_sugra_gauge()` used by `compare_hamada.py`. Not run directly.

### `compare_hamada.py`
Compares the pipeline's SUGRA catalog against Hamada's TSV reference data.
Match key:

    (T_H, LST gauge content, SUGRA gauge content, Delta)
    where Delta = H_charged - V + 29 T

```bash
# Default T_H range = 1..7
python3 compare_hamada.py

# Restrict to specific range
python3 compare_hamada.py 1 5
python3 compare_hamada.py 6 7
```

Output: a per-T_H table of `common`, `SJ-only`, and `YH-only` counts, plus the
first few example keys for any mismatches.

Externals excluded from the comparison: `hat1m1(hat1-> 1)`, `hat1m2(hat1->2)`, `su8(-2 attaching to -2 with 2 intersection number)`,
`su8mix(-2 attaching to -4 with 2 intersection number)`, `so16n2(-4 attaching to -2 with 2 intersection number)`, `so16n2mix`, `nhc_2_4`.

## Methodology notes

- The pipeline already rejects entries with true (fiber-weighted) c_ext > 16
  at the leaf stage, so all `.cat` entries pass that filter by construction.
