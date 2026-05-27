# SUGRA Base Generator Pipeline

6D supergravity base classification via incremental external curve attachment.

## Requirements

- C++17 compiler (g++ 9+)
- [Eigen3](https://eigen.tuxfamily.org/) (`brew install eigen` on macOS, `apt install libeigen3-dev` on Debian/Ubuntu)
- `unified.cat` (LST catalog, included)

## Build

```bash
make
```

If Eigen3 is in a non-standard location:

```bash
make EIGEN=/path/to/eigen3
```

## Quick Run

Use the wrapper script to run the full pipeline for a chosen tensor range:

```bash
./run_pipeline.sh T_MIN T_MAX
```

Examples:

```bash
./run_pipeline.sh 0 10    # T = 0..10  → catalogs in cat_*_t010 / r<N>{s,ns}ext_t010
./run_pipeline.sh 11 20   # T = 11..20 → catalogs in cat_*_t1120 / ...
./run_pipeline.sh 21 30   # T = 21..30 → catalogs in cat_*_t2130 / ...
```

The suffix is `t<T_MIN><T_MAX>` (T_MAX zero-padded to 2 digits). All four phases are executed:
1. Phase 1 — single-curve external attach (`gen_sugra_phase1`)
2. NHC ext — attach a full NHC chain (`gen_sugra_nhc_ext`)
3. Merge — combine Phase 1 + NHC ext outputs
4. NHC ext Phase 2 chain — multi-round external stacking (`run_chain.sh` driving `gen_sugra_nhc_ext_phase2` for up to 9 rounds)

Final outputs live in `cat_nhc_ext_nhc_ext_unified_<suffix>_r<N>{s,ns}ext/`, one directory per round.

To turn `.cat` files into a LaTeX report:

```bash
./cat2tex
```

## Modes / Flags

| Flag | Effect |
|------|--------|
| `--use-lst-T` | Canonical mode. Use the LST's `T` (number of curves in the LST base) instead of the SUGRA signature for the anomaly's tensor count. |
| `--save-nonsugra` | Also write entries that fail anomaly/uniqueness checks to a parallel `*_nonsugra/` directory, so they can feed downstream rounds (recoverable later). |
| `--no-su2` | Exclude the external that attaches an `su(2)` (i.e. a single `(-2)` curve) only to a `(-1)` LST curve. This single-curve external blows up the combinatorics, so it is skipped by default in production runs. |
| `--no-223` | Skip the `(-2)-(-2)-(-3)` NHC cluster (`nhc_2_2_3`) external attach. Useful for comparing with references that do not enumerate this cluster; also the slowest cluster, so this mode is much faster. |
| `--det-sq` | Reject IFs whose `|det|` is not a perfect square (extra geometric constraint). |
| `--out <dir>` | (Phase 1 only) Override the output directory name. |

## Per-Phase Binaries

| Binary | Stage |
|--------|-------|
| `gen_sugra_phase1` | Attach one single-curve external to each LST base. |
| `gen_sugra_nhc_ext` | Attach a whole NHC chain (`-2-3`, `-2-3-2`, `-2-2-3`, `-2-4`) as the external block. |
| `gen_sugra_phase2` | Stack a second single-curve external on Phase 1 results. |
| `gen_sugra_nhc_ext_phase2` | Stack a second external (single curve or NHC) on NHC ext results. Chains naturally via `run_chain.sh`. |
| `gen_sugra_phase3` | Glue two Phase 1 entries that share an external curve. |
| `cat2tex` | Convert `.cat` files into LaTeX tables. |

## Mathematica Interface

`SUGRACatalog.m` loads any of the produced `.cat` files for interactive analysis. Key entry points: `LoadCatalog`, `FindEntries`, `ShowEntry`, `EntryIF`, `TMin`, `Clusters`, `FiberClass`, `GlueEntries`.

- **`tutorial_en.m`** — English tutorial covering load, inspect, filter, blowdown, glue, manual IF input.
- `guide.m`, `tutorial.m` — Korean originals (kept for reference).

In Mathematica:
```mathematica
SetDirectory[NotebookDirectory[]];
<< "SUGRACatalog.m";
entries = LoadCatalog["cat_1ext_nosu2_t010"];
CatalogSummary[entries]
ShowEntry[entries[[1]]]
```
For all exported names: `? SUGRACatalog`*`.

## Analysis Scripts

The `analysis/` directory contains Python utilities for inspecting catalog
output and comparing with reference (Hamada) data:

```
analysis/
    README.md                       Notes on each script
    gen_classification_csv.py       Catalog parser (used by other scripts)
    compare_delta_T7.py             Basic SJ vs Hamada comparison
    compare_filtered.py             Same, with true c_ext ≤ 16 filter
    check_sj_only_cext.py           c_ext distribution of SJ-only entries
    dump_sj_only_179.py             Dump intersection forms of SJ-only entries
    find_yh_only_8.py               Locate YH-only entries in Hamada's TSV
```

These scripts expect Hamada's `data/T*/Ext*.tsv` files in the parent
directory (not bundled here). See `analysis/README.md` for details.

## Catalog Format

Each entry in a `.cat` file follows:

```
ENTRY <id>
SPEC <spec_id> <tag> <ext_si> <target_si> <int_num> <hat1>
BASE <catalog_id> <type> <base_T> <IF_size>
ATTACH <target_idx> <ext_curve_idx>
PHYSICS <T> <H_charged> <V> <H_neutral> ...
IF <n>
 ... n × n intersection-form matrix ...
REMAINING <count>
 ... per-curve curve_idx, self-int, available cc budget ...
END
```

The `REMAINING` block carries per-curve cc budget for downstream rounds.

## Repository Layout

```
sugra_generator.h         Shared core: NHC table, gauges, anomaly, dedup helpers
theory_sugra.h            LST loader + IF reconstruction
Tensor.{h,C}              Intersection-form matrix utilities
Topology_enhanced.*       Topology graph infrastructure
TopoLineCompact_*         Compact topology line representation
gen_sugra_phase1.cpp      Phase 1 binary
gen_sugra_phase2.cpp      Phase 2 binary
gen_sugra_phase3.cpp      Phase 3 binary
gen_sugra_nhc_ext.cpp     NHC ext binary
gen_sugra_nhc_ext_phase2.cpp
                          NHC ext Phase 2 binary
cat2tex.cpp               LaTeX exporter
run_chain.sh              Multi-round chain driver
Makefile                  Build rules
unified.cat               LST base catalog (input)
SUGRACatalog.m            Mathematica analysis package
guide.m, tutorial.m       Mathematica usage examples
```
