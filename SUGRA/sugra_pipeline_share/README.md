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

The suffix is `t<T_MIN><T_MAX>` (T_MAX zero-padded to 2 digits). All phases run automatically:
1. Phase 1 — single-curve external attach (`gen_sugra_phase1`)
2. NHC ext — attach a full NHC chain (`gen_sugra_nhc_ext`)
3. Merge — combine Phase 1 + NHC ext outputs
4. NHC ext Phase 2 chain — multi-round external stacking (`run_chain.sh` driving `gen_sugra_nhc_ext_phase2` for up to 9 rounds)
5. Collect — all output directories are moved into a single per-range folder `T_<T_MIN>_<T_MAX>/`.

So after `./run_pipeline.sh 0 10` everything lives under `T_0_10/`:
```
T_0_10/
    cat_1ext_nosu2_t010/           cat_1ext_nosu2_t010_nonsugra/
    cat_nhc_ext_t010/              cat_nhc_ext_t010_nonsugra/
    cat_phase1_nhc_merged_t010/    cat_phase1_nhc_merged_t010_nonsugra/
    cat_nhc_ext_nhc_ext_unified_t010_r<N>sext/    (Phase 2 chain, SUGRA)
    cat_nhc_ext_nhc_ext_unified_t010_r<N>nsext/   (Phase 2 chain, non-SUGRA chain)
    ... (+ _nonsugra variants)
```

### Two variants: with / without the su2 external

| Script | `su2` ((-2)→(-1)) external | Output folder |
|--------|----------------------------|---------------|
| `run_pipeline.sh`     | **excluded** (`--no-su2`) | `T_<min>_<max>/`     |
| `run_pipeline_su2.sh` | **included**              | `T_<min>_<max>_su2/` |

The default (`run_pipeline.sh`) skips the single `su2` external because it
explodes the combinatorics. `run_pipeline_su2.sh` keeps it, for a complete
enumeration at the cost of much longer runtime and a larger catalog.

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

The `analysis/` directory contains a Python utility for comparing the
pipeline catalog against a TSV reference data set:

```
analysis/
    README.md                  Notes on the script
    gen_classification_csv.py  Catalog parser (used by the comparison)
    compare_hamada.py          SJ-vs-YH comparison by (T_H, LST gauge, SUGRA gauge, Δ) key
```

Usage:
```bash
python3 analysis/compare_hamada.py         # default T_H range = 1..7
python3 analysis/compare_hamada.py 1 5     # restrict to specific range
```

The script reads TSV reference data from `<pipeline_root>/data/T*/Ext*.tsv`.
Supply your own reference data there to enable the comparison.

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
