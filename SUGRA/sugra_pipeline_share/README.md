# SUGRA Base Generator Pipeline

6D supergravity base classification via incremental external curve attachment.

## Requirements

- C++17 compiler (g++ 9+)
- [Eigen3](https://eigen.tuxfamily.org/) (`brew install eigen` on macOS, `apt install libeigen3-dev` on Debian/Ubuntu)
- `unified.cat` (LST catalog, included)
- **OpenMP** (optional but recommended — see below)

### OpenMP

`gen_sugra_nhc_ext` and `gen_sugra_nhc_ext_phase2` parallelize their main
enumeration loops with OpenMP. Speedup is roughly 3–4× on an 8-thread machine
with output bit-exactly preserved relative to the serial run.

| Platform | Default compiler | OpenMP support |
|----------|------------------|----------------|
| Linux (gcc) | `g++` | built-in `-fopenmp` |
| macOS (Homebrew gcc or LLVM) | `g++-13`/`clang++` | built-in `-fopenmp` |
| macOS (Apple Clang) | `g++` (clang shim) | needs Homebrew `libomp` (`brew install libomp`) |

The Makefile auto-detects which mode applies. If neither is present, it
prints a warning and falls back to a serial build (no functional change,
just slower).

## Build

```bash
make
```

If Eigen3 is in a non-standard location:

```bash
make EIGEN=/path/to/eigen3
```

Control the thread count at run time with `OMP_NUM_THREADS` (e.g.,
`OMP_NUM_THREADS=8 ./run_pipeline.sh 0 10`). Defaults to the number of
hardware threads.

### OpenMP-aware wrapper scripts

For convenience, two extra wrappers take the thread count as an argument:

```bash
./run_pipeline_omp.sh     T_MIN T_MAX [N_THREADS]
./run_pipeline_su2_omp.sh T_MIN T_MAX [N_THREADS]
```

Examples (useful when sharing a server):

```bash
./run_pipeline_omp.sh 0 10 8        # T=0..10, 8 threads
./run_pipeline_omp.sh 21 30 32      # T=21..30, 32 threads
./run_pipeline_su2_omp.sh 0 10 16   # same as run_pipeline_su2.sh, 16 threads
```

If `N_THREADS` is omitted, OpenMP picks the default (all hardware threads).

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

## External Specs and Intersection Numbers

### Single-curve externals (`gen_sugra_phase1.cpp`, `gen_sugra_phase2.cpp`)

The spec table is built in `build_all_specs()` near the top of
`gen_sugra_phase1.cpp`. Each line is:

```cpp
add(ext_si, target_si, int_num, is_hat1, "label", "tag");
```

- `ext_si`       — self-intersection of the external curve
- `target_si`    — self-intersection of the LST curve it attaches to
- `int_num` (k)  — **intersection number** between the external and its target
- `is_hat1`      — special "hat-1" curve flag
- `tag`          — short name used in catalog filenames

Current table:

```cpp
for (int k = 1; k <= 2; k++)
    add(-2, -1, k, false, "su(2) k=...", "su2");   // (-2)->(-1), k = 1, 2
add(-2,  -3, 1, false, "(-2)->(-3)",  "su2n3");    // k = 1
add(-2,  -4, 2, false, "su8->so16",   "su8");      // k = 2
for (int k = 1; k <= 2; k++)
    add(-3, -1, k, false, "su(3) k=...", "su3");   // (-3)->(-1), k = 1, 2
add(-3,  -2, 1, false, "(-3)->(-2)",  "su3n2");    // k = 1
add(-4,  -1, 1, false, "so(8)",       "so8");      // k = 1
add(-4,  -2, 2, false, "(-4)->(-2)",  "so16n2");   // k = 2
add(-5,  -1, 1, false, "f4",          "f4");
add(-6,  -1, 1, false, "e6",          "e6");
add(-7,  -1, 1, false, "e7'",         "e7p");
add(-8,  -1, 1, false, "e7",          "e7");
add(-12, -1, 1, false, "e8",          "e8");
add(-1, -1, 1, true,  "hat1->(-1)",   "hat1m1");
add(-1, -2, 1, true,  "hat1->(-2)",   "hat1m2");
```

**To change an intersection number** for a single-curve external, edit the
`int_num` (3rd) argument of the relevant `add(...)` line — or the loop bound
`k <= 2` to allow more values. For example, to allow `su(3)→(-1)` up to `k=3`:

```cpp
for (int k = 1; k <= 3; k++)                       // was k <= 2
    add(-3, -1, k, false, "su(3) k=...", "su3");
```

`gen_sugra_phase2.cpp` has the **same** `build_all_specs()` block — change both
files (Phase 1 and Phase 2) to keep single-curve externals consistent across rounds.

> Note: `compute_ext_cc()` in the same files maps `(ext_si, target_si, k)` to the
> central charge. If you add a genuinely new `(ext_si, target_si)` combination,
> add a matching case there so the c_ext budget check stays correct.

### NHC cluster externals (`gen_sugra_nhc_ext.cpp`, `gen_sugra_nhc_ext_phase2.cpp`)

NHC chains are defined in `build_nhc_ext_specs()` (top of `gen_sugra_nhc_ext.cpp`):

```cpp
specs.push_back({id++, "nhc_2_3",   "(-2)-(-3) NHC",     {-2,-3},    {{0,1,1}},          8,  17, 3.8});
specs.push_back({id++, "nhc_2_3_2", "(-2)-(-3)-(-2) NHC",{-2,-3,-2}, {{0,1,1},{1,2,1}}, 16,  27, 5.5});
specs.push_back({id++, "nhc_2_2_3", "(-2)-(-2)-(-3) NHC",{-3,-2,-2}, {{0,1,1},{1,2,1}},  8,  17, 3.8});
specs.push_back({id++, "nhc_2_4",   "(-2)-(-4) NHC",     {-2,-4},    {{0,1,2}},         128, 183, 27.6});
```

Fields: `{id, tag, label, self_ints, internal_ints, H, V, cc_total}`.

- `self_ints`       — the chain's curve self-intersections
- `internal_ints`   — edges **inside** the chain as `{a, b, k}` = curve `a`↔curve `b` with intersection `k`
  (e.g. `nhc_2_4`'s `{0,1,2}` means the `-2` and `-4` curves meet with k=2)
- `H`, `V`, `cc_total` — pre-computed hyper/vector counts and total central charge

**To change an internal chain intersection number**, edit the `k` in the
corresponding `internal_ints` triple.

**The intersection number between an NHC curve and the LST target** is decided
during enumeration (not in the spec table). It is set in the enumeration loop of
`gen_sugra_nhc_ext.cpp` (and the matching block in `gen_sugra_nhc_ext_phase2.cpp`):

```cpp
// single-target options
for (int t : curve_targets) {
    int kmax = restrict_223 ? 1 : (base_IF(t,t) == -1 ? 2 : 1);   // <-- max intersection number
    for (int k = 1; k <= kmax; k++) { ... }
}
// multi-target options
enumerate_subsets_v2(curve_m1_targets, 2, max_multi, [&](const auto& targets) {
    enumerate_int_nums((int)targets.size(), 2, [&](const auto& knums) { ... });   // <-- 2 = max k
});
```

- Change `kmax` to allow larger intersection numbers between an NHC curve and a
  single LST target (the `base_IF(t,t) == -1 ? 2 : 1` rule caps it at 2 only when
  the target is a `-1` curve).
- Change the `2` passed to `enumerate_int_nums` to raise the multi-target cap.
- `restrict_223` is a flag that forces the `nhc_2_2_3` cluster to single-target,
  k=1 only (it is the most expensive chain). Remove the `restrict_223 ? 1 :`
  prefix to enumerate it fully.

Increasing any of these sharply increases runtime.

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
