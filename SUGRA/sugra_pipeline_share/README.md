# SUGRA Base Generator Pipeline

6D supergravity block classification via incremental external curve attachment.

## Recent Changes (2026-06) — Performance & Memory

Speed and memory work on the high-`T` pipeline. **All output-preserving** — validated
byte-identical against the `T=0..10` baseline (full 9-round chain), so the catalogs are
unchanged, only faster and lower-memory. The `.cat` file format is unchanged, so the
Mathematica / analysis tooling is unaffected.

**Faster**
- Incremental signature via a **cached Schur complement** — reuse the base form's inertia
  instead of a fresh eigensolve per candidate attachment.
- **NHC-gauge-aware cc** pruning in the cluster enumeration (~74M → <1M candidates).
- **Faster `.cat` parsing**, **lazy + parallel** eigensolve/SVD, **structural pruning** of
  the `so7` cluster, and **deferred** full-IF construction (built only after the cheap
  signature filter passes).
- Impact: `nhc_ext` at `T=108` **~40–60 min → ~20–40 s**; a heavy chain slice **~30 s → ~12 s**.

**Lower memory**
- **Sparse intersection-form storage** (store nonzeros, reconstruct on demand) — peak RAM
  ~30× lower; e.g. the `T=108` chain load drops from ~17 GB to ~1 GB.

**Correctness**
- Fixed a high-`T` **determinant overflow** (32-bit `int` → `int64`) that corrupted the
  `|det|` / unimodular filters.

A standalone correctness check for the Schur signature is in `test_schur.cpp` (5000 random
integer forms, 0 inertia mismatches vs the full eigensolve).

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
| `--mix-k N` | Override the multi-target enumeration cap `mixed_int_max` (default `5`). Affects `(-2)`/`(-3)` `→` multiple `(-1)` enumeration and the mixed multi-target attach paths (`su2n3mix`, `su8mix`, ...). Higher `N` means broader enumeration and longer runtime. |
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

`gen_sugra_phase2.cpp` and `gen_sugra_nhc_ext_phase2.cpp` have the **same**
`build_all_specs()` block — change all three files (Phase 1 and both Phase 2 variants)
to keep single-curve externals consistent across rounds.

> Note: `compute_ext_cc()` in the same files maps `(ext_si, target_si, k)` to the
> central charge. If you add a genuinely new `(ext_si, target_si)` combination,
> add a matching case there so the c_ext budget check stays correct.

> Important: the `int_num` (3rd argument) and the `k <= 2` loop bound above only
> govern the **single-target** attach path. The multi-target paths below
> (v6 multi-target and the mixed multi-target attaches) enumerate `int_num`
> independently up to `mixed_int_max` (default 5, override with `--mix-k N`).

### Multi-target attachments (v6 multi & mix variants)

In addition to the single-target attaches described above, `gen_sugra_phase1.cpp`
and `gen_sugra_nhc_ext_phase2.cpp` enumerate multi-target attachments that don't
appear in the `build_all_specs()` table. These produce their own catalog files
distinct from the single-curve specs.

#### v6 multi-target (`(-2)` or `(-3)` ext to multiple `(-1)` LST curves)

For each `bare su2` (`ext_si=-2, target_si=-1`) and `bare su3` (`ext_si=-3, target_si=-1`)
spec, a separate enumeration attaches the external to **multiple `(-1)` LST curves
simultaneously** with independently chosen intersection numbers:

```cpp
// gen_sugra_phase1.cpp Mode 3 (v6 multi-target)
if (spec.ext_si == -2) {
    max_k = config.mixed_int_max;   // default 5 (mixed_int_max)
} else if (spec.ext_si == -3) {
    max_k = config.mixed_int_max;   // default 5
} else if (spec.ext_si == -4) {
    max_k = 2;                      // so8 capped at 2 (matches max_level_for_cc)
} else {
    max_k = 1;                      // f4 and higher
}
enumerate_subsets_v2(m1_curves, 2, max_tgt, [&](targets) {
    enumerate_int_nums(targets.size(), max_k, [&](int_nums) { ... });
});
```

These entries land in the standard `su2.cat` / `su3.cat` / `so8.cat` etc. files,
but with multiple `(-1)` connections each carrying its own `int_num`. The
single-target spec table's `k <= 2` does **not** apply here.

`gen_sugra_nhc_ext_phase2.cpp` has a matching `attach_v6_multi` block with the same
caps. Keep both in sync when changing the limits.

#### Mix variants (separate catalog files)

The mix variants attach a single external curve to a **non-`(-1)` LST target plus
additional `(-1)` LST curves simultaneously**, forming a composite NHC-like
sub-structure. Each variant produces its own `.cat` file under the variant tag:

| tag | physical attach | source |
|---|---|---|
| `su2n3mix` | `(-2)` ext → one `(-3)` target (int=1) + subset of `(-1)` curves (int=1..mixed_int_max) | `gen_sugra_phase1.cpp` Mode 4 |
| `su8mix` | `(-2)` ext → one `(-4)` target (int=2) + subset of `(-1)` curves | `gen_sugra_phase1.cpp` Mode 5 |
| `su3n2mix` | `(-3)` ext → one `(-2)` target (int=1) + subset of `(-1)` curves (int=1) | `try_mixed` lambda |
| `so16n2mix` | `(-4)` ext → one `(-2)` target (int=2) + subset of `(-1)` curves (int=1) | `try_mixed` lambda |
| `so7` | `(-3)` ext → two `(-2)` LST curves (no `(-1)` connection), forms `(-2)-(-3)-(-2) = SO(7)` chain | inline block |
| `so7mix` | `so7` + additional `(-1)` connections | inline block |

The mixed multi-target enumeration cap for the `(-1)` connections is the same
`mixed_int_max` knob (`--mix-k N`). The `(-3)`, `(-4)`, and `(-2)` non-`(-1)`
targets always use the fixed intersection number shown above (`int=1` or `int=2`).

`gen_sugra_nhc_ext_phase2.cpp` has matching attach functions in `attach_mixed_multi`
(calling `attach_mixed_generic` for the four `mix` variants and inline blocks for
`so7` / `so7mix`).

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
