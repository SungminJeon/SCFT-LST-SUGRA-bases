# SUGRA Base Generator Pipeline

6D supergravity block classification via incremental external curve attachment.

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

The pipeline produces catalogs of 6D non-Higgsable Gravity Blocks at successive rounds. A typical run for `T = 0..10` looks like:

```bash
# Phase 1: single external attach (skip su2 since NHC covers it)
./gen_sugra_phase1 unified.cat 10 0 --use-lst-T --no-su2 --save-nonsugra \
    --out cat_1ext_nosu2_t010

# NHC ext: attach a full NHC chain as the external block
./gen_sugra_nhc_ext unified.cat 10 0 --use-lst-T --save-nonsugra

# Merge Phase 1 output + NHC ext output for Phase 2 input
mkdir -p cat_phase1_nhc_merged cat_phase1_nhc_merged_nonsugra
cp cat_1ext_nosu2_t010/*.cat            cat_phase1_nhc_merged/
cp cat_nhc_ext/*.cat                     cat_phase1_nhc_merged/
cp cat_1ext_nosu2_t010_nonsugra/*.cat    cat_phase1_nhc_merged_nonsugra/ 2>/dev/null
cp cat_nhc_ext_nonsugra/*.cat            cat_phase1_nhc_merged_nonsugra/ 2>/dev/null

# Phase 2 chain: keep stacking externals until the round produces no new bases
./run_chain.sh ./gen_sugra_nhc_ext_phase2 cat_phase1_nhc_merged nhc_ext_unified 9 \
    "--use-lst-T --no-su2" cat_phase1_nhc_merged_nonsugra
```

Outputs are written to `cat_nhc_ext_nhc_ext_unified_r<N>{s,ns}ext/` for each round `N`.

To turn `.cat` files into a LaTeX report:

```bash
./cat2tex
```

## Pipeline driver scripts (`run_pipeline*.sh`)

Instead of running each binary by hand, the `run_pipeline*.sh` wrappers run the
**whole** chain (build → phase1 → nhc_ext → merge → phase2 chain → finalize) for a
tensor range `T_MIN..T_MAX` and leave the deduped catalog in `<OUTDIR>/final_byT/`.
They abort on any step error and call `finalize_sugra.sh --archive` at the end.

| Script | `su2` (single `-2`→`-1`) | `2mix` | Output dir |
|--------|--------------------------|--------|------------|
| `run_pipeline.sh T_MIN T_MAX`            | **off** (`--no-su2`) | on  | `cat_*_t<min><max>` |
| `run_pipeline_su2.sh T_MIN T_MAX`        | on                   | on  | `T_<min>_<max>_su2/` |
| `run_pipeline_su2_no2mix.sh T_MIN T_MAX` | on                   | **off** (`--no-2mix`) | `T_<min>_<max>_su2_no2mix/` |

```bash
./run_pipeline_su2.sh 0 10      # T = 0..10  →  T_0_10_su2/final_byT/
./run_pipeline_su2.sh 11 20     # T = 11..20 →  T_11_20_su2/final_byT/
```

After a run, **`<OUTDIR>/final_byT/` is the answer** — everything else is intermediate.

## Targeting specific externals

phase1 / nhc_ext enumerate many external *specs*. To restrict a run to a few:

| Flag | Effect |
|------|--------|
| `--only-tags <t1,t2,…>` | (phase1) keep only these spec tags. |
| `--only-223`            | (nhc_ext) only the `(-2)-(-2)-(-3)` (`nhc_2_2_3`) cluster. |
| `--no-su2` / `--no-2mix` / `--no-223` | drop that family. |
| `--mix-k N`             | set the mixed (-1)-leg intersection-number ceiling (`mixed_int_max`, default **5**; `2mix` uses its own default **3**). |

**Gotcha — `so7`, `so7mix`, `2`, `2mix` live under the `su3` spec.** They are
`-2`/`-3` attachments built from the su3 target enumeration and are emitted
*inline* there, so to generate them you must request **`su3`** in `--only-tags`
(not `so7`/`2`). That also produces `su3` / `su3n2mix` seeds — drop those
afterward if you only want the so7/2 family. (`2` = a gauge-less `-2` glued onto
an isolated `nhc_2_3`'s `-2`, forming a `nhc_2_2_3`; `2mix` = the same with extra
`(-1)` legs.)

## External attachment specs & multiplicities

Externals attach via several modes. **Single-curve** externals (`attach_single`,
one target curve of the listed self-intersection):

| tag | ext curve | target | int number |
|-----|-----------|--------|------------|
| `su2` | −2 | −1 | 1–2 |
| `su3` | −3 | −1 | 1–2 |
| `su2n3` | −2 | −3 | 1 |
| `su3n2` | −3 | −2 | 1 |
| `su8` | −2 | −4 | 2 |
| `so16n2` | −4 | −2 | 2 |
| `so8` | −4 | −1 | 1–2 |
| `f4` | −5 | −1 | 1 |
| `e6` | −6 | −1 | 1 |
| `e7p` | −7 | −1 | 1 |
| `e7` | −8 | −1 | 1 |
| `e8` | −12 | −1 | 1 |
| `hat1m1` | −1 | −1 | 1 |
| `hat1m2` | −1 | −2 | 1 |

**Multi-target / cluster** modes:

| mode | external | # targets | int per leg |
|------|----------|-----------|-------------|
| `attach_v6_multi` | `su2` / `su3` (−2/−3) | 2 … all `(−1)` | **1–5** |
| `attach_v6_multi` | `so8`…`e8` | 2 … (cc-limited, ≤9) | cc-limited (≤9) |
| `attach_mixed_generic` | mixed (base −2/−3 + `(−1)` legs) | 1 … many `(−1)` | 1–5 |
| `attach_nhc_cluster` | `nhc_2_3` / `nhc_2_3_2` / `nhc_2_4` | each cluster curve → 1–3 `(−1)` (`nhc_2_3` also a −2) | 1–2 |
| `attach_nhc_cluster` | `nhc_2_2_3` | 1 target / curve | 1 |

So `su2`/`su3` reach intersection number **5** via the multi-target / mixed modes
(each `(−1)` leg `1..mixed_int_max = 5`). The single-target spec is enumerated only
at int **1–2**: a higher multiplicity on a *single* `(−1)` is reached only when the
gauge is shared across `≥2` `(−1)` curves (the `v6_multi` / mixed path).

The single-target cap of **2** is a deliberate **chain-load bound**. Measured at
`T=10` (no-canonical, no-2mix): raising the single `su2`/`su3` cap to 5 adds only
**+2.9 %** distinct blocks (8,191,978 → 8,430,707) for **+95 %** wall time (≈1.95×)
— i.e. ~3 % completeness traded for a large speed-up, and the gap widens at higher
`T`. So the cap stays at 2. *(Re-measured 2026-06-30 on current code, chain MAX=12;
an earlier MAX=9 run gave +2.9 % / +79 %.)*

The **chain depth itself** is capped at **9 rounds** (`run_chain.sh` MAX=9 in
`run_pipeline*.sh`). A full `T_H<=10` re-run to convergence (chain MAX=13,
`--save-nonsugra`, clean-final lineage `--no-2mix --no-canonical`) shows the chain
produces **399 distinct blocks past round 9** (rounds 10/11/12 = 365/33/1,
converging at round 13), all at `T_H<=9`. **Every one is invalid:** they are deep
`su2`/`su3`-stacking on the `(−1)` curves of a minimal (mostly `T_H=1`) base, filled
to the c.c. budget (8 per `(−1)`; `c(su2)=1`, `c(su3)=2`). A genuine SUGRA base
attaching that many `su2`/`su3` needs a matching number of `su2n3`/`su2n3mix`
partners, which attach at most 6 — here `su2>=7` (needs >6 partners; 378 blocks) or
`su2<=6` with no `su2n3` partner (21 blocks), so **0/399 survive**. **So the 9-round
cap loses no valid blocks** — a justified practical bound, not a truncation of real
data. (Evidence dataset: `round10plus_cap9_excluded/`. `run_chain.sh` auto-stops on
an empty round, so MAX=9 is only a safety ceiling.)

## ⚠️ Intersection-number blow-up (light gauge externals)

The **light** gauge externals — `su2` (a `-2`), `su3` (a `-3`), the
`nhc_2_2_3` / `223` cluster, and especially their *mixed* variants
(`2mix`, `su3n2mix`, `so7mix`, …) — can produce an **enormous** number of bases,
because:

1. they attach to the LST's `(-1)` curves, which are **abundant** (and grow with `T`), and
2. each `(-1)` leg of a *mixed* external enumerates its **intersection number
   `1..mixed_int_max`** (`enumerate_int_nums`).

So the base count grows roughly as

```
(choices of which (-1) legs)  ×  (mixed_int_max ^ number of legs)
```

i.e. it **explodes** as you raise `mixed_int_max` (`--mix-k`) or the leg count.
Built-in guards: the number of `(-1)` legs is hard-capped (`max_m1` = 4 for most
mixes, 2 for `2mix`) and `2mix` defaults to a lower ceiling (3) — these caps are
intentional anti-blow-up limits; **do not raise them casually**.

By contrast the **heavy** externals (`e6`, `e7`, `e7p`, `e8`, `f4`) are single
curves attached with intersection number **1** (`so8` at **1–2**), so they do
*not* blow up.

**Guidance:** for production keep the default `mixed_int_max`; raise `--mix-k`
only for a **narrow** `T` window and expect super-linear growth. If a run is
unexpectedly slow or huge, the culprit is almost always the light-gauge mixed
enumeration (`2mix`, `so7mix`, …) — scope it with `--no-2mix` / `--no-su2` /
`--only-tags`, not the heavy externals.

## Final deduplicated catalog (`finalize_sugra.sh`)

A pipeline run leaves a **zoo** of directories. Only some hold final SUGRA bases:

| Directory | What it is | In the final catalog? |
|-----------|-----------|-----------------------|
| `cat_1ext_*` | Phase 1 seeds | **yes** |
| `cat_nhc_ext_*` (seed + every round `…_rN{s,ns}ext`) | NHC seed + chain rounds | **yes** |
| `*_nonsugra` | intermediates feeding the next round's non-SUGRA lane | no |
| `cat_phase1_nhc_merged_*` | redundant (= phase1 + nhc seed; Phase 2 input) | no |

**Trap:** `…_rN nsext` (no `_nonsugra` suffix) is a **SUGRA** directory — it is the
non-SUGRA→SUGRA *recovery* lane, keep it. Only literal `*_nonsugra` dirs are dropped.
The SUGRA lane (`sext`) and the recovery lane (`nsext`) reach some of the **same**
bases, so they must be deduplicated **together** per `(T, combo)`; plain
concatenation leaves cross-stream duplicates.

`finalize_sugra.sh` does exactly this and produces the single canonical catalog:

```bash
./finalize_sugra.sh [--archive] <OUT_DIR> [SRC_GLOB ...]
#   -> <OUT_DIR>/final_byT/T<n>/<combo>.cat     (THE deduped SUGRA catalog)
```

- Source dirs are auto-detected (flat `cat_1ext_*`/`cat_nhc_ext_*`, or per-T
  `byT/` + `chain/`); pass `SRC_GLOB`s explicitly for non-standard layouts
  (e.g. `chain_su2/`) so unrelated experiment variants are not merged.
- `--archive` moves the leftover zoo into `<OUT_DIR>/_intermediate/`. Use it only
  for throwaway flat-pipeline outputs; **omit** it on per-T family masters whose
  `byT/`+`chain/` are the durable, `.done`-resumable slices.
- `T` for each entry comes from its `PHYSICS` `anomaly.T` field (`consolidate_dedup.py
  --byphys`), so grouping is correct regardless of directory layout.

The `run_pipeline*.sh` scripts call this automatically as their final step, so
after a run **`<OUT_DIR>/final_byT/` is the answer** — everything else is
intermediate. `finalize_sugra.sh` is idempotent and safe to re-run by hand.

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

## Analysis — Python toolkit (`sugra_toolkit/`) — **recommended**

**This is the current, maintained way to read and analyse the catalog.** No build
step: `pip install -r sugra_toolkit/requirements.txt` (just numpy + matplotlib), then

```python
import sugra
cat = sugra.Catalog("path/to/final_byT")      # e.g. "sugra blocks clean final"
cat.tensors(); cat.combos(8)
b = cat.load(8, "e6+su3")[0]
b.IF, b.physics, b.t_min()                      # numpy IF, {Hc,V,Hn,det,sig}, T_min
cat.find(t_h=8, externals=["e6", "su3"])        # order-free multiset query
```

Catalog-wide statistics & figures (also in the toolkit):
- `sugra_toolkit/stats.py`, `gauge_analysis.py` — scan a catalog → summary CSVs.
- `sugra_toolkit/make_plots.py` — publication figures (`.png` + `.pdf`).
- **Precomputed outputs:** `gauge_out/` and `stats_final/` (CSVs + plots), and a
  curated headline summary in `sugra_toolkit/catalog_summary/`.

See **`sugra_toolkit/README.md`** for the full API, the statistics snapshot, and
the data layout.

## Mathematica interface (legacy — unmaintained)

> **Note:** `SUGRACatalog.m` is an **older** interface that is **no longer
> maintained** and does **not** reflect the current catalog format/conventions
> (su2/su2n3mix tagging, 1-tag-per-NHC-cluster, the `1`/special_m1 external, …).
> Use the Python toolkit above. Kept for reference only.

`SUGRACatalog.m` loads `.cat` files: `LoadCatalog`, `FindEntries`, `ShowEntry`,
`EntryIF`, `TMin`. See `guide.m`, `tutorial.m`.

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
run_pipeline*.sh          One-command pipeline drivers (build → … → finalize)
finalize_sugra.sh         Dedup post-step → final_byT (THE clean catalog)
Makefile                  Build rules
unified.cat               LST base catalog (input)

sugra_toolkit/            Python analysis toolkit (recommended) — loader, stats,
                          gauge analysis, plots, catalog_summary/  [see its README]
gauge_out/                Precomputed gauge-analysis CSVs + figures
stats_final/              Precomputed stats CSVs (by_combo, by_LSTGauge, …)
Frozen blocks/            Frozen gravity blocks (su8/su16/so16) — explicit
                          .cat data + analysis CSVs/plot/tex (see sugra_toolkit)

SUGRACatalog.m            Mathematica analysis package (legacy, unmaintained)
guide.m, tutorial.m       Mathematica usage examples (legacy)
```
