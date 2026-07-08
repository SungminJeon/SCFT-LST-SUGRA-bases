# LST-SUGRA-bases

This repository contains data and code used in the classification program for
six-dimensional LST, and supergravity tensor bases. It is intended as a
companion repository to the paper 

**6d Supergravity Blocks, by Yuta Hamada, Seongmin Jeon, and Hee-Cheol Kim**.
[arXiv:2607.05496](https://arxiv.org/abs/2607.05496)

The public tree is organized around two parts:

- `LST/`: tabulated little-string-theory bases that can appear inside
  supergravity blocks.
- `SUGRA/sugra_pipeline_share/`: code and auxiliary data for generating,
  deduplicating, and inspecting non-Higgsable supergravity blocks.

Large generated catalogs are archived on Zenodo under

**DOI:** [10.5281/zenodo.20979750](https://doi.org/10.5281/zenodo.20979750)

The GitHub repository contains the source tree, and pipeline code. 
Pipeline runs create directories such as `cat_*`,
`T_<min>_<max>_*`, and `final_byT/`; see the SUGRA pipeline README for details.

## Repository layout

```text
SCFT-LST-SUGRA-bases/
├── LST/
│   ├── README.md
│   └── data/
│       └── THNNN.tsv
└── SUGRA/
    └── sugra_pipeline_share/
        ├── README.md
        ├── Makefile
        ├── unified.cat
        ├── gen_sugra_*.cpp
        ├── run_pipeline*.sh
        ├── finalize_sugra.sh
        ├── cat2tex.cpp
        ├── sugra_toolkit/
        ├── stats_final/
        ├── gauge_out/
        └── Frozen blocks/
```

### `LST/`

`LST/data/` contains one TSV file `THNNN.tsv` for each value of the LST tensor
number $T^H$. The included files run from `TH001.tsv` through `TH192.tsv` and
contain 12,762 rows in total.

Each row records one LST base. The columns are documented in `LST/README.md` and
include the catalog type, catalog id, explicit base list, value of
$29T^H + H - V$, non-Higgsable gauge algebra, symmetry, and H-string weight.

### `SUGRA/sugra_pipeline_share/`

This directory contains the SUGRA block-generation pipeline. The pipeline starts
from `unified.cat`, attaches external curves or NHC chains, checks the relevant
constraints, and deduplicates the resulting block catalog.

Important files and directories include:

| Path | Purpose |
|---|---|
| `unified.cat` | Input LST catalog for the SUGRA pipeline. |
| `gen_sugra_phase1.cpp`, `gen_sugra_phase2.cpp`, `gen_sugra_phase3.cpp` | C++ generators for successive external-attachment stages. |
| `gen_sugra_nhc_ext.cpp`, `gen_sugra_nhc_ext_phase2.cpp` | Generators for attaching NHC external blocks. |
| `run_pipeline*.sh` | One-command pipeline drivers. |
| `run_chain.sh` | Multi-round chaining driver. |
| `finalize_sugra.sh` | Deduplication/finalization script; the clean output is `final_byT/`. |
| `cat2tex.cpp` | Converter from `.cat` files to LaTeX tables. |
| `sugra_toolkit/` | Python tools for loading catalogs, making summaries, and plotting. |
| `stats_final/`, `gauge_out/` | Precomputed CSV/plot summaries. |
| `Frozen blocks/` | Explicit data and summaries for frozen gravity blocks. |

For the full list of modes, flags, external tags, and output conventions, see
`SUGRA/sugra_pipeline_share/README.md`.

## Quick start

### Inspect the LST list

Open `LST/README.md` for the TSV column definitions, then read the files in
`LST/data/`. For example, `TH010.tsv` contains the LST bases with $T^H=10$.

### Build the SUGRA pipeline

The C++ pipeline requires a C++17 compiler and Eigen3.

```bash
cd SUGRA/sugra_pipeline_share
make
```

If Eigen3 is installed in a non-standard location, use:

```bash
make EIGEN=/path/to/eigen3
```

### Run a small pipeline job

```bash
./run_pipeline.sh 0 2
```

For production runs and for the `su2` / `2mix` variants, use the documented
`run_pipeline*.sh` scripts in `SUGRA/sugra_pipeline_share/README.md`.

### Read a generated catalog in Python

The maintained catalog interface is the Python toolkit in `sugra_toolkit/`.

```bash
cd SUGRA/sugra_pipeline_share/sugra_toolkit
pip install -r requirements.txt
```

```python
import sugra
cat = sugra.Catalog("../path/to/final_byT")
print(cat.tensors())
```

The toolkit can load `.cat` files, filter by tensor number or external content,
compute summary tables, and reproduce the precomputed statistics/plots. See
`SUGRA/sugra_pipeline_share/sugra_toolkit/README.md` for examples.

## Data formats

### LST TSV files

The files `LST/data/THNNN.tsv` are tab-separated tables. The column `list`
encodes the LST base as a linear or nested list: adjacent entries intersect once,
and nesting records branching. The H-string charge is recorded in the column
`H-String Weight` as coefficients in the ordered basis of curves.

### SUGRA `.cat` files

The SUGRA pipeline writes catalogs in a block format. Each entry contains the
external specification, source LST base, attachment data, anomaly/gauge summary,
intersection-form matrix, and remaining curve budgets. The clean deduplicated
catalog produced by the pipeline is placed under `final_byT/`.

## Citation

If you use this repository, please cite the companion paper. 

```bibtex
@article{Hamada:2026zta,
    author = "Hamada, Yuta and Jeon, Seongmin and Kim, Hee-Cheol",
    title = "{6d Supergravity Blocks}",
    eprint = "2607.05496",
    archivePrefix = "arXiv",
    primaryClass = "hep-th",
    reportNumber = "KEK-TH-2849, RIKEN-iTHEMS-Report-26",
    month = "7",
    year = "2026"
}
```
