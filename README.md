# Single-cell IDS

Applying the **InterDependence Score** (Radhakrishnan, Jain, Uhler & Lander, PNAS
2025, [PMC12403096](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12403096/)) to
single-cell RNA-seq: discover gene programs bottom-up from gene–gene dependence,
without imposing cell-type labels, then read the programs back against cell types.

**Design principle: IDS itself is never reimplemented.** All dependence maths goes
through the authors' `ids.compute_IDS` unchanged. Everything in `sc_ids/` is
single-cell scaffolding around it — batching, QC, stratification, program ranking,
cell-type enrichment.

## Layout

```
src/
├── ids/                  THE AUTHORS' PACKAGE — do not modify
└── sc_ids/               ours
    ├── pipeline.py       load → metadata join → QC → preprocess → IDS → modules (+ CLI)
    ├── analysis.py       end-to-end run: the above, then ranking, enrichment, figures
    ├── batched.py        gene-axis blocking (paper's d′) — exact, validated to 1e-15
    ├── streaming.py      cell-axis streaming, two-pass (paper's out-of-core setting)
    ├── programs.py       program ranking rule + cell-type enrichment
    ├── celltypes.py      granular subtype → main cell type
    ├── viz.py            PCA → Harmony → UMAP, Leiden, plots
    ├── datasets/         per-dataset particulars (see below)
    └── resources/        cell-type mapping files, one per dataset
docs/                     method decisions and results, read these first
results/<run>/            outputs of a run
tests/                    exactness tests for batching/streaming, logic tests for programs
```

### `datasets/` — where dataset-specific facts live

A `Dataset` describes what is true of one study and nothing else: how barcodes are
keyed, which column is the QC gate, which is the cell type, which columns define a
batch-free stratum, and starting parameters. The generic code takes these as
arguments and knows about no particular study.

**Adding a dataset means adding one file there.** `datasets/skin.py` is a checklist
of what must be answered before a skin dataset can be run; it deliberately
registers nothing until those answers exist.

Note that `defaults` are *starting points that worked on that dataset*, not settings
that transfer. `min_ids` in particular must be re-derived from each dataset's own
IDS histogram.

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e .
pip install scanpy anndata scipy "harmonypy==0.0.10" umap-learn leidenalg igraph networkx scikit-learn

python -m sc_ids.pipeline            # smoke test on built-in synthetic data
python tests/test_batched.py && python tests/test_streaming.py && python tests/test_programs.py
```

`sc.external.pp.harmony_integrate` is broken with harmonypy ≥ 2.0.0 — use
`viz.run_harmony()`, which calls harmonypy directly and handles both conventions.

## Running

```bash
python -m sc_ids.analysis \
  --dataset gse194315 \
  --data data/GSE194315_PBMC-01-07_processed_data_files --prefix "PBMC-01-1." \
  --metadata data/GSE194315_CellMetadata-AS_TotalCiteseq_20220711.tsv \
  --discovery-col Subject --discovery-value PSA26 \
  --out results/PSA26
```

Programs are **discovered** on one stratum and **scored** on every cell in the file,
so cell-type labels never influence discovery and the enrichment table is an
out-of-sample check. Thresholds come from `--dataset` unless overridden, and the run
prints which values it inherited.

Data (GSE194315, ~3 GB) is gitignored and kept locally under `data/`.

## Two representations, on purpose

UMAP and clustering use log1p-normalized data reduced to PCs. IDS uses normalized
counts min-max scaled to [0, 8] with **no log1p** — log1p before that scaling
collapses IDS (a real module dropped from 0.89 to 0.09). Program *scores* use a
third: `sc.tl.score_genes` on log space, because it subtracts a matched control set.
Computed separately from the same counts; do not conflate them.

## Docs

- `docs/agenda1_and_3_results.md` — the cell-type mapping, the program-selection
  rule and why it is fixed in advance, and the cell-type enrichment results
- `docs/agenda2_combining_strata.md` — how to combine stratified runs, and why
  Harmony cannot be used for IDS
- `docs/scope_skin_extension.md` — what extending to a second (skin) dataset covers,
  and the four gaps between two program lists and a shared-program claim
- `docs/architecture.html` — visual explainer of the batching design

## Upstream

- IDS paper repo: https://github.com/aradha/interdependence_scores
- Lab single-cell pipeline (QC/Harmony/UMAP reference):
  https://github.com/CutaneousBioinf/single-cell_pipeline
