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
    ├── units.py          the data as one file per lane / sample; lane-by-lane scoring (matches scanpy)
    ├── discovery.py      program discovery on those files: sampling, gene panel, IDS, programs
    ├── datasets/         per-dataset particulars (see below)
    └── resources/        cell-type mapping files, one per dataset
scripts/                  analysis scripts, one folder per pipeline stage (see "Pipeline" below)
data/                     INPUT: raw GEO downloads. Read-only, gitignored
results/                  OUTPUT: one folder per stage, one subfolder per run
deliverables/             OUTPUT for other people: workbooks sent to Dr. Patrick
docs/                     abstract, reports, method notes, data inventory (read data_inventory.md first)
tests/                    exactness tests for batching/streaming, logic tests for programs
_archive/                 superseded files, kept for one review pass (gitignored; see MANIFEST.md)
```

Paths are defined once in `src/sc_ids/paths.py`. If the data live elsewhere, set `SC_IDS_DATA`
(and `SC_IDS_PBMC_DIR` for the PBMC matrices) instead of editing scripts.

### `datasets/` — where dataset-specific facts live

A `Dataset` describes what is true of one study and nothing else: how barcodes are
keyed, which column is the QC gate, which is the cell type, which columns define a
batch-free stratum, and starting parameters. The generic code takes these as
arguments and knows about no particular study.

**Adding a dataset means adding one file there.** `datasets/psoriasis_skin.py` is the
worked example of a second dataset: its docstring lists the seven things that differed
from the blood data and the code each one needed.

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
  --out results/01_qc_labelling/PSA26
```

Programs are **discovered** on one stratum and **scored** on every cell in the file,
so cell-type labels never influence discovery and the enrichment table is an
out-of-sample check. Thresholds come from `--dataset` unless overridden, and the run
prints which values it inherited.

Data (GSE194315, ~3 GB) is gitignored and kept locally under `data/`.

## Pipeline: raw data to the program catalogue

Scripts are numbered by stage. Each stage reads the one before it and writes only to its own
`results/` folder. Run from the repo root.

| Stage | Run in this order | Writes to | What it does |
|---|---|---|---|
| `01_qc_labelling` | `qc_joint_pipeline` → `cluster_pipeline` → `theta_tuning` → `rank_markers` → `nk_percell`, `nk_deg_export`, `nk_deg_workbook` | `results/01_qc_labelling/` | QC both tissues, align blood and skin on a UMAP, label T/NK and myeloid sub-clusters, check the NK call. Blood: lane PBMC-01-1 only. |
| `02_ids_programs` | `build_units` (blood, skin, then each with `--index`) → `ids_programs_disease` (whole tissue, then `--lineage T+NK`, `Myeloid`, `Keratinocyte`; each once to see the IDS percentiles, then `--reuse --min-ids`) | `results/02_ids_programs/`, data files in `results/00_units/` | Write one QC'd file per blood lane (all 28) and per skin sample, then find IDS programs in each tissue × disease stratum. The catalogue uses the nine October lanes. |
| `03_catalogue` | `program_catalogue` → `disease_contrasts` → `distinct_programs` → `catalogue_workbook` (`keratinocyte_detection` is independent) | `results/03_catalogue/`, `deliverables/` | Permutation test, activity score, merge 194 detections into distinct programs, donor-level disease tests, the workbook. |
| `04_method_checks` | `validate_scoring`; `score_units` → `disease_tests`; `discovery_runs` → `discovery_compare` → `variant_blood_tests`; `permutation_runs` | `results/04_method_checks/` | The same pipeline with one choice changed at a time: all 28 lanes (more donors), subtype-level tests, other random samples / normalization / cap / threshold, 200 shuffles. Summary: `results/04_method_checks/README.md`. |

Which lanes are used and why: `docs/data_inventory.md`. Old paths: `docs/path_map.md`.

### Data layout: one file per lane, everywhere

Stages 02, 03 and 04 all read the same files: one QC'd file per blood lane and per skin sample, written by
`scripts/02_ids_programs/build_units.py` and read through `sc_ids.units`. No script loads a whole tissue, so
peak memory is set by the largest single lane (about 2.4 GB at the highest point, the gene-panel step) and does
not grow when lanes are added. Scoring is done lane by lane and matches scanpy on a merged object to rounding
error (`scripts/04_method_checks/validate_scoring.py`).

The switch from merged objects (November 2026) changed no result: gene panels, all 194 detections, the
catalogue, donor means and every disease test were reproduced from the per-lane files.

Stage 01 is the one exception. Harmony and UMAP need every cell in one object, so that stage still builds
its own merged, embedded object for the lanes it covers.

## Two representations, on purpose

UMAP and clustering use log1p-normalized data reduced to PCs. IDS uses normalized
counts min-max scaled to [0, 8] with **no log1p** — log1p before that scaling
collapses IDS (a real module dropped from 0.89 to 0.09). Program *scores* use a
third: `sc.tl.score_genes` on log space, because it subtracts a matched control set.
Computed separately from the same counts; do not conflate them.

## Docs

- `docs/data_inventory.md` — the PBMC lanes in use (9 of 28) and what each is for, the unused lanes, skin donors
- `docs/abstract/abstract_notes.md` — methods, limitations and numbers behind the KSID abstract
- `docs/notes/agenda1_and_3_results.md` — the cell-type mapping, the program-selection
  rule and why it is fixed in advance, and the cell-type enrichment results
- `docs/notes/agenda2_combining_strata.md` — how to combine stratified runs, and why
  Harmony cannot be used for IDS
- `docs/notes/scope_skin_extension.md` — what extending to a second (skin) dataset covers,
  and the four gaps between two program lists and a shared-program claim
- `docs/notes/architecture.html` — visual explainer of the batching design
- `docs/reports/` — QC, clustering and theta/labelling reports (stage 01)
- `docs/path_map.md` — where files moved in the October 2026 reorganisation

## Upstream

- IDS paper repo: https://github.com/aradha/interdependence_scores
- Lab single-cell pipeline (QC/Harmony/UMAP reference):
  https://github.com/CutaneousBioinf/single-cell_pipeline
