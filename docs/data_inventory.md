# Data inventory

What goes in (`data/`, read-only, gitignored) and which parts of it each analysis uses.

## Blood: GSE194315, PBMC CITE-seq

28 lanes (7 runs × 4 lanes), files named `PBMC-<run>-<lane>.{barcodes,features,matrix}`. Cell metadata:
`GSE194315_CellMetadata-AS_TotalCiteseq_20220711.tsv` (804,968 droplets; 282,705 flagged `IncludedInStudy`).

**A lane is not a sample.** Each lane pools several genetically distinct donors (demultiplexed with demuxlet), and the same
donor appears in more than one lane. Disease state is therefore a property of the donor, not of the lane, and tests average
over a donor's cells across lanes first.

### The 9 lanes in use

Counts are after our filters (source study's `IncludedInStudy`, demuxlet singlets, Healthy / PSO / PsA donors only,
min 200 genes, max 10% mitochondrial, platelets and erythrocytes removed). Defined in
`sc_ids.datasets.gse194315.LANES_USED`. Every lane is built by `scripts/02_ids_programs/build_units.py`; the catalogue reads these nine.

| Lane | Cells | Donors | Healthy | PSO | PsA | What it is for |
|---|---:|---:|---:|---:|---:|---|
| PBMC-01-1 | 12,695 | 9 | 0 | 0 | 9 | The lane the first work was done on. PsA only, so no disease contrast is possible with it alone. The only blood lane in the QC, UMAP, theta and labelling work (stage 01). |
| PBMC-03-1 | 6,862 | 11 | 7 | 3 | 1 | Healthy and PSO donors in the same lane, so healthy vs PSO is not confounded with lane. |
| PBMC-03-3 | 6,723 | 11 | 7 | 3 | 1 | Same as 03-1. |
| PBMC-05-1 | 5,790 | 16 | 9 | 6 | 1 | Same role as 03. |
| PBMC-05-2 | 5,980 | 16 | 9 | 6 | 1 | Same role as 03. |
| PBMC-06-1 | 9,378 | 18 | 0 | 6 | 12 | PSO and PsA donors in the same lane, so PSO vs PsA is not confounded with lane. No healthy donors. |
| PBMC-06-2 | 9,219 | 18 | 0 | 6 | 12 | Same as 06-1. |
| PBMC-06-3 | 8,370 | 18 | 0 | 6 | 12 | Same as 06-1. |
| PBMC-06-4 | 9,109 | 18 | 0 | 6 | 12 | Same as 06-1. |
| **Total** | **74,126** | **54** | **16** | **15** | **23** | Cells by group: 17,118 healthy, 14,338 PSO, 42,670 PsA. |

How the lanes map to the blood disease tests (`scripts/03_catalogue/disease_contrasts.py`):

| Contrast | Lanes used | Donors | Valid? |
|---|---|---|---|
| Healthy vs PSO | 03 + 05 | 16 vs 8 | Yes, both groups in the same lanes |
| PSO vs PsA | 06 only | 6 vs 12 | Yes, same lanes |
| PSO vs PsA | all nine | 14 vs 23 | Yes, but PBMC-01-1 adds PsA cells in a lane with no PSO |
| Healthy vs PsA | all nine | 16 vs 23 | No: no lane contains both groups |

### All 28 lanes after QC (built November 10)

Every lane now has its own QC'd file (`scripts/02_ids_programs/build_units.py`), with the same filters as the nine above.
Antibody rows are dropped before QC: lanes from run 02 onward carry 282 antibody features, and counting them in a cell's
total lowers its mitochondrial share and lets extra cells through. Donor columns count donors with at least one cell in the lane.

| Lane | Cells | Healthy | PSO | PsA | In the October nine |
|---|---:|---:|---:|---:|---|
| PBMC-01-1 | 12,695 | 0 | 0 | 9 | yes |
| PBMC-01-2 | 10,951 | 0 | 0 | 9 |  |
| PBMC-01-3 | 13,984 | 9 | 0 | 0 |  |
| PBMC-01-4 | 14,192 | 9 | 0 | 0 |  |
| PBMC-02-1 | 3,238 | 9 | 8 | 1 |  |
| PBMC-02-2 | 2,505 | 9 | 8 | 1 |  |
| PBMC-02-3 | 3,713 | 9 | 8 | 1 |  |
| PBMC-02-4 | 3,829 | 9 | 8 | 1 |  |
| PBMC-03-1 | 6,862 | 7 | 3 | 1 | yes |
| PBMC-03-2 | 6,238 | 8 | 4 | 1 |  |
| PBMC-03-3 | 6,723 | 7 | 3 | 1 | yes |
| PBMC-03-4 | 3,248 | 8 | 3 | 1 |  |
| PBMC-04-1 | 6,915 | 4 | 9 | 0 |  |
| PBMC-04-2 | 7,391 | 4 | 9 | 0 |  |
| PBMC-04-3 | 6,954 | 4 | 9 | 0 |  |
| PBMC-04-4 | 6,331 | 4 | 9 | 0 |  |
| PBMC-05-1 | 5,790 | 9 | 6 | 1 | yes |
| PBMC-05-2 | 5,980 | 9 | 6 | 1 | yes |
| PBMC-05-3 | 4,809 | 9 | 6 | 1 |  |
| PBMC-05-4 | 5,513 | 9 | 6 | 1 |  |
| PBMC-06-1 | 9,378 | 0 | 6 | 12 | yes |
| PBMC-06-2 | 9,219 | 0 | 6 | 12 | yes |
| PBMC-06-3 | 8,370 | 0 | 6 | 12 | yes |
| PBMC-06-4 | 9,109 | 0 | 6 | 12 | yes |
| PBMC-07-1 | 5,174 | 0 | 1 | 5 |  |
| PBMC-07-2 | 5,310 | 0 | 1 | 5 |  |
| PBMC-07-3 | 1,613 | 0 | 1 | 5 |  |
| PBMC-07-4 | 5,380 | 0 | 1 | 5 |  |
| **Total** | **191,414** | **29** | **24** | **28** | 74,126 cells |

Which runs can support which comparison: runs 02 to 05 hold healthy and PSO donors in the same lanes; runs 06 and 07 hold
PSO and PsA donors in the same lanes; no lane holds healthy and PsA donors in usable numbers (run 01 has both, in different lanes).
Comparing donors across runs is not valid: healthy donors from run 01 against healthy donors from runs 02 to 05 differ in
162 of 372 tests (`results/04_method_checks/run_effect_controls.csv`).

### The 19 lanes not in use

The code and notes do not record why these were left out. Whether to add more is an open choice; the healthy vs PSO test
has only 16 vs 8 donors, so lanes that hold both groups are the most direct way to add power.

| Run | Lanes not used | Donors in run | Groups present in the run |
|---|---|---:|---|
| PBMC-01 | 01-2, 01-3, 01-4 | 9 | 01-2: same PsA donors as 01-1 (technical replicate). 01-3, 01-4: healthy only |
| PBMC-02 | 02-1 to 02-4 | 20 | Healthy, PSO, PsA, plus AXI and PSX (not used in our analyses) |
| PBMC-03 | 03-2, 03-4 | 20 | Healthy, PSO, PsA, AXI, PSX |
| PBMC-04 | 04-1 to 04-4 | 20 | Healthy, PSO, AXI, PSX (no PsA) |
| PBMC-05 | 05-3, 05-4 | 20 | Healthy, PSO, PsA, AXI, PSX |
| PBMC-07 | 07-1 to 07-4 | 16 | PSO, PsA, AXI, PSX (no healthy) |

AXI and PSX are two further groups in the source metadata. They are dropped everywhere: the build script treats them as
different diseases from Healthy, PSO and PsA.

## Skin: GSE173706

`data/GSE173706/` holds 33 count files (`GSM52771xx_<state>-<donor>.csv.gz`), metadata in `GSE173706_psorSkinMeta.tsv`.
22 donors; states NS (normal skin), PN (psoriasis, non-lesional), PP (psoriasis, lesional). A donor can contribute a PN and a
PP sample. After QC: 67,367 cells.

## Disease labels used across both tissues

| Harmonised | Blood | Skin |
|---|---|---|
| healthy | Healthy | NS |
| psoriasis | PSO | PN, PP (kept apart as `disease_fine`) |
| PsA | PSA | no skin counterpart, so PsA vs psoriasis is a blood-only contrast |

## Cell-type labels

All catalogue results use the **source annotations** as given (blood `CellType`, skin `fullsubtype`), grouped into coarse
lineages by `sc_ids.lineage.LINEAGE`. The only labels we defined ourselves are the T/NK and myeloid sub-cluster labels from
stage 01 (theta = 4), and those cover lane PBMC-01-1 plus skin T and myeloid cells only.
