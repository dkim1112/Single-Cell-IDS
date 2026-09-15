# Joint skin + PBMC: shared cell types, aligned UMAP, shared gene programs

Three questions: (1) drop cell types that exist in only one tissue, (2) overlap the two
tissues on one UMAP and fix batch with Harmony, (3) does program overlap improve, and
which high-scoring programs appear in both.

---

## 1. Which cell types are actually shared

**Only two: T cell and Myeloid cell.**

| | skin (GSE173706) | PBMC (GSE194315) |
|---|---|---|
| **shared** | T cell 4,466 · Myeloid cell 1,435 | T cell 9,297 · Myeloid cell 5,316 |
| tissue-only | Keratinocyte 41,922 · Fibroblast 9,326 · Endothelial 3,734 · Smooth muscle 3,354 · Melanocyte 1,987 · Mast cell 563 · Eccrine gland 320 · Nerve cell 271 | NK Cell 1,700 · B Cell 1,388 · Other 157 (Eryth, Platelet, HSPC, ILC, Plasmablast) |
| **after filtering** | 67,378 → **5,901** (8.8% kept) | 17,858 → **14,613** (81.8% kept) |

Note how asymmetric this is: skin loses **91%** of its cells. The comparison is
therefore no longer about skin as a tissue — it is about the immune cells inside skin.

NK deserves a comment. Skin certainly contains NK cells biologically, but GSE173706's
annotation has no NK label (its cytotoxic cells are called Cytotoxic T / Tc17 / gdT).
Rather than force an equivalence the annotations do not support, NK was treated as
PBMC-only. The cytotoxic *program* turns out to cross anyway — see §3.

Joint object: **20,514 cells × 32,991 shared genes**.

---

## 2. Aligning the two tissues

PCA → Harmony on `tissue` → UMAP. Mixing score = fraction of a cell's 30 nearest
neighbours from a different group, divided by the value expected if fully intermingled.

| key | perfect | before | after |
|---|---|---|---|
| **tissue** | 0.410 | **0.039** | **0.322** |
| main cell type | 0.442 | 0.003 | 0.003 |
| subtype | 0.879 | 0.349 | 0.438 |
| donor | 0.910 | 0.738 | 0.768 |

Tissue mixing went from 9.5% of perfect to 78% of perfect, and **cell-type separation
did not move at all (0.003 → 0.003)** — corrected without smearing T into myeloid.

### The alignment is partial, and asymmetric

Neighbour composition after correction, which is the honest read:

| tissue | cell type | n | % neighbours from the other tissue | expected if fully mixed |
|---|---|---|---|---|
| skin | T cell | 4,466 | 42.9% | 71.2% |
| skin | Myeloid cell | 1,435 | 30.4% | 71.2% |
| PBMC | T cell | 9,297 | **3.2%** | 28.8% |
| PBMC | Myeloid cell | 5,316 | **1.0%** | 28.8% |

Skin cells move into PBMC territory; PBMC cells largely do not move to meet them. So
the skin immune cells are being absorbed onto the *edge* of the PBMC cloud rather than
interleaving with it. Partly this is size (PBMC is 2.5× larger), but it also says the
two populations are not the same cells — tissue-resident immune cells differ from
circulating ones, which is real biology rather than a correction failure.

**Do not over-read the UMAP here.** Harmony corrects the embedding, never expression,
so none of §3 depends on it.

Figures: `results/joint/umap_{before,after}_{tissue,main}.png`.

---

## 3. Programs, re-discovered inside the shared cell types

IDS was re-run per tissue on **equal-sized sets (5,901 cells each)** of T + Myeloid
cells only, gene floor applied *within* that restricted set.

### First: this fixed the gene-floor problem

The blocker identified last time was that skin is 62% keratinocyte, so a 2% detection
floor across all cells deleted T and myeloid marker genes. Restricting to the shared
cell types removes that:

| PBMC program | genes passing the 2% floor in skin: all cells → T+myeloid only |
|---|---|
| **P04_PYHIN1 (T/NK cytotoxic)** | **11/25 → 23/25** |
| P03_FCRL5 (B cell) | 3/22 → 4/22 |

GNLY 1.07% → 9.49%, PRF1 0.75% → 8.42%, NKG7 1.89% → 20.32%.

The B-cell program does **not** recover, and that is the correct answer: PAX5 0.00%,
TCL1A 0.00%, MS4A1 1.93%. Those genes are genuinely absent from skin rather than
floor-masked — which is exactly why B cells were dropped from the comparison.

### New programs

**Skin (T + myeloid), 9 programs:**

| | program | n | median IDS | genes |
|---|---|---|---|---|
| 1 | P01_MT2A | 3 | 0.782 | MT2A, MT1E, MT1X |
| 2 | P02_MARCKSL1 | 3 | 0.592 | MARCKSL1, LAMP3, FSCN1 — migratory DC |
| 3 | P03_ACTB | 3 | 0.549 | ACTB, TMSB4X, PFN1 — housekeeping |
| 4 | P04_GZMA | 3 | 0.540 | GZMA, CCL5, NKG7 |
| 5 | P05_CCL3 | 3 | 0.520 | CCL3, CCL3L3, CCL4L2 |
| 6 | P06_S100A9 | 3 | 0.505 | S100A9, S100A8, S100A7 |
| 7 | P07_HSPD1 | 9 | 0.505 | heat shock — dissociation artifact |
| 8 | P08_G0S2 | 11 | 0.479 | G0S2, IL1B, CXCL8, CXCL3, CXCL2, EREG, SOD2, BCL2A1 |
| 9 | P09_C1QA | 21 | 0.438 | C1QA/B/C, CD14, MS4A6A, AIF1, CYBB, STAB1 |

**PBMC (T + myeloid), 4 programs:**

| | program | n | median IDS | genes |
|---|---|---|---|---|
| 1 | P01_HBB | 3 | 0.945 | HBB, HBA2, HBA1 — ambient |
| 2 | P02_GNLY | 11 | 0.544 | GNLY, FGFBP2, GZMA, CTSW, PRF1, KLRD1, GZMH, GZMB, CCL5, CST7, NKG7 |
| 3 | P03_FCGR3A | 3 | 0.529 | FCGR3A, CDKN1C, TCF7L2 — CD16 mono |
| 4 | P04_CCDC50 | 3 | 0.510 | CCDC50, PLD4, TCF4 — pDC |

### Does the overlap improve? Yes, but module-list overlap is the wrong test

| | before (all cell types) | now (shared cell types) |
|---|---|---|
| best Jaccard | 0.083 (one gene) | **0.300** |
| pairs with zero overlap | 17 of 21 | 35 of 36 |

The single hit is real: skin **P04_GZMA** (GZMA, CCL5, NKG7) is entirely contained in
PBMC **P02_GNLY** — 3 of 3 genes. Jaccard reads 0.300 only because PBMC's version has
11 genes and skin's has 3; containment is 100%.

But counting module overlap asks whether two clustering runs drew the same boundary,
which is a brittle question. **The better test: take a program found in one tissue and
ask whether its genes still co-depend in the other tissue's IDS matrix.**

---

## The answer: which programs appear in both

Median IDS among a program's genes, computed in each tissue. A program counts as shared
if it clears the *other* tissue's noise floor (skin p99 = 0.169, PBMC p99 = 0.208).

| found in | program | IDS at home | IDS in the other tissue | what it is |
|---|---|---|---|---|
| skin | **P06_S100A9** | 0.505 | **0.939** | **S100A8 / A9 / A7 — calprotectin & psoriasin** |
| skin | **P04_GZMA** | 0.540 | **0.696** | **cytotoxic effector** |
| PBMC | **P02_GNLY** | 0.544 | **0.302** | **cytotoxic effector (fuller version)** |
| skin | **P09_C1QA** | 0.438 | **0.309** | **macrophage / complement** |
| skin | **P08_G0S2** | 0.479 | 0.217 | inflammatory myeloid — marginal, p99 is 0.208 |
| PBMC | P01_HBB | 0.945 | 0.988 | ambient haemoglobin — **not biology** |
| skin | P03_ACTB | 0.549 | 0.288 | cytoskeleton — **housekeeping** |
| PBMC | P04_CCDC50 | 0.510 | 0.173 | pDC — at the floor, treat as not shared |

**Did not cross:** skin's metallothionein (0.080 in PBMC), skin's heat-shock program
(0.082 — dissociation artifact, correctly tissue-specific), skin's migratory-DC program
(only 1 of 3 genes measured in PBMC), and PBMC's CD16-monocyte program (0.016 in skin).

### Four genuinely shared programs

1. **Cytotoxic effector.** Found independently in both tissues and crossing in both
   directions (skin→PBMC 0.696, PBMC→skin 0.302). The strongest and least ambiguous
   result. PBMC's version is fuller (GNLY, PRF1, GZMB, KLRD1) because blood has more
   cytotoxic cells; skin's core is GZMA/CCL5/NKG7.
2. **S100A8 / A9 / A7.** Skin 0.505 → PBMC **0.939**. Worth the most attention: S100A7
   is psoriasin, S100A8/A9 are calprotectin, and they are the canonical psoriasis
   lesional biomarkers. The same co-expression holds in circulating myeloid cells.
3. **Macrophage / complement.** C1QA/B/C with CD14 and MS4A6A. Skin 0.438 → PBMC 0.309.
4. **Inflammatory myeloid.** IL1B, CXCL8, CXCL2/3, EREG. Skin 0.479 → PBMC 0.217, only
   just above the floor — report as suggestive, not established.

---

## Caveats worth stating before this goes anywhere

- **One discovery run per tissue.** No recurrence. Everything above is n=1 on each side.
- **Skin lost 91% of its cells.** This is a statement about immune cells in skin, not
  about skin.
- **PBMC pools 9 donors, skin pools 22.** Neither side is donor-controlled here; that was
  traded away to get enough cells. Donor mixing was 0.738 before correction, so donors
  were already reasonably intermingled, but it is not a clean design.
- **Shared cell type is not the same as shared signal.** Programs 1, 3 and 4 are things
  T cells and macrophages do in any tissue. The interesting claim for the abstract is
  program 2, because S100A8/A9 is an activation state rather than a lineage identity.
- **Ambient RNA still unaddressed.** Haemoglobin is the top PBMC program again and it
  "crosses tissues" at 0.988, which is a preview of exactly the false positive SoupX
  needs to remove before any shared-program claim is published.

## Files

`results/joint/` — `programs_{skin,PBMC}.csv`, `program_crosstissue_scores.csv` (the
table above), `program_overlap_shared_celltypes.csv`, `gene_floor_rescue.csv`,
`joint_mixing.csv`, `genes_{skin,PBMC}.json`, and the four UMAPs.
