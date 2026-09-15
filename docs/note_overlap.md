# [Do the programs overlap now? — and which ones appear in both]

## Setup

- IDS re-run separately on each tissue, using **only T + Myeloid cells**.
- Equal size on both sides: **5,901 cells each**, so neither tissue wins by having more cells.
- The 2% gene floor is now applied **inside** that restricted set, not across all cells.

---

## (a) First — this fixed the gene floor problem

Before: skin is 62% keratinocyte, so a floor set across all skin cells deleted T-cell and myeloid marker genes. They looked "absent from skin" when they were just rare in a keratinocyte-dominated pool.

| PBMC program | genes clearing the 2% floor in skin (all cells → T+myeloid only) |
|---|---|
| **P04_PYHIN1 (T/NK cytotoxic)** | **11/25 → 23/25** |
| P03_FCRL5 (B cell) | 3/22 → 4/22 |

- GNLY 1.07% → 9.49%, PRF1 0.75% → 8.42%, NKG7 1.89% → 20.32%.
- The T/NK genes were never missing. They were hidden by the wrong denominator.
- The B-cell program does **not** recover — PAX5 0.00%, TCL1A 0.00%, MS4A1 1.93%.
- That is the right answer, not a failure: those genes really are absent from skin. It also confirms dropping B cells was correct.

→ The floor bias is a **denominator problem**, and restricting to shared cell types is the fix we actually applied. (A per-cell-type floor inside the code is still not implemented.)

---

## (b) Programs found after the restriction

**Skin (T + myeloid) — 9 programs**

| # | program | n | med IDS | reading |
|---|---|---|---|---|
| 1 | P01_MT2A | 3 | 0.782 | metallothionein — stress |
| 2 | P02_MARCKSL1 | 3 | 0.592 | migratory DC (LAMP3, FSCN1) |
| 3 | P03_ACTB | 3 | 0.549 | housekeeping |
| 4 | P04_GZMA | 3 | 0.540 | cytotoxic (GZMA, CCL5, NKG7) |
| 5 | P05_CCL3 | 3 | 0.520 | chemokine |
| 6 | P06_S100A9 | 3 | 0.505 | S100A9, S100A8, S100A7 |
| 7 | P07_HSPD1 | 9 | 0.505 | heat shock — dissociation artifact |
| 8 | P08_G0S2 | 11 | 0.479 | inflammatory myeloid (IL1B, CXCL8) |
| 9 | P09_C1QA | 21 | 0.438 | macrophage (C1QA/B/C, CD14) |

**PBMC (T + myeloid) — 4 programs**

| # | program | n | med IDS | reading |
|---|---|---|---|---|
| 1 | P01_HBB | 3 | 0.945 | ambient haemoglobin — ARTIFACT |
| 2 | P02_GNLY | 11 | 0.544 | cytotoxic (GNLY, PRF1, GZMB, KLRD1…) |
| 3 | P03_FCGR3A | 3 | 0.529 | CD16 monocyte |
| 4 | P04_CCDC50 | 3 | 0.510 | pDC |

---

## (c) Overlap improved — but this is the wrong test

| | before (all cell types) | now (shared only) |
|---|---|---|
| best Jaccard | 0.083 (one gene) | **0.300** |
| pairs with zero overlap | 17 of 21 | 35 of 36 |

- The one hit is real: skin **P04_GZMA** (GZMA, CCL5, NKG7) sits entirely inside PBMC **P02_GNLY** — 3 of 3 genes.
- Jaccard only reads 0.300 because PBMC's version has 11 genes and skin's has 3. **Containment is 100%.**
- Problem with this metric: comparing gene lists asks *"did two clustering runs draw the same boundary?"* Two runs can find the same biology and still get a low Jaccard just because one drew a bigger circle.

→ Better test: **take a program found in one tissue, and ask whether its genes still co-depend in the OTHER tissue's IDS matrix.** No clustering involved — the program is just carried over and re-scored.

---

## (d) The answer — which programs appear in both

A program counts as shared if its median IDS clears the other tissue's noise floor.
Noise floor = 99th percentile of IDS across all gene pairs: **skin 0.169, PBMC 0.208.**

| found in | program | IDS at home | IDS in other tissue | what it is |
|---|---|---|---|---|
| skin | **P06_S100A9** | 0.505 | **0.939** | S100A8/A9/A7 — calprotectin & psoriasin |
| skin | **P04_GZMA** | 0.540 | **0.696** | cytotoxic effector |
| PBMC | **P02_GNLY** | 0.544 | **0.302** | cytotoxic effector (fuller) |
| skin | **P09_C1QA** | 0.438 | **0.309** | macrophage / complement |
| skin | P08_G0S2 | 0.479 | 0.217 | inflammatory myeloid — marginal (floor 0.208) |
| PBMC | P01_HBB | 0.945 | 0.988 | ambient haemoglobin — **not biology** |
| skin | P03_ACTB | 0.549 | 0.288 | housekeeping — **not informative** |
| PBMC | P04_CCDC50 | 0.510 | 0.173 | pDC — at floor, not shared |

**Did not cross:** skin metallothionein (0.080), skin heat-shock (0.082 — correctly tissue-specific, it's a dissociation artifact), skin migratory DC (only 1 of 3 genes measured in PBMC), PBMC CD16 monocyte (0.016).

### Four genuinely shared programs

1. **Cytotoxic effector** — the strongest result. Found *independently* in both tissues and crosses in *both* directions (skin→PBMC 0.696, PBMC→skin 0.302). PBMC's version is fuller (GNLY, PRF1, GZMB, KLRD1) because blood has more cytotoxic cells; skin's core is GZMA/CCL5/NKG7.
2. **S100A8 / A9 / A7** — skin 0.505 → PBMC **0.939**. The one worth the abstract: S100A7 is psoriasin, S100A8/A9 are calprotectin, and these are the canonical psoriasis lesional biomarkers. The same co-expression holds in circulating myeloid cells.
3. **Macrophage / complement** — C1QA/B/C with CD14, MS4A6A. Skin 0.438 → PBMC 0.309.
4. **Inflammatory myeloid** — IL1B, CXCL8, CXCL2/3, EREG. Skin 0.479 → PBMC 0.217, only just over the floor. Report as suggestive, not established.

---

## Caveats — state these before the abstract

- **n = 1 per tissue.** One discovery run each, no recurrence. Nothing above is replicated.
- **Skin lost 91% of its cells.** This is a claim about immune cells *in* skin, not about skin.
- **Not donor-controlled.** PBMC pools 9 donors, skin pools 22. Traded away to get enough cells.
- **Shared cell type ≠ shared signal.** Programs 1, 3, 4 are things T cells and macrophages do in any tissue. Program 2 is the interesting one because S100A8/A9 is an **activation state**, not a lineage identity — that's the part that actually says "skin inflammation shows up in blood."
- **Ambient RNA still unhandled.** Haemoglobin is PBMC's top program again and it "crosses tissues" at 0.988. That is a preview of exactly the false positive SoupX has to remove before any shared-program claim is publishable.
- **No script in the repo.** This run was ad-hoc; outputs are in `results/joint/` but the code that made them was never committed.

**Files:** `results/joint/` — `programs_{skin,PBMC}.csv`, `program_crosstissue_scores.csv`, `program_overlap_shared_celltypes.csv`, `gene_floor_rescue.csv`, `joint_mixing.csv`, `genes_{skin,PBMC}.json`, 4 UMAPs.
