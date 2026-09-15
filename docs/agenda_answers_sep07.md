# Answers — Agendas 1, 2, 3 (Sep 7)

---

# Agenda 3 — the split keratinocyte / endothelial clusters

**Your read was half right, and the half that was wrong is the more interesting half.**

## The two blobs you pointed at are biology, not batch

I split each lineage's UMAP cloud into its two visible blobs and cross-tabulated.

**Endothelial.** The two blobs have almost identical disease and chemistry
composition (49.6% vs 56.2% PP; 62% vs 70% v2_lib), so neither explains the split.
What does explain it is our own program **P03_PROX1**:

| | blob 1 (472 cells) | blob 2 (3,262 cells) |
|---|---|---|
| P03_PROX1 mean score | **6.63** | **0.41** |

PROX1 / MMRN1 / CCL21 are lymphatic endothelial markers, and PROX1 is the master
regulator of lymphatic endothelium. So blob 1 is **lymphatic** endothelium and
blob 2 is **vascular** endothelium.

The published annotation has only one label, `Endothelial`. **Our IDS program
resolved a subtype the annotation does not contain.** That is a point in favour of
the method, and it is the strongest single result we have so far.

**Keratinocyte.** Here the split is the epidermal differentiation axis:

| | blob 1 | blob 2 |
|---|---|---|
| Basal | 53.5% | 1.3% |
| Spinous | 46.4% | 78.4% |
| Supraspinous | 0.1% | 20.3% |

blob 1 is basal, blob 2 is differentiated. The disease skew you would have seen
(blob 1 is 74% PP) follows from that: lesional psoriatic epidermis is
hyperproliferative and basal-expanded. It is the disease, expressed through the
differentiation axis — not a technical artifact.

Both blobs draw on 20 donors and 29 samples each, so neither is one sample forming
its own island, which is what a pure batch effect looks like.

## But you were right that there IS a batch effect — it is chemistry

Mixing score = fraction of a cell's 30 nearest neighbours from a different group,
divided by the value expected if the groups were fully intermingled. 1.0 = mixed,
0 = completely separated.

| key | before Harmony | after Harmony |
|---|---|---|
| **chemistry (v2/v3)** | **0.393** | **0.667** |
| orig.ident | 0.606 | 0.671 |
| donor | 0.609 | 0.679 |
| disease | 0.285 | 0.366 |
| celltype | 0.003 | 0.004 |
| fullsubtype | 0.105 | 0.105 |

To rule out the disease confound I held **both cell type and disease fixed** and
asked whether v2 and v3 cells still sit apart. They do, everywhere:

| cell type | disease | v2/v3 | mixing ratio |
|---|---|---|---|
| Spinous Keratinocyte | PP | 7797/2749 | 0.303 |
| Basal Keratinocyte | PP | 5347/1071 | 0.294 |
| SFRP2+ Fibroblast | PP | 2705/1176 | 0.271 |
| Endothelial | PP | 1595/471 | 0.383 |
| Naive T cell | PP | 490/205 | 0.600 |

Same cell type, same disease, still separated by chemistry. That is unambiguous
technical batch.

## Correction, and its limits

Harmony on `chemistry`: chemistry mixing 0.393 -> 0.667, and **cell type separation
is untouched (0.003 -> 0.004)** — corrected without smearing the biology. Both
biological splits above survive.

Three things to be clear about:

1. **Why chemistry and not sample.** `orig.ident` would be the finer batch unit,
   but one sample = one donor + one phenotype, so correcting on it would erase the
   PN-vs-PP difference — the biology we are after. Chemistry is technical and, for
   the 11 paired donors, PN and PP always share the same chemistry, so the PN/PP
   contrast is chemistry-matched and survives correction.
2. **The NS comparison is still compromised.** NS is 100% v3_lib, so disease mixing
   moved 0.285 -> 0.366: correcting chemistry necessarily pulls NS toward the rest.
   Donor 30696 appears as both `30696` (v2) and `30696V3` (v3) and is the control
   for measuring how much of an NS-vs-PP difference is really chemistry.
3. **This does not change any IDS result.** Harmony corrects the embedding, never
   expression. IDS ran inside one sample (PP_8940), so it never saw cross-sample
   batch. The programs, the enrichment table and the rankings are all unchanged.
   What improved is the picture.

Files: `results/skin-batch/` (before/after UMAPs, `batch_mixing.csv`).

**On combining skin with PBMC afterwards** — see the last section of Agenda 2. The
gene floor, not batch, is currently the blocking problem for that.

---

# Agenda 1 — programs per cell type

`programs_by_celltype()` transposes the enrichment table: instead of "given this
program, which cell type", it answers **"given this cell type, which programs"**.

**Ranking by mean score alone does not work.** The mitochondrial program scored ~5
in all 26 skin cell types and won 24 of them — strongest everywhere, therefore
informative about nothing. So the ranking key is

> **score_gap = this cell type's mean score − the median mean score across all cell types**

i.e. how far above the *typical* cell type this one sits. A program every type
expresses equally has a gap near zero whatever its absolute level. Subtraction
rather than a ratio, because mean scores go negative.

This also fixes the mast-cell false negative from last time: P06_FCER1G has
enrichment 1.04 (looks non-specific) but score_gap **+4.79** in mast cells.

## Main program per cell type — skin

| cell type | n | top program | gap | enrich |
|---|---|---|---|---|
| Melanocyte | 1,987 | P10_QPCT | +5.10 | 4.07 |
| Mast cell | 563 | P06_FCER1G | +4.79 | 1.04 |
| Eccrine gland | 316 | P08_CHI3L1 | +2.96 | 2.56 |
| Treg | 662 | P09_TNFRSF18 | +2.80 | 9.44 |
| Tc17 cell | 459 | P09_TNFRSF18 | +2.79 | 9.44 |
| IFN T cell | 70 | P09_TNFRSF18 | +2.62 | 9.44 |
| COL11A1+ Fibroblast | 1,097 | P07_COL11A1 | +2.29 | 8.17 |
| Naive / Cytotoxic / Stressed T | 1,154 / 527 / 1,428 | P09_TNFRSF18 | +2.3 / +2.27 / +2.24 | 9.44 |
| Smooth muscle | 3,356 | P02_MT2A | +1.75 | 1.59 |
| RAMP1+ Fibroblast | 131 | P07_COL11A1 | +1.50 | 8.33 |
| Endothelial | 3,734 | P03_PROX1 | +1.30 | 3.85 |
| Supraspinous Keratinocyte | 4,190 | P05_PNPLA1 | +1.09 | 4.10 |
| SFRP4+ Fibroblast | 217 | P07_COL11A1 | +0.84 | 7.71 |
| Basal Keratinocyte | 11,715 | P02_MT2A | +0.27 | 1.49 |
| Spinous Keratinocyte | 26,018 | P08_CHI3L1 | +0.27 | 1.55 |

Every T-cell subtype converges on the same program, every fibroblast subtype on
the same program — which is the expected behaviour and a sanity check.

**The two keratinocyte populations are where we have the least.** Spinous
keratinocytes are 26,018 cells — 39% of the dataset — and their best program is a
+0.27 gap. That is the 279-gene cornified-envelope community the size filter drops.
Raising `--min-ids` is the fix and is still the top open item.

## Main program per cell type — PBMC

| cell type | top program | gap |
|---|---|---|
| B naive / intermediate / memory | P03_FCRL5 | +3.91 / +3.37 / +2.56 |
| NK | P04_PYHIN1 | +3.84 |
| CD16 Mono | P02_HES4 | +4.28 |
| CD8 TEM / CD4 CTL / gdT / MAIT | P04_PYHIN1 | +2.63 / +2.69 / +2.26 / +1.64 |
| CD14 Mono | P02_HES4 | +0.40 |
| CD4 TCM / CD4 Naive / Treg | (nothing above +0.15) | — |

PBMC's four programs cover B, NK/cytotoxic and CD16 monocytes and leave the CD4
compartment — the largest cell type in the dataset — with no program at all.

Files: `results/by-celltype/`.

---

# Agenda 2 — clarification questions

## Q1. Did you compare our filtering with the pipeline's?

**No, I had not. You are right to ask.** I inherited GSE194315's own
`IncludedInStudy` flag instead of running the lab pipeline's QC. I have now done
the comparison.

Lab pipeline thresholds, from `config/config_sample.json` in
CutaneousBioinf/single-cell_pipeline: `min_gene_per_cell: 200`,
`max_mt_percentage: 10.0`, `max_doublet_score: 0.15`.

On PBMC-01-1 (27,299 cells):

| filter | keeps | % |
|---|---|---|
| min_gene >= 200 | 27,206 | 99.7% |
| pct_mt <= 10 | 18,727 | 68.6% |
| singlet (demuxlet, standing in for doublet score) | 20,517 | 75.2% |
| **lab, all three** | **14,186** | **52.0%** |
| **study `IncludedInStudy`** | **17,858** | **65.4%** |

They agree on **76.1%** of cells. The disagreement:

- **study keeps, lab drops: 5,094 cells — every single one fails `pct_mt <= 10`,
  and none fail on gene count or doublets.** The whole disagreement in that
  direction is one threshold: the lab's 10% mitochondrial cut is stricter than what
  GSE194315's authors used.
- lab keeps, study drops: 1,422 cells. All singlets, median 1,745 genes, median
  8.0% mito — i.e. healthy-looking cells the study dropped for reasons its metadata
  does not expose.

**Two pipeline steps we have never run**, and one of them matters a lot:

- `run_SoupX.py` — ambient RNA removal. This is the fix for the haemoglobin program
  that ranked #1 in PBMC and the mitochondrial program that ranks #4 in skin. I had
  proposed building a gene-family exclusion list; **SoupX is the principled version
  of that and the lab already has it.** Use it instead.
- `run_Doublet.py` — a doublet score. For PBMC, demuxlet genotype demultiplexing is
  arguably stronger evidence than a score, but skin is not multiplexed and has had
  no doublet filtering at all.

Files: `results/qc-compare/`.

## Q2. What do "paired donors" and "rerun" mean?

**Paired donors** — yes, exactly what you assumed: the same patient contributed
**both a non-lesional (PN) and a lesional (PP) biopsy**. 11 of the 22 donors: 929,
8940, 31170, 31277, 369PC, 5851, 7802ED, 8659ED, 9709PC, 30696, 30696V3. Nothing
to do with skin-vs-blood; the two cohorts are unpaired there as far as we know.
It matters because a PP-vs-PN comparison within one patient removes donor genetics,
and those pairs are always the same chemistry, so it removes batch too.

**"Rerun"** — my wording was unclear. I ran the identical analysis twice, changing
only which sample IDS was computed on. First on **PP_929**, which turned out to be
96% keratinocytes and 32nd of 33 samples for cell-type diversity; the result was
mostly noise. Then on **PP_8940** (36% keratinocyte, 6 cell types with >=50 cells),
which gave the 10 programs. Same code, same parameters, different discovery sample.
The PP_929 output was not kept.

## Q3. What do `agree` and `specific` mean?

Both come from `specificity_summary.csv`, one row per program.

Each program's cell types are ranked twice — once by **proportion** (in how many
cells is the program on at all, score > 0) and once by **mean score** (how strongly).

- **`agree`** — do the two rankings name the same cell type? `False` means they
  disagree and you should look at the genes yourself. For P03 proportion said
  fibroblast and mean score said endothelial; PROX1/MMRN1/CCL21 settle it in favour
  of endothelial, so proportion was wrong.
- **`specific`** — the automatic verdict on "is this program tied to one cell type",
  computed from enrichment and the gap to the typical cell type.

**`specific` is not reliable in skin and should be overruled.** It reads False for
the mast-cell program, because with a `score > 0` cutoff that program is "on" in
96.4% of all skin cells, so every cell type reaches proportion ~1.00 and the
enrichment ratio collapses to 1.04. The biology is not in doubt — CPA3, TPSAB1,
TPSB2, LTC4S, CTSG. Use `score_gap` (Agenda 1) instead: mast cells sit +4.79 above
the typical cell type there. score_gap does not saturate.

## Q4. How do the skin programs compare with PBMC, for the same cell types?

**They barely overlap at all — and the reason is a bug in our own setup, not biology.**

Of 21 program pairs, 17 share no genes. The best is a Jaccard of 0.083, which is one
gene. So on the face of it, skin and blood share nothing.

That conclusion is not safe, because the comparison is crippled by the per-dataset
gene floor. Checking every gene of PBMC's two clearest programs against the skin
matrix:

| PBMC program | genes | absent from skin data | present but below the 2% floor | survived |
|---|---|---|---|---|
| P04_PYHIN1 (T/NK cytotoxic) | 25 | **0** | **13** | 12 |
| P03_FCRL5 (B cell) | 22 | 1 | **17** | 4 |

Nothing is missing from the skin data. GNLY is detected in 1.07% of skin cells,
PRF1 0.75%, KLRF1 0.54%, MS4A1 0.19%, PAX5 0.02% — all real, all just under the 2%
floor. Skin is 62% keratinocyte and 6.6% T cell, so a floor calibrated where T and
B cells are the majority deletes their marker genes in a tissue where they are not.

This is exactly the failure predicted in `docs/scope_skin_extension.md` §3, now
measured rather than argued.

**The fix, and it is a prerequisite for the whole skin<->blood comparison:** the
gene floor must be applied **within cell type**, not across all cells. A T-cell gene
should be judged by its detection among T cells. Until that changes, any statement
about shared or unshared programs across the two tissues is uninterpretable.

There is probably also a real biological difference underneath — PBMC's T/NK program
is cytotoxic effector (GNLY, PRF1, granzymes) while skin's is T identity/regulatory
(TNFRSF18, TNFRSF4, CD2, CD247) — but we cannot separate that from the floor
artifact until the floor is fixed.

Files: `results/cross-tissue/program_overlap_skin_vs_pbmc.csv`.

---

# What this changes about the plan

1. **Per-cell-type gene floor.** Now the top blocker — it silently deletes the
   marker genes of any minority cell type, and it makes the cross-tissue comparison
   meaningless. Fix before anything else.
2. **SoupX instead of a hand-written exclusion list.** The lab pipeline already
   solves the ambient-RNA problem I proposed building around.
3. **Decide the mitochondrial threshold.** The lab uses 10%; GSE194315's authors
   were looser. Whichever we pick, it should be the same for both tissues.
4. Raise `--min-ids` (the 279-gene keratinocyte community) — unchanged from before.
5. Agenda 2 consensus across the 33 skin samples — unchanged from before.
