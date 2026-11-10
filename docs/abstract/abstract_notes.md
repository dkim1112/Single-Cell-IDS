# KSID abstract — motivation and methodology notes

Working notes, written alongside the first IDS run. Numbers marked `[TBD]` are
waiting on that run. Claims marked `[CITE]` need a reference before they go in
the abstract; they are things I believe are standard but have not verified
against a paper in this session.

---

## 1. Motivation

### The clinical question

Psoriasis is not confined to the skin. A substantial minority of patients
develop psoriatic arthritis, and there is no reliable way to tell in advance
which ones will `[CITE — need the progression rate and a source]`. If the
systemic component of psoriasis leaves a trace in circulating immune cells,
then blood is the accessible tissue in which to look for it, and skin is where
the disease is actually visible. The question that connects them is whether the
same molecular processes are running in both places.

### Why the usual approach does not answer it

Cross-tissue single-cell comparisons are normally done one of two ways, and both
inherit a problem from the cell-type labels:

- **Compare cell-type proportions.** Requires that a "CD4 TCM" in the blood
  dataset and a "CD4 TCM" in the skin dataset mean the same thing. They were
  annotated by different groups, with different panels, at different
  resolutions. In our own data the skin and blood annotations do not even use
  the same vocabulary for the same lineages.
- **Compare differential expression within each matched cell type.** Same
  dependence on the labels, plus it needs the matching to be one-to-one. Our
  myeloid result shows it is not: blood myeloid cells are monocytes and skin
  myeloid cells are Langerhans and dermal dendritic cells, so there is no
  correspondence to compute DE across.

### Why a gene-program framing helps

A gene program is discovered from gene–gene statistical dependence across
cells, **without using the cell-type labels at all**. That changes the question
from "is this blood cell type the same as that skin cell type" (which needs the
labels to be trustworthy) to "is this set of co-dependent genes the same set in
both tissues" (which needs only the gene names). Three consequences:

1. The comparison survives our annotation problems. We can run it now, on the
   current provisional labels, because the labels are used only *after*
   discovery, to describe where a program is expressed.
2. A program can be shared across tissues while living in *different* cell
   types in each — which is a result the cell-type-matching approach cannot
   express, and is biologically plausible for an interferon or stress response.
3. Both outcomes are informative. A program present in both tissues is a
   candidate systemic axis of the disease; a program confined to one tissue
   tells us what the other tissue does not report.

### Why IDS rather than correlation or NMF

- **Against Pearson correlation:** IDS measures general statistical dependence,
  not linear association. In single-cell data a program is frequently *off* in
  most cells and sharply *on* in a subset, giving a threshold-like or
  non-monotone relationship between its member genes that a correlation
  coefficient understates.
- **Against NMF / cNMF and similar factorisations:** those impose a chosen
  number of factors and give every gene a loading on every factor. IDS gives a
  gene–gene dependence matrix, from which most genes end up in **no** program —
  which is the honest outcome, since most genes are not part of a coordinated
  module. The number of programs is read off the data rather than set in
  advance.
- **Practical point for this project:** the method is the lab's existing
  direction of work, and the implementation is the authors' own
  (`ids.compute_IDS`), called unmodified.

### The gap this fills

To my knowledge no study has asked which gene programs are shared between
psoriatic blood and psoriatic skin, nor whether psoriatic arthritis carries
blood programs that cutaneous-only psoriasis does not `[CITE — needs a
literature check before claiming novelty]`.

---

## 2. What this design can and cannot answer

This is the part to be most careful about, because the three disease states are
not available in both tissues.

| contrast | blood | skin | answerable? |
|---|---|---|---|
| healthy vs psoriasis | yes | yes | **yes, in both tissues** |
| PsA vs cutaneous psoriasis | yes | no skin PsA samples | **blood only** |
| lesional vs non-lesional | not applicable | yes (PP vs PN) | **skin only** |
| programs shared across tissues | — | — | yes, within healthy and within psoriasis |

So the abstract can claim a cross-tissue comparison of psoriasis, and a
blood-based comparison of PsA against cutaneous psoriasis. It cannot claim a
cross-tissue comparison *of PsA*. Any sentence implying skin PsA data would be
false.

Disease-state counts after QC (donors in brackets):

| tissue | healthy | psoriasis | PsA |
|---|---|---|---|
| blood | 17,118 cells (16) | 14,338 (15; 14 in donor-level tests) | 42,670 (23) |
| skin | 10,747 (8) | 56,620 (25: 11 non-lesional, 14 lesional) | — |

**And one structural fact that changes the design.** No lane in GSE194315 contains
both PsA and healthy subjects. Nine of the eleven PsA donors we have come only
from PBMC-01-1, which has no healthy and no psoriasis subjects at all, so removing
it leaves two PsA donors and the contrast cannot be tested. The PsA comparison in
the current data is therefore not separable from a lane effect — see the results
note. The fix is specific: lanes PBMC-06-1 to 06-4 each carry **12 PsA and 6
psoriasis donors in the same lane**, which gives a fully within-lane
PsA-versus-psoriasis comparison. PsA-versus-healthy is not available within any lane in this dataset and will always be cross-lane.

*Update, Oct 22: the PBMC-06 lanes have since been added (nine lanes in total), so the within-lane PsA-versus-psoriasis comparison now exists (6 psoriasis vs 12 PsA donors in the PBMC-06 lanes; 14 vs 23 across all lanes). The counts in this section's earlier paragraph and table describe the five-lane data.*

---

## 3. Methodology

### 3.1 Data

- **Blood.** GSE194315, PBMC CITE-seq. Nine 10x lanes (PBMC-01-1, 03-1, 03-3, 05-1, 05-2, 06-1 to 06-4; the 06 lanes were added later). Each lane pools nine genetically distinct subjects, demultiplexed
  with demuxlet, so a lane is not a sample. Kept: the source study's own
  inclusion verdict (`IncludedInStudy`) and demuxlet singlets
  (`DemuxletDropletType == SNG`), then subjects with `Status` in Healthy, PSO or
  PSA. AXI and PSX subjects dropped as different diseases.
- **Skin.** GSE173706, 22 donors, with the source study's own filtering applied.
  Disease states PP (psoriasis lesional), PN (psoriasis non-lesional), NS
  (normal skin).
- **Lane choice is deliberate and load-bearing.** Lanes 03 and 05 each contain
  Healthy, PSO and PSA subjects *within the same lane*, so disease state is not
  confounded with the sequencing run. This is why four lanes were added rather
  than one: the lane we started with (PBMC-01-1) is 100% PsA, and on its own it
  cannot distinguish a disease effect from a batch effect.

### 3.2 Cell QC

The lab pipeline's thresholds, applied identically to both tissues: at least 200
genes detected per cell, at most 10% mitochondrial reads. No doublet-score
filter — the skin data has none available, and for blood the demuxlet singlet
gate already removes cross-donor doublets. Blood 74,126 cells after QC across nine lanes (the before-QC count was logged only for the first five lanes: 47,186 -> 38,182); skin 67,378 -> 67,367 (the skin study had already
applied a 10% mitochondrial cut, so almost nothing is removed a second time).

### 3.3 Representation for IDS (three representations, on purpose)

IDS runs on library-size-normalised counts (`normalize_total`, target 1e4),
per-gene min–max scaled to [0, 8], **with no log1p step**. The omission is
deliberate and was established empirically in this project: the IDS Gaussian
feature map carries an `exp(-B·x²)` envelope with B = 1/2, so values near 8
contribute essentially nothing. Log-transforming before the [0,8] scaling
spreads cells evenly across that range and pushes much of the data into the
saturated region — measured as a strong module's mean IDS falling from about
0.89 to about 0.09. Raw and library-normalised counts are right-skewed, with most
mass near zero where the features are informative, which is what the [0,8]
scaling is calibrated for.

Per-cell program *scores* are computed on the standard log1p representation
instead, because `sc.tl.score_genes` subtracts the mean of an
expression-matched control gene set and that subtraction is meaningless in [0,8]
space. **Programs are discovered in IDS space and scored in log space.** The
distinction is worth one sentence in the methods, since it looks like an
inconsistency until it is explained.

### 3.4 Gene panel

IDS forms a (k·d) × (k·d) correlation tensor internally, so the number of genes
d is the binding constraint. Two choices here:

- A **single fixed gene panel is used for every stratum**, rather than
  per-stratum highly variable genes. Program comparison across tissues and
  diseases is a gene-set comparison, and it is only meaningful if every stratum
  had the opportunity to use the same genes. Per-stratum HVG selection would
  make a program look tissue-specific simply because its genes were never
  offered to the other tissue.
- The panel is the **union of highly variable genes across the six strata**,
  plus a curated marker set added unconditionally. The curated addition guards a
  known failure mode of mean–variance HVG selection: a program that is off in
  most cells has a low mean and may not rank as highly variable, so exactly the
  bimodal programs we are looking for can be dropped before IDS sees them.

Panel: 4,240 genes (union of per-stratum highly variable genes, plus 71 curated genes not already in it). [Oct 22: the current panel file for the nine-lane run has 4,273 genes, 4,249 highly variable plus 24 curated.]

**A detection floor is then applied within each stratum** (the registered default,
2% of cells), and it is load-bearing. Without it the ranking is topped by
three-gene modules of lncRNAs and pseudogenes at median IDS 1.00 — genes detected
in a handful of cells whose few non-zero cells coincide, so IDS saturates. With it,
the panel drops to 831–884 genes in blood and 1,154–1,345 in skin, and the
programs become interpretable. The floor is applied per stratum rather than to the
shared panel, so every stratum is offered the same genes and a program can be
tissue-specific for a biological reason rather than because panel construction
never showed the gene to the other tissue.

**What the floor costs, stated plainly.** IL17A, IL4 and IL13 fall below it in
every stratum, so the IL-17 axis — the pathway psoriasis is best known for — is
invisible to this whole-tissue analysis. GZMA and NKG7 fall below it in healthy and
non-lesional skin but not lesional skin. Both are consequences of measuring
detection across all cells of a tissue; the remedy is to run IDS within lineage,
where a T-cell cytokine is no longer diluted by keratinocytes.

### 3.5 IDS computation

The authors' `ids.compute_IDS` is called unmodified, with `num_terms=6`,
`p_norm='max'`, `bandwidth_term=0.5`. The gene–gene matrix is assembled
block-wise over the gene axis (diagonal blocks as `compute_IDS(X_i)`,
off-diagonal as `compute_IDS(X_i, Y=X_j)`, lower triangle by symmetry), which
reproduces a single full call to floating-point precision while keeping peak
memory bounded. IDS is never reimplemented anywhere in this project.

Permutation p-values (the authors' per-gene-pair test, 60 shuffles) are computed after programs are found, on the union of each condition's program genes, not inside the IDS matrix itself. Edge selection when building programs still uses an effect-size threshold only. Programs are then required to have at least 90% of their gene pairs beat all 60 shuffles; that program-level rule is ours, not the authors'. (Earlier version of this note: permutation was not computed at this scale.)

### 3.6 From the IDS matrix to programs

An edge is kept between two genes when their IDS is at or above `min_ids` and
the edge is mutual. `min_ids` is **read off each run's own off-diagonal IDS
distribution** — the value where the background bulk ends and the program mode
begins — and is not carried over from another dataset. The reason an
effect-size floor is needed at all: with thousands of cells, library-size
normalisation induces weak but statistically significant dependence across many
genes, so a significance filter alone lets background edges bridge real programs
into one giant component.

Programs are then the communities of that graph under greedy modularity
(weighted by IDS), subject to a **selection rule fixed before any biology was
looked at**: keep programs of 3 to 50 genes; rank by recurrence across strata
first, then median within-program IDS; report all that pass rather than a
hand-picked subset. Median rather than mean, so a program cannot be inflated by
two near-duplicate genes. Programs excluded by the size filter are written out
every run, so the boundary is auditable.

### 3.7 Where each program is expressed

Every cell receives a score for every program (`sc.tl.score_genes`, control set
matched on expression bin). "Expressing" is `score > 0`, fixed in advance and
not tuned: the score is already background-subtracted, so 0 is the value
expected when a program's genes are no higher than equivalently expressed genes
elsewhere.

Programs are then summarised over three axes — cell type, tissue, disease state.

### 3.8 Statistics — the part most likely to be challenged

- **The unit of replication is the donor, not the cell.** Program scores are
  averaged within (donor × cell type) before any disease comparison, and the
  comparisons are made across donors. Testing across cells would treat 2,000
  cells from one person as 2,000 independent observations and would make almost
  every contrast significant.
- **Cell-count imbalance.** The strata differ in size by more than an order of
  magnitude (skin keratinocytes dominate the skin data; PBMC-01-1 dominates the
  PsA cells). Donor-level aggregation handles this for the disease contrasts;
  for program discovery, strata are capped at a common cell count so that one
  large stratum does not set the gene panel on its own.
- **One residual confound, stated rather than buried.** PBMC-01-1 contributes
  about 17,900 PsA cells while lanes 03 and 05 contribute roughly 3,400 between
  them, so the PsA group is dominated by a single lane even though PsA also
  appears inside the balanced lanes. The PsA-versus-psoriasis contrast is
  therefore reported **twice**: once with all lanes, and once with PBMC-01-1
  excluded, using only the within-lane PsA subjects. A result that holds in both
  is a disease effect; one that appears only with PBMC-01-1 included is a lane effect wearing a disease label. [Oct 22: the current analysis uses the PBMC-06 lanes instead (6 psoriasis vs 12 PsA donors in the same lanes) and, separately, all lanes (14 vs 23); no blood comparison was significant.]
- **Skin lesional and non-lesional samples are paired within donor**, so that
  contrast is paired and should not be treated as two independent groups.
- **Multiple testing** across programs × cell types × contrasts is controlled
  with Benjamini–Hochberg, reported alongside effect sizes rather than instead
  of them.

### 3.9 Matching programs across tissues

Two independent criteria, reported together:

1. **Gene overlap** — Jaccard index between a blood program's gene set and a
   skin program's, with a null from random gene sets drawn from the same panel.
2. **Shared usage pattern** — correlation of the two programs' scores across the
   cell types present in both tissues.

A program pair that matches on genes but not on usage is the interesting case:
the same genes running in a different cell type in the other tissue.

---

## 4. Open items

- `[CITE]` the PsA progression rate, and a novelty check on cross-tissue
  psoriasis program comparisons.
- Final cell-type labels. The run uses the current provisional annotations; the
  NK/γδ cluster is awaiting outside confirmation. Programs are discovered without
  labels, so relabelling changes the *description* of a program's location, not
  the program.
- Whether the gene filter question affects any of this. It should not: the panel
  here is built for IDS and is unrelated to the 1% filter used for the Harmony
  embedding, but worth stating explicitly so the two are not conflated.
