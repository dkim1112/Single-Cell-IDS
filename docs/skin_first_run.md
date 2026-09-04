# GSE173706 psoriasis skin — first run

```bash
# once: 33 gene-x-cell CSVs + metadata -> one sparse .h5ad (~10 min, 1.3 GB)
python -m sc_ids.datasets.psoriasis_skin --scratch /tmp

python -m sc_ids.analysis --dataset psoriasis_skin \
  --data data/GSE173706_skin.h5ad \
  --discovery-col orig.ident --discovery-value PP_8940 \
  --min-ids 0.5 --out results/skin-PP8940
```

Discovery: 2,764 cells (PP_8940, lesional) x 11,998 genes. Evaluation: all 67,378 cells.

## The dataset, as configured

`datasets/psoriasis_skin.py` records the four things that differ from GSE194315.
All were read off the metadata, not assumed:

- **Barcode key.** The metadata is one merged Seurat object, so cells are
  `<barcode>-<n>` with n the sample's index in the merge (1..33), not the `-1`
  each per-sample file carries. Suffix<->sample is exactly 1:1. Getting this wrong
  does not error, it silently matches ~1/33 of the cells.
- **QC gate = membership.** No `IncludedInStudy` column; only filtered cells appear
  in the metadata at all. **96,057 raw cells -> 67,378, matching the metadata
  exactly.** That was the stated test and it passes.
- **Cell-type mapping came free.** `fullsubtype` (26) -> `celltype` (10) are both in
  the metadata; no file needed from the lab. Seven granular labels carry 1-4 cells
  under a different main label (annotation noise) and are resolved by majority vote.
- **Ensembl IDs -> symbols**, using GSE194315's own `features.tsv.gz` (32,991/33,538
  mapped). This is not cosmetic: it puts skin and blood programs in ONE symbol
  space, without which no cross-tissue comparison is possible.

## Threshold

Skin's IDS background is **lower** than PBMC's (p99 0.158 vs 0.195; p99.9 0.261),
and neither distribution has a real gap — both decay smoothly, so PBMC's "gap at
0.5" was always a defensible cut rather than a discovered boundary. `min_ids = 0.5`
was kept for skin, deliberately: it sits further into skin's tail than into PBMC's
(0.007% vs 0.02% of pairs), and one absolute threshold in both tissues is what
makes a program edge mean the same thing in each — a precondition for comparing them.

## The stratum choice matters more than anything else here

The first run used **PP_929**, picked for cell count and paired-donor status. It is
**96% keratinocytes** (3,163 of 3,292 cells; zero fibroblasts, zero endothelium)
and ranks **32nd of 33 samples** for cell-type evenness. It produced one real
program (melanocyte, from 16 cells) and a lot of noise.

Re-running on **PP_8940** — same settings, 2,764 cells, 6 cell types with >=50
cells, 36% keratinocyte — changed the result completely.

This is the handoff's own rule biting in a new place: *stratification level
determines what IDS can see*. In GSE194315 a file pooled nine donors and every cell
type, so "one sample" was a rich stratum. **A skin sample is one biopsy dominated by
one lineage.** Cell-type evenness, not cell count, is the criterion for choosing a
skin stratum, and it should be checked before every run.

## Result: 10 programs, 8 interpretable

| rank | program | n | med IDS | top cell type | enrich | reading |
|---|---|---|---|---|---|---|
| 1 | P01_HLA-A | 5 | 0.629 | Treg | 1.02 | **MHC class I** (HLA-A/B/C/E, B2M) — every nucleated cell. Artifact, flagged diffuse |
| 2 | P02_MT2A | 3 | 0.562 | Melanocyte | 1.65 | metallothioneins — a stress response, density 1.00 |
| 3 | P03_PROX1 | 5 | 0.561 | *see below* | 4.98 | **lymphatic endothelium**: PROX1, MMRN1, CCL21, TFF3, TFPI |
| 4 | P04_MT-ND1 | 9 | 0.527 | PVM | 1.01 | **mitochondrial genome**. Artifact, flagged diffuse |
| 5 | P05_PNPLA1 | 3 | 0.501 | Supraspinous Keratinocyte | 4.10 | cornification (PNPLA1, HOXC13) |
| 6 | P06_FCER1G | 16 | 0.477 | Mast cell | 1.04 | **mast cell**: CPA3, TPSAB1/B2, LTC4S, HPGDS, CTSG |
| 7 | P07_COL11A1 | 9 | 0.441 | RAMP1+ Fibroblast | 8.33 | fibroblast ECM: COL11A1, HAPLN1, TNMD, MFAP5 |
| 8 | P08_CHI3L1 | 10 | 0.437 | Eccrine gland | 2.56 | **eccrine gland**: DCD, SCGB2A2, MUCL1, KRT7 |
| 9 | P09_TNFRSF18 | 46 | 0.396 | Treg | **9.44** | **T cell / Treg**: TNFRSF18 (GITR), TNFRSF4 (OX40), CD2, CD247, PTPRC |
| 10 | P10_QPCT | 10 | 0.362 | Melanocyte | 4.07 | **melanocyte**: MITF, PMEL, OCA2, GPR143 |

Eight land on a coherent, textbook cell population, discovered bottom-up with no
label ever entering the computation. P09's enrichment of 9.44 is the highest seen in
either tissue. P10 **replicates** the melanocyte program independently found in
PP_929 (MITF shared) — informal evidence of the recurrence that Agenda 2 formalises.

### Three things the run exposes

**1. `score > 0` saturates far more in skin than in blood.** P06 is unmistakably
mast-cell (CPA3, tryptases, LTC4S) and its mean score in mast cells is 5.64 against
1.64 in the next cell type — but *every* cell type reaches proportion ~1.00, so
enrichment is 1.04 and `specific` reads False. That is a **false negative caused by
the cutoff, not by the biology**. In skin, `mean_score` is doing the discriminating
and proportion is nearly uninformative. The cutoff should be revisited before the
cross-tissue comparison, and changed for both tissues together or not at all.

**2. `agree = False` earned its place.** For P03, proportion ranks SFRP4+ Fibroblast
first; mean score ranks Endothelial first. PROX1/MMRN1/CCL21 are lymphatic
endothelial markers, so **the mean-score answer is the biologically correct one and
the proportion answer is wrong.** Quoting the proportion ranking alone would have
mislabelled this program.

**3. The size filter is costlier here than in PBMC.** It dropped communities of 52,
99, 121, 211 and **279** genes. The 279-gene one contains IVL, SPRR1A/1B/2A-2G and
S100A7/A8/A9 — involucrin, small proline-rich proteins and psoriasin/calprotectin,
i.e. **the cornified-envelope programme that is the central biology of psoriatic
skin.** Its density is 0.13, so it is a real core with a large loose fringe, not a
279-gene program. The fix is the same one PBMC pointed to: raise `--min-ids` so
communities form tightly, rather than raising the size cap after seeing what it cut.

## Next, in order

1. **Gene-family exclusion (MT-, RP-, HLA-).** Two of ten programs here are
   artifacts, and a mitochondrial or MHC program will appear in *both* tissues and
   register as "shared" — the worst false positive available for the abstract.
   Decide the list before looking at any overlap.
2. **Raise `--min-ids` and re-run**, to test whether the 279-gene keratinocyte
   community resolves into a reportable cornification program.
3. **Dissociation-stress genes.** PP_929 produced a JUN/FOS/EGR1/HSPA1A program, the
   classic solid-tissue dissociation artifact. It did not survive in PP_8940, but it
   will recur; treat it like the family exclusions.
4. **Agenda 2 across the 33 samples.** 25 samples clear 1,000 cells, 32 clear 500;
   only NS_AR008 (42 cells) is unusable. Recurrence is what makes the composition
   problem in section 3 stop mattering.
5. **PN vs PP within donor** — 11 donors have both. Not built yet; see
   `docs/scope_skin_extension.md`.

## Watch out: chemistry is confounded with phenotype

NS is 100% `v3_lib`; PN and PP are ~2/3 `v2_lib`. Any NS-vs-PP comparison is partly
a v2-vs-v3 comparison. Donor 30696 appears as both `30696` (v2) and `30696V3` (v3)
and is the control for measuring how big that effect is. Stratifying on
`orig.ident` holds chemistry constant *within* a run; only across-phenotype
comparisons are exposed.
