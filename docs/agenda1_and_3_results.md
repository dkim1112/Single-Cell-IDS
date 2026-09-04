# Agendas 1 and 3 — what was built, and what the first run says

Run reproduced with:

```bash
python -m sc_ids.analysis \
  --data data/GSE194315_PBMC-01-07_processed_data_files --prefix "PBMC-01-1." \
  --metadata data/GSE194315_CellMetadata-AS_TotalCiteseq_20220711.tsv \
  --discovery-col Subject --discovery-value PSA26 \
  --min-cells-frac 0.02 --min-ids 0.5 --out results/PSA26
```

Discovery: 4,259 cells (PSA26, lane PBMC-01-1) × 9,542 genes.
Evaluation: all 17,858 `IncludedInStudy` cells in the lane, 9 subjects.
IDS off-diagonal: median 0.050, p99 0.195, max 0.974.

Programs are **discovered** on one subject and **scored** on the whole lane, so the
cell-type labels never touch discovery — the enrichment table below is an
out-of-sample check, not a restatement of the input.

---

## Agenda 1a — CellTypeMain

`celltype_mapping.txt` lives at `src/sc_ids/resources/celltype_mapping.txt` and is
loaded by `sc_ids.celltypes`. Coverage is exact in both directions: **30 granular
labels in the data, 30 in the map, zero unmapped cells.** (The handoff said 29; it
is 30.)

**One defect in the file as supplied.** `NK_CD56bright` and `NK Proliferating` map
to `"K Cell"`, not `"NK Cell"`. Left uncorrected, that splits the NK lineage into
two colours in every figure — `NK` in one, its two subtypes in a phantom sixth
category. `load_celltype_map` repairs it on load, prints a warning naming the two
labels, and leaves the file byte-identical so a re-send from the lab diffs cleanly.
**Worth telling whoever produced the file.**

Result: `results/PSA26/umap_celltype_main.png` — five clean lineages, legible
legend, against `umap_celltype_granular.png` with 30.

## Agenda 1b — the selection rule, fixed in advance

Stated before looking at anything: **size 3–50 genes; rank by recurrence, then
median IDS; report all that pass.** Rationale in `src/sc_ids/programs.py`.

Median rather than mean because mean rewards a large module carried by two
near-duplicate genes; `tests/test_programs.py` contains that exact case.

Recurrence needs more than one stratum, so this run ranked on median IDS alone and
said so in its own output rather than implying the primary key had been applied.

**Four programs passed, and the rule cost something real.** Excluded, recorded in
`results/PSA26/programs_dropped_by_size_filter.csv`:

| reason | n genes | first genes |
|---|---|---|
| too small | 2 | HLA-DQA1, HLA-DQB1 |
| too large | 57 | NEXN, SSX2IP, LMNA, TAGLN2, RGS18, … PF4 |
| too large | 92 | RPL22, RPL11, RPS8, RPL5, … (ribosomal) |
| too large | 207 | RBP7, CSF3R, S100A8, S100A9, S100A12, … (myeloid) |

The 92-gene ribosomal and 207-gene myeloid communities are exactly what the cap is
for — compartments, not programs. **The 57-gene one is not.** It contains
TAGLN2/RGS18/PF4, i.e. the platelet program the previous session reported as a
genuine finding, now buried in a 57-gene community and sitting seven genes above the
cap.

The wrong response is to move the cap to 60 after seeing this — that is the
cherry-picking the rule exists to stop, and 60 is where the previous session's
undisclosed cap already was. Two defensible responses:

1. Change the cap in advance for all future runs and re-report everything under it.
2. Better: treat the size as a symptom. Greedy modularity at `min_ids=0.5`
   over-agglomerates; raising `--min-ids` shrinks communities on a principled axis
   (edge strength) instead of trimming them on an arbitrary one (member count).
   Density does not separate these cases cleanly here (0.43 for the 57-gene platelet
   community vs 0.50 for the ribosomal one), so density is not a drop-in substitute.

## The four programs

| rank | program | n | median IDS | density | genes |
|---|---|---|---|---|---|
| 1 | P01_HBB | 3 | 0.902 | 1.00 | HBB, HBA2, HBA1 |
| 2 | P02_HES4 | 7 | 0.485 | 0.43 | HES4, FCGR3A, IFITM3, CDKN1C, TCF7L2, NAP1L1, LINC01272 |
| 3 | P03_FCRL5 | 22 | 0.480 | 0.42 | FCRL5, RALGPS2, IGKC, BANK1, PAX5, MS4A1, TCL1A, IGHD, IGHM, … |
| 4 | P04_PYHIN1 | 25 | 0.431 | 0.32 | PYHIN1, GNLY, CD8A, SPON2, FGFBP2, HOPX, GZMA, CTSW, PRF1, KLRF1, … |

P02 is textbook CD16⁺ monocyte (FCGR3A *is* CD16; HES4, CDKN1C, TCF7L2 are its
canonical partners). P03 is B-cell identity. P04 reproduces the cytotoxic program
from the handoff. **P01 is haemoglobin — ambient RNA, not biology — and the
pre-registered rule ranks it first**, because it is a tight three-gene clique with
the highest median IDS in the run. The rule is honest, not smart: it has no notion
of "interesting". This is handoff open item #3 (a gene-family exclusion flag) and
this run is the argument for building it.

---

## Agenda 3 — cell-type enrichment

`results/PSA26/enrichment_granular.csv` (all 30 subtypes), with
`specificity_summary.csv` on top. Granular only — a 5-main-type version would answer
the question the UMAP already answers. "Expressing" is **`sc.tl.score_genes` > 0**,
fixed in advance: the score is already background-subtracted against an
expression-matched control set, so 0 is the principled null and there is no knob.

| program | top type by proportion | prop | mean | enrichment | note |
|---|---|---|---|---|---|
| P03_FCRL5 | B naive | 1.00 | 3.54 | **7.46** | saturated: 3 types tie at 1.00 |
| P04_PYHIN1 | NK | 1.00 | 4.06 | 2.44 | saturated: 4 types tie at 1.00 |
| P02_HES4 | CD16 Mono | 1.00 | 4.29 | 1.57 | 2nd is Platelet at 0.88 |
| P01_HBB | cDC2 | 0.23 | 0.34 | 1.47 | **disagrees**: by mean score → Platelet |

It works — every identity program lands on the right lineage without ever having
seen a label. But the run also broke the metric the agenda proposed, twice, and both
breaks are now instrumented rather than papered over.

**1. Proportion saturates.** With `score > 0`, a strong program hits `prop = 1.000`
in *every* type that expresses it at all — three types for P03, four for P04. The
sort key then cannot order the head of its own list, and the "most associated" cell
type was whichever one pandas met first. `celltype_enrichment` now breaks the tie by
mean score (which moves P03 from B intermediate to **B naive**, mean 3.54 vs 3.00)
and `specificity_summary` reports `n_tied_at_top_prop` and a `saturated` flag so a
saturated ranking is never quoted as if it were resolved.

**2. Proportion and strength can point at different cell types.** P01_HBB has its
*highest mean score in Eryth* (0.91) while ranking near the **bottom** of the same
list by proportion (0.14) — haemoglobin is enormous in a handful of erythroid cells
and absent in the rest. Sorting on proportion alone gets that program exactly wrong.
Both rankings are now columns (`rank_in_program`, `rank_by_mean_score`) and
`specificity_summary` carries an `agree` flag. **Report proportion for "how
widespread" and mean score for "how strong", and treat disagreement as the finding
it is.**

Two further columns exist for the same reason:

- **Wilson 95% CI** on every proportion. Cell types here span 2 to 4,705 cells; a
  proportion from 3 ASDC cells is not comparable to one from 4,705 monocytes, and
  without the interval a rare subtype tops the list on sampling noise alone. (It did
  in the first run: P01's leader was CD4 Proliferating, n=12.) `specificity_summary`
  therefore restricts its "top" to types with ≥100 cells; the full table keeps
  everything.
- **`enrichment` = prop / prop_overall.** P02_HES4 is "expressed" by 60–77% of every
  cell type in the lane. Its proportion of 1.00 in CD16 Mono looks decisive until you
  see the enrichment is 1.57 against P03's 7.46. Proportion alone cannot tell a
  cell-identity program from a housekeeping one.

**Agenda 3's premise holds.** The granular heatmap resolves what the UMAP cannot:
P04 is visibly a *program*, not a type — NK, NK Proliferating, NK_CD56bright, MAIT,
CD4 CTL and CD8 TEM all at proportion ≈ 1.00, which on a UMAP is one diffuse smear
across the T/NK continuum. Columns are grouped by lineage using the Agenda 1
mapping, so the block structure is readable.

---

## Not committed

`ids_matrix.npy` (364 MB) and `embedded.h5ad` (386 MB) exceed the transfer limit and
stayed in the compute container. Both regenerate from the command at the top in
about 12 minutes on 2 cores.
