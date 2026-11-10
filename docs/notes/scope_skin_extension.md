# Extending to skin: what "the same thing" actually covers

Written before touching a skin dataset, so the scope is agreed rather than discovered.

**Assumption made in the absence of an answer: the skin and blood cohorts are
UNPAIRED (different patients).** If they are paired, section 5 changes and
directionality becomes testable. This is the first question for Matthew.

---

## 1. What the PBMC work actually is, step by step

| # | step | reusable on skin? | what has to change |
|---|---|---|---|
| 1 | load 10x, join per-cell metadata, apply the study's QC gate | **no** | barcode key, QC column, doublet handling are all study-specific. Now isolated in `sc_ids/datasets/` |
| 2 | collapse granular labels → `CellTypeMain` | machinery yes, **mapping no** | skin has keratinocytes, fibroblasts, endothelium, melanocytes. Needs a new mapping file **from the lab** |
| 3 | pick a stratum, run IDS inside it | **no** | (lane × subject) assumes demuxlet-pooled lanes. If each skin file is one donor, there is no subject axis |
| 4 | threshold the IDS matrix, extract modules | code yes, **thresholds no** | `min_ids=0.5` came from PBMC's own IDS histogram. Re-derive from skin's |
| 5 | rank programs (3–50 genes, median IDS) | **yes, as is** | — |
| 6 | cell-type enrichment (proportion expressing) | **yes, as is** | — |

So "the same thing" reuses steps 4–6 and rebuilds 1–3. That part is genuinely
straightforward, and the `datasets/` split now means step 1 is one new file.

**But steps 1–6 produce two independent program lists. That is not a shared-program
result, and the abstract needs one.**

---

## 2. What "the same thing" does not cover

Six things stand between two program lists and the abstract's claim. Four are
unbuilt; one is unanswerable as currently designed.

| gap | status | why it matters |
|---|---|---|
| **Ambient-RNA / gene-family exclusion** | not built | PBMC ranked haemoglobin (HBB/HBA1/HBA2) **first**. Skin dissociation is worse — expect ambient keratin (KRT\*) and collagen (COL\*). Ambient programs will appear in *both* tissues and register as "shared", which is the worst possible false positive for this paper |
| **Consensus across strata (Agenda 2)** | designed, not built | Single-run program lists are unstable. An overlap computed between two unstable lists is doubly unstable. Also turns `recurrence` on, which is the primary ranking key |
| **Cross-tissue program matching** | does not exist | There is currently no definition of "the same program in skin and blood", and no null to test overlap against |
| **The shared-cell-type confound** | not addressed | See §3. This is the one that can invalidate the result |
| **Directionality** | not testable if unpaired | See §5 |
| **IDS vs Pearson** | never run | The abstract's novelty is the method. A reviewer will ask what IDS found that correlation did not |

Rough split: "the same thing" is about 40% of what the abstract needs.

---

## 3. The confound that decides whether this works

**Blood and skin both contain immune cells.** Run the pipeline on skin and you will
get a B-cell identity program (MS4A1, CD79A, PAX5) much like PBMC's P03_FCRL5, and a
myeloid one, and a T-cell one. Their gene overlap with the blood programs will be
high.

Concluding "skin and blood share a B-cell program" would be **true and worthless**.
It says both tissues contain B cells. It says nothing about immune signalling
crossing between them.

The interesting claim needs programs shared **beyond what cell-type composition
explains** — a program that is a cell *state* rather than a cell *identity*. An
interferon-response program appearing in skin keratinocytes *and* blood monocytes
would qualify. A B-cell identity program appearing in both B-cell populations would not.

**The discriminator already exists in what we built.** `specificity_summary.csv`
separates these:

- **Identity program**: `specific = True` with high `top_enrichment`. P03_FCRL5 sits
  at 7.46, confined to B cells. Trivially shared across tissues.
- **State / cross-type program**: high recurrence but low `top_enrichment`, several
  cell types tied at the top. P04_PYHIN1 is halfway there — enrichment 2.44, four
  types at proportion 1.00, spanning NK, MAIT, CD4 CTL and CD8 TEM.

So Agenda 3 was not a side quest. It is the filter that decides which shared
programs are worth putting in the abstract.

**Practical consequence:** the headline candidates are programs that are recurrent
across strata, shared across tissues, and **not** cell-type-specific. Plan the
comparison around that column, not around the top of the ranked list.

A second, smaller version of the same problem: mitochondrial and ribosomal programs
appear in every dataset ever sequenced. They will "share" perfectly. That is what
the gene-family exclusion in §2 is for, and it has to be decided **before** looking
at the overlap.

---

## 4. Defining "the same program" across tissues

Nothing here exists yet; these are the decisions to make.

1. **Matching rule.** Gene-set overlap (Jaccard, or hypergeometric on the shared
   gene universe) between a skin program and a blood program. Must be computed on
   the **intersection of genes measured in both**, not on all genes — otherwise
   tissue-specific genes inflate the denominator.
2. **The null.** Overlap of what size is surprising? Permute program membership
   within each tissue and rebuild the overlap distribution. Without this, "these two
   programs share 8 genes" is not a result.
3. **Granularity.** Programs are not atomic (a 57-gene community is a core plus a
   fringe — see `docs/agenda1_and_3_results.md`). Matching at the **gene-pair** level
   rather than the program level sidesteps this, exactly as Agenda 2 argues for
   combining strata. Prefer that if the compute allows.
4. **Cell-type control.** Report the overlap twice: once raw, once after excluding
   programs flagged `specific = True`. The difference is the answer to §3.

---

## 5. Directionality — be honest about this early

Under the unpaired assumption, **directionality is not recoverable.** Cross-sectional
expression from two separate cohorts contains no information about which tissue's
program preceded the other's. No statistical treatment of that data will produce a
direction; anything that appears to is an artifact of the modelling choice.

Three ways to earn a directional claim, in descending order of feasibility:

- **Paired samples** (same patient, skin + blood). Then a shared program's score can
  be correlated *within* patient across tissues, and its association with
  skin-side severity (PASI, lesional vs non-lesional) versus systemic outcome
  (arthritis status) can be compared. Still not causal, but it is evidence.
- **An external anchor**: disease duration, treatment response, or skin-only vs
  skin+joint patients. Direction becomes an ordering over patient groups.
- **Reframe.** Claim shared programs; state direction as the hypothesis the work
  motivates. Weaker abstract, but defensible, and it is what the data supports.

Deciding this now matters because it changes what data you need to ask for.

---

## 6. Recommended order (no deadline, so dependency order)

1. **Ambient / gene-family exclusion.** Cheap, fixes a known PBMC defect, and is a
   prerequisite for any cross-tissue comparison. Decide the exclusion list in advance.
2. **Agenda 2 consensus on PBMC.** Gives a stable blood program list and turns
   recurrence on as the primary ranking key.
3. **IDS vs Pearson on that consensus.** Answers the novelty question while there is
   still only one dataset to reason about. Can run in parallel with 4.
4. **Skin.** New `datasets/skin.py`, mapping file from the lab, stratum definition,
   thresholds re-derived from skin's own IDS histogram. Then steps 1–6 of §1.
5. **Cross-tissue matching** (§4), reported with and without the identity control.
6. **Directionality**, only if paired data exists (§5).

Steps 1–3 are all on data already in hand and none of them are blocked.

---

## 7. Open questions for Matthew

1. **Are the skin and blood cohorts paired?** Decides whether directionality is a
   result or a hypothesis (§5).
2. **Which skin dataset**, and does it publish a per-cell QC verdict and a cell-type
   annotation, like GSE194315's `IncludedInStudy` and `CellType`?
3. **Is there a skin granular→main cell-type mapping?** Same two-column format as
   `celltype_mapping.txt`. We should not invent one.
4. **Is the abstract's claim "shared programs" or "skin → systemic direction"?**
   The second one requires data we may not have.
5. Also worth passing on: the mapping file already sent has a typo —
   `NK_CD56bright` and `NK Proliferating` map to `"K Cell"` rather than `"NK Cell"`.
