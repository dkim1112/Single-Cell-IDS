# Agenda 2 — How to combine results across stratified IDS runs

**Status:** design answer, nothing run yet. Written to be argued with before it is built.

---

## 0. The premise, stated precisely

The starting observation is right but worth sharpening, because the exact reason
matters for what follows.

Harmony does not produce corrected expression. It corrects the **PCA embedding**:
`obsm['X_pca'] -> obsm['X_pca_harmony']`, with `adata.X` untouched
(measured: `max|ΔX| = 0.000e+00`). IDS consumes gene expression. So Harmony is not
merely "not used" by IDS — it is *unusable*, and this was verified rather than
assumed: IDS before vs after Harmony is bit-identical (`r = 1.0000`, `max|ΔIDS| = 0`).

The obvious workaround does not work either. Rebuilding gene-space expression from
the 20 Harmony PCs gives `r = 0.26` against the true IDS — but rebuilding from the
*uncorrected* PCs is equally bad (`r = 0.26`, and the two agree at `r = 0.9995`).
The damage is the 20-PC bottleneck, not Harmony. There is no version of "project it
back" that survives.

The one method that *does* edit expression, ComBat, changes IDS at `r = 0.82`
**even when batch labels are assigned at random** — it manufactures structure. It
is not a safe substitute.

**So: correction is off the table, and stratification is the only batch control
available.** One further trap, already noted in the handoff and easy to fall into:
Harmony-corrected Leiden clusters and the study's `CellType` labels both *span
lanes*. Grouping by them controls cell-type heterogeneity, **not batch**. Only
restricting to a single file does that.

---

## 1. Combine at the gene-pair level, not the module level

Two arguments, one practical and one statistical.

**Practical.** "The same module in two runs" is not well defined. Modules come out
of a hard threshold followed by greedy modularity community detection — both
discontinuous. A single edge crossing `min_ids` can split one community into two or
fuse two into one. Matching module *sets* across 18 runs is then a set-matching
problem with no canonical solution; Jaccard-based matching is a choice, and an
unstable one.

**Statistical.** A gene pair is atomic and identically indexed in every stratum
(given a shared gene list), so `IDS[i,j]` is directly comparable run to run.
Clustering per stratum and then merging means making ~d²/2 threshold decisions
*per stratum* and then a matching decision on top; the errors compound. Combining
first and clustering **once** on the consensus makes those decisions a single time,
on a less noisy input.

**Therefore:** combine the matrices, then run `extract_modules(..., method="greedy")`
exactly once, on the consensus.

---

## 2. What a stratum should be

**Stratum = (file × subject).** All cell types kept.

- **File** = one 10x lane = one machine run. This is the only thing that actually
  removes batch, and it replaces what Harmony would have done.
- **Subject** is not optional. The measured effects say donor biology is *larger*
  than batch: same donor across lanes `r = 0.43`; different donors within one lane
  `r = 0.23`. Pooling donors inside a stratum leaves the bigger confound in place
  while carefully removing the smaller one.
- **Do not stratify by cell type.** This is deliberate and it is the one place the
  intuition points the wrong way. Inside a single T-cell cluster, IDS returned
  mitochondrial, ambient-platelet and myeloid genes — contamination, not T-cell
  biology. When every cell in a stratum is one type, the identity genes stop
  varying and IDS cannot see them. Keeping all cell types is what makes
  cell-identity programs detectable at all.

For the PSA arm this is 2 lanes × 9 subjects = **18 runs** (the HC arm gives another
18 from lanes 01-3/01-4).

### The correction the handoff's plan needs

Those 18 strata are **not 18 independent samples**. PBMC-01-1 and PBMC-01-2 carry
*the same nine subjects* on two lanes — they are 9 donors × 2 technical replicates.
A single "fraction of 18 strata" recurrence conflates technical reproducibility with
biological generalization, and inflates it: a donor-specific artifact that
reproduces on both lanes of one donor contributes 2/18 exactly as two different
donors would.

**Compute recurrence hierarchically instead:**

1. `technical_recurrence(pair, donor)` = fraction of that donor's 2 lanes where
   `IDS ≥ t`. Values 0, 0.5, 1. This measures reproducibility.
2. `biological_recurrence(pair)` = fraction of the **9 donors** with
   `technical_recurrence = 1`. This measures generalization.

Rank on **biological recurrence**; use technical recurrence as a QC read-out. A pair
at biological recurrence 8/9 is a real finding. A pair at 18/18 flat is the same
number reported in a way that hides which of the two things it demonstrates.

---

## 3. The shared gene set — do not take a hard intersection

The handoff's step 2 says "intersect genes passing the floor in every stratum."
That has a bias, not just a cost: the intersection is systematically enriched for
*ubiquitously expressed* genes. Any gene specific to a cell type that happens to be
rare in one donor drops out of all 18 runs because of that one donor. The genes most
likely to be lost are exactly the cell-identity markers the pooled-run design is
built to find.

**Better:**

1. Fix the gene list **once**, from the pooled included cells of both lanes, using
   the usual floor (`--min-cells-frac 0.02`). This is the index for every stratum,
   so all 18 matrices are aligned.
2. Compute IDS on that same list in every stratum.
3. Keep a per-stratum **validity mask**: a pair is valid in a stratum only if both
   genes clear the floor *in that stratum*.
4. Take the median and the recurrence over **valid strata only**, and require a
   minimum number of valid strata (say ≥ 6 of 9 donors) for a pair to be eligible.

This keeps a gene that is well-detected in eight donors and absent in one, which the
hard intersection throws away, while still refusing to score pairs on strata where
the measurement is not there.

---

## 4. Equalize cell counts across strata — this is a bias, not just noise

The handoff notes that the median "weights strata equally regardless of size."
The problem is worse than unequal weighting.

**IDS is a maximum of absolute correlations** over 6×6 feature pairs. The maximum of
several noisy correlations is **biased upward at small n**. A small stratum does not
merely produce a noisier IDS — it produces a *systematically higher* one. Feeding
600-cell and 4,000-cell strata into the same median therefore lets the small strata
push the consensus up, and the pairs that benefit most are the weak ones, which is
precisely where the threshold sits.

(The handoff's own observation is consistent with this: a 482-cell subset agreed
with the full run at only `r = 0.23`.)

**Fix: subsample every stratum to the same n before computing IDS.** IDS is linear
in cells, so this is nearly free. Take n = the smallest stratum that clears the
~1,000-cell floor (drop strata below it rather than letting them vote), and repeat
the draw 3–5 times per stratum, averaging within stratum, to avoid making the answer
depend on one random subsample. Then every estimate carries the same bias and the
same variance, and the median across strata means what it looks like it means.

Report the chosen n and how many strata were dropped.

---

## 5. Thresholds: derive them, don't pick them

The handoff proposes `IDS ≥ threshold` per stratum and recurrence `≥ ~0.7`. Both
numbers are currently guesses. Both can be derived.

**The per-stratum threshold `t`.** Take one stratum, permute the cell order of one
gene, recompute IDS. That gives the distribution of IDS reachable when there is no
dependence, at that n and that gene count. Set `t` at, say, the 99.9th percentile of
that null. (Note the practical blocker: p-values are *not currently supported under
gene blocking* — a shared cell permutation across all blocks is needed first. Open
item #2 in the handoff. That work is a prerequisite for this step, or `t` must come
from the off-diagonal histogram's background/module gap as it does today.)

**The recurrence cutoff.** Once `t` gives a per-stratum false-positive rate `p`, the
number of strata a null pair passes is roughly `Binomial(9, p)`. With `p = 0.01` and
9 donors, `P(≥5) ≈ 1e-8`; against ~4.5×10⁷ gene pairs that is well under one
expected false positive, so **5/9 is already sufficient** and 0.7 (≈7/9) is stricter
than the error budget requires. Being stricter than necessary is not free here: it
preferentially deletes programs that are real but donor-variable, which in a
psoriatic-arthritis cohort may be the interesting ones.

Caveat, stated rather than buried: the binomial assumes independent strata and
independent pairs. Neither holds — genes are correlated and donors share biology.
Treat this as an order-of-magnitude guide for choosing the cutoff, not as an FDR.

**Consensus matrix** = median IDS over valid strata, set to 0 wherever biological
recurrence < the derived cutoff. Then `extract_modules(consensus, method="greedy")`,
once, and rank with the Agenda 1b rule (size filter 3–50, then recurrence, then
median IDS — recurrence is now available, so the primary key is finally usable).

---

## 6. What not to do

| Don't | Why |
|---|---|
| Pool all cells and run IDS once | A two-component mixture creates dependence between *any* two genes whose means differ between components. Lane and donor both act as hidden third variables. This is the failure mode stratification exists to prevent. |
| Run Harmony first | Verified no-op on IDS (`r = 1.0000`). |
| Reconstruct expression from Harmony PCs | `r = 0.26`; the 20-PC bottleneck destroys gene-level structure. |
| Use ComBat | Changes IDS at `r = 0.82` under *random* batch labels — it invents structure. |
| Take the **mean** across strata | One stratum with an artifact drags the mean; the median doesn't. |
| Combine modules instead of pairs | §1. |
| Report recurrence over 18 strata | §2 — conflates technical replication with biological generalization. |

---

## 7. Validation the consensus should carry

1. **Split-half.** Build the consensus from 5 donors and again from the other 4.
   Correlate the two consensus matrices; compare the module overlap. A consensus
   that does not replicate across donor halves is not a consensus.
2. **Negative control.** Shuffle donor labels across strata and rebuild. Recurrence
   should collapse. If it doesn't, the recurrence statistic is measuring something
   structural rather than biological.
3. **The pooled-run contrast.** Count how many programs from a naive pooled run
   survive into the consensus. **The difference is the batch/donor artifact rate** —
   a single number that quantifies whether stratification was worth doing, and the
   most direct empirical answer to the question that motivated this agenda.
4. **IDS vs Pearson (handoff open item #6).** Run this identical
   stratify-and-combine pipeline with Pearson substituted for IDS and compare the
   consensus programs. This is a far more convincing venue for that comparison than a
   single run: a nonlinear program that survives 9 donors and that Pearson never
   finds is the argument for IDS's complexity. If the two consensuses agree, that is
   worth knowing early.

---

## 8. Cost, and how to keep it bounded

18 runs at roughly the cost of the single run already measured
(4,259 cells × 9,542 genes ≈ 10 min on 2 cores), less after subsampling to a common
n — call it a few hours, embarrassingly parallel across strata.

Storage is the real constraint: 18 × 9,542² float32 = **6.6 GB** if every matrix is
held in memory. Two ways out:

- **Two-pass, mirroring the paper's own structure.** Pass 1 streams the strata and
  accumulates only a `uint8` recurrence count per pair (~45 MB for 4.5×10⁷ pairs)
  plus a running sum. Pass 2 revisits only the pairs that cleared recurrence — a
  small fraction — to compute their exact median.
- **Or simply write each stratum's upper triangle to disk as float16** (~182 MB per
  stratum, ~3.3 GB total) and take the median in gene-block chunks. Simpler, and
  3.3 GB on disk is not a problem.

Prefer the second unless memory forces the first; the exact median for every pair is
worth 3.3 GB of disk.

---

## 9. Summary — the recommended procedure

1. Fix one gene list from the pooled included cells of both lanes (`floor 0.02`).
2. Define 18 strata as (lane × subject); drop any below ~1,000 cells.
3. Subsample every stratum to a common n; repeat 3–5× and average within stratum.
4. Compute IDS per stratum on the fixed gene list; keep a per-stratum validity mask.
5. Derive `t` from a permutation null (needs the shared-permutation work first).
6. Per pair: median IDS over valid strata; technical recurrence within donor;
   biological recurrence across the 9 donors.
7. Consensus = median, zeroed where biological recurrence < the derived cutoff
   (expected ≈ 5/9, not 0.7).
8. `extract_modules(consensus, method="greedy")` **once**; rank with the Agenda 1b
   rule, now with recurrence as the primary key.
9. Run the four validations in §7 before believing any of it.

**The three changes to the plan as previously written:** hierarchical recurrence
instead of flat 18-way (§2), a validity mask instead of a hard gene intersection
(§3), and equal-n subsampling to remove the small-stratum upward bias in IDS (§4).
