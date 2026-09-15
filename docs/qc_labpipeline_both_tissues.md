# Lab-pipeline QC applied to both datasets

**Question:** what proportion of the cells marked as passing QC pass under our lab's procedures — in PBMC, and in skin (where having a cell-type annotation is the QC pass mark).

Run on **all 28 PBMC lanes** and **all 33 skin samples**, not a single sample. Lab thresholds from `config/config_sample.json`: `min_gene_per_cell 200`, `max_mt_percentage 10.0`, `max_doublet_score 0.15`.

---

## Headline

| | study kept | pass lab QC | |
|---|---|---|---|
| **PBMC** (GSE194315) | 282,705 | **246,474** | **87.2%** |
| **skin** (GSE173706) | 67,378 | **67,367** | **99.98%** |

---

## PBMC — 87.2%, and every failure is mitochondrial

| lab filter | passes | % of study-kept |
|---|---|---|
| min_gene ≥ 200 | 282,705 | **100.0%** |
| demuxlet singlet | 282,705 | **100.0%** |
| pct_mt ≤ 10 | 246,474 | **87.2%** |
| **all three** | **246,474** | **87.2%** |

- Of the 36,231 cells we would drop, **36,231 fail on mitochondrial percentage and nothing else.** Zero fail on gene count, zero on doublets.
- The study's own thresholds are recoverable from the data: `n_genes ≥ 200` exactly, `total_counts ≥ 500` exactly, and **no mitochondrial cut at all** — the highest retained cell is 87.7% mitochondrial, and 446 retained cells are over 50%.
- So our pipeline is not stricter in three ways. It is stricter in exactly one, and the study is effectively unfiltered there.
- The doublet criterion is redundant here: the study already required a demuxlet singlet, so every kept cell is SNG.

**Per-lane spread is wide — 70.6% to 99.2%:**

| worst | | best | |
|---|---|---|---|
| PBMC-01-2 | 70.6% | PBMC-03-1 | 99.2% |
| PBMC-01-1 | 71.5% | PBMC-03-2 | 98.8% |
| PBMC-01-4 | 75.2% | PBMC-03-3 | 98.8% |

Median mitochondrial % tracks this exactly — 8.4% in the PBMC-01 lanes vs 3.3–3.7% in PBMC-03. This is a batch-level property of the run, not scattered cell damage.

> **Correction to my earlier answer.** I previously reported 52.0% lab-keep and 76.1% agreement. That was **PBMC-01-1 alone**, which turns out to be the second-worst of 28 lanes. The cohort-wide figure is 87.2%.

---

## Skin — 99.98%, because the study already used our mitochondrial threshold

| lab filter | passes | % of annotated cells |
|---|---|---|
| pct_mt ≤ 10 | 67,378 | **100.0%** |
| min_gene ≥ 200 | 67,367 | 99.98% |
| **both** | **67,367** | **99.98%** |

Only **11 cells** in the whole dataset would be dropped, all on gene count (minimum 101 genes).

The reason for the near-perfect agreement is worth stating plainly: **the maximum mitochondrial percentage among annotated skin cells is exactly 10.00%.** That is not a coincidence — GSE173706 applied a ≤10% cut, the same one we use. Their minimum total counts is 501, so they used a ~500 count floor too. The only thing they did not impose is a 200-gene floor.

Skin is also just cleaner: median mitochondrial 3.4% vs 6.1% in PBMC (p99 9.5% vs 29.4%).

**Doublets are the open gap.** Skin is not multiplexed, so there is no demuxlet call and `run_Doublet.py` has never been run on it. The 99.98% figure covers two of the three lab filters, not three. This is the one place where our QC is genuinely untested on skin.

---

## The cells each study threw away

| | study dropped | would pass our genes+mito | |
|---|---|---|---|
| PBMC | 229,015 | 57,661 | 25.2% |
| skin | 28,710 | 15,800 | 55.0% |

Skin's discards are the more interesting half: **55% of the cells GSE173706 removed look fine on count and mitochondrial grounds.** Their filtering removed those cells for a reason the deposited metadata does not expose — most likely doublet or cluster-level QC during annotation. We are inheriting that decision without being able to audit it.

---

## What this changes

1. **The mitochondrial threshold is the only real disagreement, and it only exists in PBMC.** Skin already agrees with us. If we adopt 10% for both tissues, skin is unchanged and PBMC loses 12.8% of its cells — concentrated in the PBMC-01 and PBMC-07 lanes.
2. **This cuts against a batch confound rather than creating one.** The lanes we would prune hardest are the same lanes with the highest median mitochondrial content, which is exactly the technical variation we do not want IDS reading as biology.
3. **Doublets remain unaddressed in skin**, and are the only lab filter we cannot currently evaluate there.

---

## Files

- `results/qc-compare/labqc_vs_study_both_tissues.csv` — the headline table
- `results/qc-compare/pbmc_percell_qc.csv.gz` — 511,726 cells, per-cell n_genes / counts / pct_mt / pass flags
- `results/qc-compare/skin_percell_qc.csv.gz` — 96,088 cells, same columns
- `results/qc-compare/pbmc_qc_by_lane.csv`, `skin_qc_by_sample.csv` — per-batch breakdown
- `scripts/qc_pbmc_percell.sh`, `scripts/qc_skin_percell.py`, `scripts/qc_compare_summarize.py` — reproduces all of it from the GEO files

*Minor correction: the raw skin cell count is **96,088**, not the 96,057 I reported earlier. Verified two ways — summing the 33 per-sample column counts, and counting header fields directly.*
