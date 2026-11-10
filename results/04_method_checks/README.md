# Method checks (November 10, 2026)

**Question.** If we choose differently (more donors, finer cell types, a different random sample, a different
normalization, more shuffles), do disease differences become detectable, especially in blood? And are the 57 programs stable?

**Answer.** No choice made a blood difference detectable. The blood result held in every valid comparison.
The program list is partly choice-dependent: 43 of the 57 programs come back in every other random sample, and 22 of 57
come back under every change tried.

Everything below was produced by `scripts/04_method_checks/`, on the per-lane data files from `scripts/02_ids_programs/build_units.py`. The 4,000-cell cap was not changed.

## First: the rebuild is exact

| Check | Result |
|---|---|
| Cells, nine October lanes | 74,126, identical to October cell for cell and in the same order |
| Cells, skin | 67,367, identical to October |
| Lane-by-lane scoring vs scanpy on a merged object | largest difference 7e-07 |
| October donor means (46,314 blood rows) | reproduced to the 5th decimal |
| October test counts | reproduced exactly: blood 0, 0, 0 and 86; skin 103, 274 and 0 |
| October detections | all 194 reproduced gene for gene |
| Edge thresholds (never recorded before) | whole-tissue 0.50, T/NK 0.35, myeloid 0.40, keratinocyte 0.50 |
| Peak memory (measured) | build 1.5 GB; scoring all 28 lanes 0.9 GB; discovery 2.3 GB at 8,000 cells; permutation 0.5 GB |

## Test 1: more donors (all 28 lanes built, 191,414 blood cells)

| Comparison | October | With more lanes |
|---|---|---|
| Healthy vs PSO, same lanes | 16 vs 8 donors: 0 of 372 | 20 vs 17 donors (runs 02-05): 0 of 372 |
| PSO vs PsA, same lanes | 6 vs 12 donors: 0 of 372 | 7 vs 17 donors (runs 06+07): 0 of 372 |
| Healthy vs PsA | different lanes: 86 of 372, not interpretable | run 01 only, same run but different lanes, 9 vs 9: 2 of 372, both under 0.02 score units |

Largest change seen in a valid comparison: 0.22 score units (healthy vs PSO), 0.18 (PSO vs PsA). Smallest adjusted p: 0.31 and 0.81.

**Comparing across runs manufactures differences.** Same disease, different runs (`run_effect_controls.csv`):

| Same-disease control | Donors | Significant of 372 |
|---|---|---|
| Healthy: run 01 vs runs 02-05 | 9 vs 20 | 162 |
| PsA: run 01 vs run 06 | 9 vs 12 | 163 |
| PSO: runs 02-05 vs run 06 | 17 vs 6 | 54 |
| PsA: run 06 vs run 07 | 12 vs 5 | 62 |
| Healthy: runs 02+04 vs run 05 | 11 vs 9 | 66 |

So pooling all 28 lanes is not a way to add power: PSO vs PsA across all lanes gives 29 "significant" and healthy vs PsA 71, and
these are run differences. It also shows the tests are not blind: they detect run differences of 0.3 to 0.9 score units easily.
About three quarters of the blood tests could have detected a change of 0.2 (80% power, Bonferroni level; normal approximation).

## Test 2: finer cell types (source subtypes instead of 13 lineages)

| Comparison | Tests | Significant |
|---|---|---|
| Healthy vs PSO, 20 vs 17 donors | 837 | 0 (smallest adjusted p 0.15) |
| PSO vs PsA, 7 vs 17 donors | 651 | 0 |
| Skin, normal vs lesional | 1,023 | 361 (lineage level: 274 of 651) |
| Skin, normal vs non-lesional | 1,023 | 82 (lineage level: 103 of 651) |

Skin: 20 distinct programs change at either level. GP47 and GP57 appear only at subtype level, GP11 and GP12 only at lineage level.

## Test 3: discovery consistency (one change at a time, nine October lanes)

| Change | Programs found | Of the 57 distinct programs, recovered | Not in October |
|---|---|---|---|
| None (seed 0) | 194 | 57 | 0 |
| Another random 4,000 cells (seeds 1 to 4) | 198 to 203 | 47 to 50 | 22 to 24 |
| Library size over all genes | 185 | 43 | 53 |
| 8,000 cells instead of 4,000 | 192 | 47 | 20 |
| Edge threshold -0.05 | 155 | 37 | 48 |
| Edge threshold +0.05 | 214 | 37 | 76 |

Recovered = some program in the rerun overlaps a member detection by at least half of the genes (the catalogue's merge rule).

- 43 of 57 are recovered in all four other samples; 4 in three, 4 in two, 2 in one, 4 in none (`discovery_program_stability.csv`).
- All 16 programs detected more than once in October are recovered in every sample. The 14 unstable ones were each detected once.
- Recovered under every change: GP01 (AP-1), GP22 (antimicrobial), GP24 (MIP-1), GP26, GP29, GP33 (cytotoxic), GP46 (interferon), GP48 (KRT6/16).
- Not recovered under at least one change: GP32 (IFIT1/2/3; threshold -0.05), GP45 (interferon; +0.05), GP20 (neutrophil chemokines; both
  thresholds), GP53 (desmosome; total normalization, +0.05), GP52 (IL-36; 2 of 4 samples), GP39 (1 of 4 samples).
- Doubling the cap to 8,000 changes the list about as much as redrawing 4,000 does.

**Do the extra programs differ in blood?** The 268 programs found only under a change were scored in all 28 lanes:
0 significant, for healthy vs PSO and PSO vs PsA, at lineage level (1,072 tests each) and subtype level (2,412 and 1,876).

## Test 4: more shuffles in the permutation test

| Rule | Pass | Kept detections | Distinct programs |
|---|---|---|---|
| October: 60 shuffles, 90% of pairs beat all 60 | 90 | 88 | 57 |
| 200 shuffles, 90% of pairs with p < 0.01 | 91 | 88 | 56 |
| 200 shuffles, 90% of pairs beat all 200 | 85 | 84 | 53 |
| 200 shuffles, every pair with p < 0.01 | 75 | 73 | 45 |

One three-gene set leaves (GP43: FLG, TPSAB1, TPSB2) and one enters (IFIT1, IFIT3, VPS33B, which merges into an existing
interferon program). Both come from the 901-cell skin myeloid search. With 200 shuffles the filter supports p < 0.01.

## Limits

- Healthy vs PSO pools four runs (02 to 05) without stratifying by run. The run mix differs (run 04 has 4 healthy and 9 PSO donors).
- Discovery was varied one change at a time and on the nine October lanes only.
- "Could have detected 0.2" is an approximation from the between-donor spread, not an exact power calculation.
- The 268 extra programs were not put through the permutation or activity filter.
- Healthy vs PsA still has no valid design in this dataset.

## Files

| File | What it holds |
|---|---|
| `tests_summary.csv`, `tests_detail.csv.gz` | every comparison at both cell-type levels (tests 1 and 2) |
| `run_effect_controls.csv` | same disease, different runs |
| `discovery_variants_summary.csv`, `discovery_program_stability.csv` | test 3 |
| `discovery_variant_programs.csv`, `variant_blood_tests_*.csv` | the 268 extra programs and their blood tests |
| `permutation_200_shuffles.csv`, `permutation_chain.csv` | test 4 |
| `scores_*.csv.gz` | donor x subtype mean scores (gitignored; rebuilt in about 2 minutes) |

## Rerun

```bash
export SC_IDS_WORK=<a local scratch folder>          # default: results/00_units (gitignored)
python scripts/02_ids_programs/build_units.py --tissue blood && python scripts/02_ids_programs/build_units.py --tissue blood --index
python scripts/02_ids_programs/build_units.py --tissue skin  && python scripts/02_ids_programs/build_units.py --tissue skin --index
python scripts/04_method_checks/validate_scoring.py --units PBMC-03-1,PBMC-03-3
for s in "blood nine" "blood all28" "skin all"; do set -- $s; python scripts/04_method_checks/score_units.py --tissue $1 --set $2; done
python scripts/04_method_checks/disease_tests.py
for a in whole-tissue T+NK Myeloid Keratinocyte; do
  python scripts/04_method_checks/discovery_runs.py --analysis $a --seeds 0,1,2,3,4
  python scripts/04_method_checks/discovery_runs.py --analysis $a --norm total
  python scripts/04_method_checks/discovery_runs.py --analysis $a --cap 8000
done
python scripts/04_method_checks/discovery_compare.py
python scripts/04_method_checks/score_units.py --tissue blood --set all28 --programs results/04_method_checks/discovery_variant_programs.csv --tag variants
python scripts/04_method_checks/variant_blood_tests.py
python scripts/04_method_checks/permutation_runs.py --tests 200 && python scripts/04_method_checks/permutation_runs.py --tests 200 --chain
```

Each script stops after about 140 seconds and resumes where it left off, so re-run a line until it prints `done`.
Needs: scanpy, anndata, networkx, tqdm; `fast_matrix_market` is optional (faster lane reading).
