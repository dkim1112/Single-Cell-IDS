"""
variant_blood_tests.py — do programs found only under a DIFFERENT discovery choice differ in blood?

Takes every program that a rerun found and October did not (other random samples, 8,000 cells,
total normalization, threshold +/- 0.05; discovery_variant_programs.csv), scores them in all 28
lanes and runs the valid blood comparisons. No permutation or activity filter is applied, which
makes this more permissive than the catalogue; technical gene sets are flagged with the catalogue's rules.

Writes variant_blood_tests_summary.csv and variant_blood_tests_detail.csv.gz
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "03_catalogue"))
from sc_ids import paths                                                      # noqa: E402
import disease_tests as T                                                     # noqa: E402
from program_catalogue import classify                                        # noqa: E402

OUT = T.OUT
V = pd.read_csv(f"{OUT}/discovery_variant_programs.csv")
gs = V.genes.map(lambda g: [x.strip() for x in g.split(",")])
V["technical"] = [bool(classify(g)[0]) or (sum(x.startswith(("RPL", "RPS", "MT-", "MTRNR")) for x in g) / len(g) >= 0.25)
                  or ("XIST" in g and len(g) <= 6) for g in gs]
V.to_csv(f"{OUT}/discovery_variant_programs.csv", index=False)
S = pd.read_csv(f"{OUT}/scores_blood_all28_variants.csv.gz", dtype={"donor": str})
plan = [("more lanes: all 16 lanes of runs 02-05", T.lanes(2, 3, 4, 5), "Healthy", "PSO"),
        ("more lanes: runs 06 + 07 (8 lanes)", T.lanes(6, 7), "PSO", "PSA")]
summ, det = [], []
for which, ids in (("all 268 variant programs", list(V.program_id)), ("non-technical only", list(V[~V.technical].program_id))):
    for level in ("lineage", "subtype"):
        for label, keep, a_, b_ in plan:
            R = T.contrast(T.donor_table(S, ids, level, keep), ids, a_, b_)
            summ.append(dict(programs=which, n_programs=len(ids), cell_type_level=level, comparison=f"{a_} vs {b_}", lanes=label, **T.summarise(R, False)))
            if which.startswith("all"):
                R.insert(0, "comparison", f"{a_} vs {b_}"); R.insert(0, "cell_type_level", level); det.append(R)
pd.DataFrame(summ).to_csv(f"{OUT}/variant_blood_tests_summary.csv", index=False)
D = pd.concat(det, ignore_index=True).merge(V[["program_id", "variants", "technical", "genes"]], on="program_id")
D.round(6).to_csv(f"{OUT}/variant_blood_tests_detail.csv.gz", index=False, compression="gzip")
pd.set_option("display.width", 250)
print(pd.DataFrame(summ)[["programs", "n_programs", "cell_type_level", "comparison", "donors_a", "donors_b", "tests", "raw_p_below_0_05", "expected_by_chance", "significant", "significant_and_change_0_2", "largest_abs_change"]].to_string(index=False))
sig = D[D.p_adj < 0.05].sort_values("p")
print(sig[["cell_type_level", "comparison", "program_id", "group", "n_a", "n_b", "mean_a", "mean_b", "delta", "p_adj", "technical", "variants", "genes"]].head(20).to_string(index=False))
