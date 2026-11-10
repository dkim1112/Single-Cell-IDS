"""
discovery_compare.py — how much do the programs depend on discovery choices?

Compares the October detections with reruns from discovery_runs.py that change ONE thing:
  seed 1-4            a different random 4,000 cells
  total normalization library size over all genes instead of the panel genes
  cap 8000            twice as many cells per stratum
  threshold -/+ 0.05  the edge threshold moved by 0.05

"Recovered" uses the catalogue's own merge rule: gene sets overlapping by at least half (Jaccard >= 0.5).
  detection level   the rerun of the SAME search (analysis x stratum) contains a matching program
  program level     ANY search in the rerun contains a program matching any member detection

Writes discovery_variants_summary.csv, discovery_program_stability.csv, discovery_novel_programs.csv
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))
from sc_ids import paths                                                      # noqa: E402
import discovery_runs as D                                                    # noqa: E402

OUT = os.path.join(paths.RESULTS, "04_method_checks")
VARIANTS = [("October rerun (seed 0)", dict(seed=0)),
            *[(f"seed {k}", dict(seed=k)) for k in (1, 2, 3, 4)],
            ("total normalization", dict(norm="total")),
            ("cap 8000", dict(cap=8000)),
            ("threshold -0.05", dict(dthr=-0.05)), ("threshold +0.05", dict(dthr=0.05))]


def jac(a, b):
    return len(a & b) / len(a | b)


def variant_programs(seed=0, norm="panel", cap=4000, dthr=0.0):
    rows = []
    for an in D.RUNS:
        for s in D.strata_of(an):
            C, info = D.ids_matrix(an, s, seed, norm, cap=cap)
            if C is None:
                continue
            t = D.programs_from(C, info["genes"], D.MIN_IDS[an] + dthr)
            for g in t.genes:
                rows.append((an, s, frozenset(g.split(", "))))
    return rows


def main():
    det = pd.read_csv(os.path.join(paths.CATALOGUE, "detections_status.csv"))
    det["set"] = det.genes.map(lambda g: frozenset(x.strip() for x in g.split(",")))
    kept = det[det.status == "kept"]
    gps = sorted(kept.program_no.unique())
    summ, stab = [], pd.DataFrame(index=gps)
    novel = []
    for name, kw in VARIANTS:
        V = variant_programs(**kw)
        by = {}
        for an, s, g in V:
            by.setdefault((an, s), []).append(g)
        best = lambda r: max([jac(r.set, g) for g in by.get((r.analysis, r.stratum), [])] or [0.0])
        b_all, b_kept = det.apply(best, axis=1), kept.apply(best, axis=1)
        allsets = [g for _, _, g in V]
        gp_ok = {gp: any(jac(m, g) >= 0.5 for m in kept[kept.program_no == gp].set for g in allsets) for gp in gps}
        stab[name] = pd.Series(gp_ok)
        nov = [(an, s, g) for an, s, g in V if max(jac(g, o) for o in det.set) < 0.5]
        novel += [(name, an, s, len(g), ", ".join(sorted(g))) for an, s, g in nov]
        summ.append(dict(variant=name, programs_found=len(V),
                         identical_gene_sets=int(sum(any(g == r.set for g in by.get((r.analysis, r.stratum), [])) for r in det.itertuples())),
                         detections_recovered_of_194=int((b_all >= 0.5).sum()),
                         kept_detections_recovered_of_88=int((b_kept >= 0.5).sum()),
                         distinct_programs_recovered_of_57=int(sum(gp_ok.values())),
                         programs_not_in_october=len(nov)))
        print(summ[-1], flush=True)
    S = pd.DataFrame(summ); S.to_csv(f"{OUT}/discovery_variants_summary.csv", index=False)
    seeds = [f"seed {k}" for k in (1, 2, 3, 4)]
    stab.insert(0, "other_seeds_recovered_of_4", stab[seeds].sum(axis=1))
    dp = pd.read_csv(os.path.join(paths.CATALOGUE, "distinct_programs.csv")).set_index("program_no")
    stab.insert(0, "category", dp.category.reindex(stab.index)); stab.insert(0, "label", dp.label_rule.reindex(stab.index))
    stab.insert(0, "n_genes", dp.n_genes.reindex(stab.index)); stab.insert(0, "n_detections", dp.n_detections.reindex(stab.index))
    stab.index.name = "program_no"; stab.to_csv(f"{OUT}/discovery_program_stability.csv")
    pd.DataFrame(novel, columns=["variant", "analysis", "stratum", "n_genes", "genes"]).to_csv(f"{OUT}/discovery_novel_programs.csv", index=False)
    print(stab.other_seeds_recovered_of_4.value_counts().sort_index(ascending=False).to_dict())


if __name__ == "__main__":
    main()
