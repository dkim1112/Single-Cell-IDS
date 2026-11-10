"""
permutation_runs.py — the authors' permutation test with more shuffles.

Same procedure as scripts/03_catalogue/program_catalogue.py (authors' ids.compute_IDS(p_val=True) on
the union of each search's program genes, up to 4,000 cells drawn with the same seed, library
size over all genes, [0,8] scaling), with the number of shuffles as a parameter. October used 60;
the authors' default is 100.

  python scripts/04_method_checks/permutation_runs.py --tests 200      # checkpointed per search
  python scripts/04_method_checks/permutation_runs.py --tests 200 --chain

--chain recomputes the selection chain (technical -> permutation -> activity -> merge) under each rule.
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))
import ids                                                                    # noqa: E402
from sc_ids import paths, units                                               # noqa: E402
import discovery_runs as D                                                    # noqa: E402

OUT = os.path.join(paths.RESULTS, "04_method_checks")


def run(tests, budget):
    t0 = time.time()
    part = os.path.join(units.WORK, f"perm{tests}"); os.makedirs(part, exist_ok=True)
    det = pd.read_csv(os.path.join(paths.CATALOGUE, "detections_status.csv"))
    for (an, st), rs in det.groupby(["analysis", "stratum"], sort=True):
        f = os.path.join(part, f"{an.replace('+','-')}__{st}.csv")
        if os.path.exists(f):
            continue
        if time.time() - t0 > budget:
            print("[budget] stopping; re-run to continue"); return False
        t = time.time()
        tissue, c = D.stratum_cells(an, st)
        pos = np.arange(len(c))
        if len(pos) > 4000:
            pos = np.random.default_rng(0).choice(pos, 4000, replace=False)
        sets = {r.program_id: [g.strip() for g in r.genes.split(",")] for r in rs.itertuples()}
        allg = set(units.genes(tissue))
        genes = sorted({g for v in sets.values() for g in v if g in allg})
        X = np.asarray(units.read_cells(tissue, c.iloc[pos], genes, norm="total").todense(), dtype=np.float64)
        lo, span = X.min(0, keepdims=True), X.max(0, keepdims=True) - X.min(0, keepdims=True)
        X = (X - lo) / np.where(span == 0, 1, span) * 8.0
        C, Pv = ids.compute_IDS(X, num_terms=6, p_norm="max", p_val=True, num_tests=tests, bandwidth_term=0.5)
        gi = {g: i for i, g in enumerate(genes)}
        rows = []
        for pid, gl in sets.items():
            ix = [gi[g] for g in gl if g in gi]
            if len(ix) < 2:
                continue
            v = Pv[np.ix_(ix, ix)][np.triu_indices(len(ix), k=1)]
            rows.append(dict(program_id=pid, n_pairs=len(v), frac_beat_all=float((v == 0).mean()),
                             frac_p_below_0_01=float((v < 0.01).mean()), frac_p_below_0_05=float((v < 0.05).mean()),
                             p_median=float(np.median(v)), p_max=float(v.max())))
        pd.DataFrame(rows).to_csv(f, index=False)
        print(f"{an} {st}: {X.shape[0]} cells x {len(genes)} genes, {len(rows)} detections ({time.time()-t:.0f}s)", flush=True)
    return True


def chain(tests):
    det = pd.read_csv(os.path.join(paths.CATALOGUE, "detections_status.csv"))
    part = os.path.join(units.WORK, f"perm{tests}")
    P = pd.concat([pd.read_csv(os.path.join(part, f)) for f in sorted(os.listdir(part))], ignore_index=True)
    P.to_csv(f"{OUT}/permutation_{tests}_shuffles.csv", index=False)
    d = det.merge(P, on="program_id", how="left")
    sets = d.genes.map(lambda g: frozenset(x.strip() for x in g.split(",")))

    def distinct(keep):
        idx = list(np.where(keep)[0]); par = {i: i for i in idx}
        def f(i):
            while par[i] != i:
                par[i] = par[par[i]]; i = par[i]
            return i
        for a_, i in enumerate(idx):
            for j in idx[a_ + 1:]:
                if len(sets[i] & sets[j]) / len(sets[i] | sets[j]) >= 0.5:
                    par[f(i)] = f(j)
        return len({f(i) for i in idx})

    nontech = ~d.technical.astype(bool)
    rules = [("October: 60 shuffles, >= 90% of pairs beat all 60", d["frac_pairs_beat_all_shuffles"] >= 0.9),
             (f"{tests} shuffles, >= 90% of pairs with p < 0.01", d["frac_p_below_0_01"] >= 0.9),
             (f"{tests} shuffles, >= 90% of pairs beat all {tests}", d["frac_beat_all"] >= 0.9),
             (f"{tests} shuffles, every pair with p < 0.01", d["frac_p_below_0_01"] >= 1.0),
             (f"{tests} shuffles, >= 90% of pairs with p < 0.05", d["frac_p_below_0_05"] >= 0.9)]
    base = nontech & (d["frac_pairs_beat_all_shuffles"] >= 0.9) & d.active.astype(bool)
    rows = []
    for name, ok in rules:
        keep = nontech & ok & d.active.astype(bool)
        rows.append(dict(rule=name, non_technical=int(nontech.sum()), pass_permutation=int((nontech & ok).sum()),
                         kept_detections=int(keep.sum()), distinct_programs=distinct(keep.values),
                         kept_in_october_but_not_here=int((base & ~keep).sum()),
                         kept_here_but_not_in_october=int((keep & ~base).sum())))
    R = pd.DataFrame(rows); R.to_csv(f"{OUT}/permutation_chain.csv", index=False)
    print(R.to_string(index=False))
    x = d[nontech]
    print("non-technical detections, agreement of the 90% rule, 60 vs", tests, "shuffles (p<0.01):")
    print(pd.crosstab(x["frac_pairs_beat_all_shuffles"] >= 0.9, x["frac_p_below_0_01"] >= 0.9).to_string())
    print("correlation of the per-detection fractions (60 beat-all vs %d beat-all): %.3f" % (tests, np.corrcoef(x["frac_pairs_beat_all_shuffles"], x["frac_beat_all"])[0, 1]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tests", type=int, default=200); ap.add_argument("--budget", type=float, default=140)
    ap.add_argument("--chain", action="store_true")
    a = ap.parse_args()
    if a.chain:
        chain(a.tests)
    elif run(a.tests, a.budget):
        print("done")
