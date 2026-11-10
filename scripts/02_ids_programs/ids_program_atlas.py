"""
ids_program_atlas.py — where is each IDS gene program expressed, and which
programs are shared between blood and skin?

Takes the per-stratum programs from ids_programs_disease.py and answers the
three questions the abstract needs:

  1. WHICH CELL TYPES. Every cell in both tissues is scored for every program,
     then summarised by coarse lineage. Scoring uses `sc.tl.score_genes` on the
     log1p representation -- NOT the [0,8] IDS representation -- because
     score_genes subtracts an expression-matched control set and that
     subtraction is meaningless after min-max scaling. Programs are DISCOVERED
     in IDS space and SCORED in log space, on purpose.

  2. WHICH DISEASE STATES. Scores are averaged within (donor x lineage) FIRST,
     and the disease comparison is made across donors. The unit of replication
     is the donor: testing across cells would treat 2,000 cells from one person
     as 2,000 independent observations and make nearly everything significant.
     The PsA-vs-psoriasis contrast is additionally run with lane PBMC-01-1
     excluded, because that lane contributes most of the PsA cells and a result
     that survives only with it included is a lane effect, not a disease effect.

  3. SHARED OR TISSUE-SPECIFIC. Programs from different strata are matched on
     gene overlap (Jaccard, against a null from size-matched random gene sets
     drawn from the same panel) and, separately, on whether their scores track
     each other across the lineages present in both tissues. A pair that matches
     on genes but not on usage is the interesting case: the same genes running in
     a different cell type in the other tissue.

Usage:
  python scripts/02_ids_programs/ids_program_atlas.py --dir results/02_ids_programs/ids-disease --out results/02_ids_programs/ids-atlas
"""
import argparse, glob, json, os, sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from sc_ids.lineage import LINEAGE                                           # noqa: E402
EXPRESSING = 0.0        # sc.tl.score_genes is background-subtracted; 0 is the null
MIN_DONOR_CELLS = 20    # a donor x lineage mean below this is not reported
EXCLUDE_LANE = "PBMC-01-1"


def load_programs(d):
    rows = []
    for f in sorted(glob.glob(os.path.join(d, "programs_*.csv"))):
        if "dropped" in f:
            continue
        t = pd.read_csv(f)
        if t.empty:
            continue
        for _, r in t.iterrows():
            rows.append({"program": f"{r['stratum']}_p{int(r['rank']):02d}",
                         "stratum": r["stratum"],
                         "tissue": "PBMC" if r["stratum"].startswith("PBMC") else "skin",
                         "disease_fine": r["stratum"].split("_", 1)[1],
                         "n_genes": int(r["n_genes"]),
                         "median_ids": float(r["median_ids"]),
                         "density": float(r["density"]),
                         "support_cells": int(r.get("cells_half_genes", -1)),
                         "genes": [g.strip() for g in r["genes"].split(",")]})
    return rows


def score(adata, progs, tag):
    """One score column per program, on the log1p representation."""
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    names = []
    for p in progs:
        g = [x for x in p["genes"] if x in adata.var_names]
        if len(g) < 2:
            continue
        sc.tl.score_genes(adata, g, score_name=p["program"], ctrl_size=50,
                          random_state=0)
        names.append(p["program"])
    print(f"  [score] {tag}: {len(names)} of {len(progs)} programs scorable", flush=True)
    return names


def summarise(adata, names, tag):
    o = adata.obs
    lin = o["subtype"].astype(str).map(LINEAGE).fillna("Other")
    out = []
    for p in names:
        v = adata.obs[p].values
        for key, m in o.groupby([lin, o["disease"].astype(str)], observed=True).groups.items():
            idx = o.index.get_indexer(m)
            if len(idx) < MIN_DONOR_CELLS:
                continue
            out.append({"program": p, "tissue": tag, "lineage": key[0],
                        "disease": key[1], "n_cells": len(idx),
                        "mean_score": float(np.mean(v[idx])),
                        "prop_expressing": float(np.mean(v[idx] > EXPRESSING))})
    return pd.DataFrame(out)


def donor_means(adata, names, tag):
    o = adata.obs
    lin = o["subtype"].astype(str).map(LINEAGE).fillna("Other")
    key = pd.DataFrame({"donor": o["donor"].astype(str).values,
                        "lineage": lin.values,
                        "disease": o["disease"].astype(str).values,
                        "lane": o["lane"].astype(str).values if "lane" in o else tag})
    rows = []
    for p in names:
        key["s"] = adata.obs[p].values
        g = key.groupby(["donor", "lineage", "disease", "lane"], observed=True)["s"]
        agg = g.agg(["mean", "size"]).reset_index()
        agg = agg[agg["size"] >= MIN_DONOR_CELLS]
        agg["program"], agg["tissue"] = p, tag
        rows.append(agg.rename(columns={"mean": "mean_score", "size": "n_cells"}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def contrasts(dm, out):
    """Disease comparisons across DONORS, within tissue and lineage."""
    rows = []
    pairs = [("PBMC", "healthy", "psoriasis"), ("PBMC", "healthy", "PsA"),
             ("PBMC", "psoriasis", "PsA"), ("skin", "healthy", "psoriasis")]
    for tissue, a, b in pairs:
        for excl in ([False, True] if (tissue == "PBMC" and b == "PsA") else [False]):
            d = dm[dm["tissue"] == tissue]
            if excl:
                d = d[d["lane"] != EXCLUDE_LANE]
            for (prog, lin), sub in d.groupby(["program", "lineage"], observed=True):
                x = sub.loc[sub["disease"] == a, "mean_score"].values
                y = sub.loc[sub["disease"] == b, "mean_score"].values
                if len(x) < 3 or len(y) < 3:
                    continue
                u, pv = stats.mannwhitneyu(x, y, alternative="two-sided")
                rows.append({"tissue": tissue, "lineage": lin, "program": prog,
                             "group_a": a, "group_b": b,
                             "excl_" + EXCLUDE_LANE: excl,
                             "n_donors_a": len(x), "n_donors_b": len(y),
                             "mean_a": float(x.mean()), "mean_b": float(y.mean()),
                             "delta": float(y.mean() - x.mean()), "p": float(pv)})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    # BH within each (tissue, contrast, exclusion) family
    df["p_adj"] = np.nan
    for k, idx in df.groupby(["tissue", "group_a", "group_b",
                              "excl_" + EXCLUDE_LANE], observed=True).groups.items():
        p = df.loc[idx, "p"].values
        o = np.argsort(p)
        adj = np.empty_like(p)
        adj[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
        df.loc[idx, "p_adj"] = np.clip(adj, 0, 1)
    df = df.sort_values(["tissue", "group_a", "group_b", "p"])
    df.to_csv(os.path.join(out, "disease_contrasts.csv"), index=False)
    return df


def match_across(progs, panel, out, seed=0):
    """Jaccard between every pair of programs from DIFFERENT strata, with a null
    from size-matched random gene sets drawn from the same panel."""
    rng = np.random.default_rng(seed)
    P = list(panel)
    rows = []
    for i, a in enumerate(progs):
        for b in progs[i + 1:]:
            if a["stratum"] == b["stratum"]:
                continue
            A, B = set(a["genes"]), set(b["genes"])
            inter = len(A & B)
            j = inter / len(A | B)
            if inter == 0:
                continue
            null = [len(set(rng.choice(P, len(A), replace=False))
                        & set(rng.choice(P, len(B), replace=False))) / len(A | B)
                    for _ in range(200)]
            rows.append({"program_a": a["program"], "program_b": b["program"],
                         "stratum_a": a["stratum"], "stratum_b": b["stratum"],
                         "cross_tissue": a["tissue"] != b["tissue"],
                         "n_a": len(A), "n_b": len(B), "n_shared": inter,
                         "jaccard": round(j, 4),
                         "jaccard_null_mean": round(float(np.mean(null)), 4),
                         "p_empirical": round(float((np.array(null) >= j).mean()), 4),
                         "shared_genes": ", ".join(sorted(A & B))})
    df = pd.DataFrame(rows).sort_values("jaccard", ascending=False)
    df.to_csv(os.path.join(out, "program_matches.csv"), index=False)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    O = lambda *p: os.path.join(args.out, *p)

    progs = load_programs(args.dir)
    panel = json.load(open(os.path.join(args.dir, "gene_panel.json")))["panel"]
    print(f"[in] {len(progs)} programs from "
          f"{len({p['stratum'] for p in progs})} strata; panel {len(panel):,} genes")
    pd.DataFrame([{**{k: v for k, v in p.items() if k != 'genes'},
                   "genes": ", ".join(p["genes"])} for p in progs]
                 ).to_csv(O("program_catalogue.csv"), index=False)

    sums, dms = [], []
    for tag, f in (("PBMC", "blood_qc.h5ad"), ("skin", "skin_qc.h5ad")):
        a = ad.read_h5ad(os.path.join(args.dir, f))
        names = score(a, progs, tag)
        sums.append(summarise(a, names, tag))
        dms.append(donor_means(a, names, tag))
        a.obs[["tissue", "disease", "donor", "subtype"] + names].to_csv(
            O(f"cell_scores_{tag}.csv.gz"), compression="gzip")
        del a

    S = pd.concat(sums, ignore_index=True)
    S.to_csv(O("program_by_lineage_disease.csv"), index=False)
    DM = pd.concat(dms, ignore_index=True)
    DM.to_csv(O("donor_lineage_means.csv.gz"), index=False, compression="gzip")
    print(f"[out] summary rows {len(S):,}  donor-level rows {len(DM):,}")

    C = contrasts(DM, args.out)
    print(f"[out] contrasts {len(C):,}; significant at p_adj<0.05: "
          f"{int((C['p_adj'] < 0.05).sum()) if len(C) else 0}")

    M = match_across(progs, panel, args.out)
    xt = M[M["cross_tissue"]] if len(M) else M
    print(f"[out] program pairs with any gene overlap {len(M):,}; cross-tissue {len(xt):,}")
    if len(xt):
        print("\nTop cross-tissue matches:")
        for _, r in xt.head(12).iterrows():
            print(f"  {r['program_a']:<20} {r['program_b']:<20} "
                  f"J={r['jaccard']:.3f} (null {r['jaccard_null_mean']:.3f}, "
                  f"p={r['p_empirical']:.3f})  shared: {r['shared_genes'][:70]}")


if __name__ == "__main__":
    main()
