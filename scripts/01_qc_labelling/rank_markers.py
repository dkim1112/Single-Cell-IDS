"""
rank_markers.py — Agenda 2 / Step 7: which genes are up in each sub-cluster, and
is the NK sub-cluster NK cells or gamma-delta T cells?

Differential expression follows the lab pipeline (run_Ranking.py): scanpy
`rank_genes_groups`, method 'wilcoxon', one cluster against all the others,
then filtered at the lab's config thresholds (p < 0.01, log2FC > 1).

The NK-vs-gamma-delta question is settled with the panel the lab suggested —
NK = FCGR3A, NKG7, GZMB, PRF1 present and CD3D ABSENT — plus the gamma-delta
T-cell receptor genes (TRDC, TRGC1, TRGC2), which are positive evidence for
gamma-delta rather than evidence against NK. The panel is applied to each
sub-cluster AND to each tissue within it, because the question is specifically
whether the skin cells sitting in the NK cluster are the same thing as the blood
cells there.

IMPORTANT — run this on the UNFILTERED gene set. The 1% gene filter removes
FCGR3A, PRF1, KLRD1, TRDC and TRGC1/2 outright, so the question cannot be
answered on that object at all.

Usage:
  python scripts/01_qc_labelling/rank_markers.py --lineage results/01_qc_labelling/clustered-nofilter/step5_t_embedded.h5ad \
      --tag t --cluster-col leiden_res_0.40 --out results/01_qc_labelling/markers
"""
import argparse, os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad

PVAL, LOG2FC = 0.01, 1.0                      # lab config_sample.json
PANEL = {
    "NK (should be present)": ["FCGR3A", "NKG7", "GZMB", "PRF1"],
    "NK (should be absent)": ["CD3D"],
    "T cell": ["CD3E", "TRAC", "TRBC2", "IL7R"],
    "gamma-delta receptor": ["TRDC", "TRGC1", "TRGC2", "TRDV1"],
    "other NK": ["KLRD1", "KLRF1", "NCAM1", "GNLY", "TYROBP"],
}


def detection(adata, genes, mask):
    """Share of cells with at least one read, and mean normalised expression."""
    X = sp.csr_matrix(adata.X)
    out = {}
    for g in genes:
        if g not in adata.var_names:
            out[g] = (np.nan, np.nan)
            continue
        v = X[np.where(mask)[0], adata.var_names.get_loc(g)].toarray().ravel()
        out[g] = (round(100 * float((v > 0).mean()), 1), round(float(v.mean()), 2))
    return out


def panel_table(adata, cluster_col, split_tissue=True):
    rows = []
    cl = adata.obs[cluster_col].astype(str).values
    tis = adata.obs["tissue"].astype(str).values
    groups = ([(c, t) for c in sorted(pd.unique(cl)) for t in ("PBMC", "skin")]
              if split_tissue else [(c, None) for c in sorted(pd.unique(cl))])
    for c, t in groups:
        mask = (cl == c) if t is None else ((cl == c) & (tis == t))
        if mask.sum() < 20:
            continue
        row = {"subcluster": c, "tissue": t or "both", "n": int(mask.sum())}
        for family, genes in PANEL.items():
            for g, (pct, mean) in detection(adata, genes, mask).items():
                row[f"{g} %det"] = pct
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lineage", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--cluster-col", default="leiden_res_0.40")
    ap.add_argument("--nk-cluster", default=None,
                    help="cluster id to split skin vs blood inside (the NK one)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    O = lambda *p: os.path.join(args.out, *p)

    e = ad.read_h5ad(args.lineage)
    print(f"[in] {args.tag}: {e.n_obs:,} cells x {e.n_vars:,} genes", flush=True)

    # ---- lab ranking step: each sub-cluster vs all the others
    sc.tl.rank_genes_groups(e, groupby=args.cluster_col, method="wilcoxon",
                            key_added="rank")
    df = sc.get.rank_genes_groups_df(e, key="rank", group=None,
                                     pval_cutoff=PVAL, log2fc_min=LOG2FC)
    df.to_csv(O(f"marker_genes_{args.tag}.csv"), index=False)
    top = (df.sort_values(["group", "scores"], ascending=[True, False])
             .groupby("group").head(20))
    top.to_csv(O(f"marker_genes_{args.tag}_top20.csv"), index=False)
    for g, sub in top.groupby("group"):
        print(f"  cluster {g}: " + ", ".join(sub.names.head(15)), flush=True)

    # ---- the marker panel, per cluster and per tissue within cluster
    panel_table(e, args.cluster_col, split_tissue=False).to_csv(
        O(f"panel_{args.tag}_by_cluster.csv"), index=False)
    pt = panel_table(e, args.cluster_col, split_tissue=True)
    pt.to_csv(O(f"panel_{args.tag}_by_cluster_tissue.csv"), index=False)
    print(pt.to_string(index=False), flush=True)

    # ---- inside the NK cluster: are the skin cells the same thing as the blood cells?
    if args.nk_cluster is not None:
        sub = e[e.obs[args.cluster_col].astype(str) == str(args.nk_cluster)].copy()
        if sub.obs["tissue"].nunique() > 1:
            sc.tl.rank_genes_groups(sub, groupby="tissue", method="wilcoxon",
                                    key_added="rank_tissue")
            d2 = sc.get.rank_genes_groups_df(sub, key="rank_tissue", group=None,
                                             pval_cutoff=PVAL, log2fc_min=LOG2FC)
            d2.to_csv(O(f"nk_cluster_skin_vs_blood_{args.tag}.csv"), index=False)
            for g, s in (d2.sort_values(["group", "scores"], ascending=[True, False])
                           .groupby("group")):
                print(f"  [NK cluster] up in {g}: " + ", ".join(s.names.head(15)), flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
