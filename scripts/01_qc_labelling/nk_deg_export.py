"""
nk_deg_export.py — assemble the NK-cluster DEG table for outside review.

Pulls together everything someone needs to confirm (or reject) the NK call on
T/NK sub-cluster 0, computed on the settled parameters (Harmony theta = 4):

  1. the up-regulated genes for the cluster vs all other T/NK sub-clusters,
     with log2 fold change, p-value and adjusted p-value (scanpy
     rank_genes_groups, Wilcoxon — the lab's run_Ranking.py method);
  2. for each of those genes, the share of cells it is detected in, split by
     tissue, because the blood and skin halves of this cluster differ;
  3. whether the gene survives the 1% gene filter, so the filter question can be
     answered on evidence rather than in principle;
  4. differential expression between the skin and blood cells INSIDE the cluster.

Usage:
  python scripts/01_qc_labelling/nk_deg_export.py --lineage results/01_qc_labelling/clustered-theta4/step5_t_embedded.h5ad \
      --cluster 0 --genes results/01_qc_labelling/qc-joint/gene_lists.json --out results/01_qc_labelling/markers-theta4
"""
import argparse, json, os
import numpy as np
import scipy.sparse as sp
import scanpy as sc
import anndata as ad

PVAL, LOG2FC = 0.01, 1.0          # lab config_sample.json thresholds


def pct_detected(a, genes, mask):
    X = sp.csr_matrix(a.X)
    idx = {g: a.var_names.get_loc(g) for g in genes if g in a.var_names}
    rows = np.where(mask)[0]
    out = {}
    for g, j in idx.items():
        v = X[rows, j].toarray().ravel()
        out[g] = round(100 * float((v > 0).mean()), 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lineage", required=True)
    ap.add_argument("--cluster", default="0")
    ap.add_argument("--cluster-col", default="leiden_res_0.40")
    ap.add_argument("--genes", required=True, help="gene_lists.json from the QC step")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    O = lambda *p: os.path.join(args.out, *p)

    e = ad.read_h5ad(args.lineage)
    cl = e.obs[args.cluster_col].astype(str).values
    tis = e.obs["tissue"].astype(str).values
    kept = set(json.load(open(args.genes))["kept"])
    print(f"[in] {e.n_obs:,} cells x {e.n_vars:,} genes; 1% filter keeps {len(kept):,}", flush=True)

    # ---------- cluster 0 vs all other sub-clusters
    sc.tl.rank_genes_groups(e, groupby=args.cluster_col, method="wilcoxon", key_added="rank")
    df = sc.get.rank_genes_groups_df(e, key="rank", group=str(args.cluster))
    up = df[(df.pvals_adj < PVAL) & (df.logfoldchanges > LOG2FC)].copy()
    up = up.sort_values("scores", ascending=False).reset_index(drop=True)

    m0 = cl == str(args.cluster)
    g = list(up.names)
    pb = pct_detected(e, g, m0 & (tis == "PBMC"))
    sk = pct_detected(e, g, m0 & (tis == "skin"))
    up["pct_detected_blood"] = [pb.get(x, np.nan) for x in up.names]
    up["pct_detected_skin"] = [sk.get(x, np.nan) for x in up.names]
    up["survives_1pct_filter"] = ["yes" if x in kept else "NO" for x in up.names]
    up["rank"] = np.arange(1, len(up) + 1)
    up = up[["rank", "names", "logfoldchanges", "pvals", "pvals_adj", "scores",
             "pct_detected_blood", "pct_detected_skin", "survives_1pct_filter"]]
    up.columns = ["rank", "gene", "log2FC", "p_value", "p_adj", "wilcoxon_score",
                  "% cells detected (blood half)", "% cells detected (skin half)",
                  "survives 1% gene filter"]
    up.to_csv(O(f"NK_cluster{args.cluster}_upregulated.csv"), index=False)
    print(f"  up-regulated genes passing p_adj<{PVAL}, log2FC>{LOG2FC}: {len(up):,}", flush=True)
    lost = up[up["survives 1% gene filter"] == "NO"]
    print(f"  of those, removed by the 1% filter: {len(lost):,} "
          f"({', '.join(lost.gene.head(25))})", flush=True)
    print(up.head(25).to_string(index=False), flush=True)

    # ---------- skin vs blood inside the cluster
    sub = e[m0].copy()
    sc.tl.rank_genes_groups(sub, groupby="tissue", method="wilcoxon", key_added="rt")
    d2 = sc.get.rank_genes_groups_df(sub, key="rt", group=None)
    d2 = d2[(d2.pvals_adj < PVAL) & (d2.logfoldchanges > LOG2FC)].copy()
    d2 = d2.sort_values(["group", "scores"], ascending=[True, False])
    d2["survives_1pct_filter"] = ["yes" if x in kept else "NO" for x in d2.names]
    d2.to_csv(O(f"NK_cluster{args.cluster}_skin_vs_blood.csv"), index=False)
    print(f"  skin-vs-blood inside cluster: {len(d2):,} genes", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
