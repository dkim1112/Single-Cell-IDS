"""
theta_tuning.py — how far should Harmony's theta go?

The goal is not maximum mixing. We want three things at once:
  1. skin and blood aligned (cells of one tissue have neighbours from the other),
  2. sub-clusters that are still DISTINCT (not one smeared blob),
  3. each sub-cluster containing BOTH tissues, not one tissue each.

Metrics per theta, on one lineage:
  skin_nbrs_from_PBMC   mean number of a skin cell's 30 nearest neighbours that
                        are blood cells (and the reverse). Higher = more aligned.
  n_clusters            Leiden clusters at the working resolution.
  clusters_with_both    clusters where the smaller tissue is >= 10% and >= 20 cells.
  min_tissue_share      smallest minority-tissue share across clusters.
  silhouette            how separated the clusters are in the corrected space.
                        Falls when correction smears everything together.
  ARI_vs_annotation     agreement between the clusters and the ORIGINAL per-tissue
                        annotations, computed within each tissue. This is the guard
                        against over-correction: if pushing theta destroys real
                        cell-type structure, agreement with the labels the source
                        studies assigned drops.

Usage:
  python scripts/01_qc_labelling/theta_tuning.py --lineage results/01_qc_labelling/clustered-nofilter/step5_t_embedded.h5ad \
      --tag t --resolution 0.40 --out results/01_qc_labelling/theta
"""
import argparse, os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
import harmonypy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import adjusted_rand_score, silhouette_score

THETAS = [None, 1, 2, 4, 8, 16, 32]      # None = no correction
TISSUE_COLORS = {"PBMC": "#1f77b4", "skin": "#ff7f0e"}


def harmony(Z, obs, theta):
    ho = harmonypy.run_harmony(Z, obs, ["tissue"], theta=theta, max_iter_harmony=20)
    H = np.asarray(ho.Z_corr)
    return H if H.shape[0] == Z.shape[0] else H.T


def metrics(e, Z, cl):
    t = e.obs["tissue"].astype(str).values
    s = e.obs["subtype"].astype(str).values
    idx = NearestNeighbors(n_neighbors=31).fit(Z).kneighbors(Z, return_distance=False)[:, 1:]
    other = (t[idx] != t[:, None]).sum(1)
    out = {"skin_nbrs_from_PBMC": round(float(other[t == "skin"].mean()), 2),
           "PBMC_nbrs_from_skin": round(float(other[t == "PBMC"].mean()), 2)}
    # within-tissue subtype purity, as before
    same_t = t[idx] == t[:, None]
    pure = ((s[idx] == s[:, None]) & same_t).sum(1) / np.maximum(same_t.sum(1), 1)
    out["subtype_purity"] = round(float(pure.mean()), 3)

    tab = pd.crosstab(cl, t)
    for c in ("PBMC", "skin"):
        if c not in tab:
            tab[c] = 0
    share = tab.min(axis=1) / tab.sum(axis=1)
    out["n_clusters"] = int(len(tab))
    out["clusters_with_both"] = int(((share >= 0.10) & (tab.min(axis=1) >= 20)).sum())
    out["min_tissue_share"] = round(float(share.min()), 3)
    out["median_tissue_share"] = round(float(share.median()), 3)

    if len(set(cl)) > 1:
        sub = np.random.default_rng(0).choice(len(cl), size=min(5000, len(cl)), replace=False)
        out["silhouette"] = round(float(silhouette_score(Z[sub], np.asarray(cl)[sub])), 3)
    else:
        out["silhouette"] = np.nan
    for c in ("PBMC", "skin"):
        m = t == c
        out[f"ARI_vs_annotation_{c}"] = round(float(adjusted_rand_score(s[m], np.asarray(cl)[m])), 3)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lineage", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--resolution", type=float, default=0.40)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    O = lambda *p: os.path.join(args.out, *p)

    e = ad.read_h5ad(args.lineage)
    base = e.obsm["X_pca"]
    key = f"leiden_res_{args.resolution:4.2f}".replace(" ", "")
    rows, embeds = [], {}
    for th in THETAS:
        name = "none" if th is None else str(th)
        Z = base if th is None else harmony(base, e.obs, th)
        w = ad.AnnData(X=sp.csr_matrix((e.n_obs, 1)), obs=e.obs.copy())
        w.obsm["Z"] = Z
        sc.pp.neighbors(w, use_rep="Z", key_added="g")
        sc.tl.leiden(w, resolution=args.resolution, key_added=key, flavor="igraph",
                     n_iterations=2, neighbors_key="g", directed=False)
        sc.tl.umap(w, neighbors_key="g")
        cl = w.obs[key].astype(str).values
        r = {"theta": name}
        r.update(metrics(e, Z, cl))
        rows.append(r)
        embeds[name] = (w.obsm["X_umap"], cl)
        print(f"[theta {name}] " + "  ".join(f"{k}={v}" for k, v in r.items() if k != "theta"), flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(O(f"theta_tuning_{args.tag}.csv"), index=False)
    print(df.to_string(index=False), flush=True)

    # one row of UMAPs by tissue, one by cluster
    t = e.obs["tissue"].astype(str).values
    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(2, len(THETAS), figsize=(3.1 * len(THETAS), 6.6))
    for j, th in enumerate(THETAS):
        name = "none" if th is None else str(th)
        U, cl = embeds[name]
        o = rng.permutation(len(U))
        axes[0, j].scatter(U[o, 0], U[o, 1], s=0.8, c=[TISSUE_COLORS[x] for x in t[o]],
                           linewidths=0, rasterized=True)
        cats = sorted(set(cl), key=int)
        cmap = plt.get_cmap("tab20")
        col = {c: matplotlib.colors.to_hex(cmap(i % 20)) for i, c in enumerate(cats)}
        axes[1, j].scatter(U[o, 0], U[o, 1], s=0.8, c=[col[x] for x in cl[o]],
                           linewidths=0, rasterized=True)
        axes[0, j].set_title(f"theta = {name}", fontsize=10)
        axes[1, j].set_title(f"{len(cats)} clusters", fontsize=9, color="#555")
        for i in (0, 1):
            axes[i, j].set_xticks([]); axes[i, j].set_yticks([])
            for s_ in axes[i, j].spines.values():
                s_.set_color("#BBBBBB")
    axes[0, 0].set_ylabel("coloured by tissue", fontsize=10)
    axes[1, 0].set_ylabel("coloured by sub-cluster", fontsize=10)
    fig.suptitle(f"{args.tag}: effect of Harmony strength (theta)", fontsize=12)
    fig.tight_layout()
    fig.savefig(O(f"theta_umaps_{args.tag}.png"), dpi=130, bbox_inches="tight")
    print("done", flush=True)


if __name__ == "__main__":
    main()
