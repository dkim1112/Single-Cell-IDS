"""
cluster_pipeline.py — Agenda 2 / Steps 2-6, on top of the QC'd joint object.

  Step 2  Leiden on the HARMONY-corrected neighbour graph, at the lab pipeline's
          four resolutions (0.20/0.40/0.60/0.80). Every cell is then relabelled
          with the majority annotation of its cluster, which is what absorbs the
          misplaced cells from the per-dataset annotations.
  Step 3  Drop clusters whose majority label lives in only one tissue
          (keratinocyte, fibroblast, mast, B, ... ). Numbers reported.
  Step 4  Pull the T/NK cluster and the myeloid cluster out as separate objects.
  Step 5  Re-embed each ON ITS OWN — new HVGs, new PCA, new Harmony, new
          clustering. With one lineage only, the big between-lineage differences
          are gone and the smaller within-lineage ones (NK vs T) can surface.
  Step 6  Label each sub-cluster by the majority of the original annotations.

Clustering must run on the Harmony graph, not the raw one: on the raw graph the
clusters would just be "skin" and "PBMC".

Usage:
  python scripts/01_qc_labelling/cluster_pipeline.py --joint results/01_qc_labelling/qc-joint/qc_joint_embedded.h5ad \
      --resolution 0.20 --out results/01_qc_labelling/clustered
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
import matplotlib
import matplotlib.patheffects
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

RESOLUTIONS = [0.20, 0.40, 0.60, 0.80]          # lab config_sample.json
TISSUE_COLORS = {"PBMC": "#1f77b4", "skin": "#ff7f0e"}
# Labels that exist in both tissues. NK is PBMC-only as a LABEL, but it clusters
# with T cells and Matthew wants the T/NK cluster kept, so it rides along inside
# that cluster rather than being kept as its own type.
SHARED = {"T cell", "Myeloid cell"}


def majority_label(obs, cluster_col, label_col="main"):
    t = pd.crosstab(obs[cluster_col], obs[label_col])
    return t.idxmax(axis=1), t


def cluster(e, resolutions, key_prefix="leiden", neighbors_key="after_harmony"):
    for r in resolutions:
        sc.tl.leiden(e, resolution=r, key_added=f"{key_prefix}_res_{r:4.2f}".replace(" ", ""),
                     flavor="igraph", n_iterations=2, neighbors_key=neighbors_key,
                     directed=False)
    return e


def mixing_counts(e, rep="X_pca_harmony", n=30, group="tissue", within=None):
    """Average number of a cell's 30 nearest neighbours that come from the other
    tissue — counts, because fractions of a fraction get confusing."""
    idx = NearestNeighbors(n_neighbors=n + 1).fit(e.obsm[rep]).kneighbors(
        e.obsm[rep], return_distance=False)[:, 1:]
    g = e.obs[group].astype(str).values
    k = (g[idx] != g[:, None]).sum(1)
    df = pd.DataFrame({"group": g, "k": k})
    if within is not None:
        df[within] = e.obs[within].astype(str).values
        return df.groupby(["group", within]).k.agg(["size", "mean"]).round(2).reset_index()
    return df.groupby("group").k.agg(["size", "mean"]).round(2).reset_index()


def plot(e, keys, out_png, title, basis="X_umap", legend_on_data=False):
    rng = np.random.default_rng(0)
    order = rng.permutation(e.n_obs)
    U = e.obsm[basis][order]
    fig, axes = plt.subplots(1, len(keys), figsize=(6.4 * len(keys), 5.2))
    axes = np.atleast_1d(axes)
    for ax, key in zip(axes, keys):
        lab = e.obs[key].astype(str).values[order]
        cats = sorted(pd.unique(lab), key=lambda x: (len(x), x))
        cmap = plt.get_cmap("tab20")
        colors = ({c: TISSUE_COLORS[c] for c in cats} if key == "tissue"
                  else {c: matplotlib.colors.to_hex(cmap(i % 20)) for i, c in enumerate(cats)})
        ax.scatter(U[:, 0], U[:, 1], s=1.2, c=[colors[x] for x in lab], linewidths=0, rasterized=True)
        if legend_on_data or key.startswith("leiden"):
            for cth in cats:
                m = e.obs[key].astype(str).values == cth
                ax.text(np.median(e.obsm[basis][m, 0]), np.median(e.obsm[basis][m, 1]), cth,
                        fontsize=8, ha="center", va="center",
                        path_effects=[matplotlib.patheffects.withStroke(linewidth=2.4, foreground="white")])
        else:
            counts = pd.Series(lab).value_counts()
            ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=colors[g], ms=6,
                                          label=f"{g} ({counts[g]:,})") for g in counts.index],
                      loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=8)
        ax.set_title(f"{title} — {key}", fontsize=11)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_color("#BBBBBB")
    fig.tight_layout(); fig.savefig(out_png, dpi=140, bbox_inches="tight"); plt.close(fig)


def embed_lineage(sub, theta, n_pcs=20, n_hvg=2000):
    """Same recipe as the lab pipeline, but with Harmony's strength exposed.

    Harmony's `theta` controls how hard it pushes the batches together. The
    default (2) is tuned for batches of the same tissue; within a single lineage
    it leaves skin and blood almost fully separate. The sweep in theta_sweep.csv
    shows a higher theta fixes that without disturbing within-tissue structure.
    """
    import harmonypy
    sc.pp.normalize_total(sub, target_sum=1e6, inplace=True)
    sc.pp.log1p(sub, copy=False)
    sc.pp.highly_variable_genes(sub, n_top_genes=min(n_hvg, sub.n_vars - 1))
    sc.tl.pca(sub)
    sub.obsm["X_pca"] = sub.obsm["X_pca"][:, :n_pcs]
    sc.pp.neighbors(sub, key_added="before_harmony", use_rep="X_pca")
    sc.tl.umap(sub, neighbors_key="before_harmony")
    sub.obsm["X_umap_before"] = sub.obsm["X_umap"].copy()
    ho = harmonypy.run_harmony(sub.obsm["X_pca"], sub.obs, ["tissue"], theta=theta,
                               max_iter_harmony=20)
    Z = np.asarray(ho.Z_corr)
    sub.obsm["X_pca_harmony"] = Z if Z.shape[0] == sub.n_obs else Z.T
    sc.pp.neighbors(sub, key_added="after_harmony", use_rep="X_pca_harmony")
    sc.tl.umap(sub, neighbors_key="after_harmony")
    return sub


def neighbour_stats(e, rep="X_pca_harmony", n=30):
    """(mean other-tissue neighbours per skin cell, per PBMC cell, subtype purity
    measured only among same-tissue neighbours)."""
    idx = NearestNeighbors(n_neighbors=n + 1).fit(e.obsm[rep]).kneighbors(
        e.obsm[rep], return_distance=False)[:, 1:]
    t = e.obs["tissue"].astype(str).values
    s = e.obs["subtype"].astype(str).values
    other = t[idx] != t[:, None]
    same_t = ~other
    pure = ((s[idx] == s[:, None]) & same_t).sum(1) / np.maximum(same_t.sum(1), 1)
    return (round(float(other[t == "skin"].sum(1).mean()), 2),
            round(float(other[t == "PBMC"].sum(1).mean()), 2),
            round(float(pure.mean()), 3))


MARKERS = {
    "monocyte": ["CD14", "S100A8", "S100A9", "VCAN", "FCN1", "LYZ"],
    "dendritic": ["FCER1A", "CD1C", "CLEC10A", "CD207", "CD1A"],
    "T": ["CD3E", "CD3D", "TRAC", "IL7R"],
    "NK": ["GNLY", "NKG7", "KLRD1", "PRF1"],
}


def marker_check(sub, cluster_col, out_csv):
    """Per sub-cluster and tissue, the share of cells expressing each marker.

    A sub-cluster that mixes the tissues is only a real match if BOTH halves
    carry the same markers; two different cell types can share a cluster and
    still look well mixed. "Detected" = at least one read, which is unaffected
    by normalisation (zeros stay zeros), so this reads the same off counts or
    normalised values.
    """
    X = sp.csr_matrix(sub.layers["counts"]) if "counts" in sub.layers else sp.csr_matrix(sub.X)
    rows = []
    cl = sub.obs[cluster_col].astype(str).values
    tis = sub.obs["tissue"].astype(str).values
    for c in sorted(pd.unique(cl)):
        for t in ("PBMC", "skin"):
            m = (cl == c) & (tis == t)
            if m.sum() < 20:
                continue
            row = {"subcluster": c, "tissue": t, "n": int(m.sum())}
            for family, genes in MARKERS.items():
                g = [x for x in genes if x in sub.var_names]
                if not g:
                    continue
                det = (X[np.where(m)[0]][:, [sub.var_names.get_loc(x) for x in g]] > 0
                       ).toarray().mean(0)
                for name, d in zip(g, det):
                    row[f"{name} ({family})"] = round(100 * float(d), 1)
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--joint", required=True)
    ap.add_argument("--resolution", type=float, default=0.20)
    ap.add_argument("--sub-resolution", type=float, default=0.40)
    ap.add_argument("--theta", type=float, default=8.0,
                    help="Harmony strength for the per-lineage re-embedding (Step 5)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    O = lambda *p: os.path.join(args.out, *p)
    log = {"resolution": args.resolution, "sub_resolution": args.sub_resolution}

    e = ad.read_h5ad(args.joint)
    print(f"[in] {e.n_obs:,} cells x {e.n_vars:,} genes", flush=True)

    # ---------------------------------------------------------------- Step 2
    e = cluster(e, RESOLUTIONS)
    key = f"leiden_res_{args.resolution:4.2f}".replace(" ", "")
    n_by_res = {f"{r:.2f}": int(e.obs[f'leiden_res_{r:4.2f}'.replace(" ", "")].nunique())
                for r in RESOLUTIONS}
    log["clusters_by_resolution"] = n_by_res
    print("[step2] clusters per resolution:", n_by_res, flush=True)

    maj, xt = majority_label(e.obs, key)
    xt.to_csv(O("step2_cluster_by_celltype.csv"))
    pd.crosstab(e.obs[key], e.obs["tissue"]).to_csv(O("step2_cluster_by_tissue.csv"))
    e.obs["cluster_label"] = e.obs[key].map(maj).astype(str)
    agree = (e.obs["cluster_label"] == e.obs["main"].astype(str))
    log["relabelled_cells"] = int((~agree).sum())
    log["relabelled_pct"] = round(100 * float((~agree).mean()), 2)
    (e.obs.loc[~agree].groupby(["tissue", "main", "cluster_label"], observed=True).size()
       .rename("n").reset_index().sort_values("n", ascending=False)
       .to_csv(O("step2_relabelled_cells.csv"), index=False))
    plot(e, ["tissue", key, "cluster_label"], O("step2_umap.png"), "Step 2")
    print(f"[step2] {log['relabelled_cells']:,} cells relabelled "
          f"({log['relabelled_pct']}%)", flush=True)

    # ---------------------------------------------------------------- Step 3
    per_tissue = pd.crosstab(e.obs["cluster_label"], e.obs["tissue"])
    keep_labels = sorted(SHARED & set(per_tissue.index))
    dropped = per_tissue.drop(index=keep_labels)
    dropped.to_csv(O("step3_dropped_clusters.csv"))
    before = e.obs.groupby("tissue", observed=True).size()
    e2 = e[e.obs["cluster_label"].isin(keep_labels)].copy()
    after = e2.obs.groupby("tissue", observed=True).size()
    log["step3"] = {"kept_labels": keep_labels,
                    "cells_before": {k: int(v) for k, v in before.items()},
                    "cells_after": {k: int(v) for k, v in after.items()}}
    print("[step3] kept", keep_labels, dict(after), flush=True)
    (pd.crosstab(e2.obs["cluster_label"], e2.obs["tissue"])
       .to_csv(O("step3_kept_clusters.csv")))

    # ------------------------------------------------------------- Steps 4-6
    sub_summary = []
    for lineage in keep_labels:
        tag = lineage.split()[0].lower()
        sub = e2[e2.obs["cluster_label"] == lineage].copy()
        sub = ad.AnnData(X=sub.layers["counts"] if "counts" in sub.layers else sub.raw.X
                         if sub.raw is not None else sub.X,
                         obs=sub.obs.copy(), var=pd.DataFrame(index=sub.var_names))
        for c in ("tissue", "main", "subtype", "donor"):
            sub.obs[c] = sub.obs[c].astype(str).astype("category")
        print(f"[step4] {lineage}: {sub.n_obs:,} cells", flush=True)
        sub = embed_lineage(sub, theta=args.theta)
        sub = cluster(sub, RESOLUTIONS)
        skey = f"leiden_res_{args.sub_resolution:4.2f}".replace(" ", "")
        smaj, sxt = majority_label(sub.obs, skey, "subtype")
        sxt.to_csv(O(f"step5_{tag}_subcluster_by_subtype.csv"))
        pd.crosstab(sub.obs[skey], sub.obs["tissue"]).to_csv(O(f"step5_{tag}_subcluster_by_tissue.csv"))
        sub.obs["subcluster_label"] = sub.obs[skey].map(smaj).astype(str)
        # Step 6, read per tissue as well: a single majority label collapses to the
        # bigger tissue's vocabulary (PBMC here), which hides what the skin cells in
        # a cluster actually are.
        lab_rows = []
        for cl, g in sub.obs.groupby(skey, observed=True):
            row = {"subcluster": cl, "n": len(g), "label_overall": smaj[cl]}
            for tis in ("PBMC", "skin"):
                gg = g[g["tissue"].astype(str) == tis]
                row[f"n_{tis}"] = len(gg)
                row[f"top_{tis}"] = (gg["subtype"].astype(str).value_counts().idxmax()
                                     if len(gg) else "")
                row[f"top_{tis}_share"] = (round(100 * gg["subtype"].astype(str)
                                                 .value_counts(normalize=True).iloc[0], 1)
                                           if len(gg) else np.nan)
            lab_rows.append(row)
        pd.DataFrame(lab_rows).to_csv(O(f"step6_{tag}_labels.csv"), index=False)
        marker_check(sub, skey, O(f"step6_{tag}_marker_check.csv"))
        plot(sub, ["tissue", skey, "subcluster_label"], O(f"step5_{tag}_umap.png"), lineage)
        plot(sub, ["tissue"], O(f"step5_{tag}_umap_before.png"), f"{lineage} before Harmony",
             basis="X_umap_before")
        mc = mixing_counts(sub); mc.insert(0, "lineage", lineage)
        mc.to_csv(O(f"step5_{tag}_mixing.csv"), index=False)
        sk, pb, pu = neighbour_stats(sub)
        sk0, pb0, pu0 = neighbour_stats(sub, rep="X_pca")
        sub.write_h5ad(O(f"step5_{tag}_embedded.h5ad"))
        sub_summary.append(dict(lineage=lineage, cells=int(sub.n_obs), theta=args.theta,
                                clusters=int(sub.obs[skey].nunique()),
                                skin_nbrs_from_PBMC_before=sk0, skin_nbrs_from_PBMC_after=sk,
                                PBMC_nbrs_from_skin_after=pb,
                                subtype_purity_before=pu0, subtype_purity_after=pu,
                                labels=", ".join(sorted(set(sub.obs["subcluster_label"])))))
        print(f"[step5] {lineage}: {sub.obs[skey].nunique()} sub-clusters -> "
              f"{sorted(set(sub.obs['subcluster_label']))}", flush=True)
    pd.DataFrame(sub_summary).to_csv(O("step5_summary.csv"), index=False)

    mixing_counts(e, within="cluster_label").to_csv(O("step2_mixing_counts.csv"), index=False)
    e.write_h5ad(O("clustered.h5ad"))
    json.dump(log, open(O("run_log.json"), "w"), indent=2)
    print(json.dumps(log, indent=2), flush=True)


if __name__ == "__main__":
    main()
