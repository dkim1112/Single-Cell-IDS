"""
viz.py — PCA -> Harmony -> neighbors -> UMAP, and plotting IDS modules on the UMAP.

Follows the lab's pipeline (CutaneousBioinf/single-cell_pipeline):
    normalize_total(1e6) -> log1p -> highly_variable_genes(2000) -> sc.tl.pca
    -> harmony_integrate(key) -> sc.pp.neighbors(use_rep=...) -> sc.tl.umap

TWO REPRESENTATIONS, ON PURPOSE — do not confuse them:
  * UMAP/clustering use log1p-normalized data reduced to PCs. That is the standard
    single-cell representation and what Harmony corrects.
  * IDS uses normalized counts min-max scaled to [0, 8] with NO log1p (the paper's
    recipe; log1p demonstrably collapses IDS — see pipeline.preprocess docstring).
  They serve different purposes and are computed separately from the same counts.

WHAT HARMONY DOES (verified empirically, see harmony_ids_test.py):
  Harmony corrects the PCA EMBEDDING (obsm['X_pca'] -> obsm['X_pca_harmony']).
  It does NOT modify adata.X. Therefore it cannot distort IDS — and equally, it
  cannot batch-correct the input IDS reads. Batch is handled for IDS by
  stratification (--subset-col), not by correction.
"""
import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def run_harmony(adata, batch_key, use_rep="X_pca", adjusted="X_pca_harmony"):
    """Batch-correct the PCA embedding with Harmony.

    Calls harmonypy directly: sc.external.pp.harmony_integrate is broken with
    harmonypy>=2.0.0 (its wrapper transposes Z_corr assuming the old PCs-x-cells
    convention, while 2.0.0 already returns cells-x-PCs). The lab pins 0.0.10.
    """
    import harmonypy
    ho = harmonypy.run_harmony(adata.obsm[use_rep], adata.obs, [batch_key])
    Z = np.asarray(ho.Z_corr)
    adata.obsm[adjusted] = Z if Z.shape[0] == adata.n_obs else Z.T
    return adata


def compute_embedding(adata, batch_key=None, n_pcs=20, n_hvg=2000,
                      target_sum=1e6, verbose=True, with_umap=True,
                      store_raw=True, copy=True):
    """normalize -> log1p -> HVG -> PCA -> (Harmony) -> neighbors -> UMAP.

    If batch_key is given, computes UMAP twice: 'X_umap_before_harmony' (on X_pca)
    and 'X_umap' (on X_pca_harmony), so correction can be verified visually.
    Returns adata with obsm populated.
    """
    # `copy=False` mutates the caller's object. That is the wrong default for
    # interactive use, but on a large cohort the defensive copy is a second full
    # matrix in memory, so a caller that has no further use for the input passes
    # False rather than being OOM-killed for tidiness.
    adata = adata.copy() if copy else adata
    # `store_raw` keeps scanpy's usual conveniences — a counts layer and .raw —
    # but each is another full copy of the matrix. On a large cohort that is the
    # difference between fitting in memory and being OOM-killed, and neither is
    # read by anything in this project, so callers working at scale pass False.
    if store_raw:
        adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum, inplace=True)
    sc.pp.log1p(adata, copy=False)
    if store_raw:
        adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=min(n_hvg, adata.n_vars - 1))
    sc.tl.pca(adata)
    adata.obsm["X_pca"] = adata.obsm["X_pca"][:, :n_pcs]
    if verbose:
        print(f"  [viz] PCA: {adata.obsm['X_pca'].shape[1]} PCs on {adata.n_obs} cells")

    if batch_key is not None:
        # neighbours (and optionally UMAP) BEFORE correction
        sc.pp.neighbors(adata, key_added="before_harmony", use_rep="X_pca")
        if with_umap:
            sc.tl.umap(adata, neighbors_key="before_harmony")
            adata.obsm["X_umap_before_harmony"] = adata.obsm["X_umap"].copy()
        # correct, then neighbours (and optionally UMAP) AFTER
        run_harmony(adata, batch_key)
        sc.pp.neighbors(adata, key_added="after_harmony", use_rep="X_pca_harmony")
        if with_umap:
            sc.tl.umap(adata, neighbors_key="after_harmony")
        if verbose:
            print(f"  [viz] Harmony on '{batch_key}' -> X_pca_harmony"
                  f"{'; UMAP before & after' if with_umap else ' (no UMAP)'}")
    else:
        sc.pp.neighbors(adata, use_rep="X_pca")
        if with_umap:
            sc.tl.umap(adata)
    return adata


def cluster_cells(adata, batch_key=None, resolution=0.6, n_pcs=20, n_hvg=2000,
                  return_adata=False, verbose=True):
    """Data-derived cell groups the way the lab pipeline makes them.

    Mirrors CutaneousBioinf run_PCA -> run_Harmony -> run_UMAP -> run_Clustering:
        normalize_total(1e6) -> log1p -> HVG(2000) -> PCA(n_pcs)
        -> Harmony(batch_key) -> neighbors(use_rep=X_pca_harmony) -> leiden(resolution)

    WHY THIS EXISTS: the lab makes its GENE-LEVEL step (rank_genes_groups) batch-robust
    not by correcting expression — Harmony never produces corrected expression — but by
    GROUPING on Harmony-corrected clusters and reading raw expression within each group.
    IDS sits in exactly that slot, so we replicate the pattern: cluster here, then run
    IDS inside each cluster on uncorrected counts.

    Operates on a copy; `adata` is NOT modified. Returns a pandas Series of cluster
    labels aligned to adata.obs_names (or the embedded AnnData if return_adata=True).
    """
    work = compute_embedding(adata, batch_key=batch_key, n_pcs=n_pcs, n_hvg=n_hvg,
                             verbose=verbose, with_umap=return_adata)
    nkey = "after_harmony" if batch_key is not None else None
    key = f"leiden_res_{resolution:4.2f}".replace(" ", "")
    kw = dict(resolution=resolution, key_added=key, flavor="igraph", n_iterations=2)
    if nkey:
        kw["neighbors_key"] = nkey
    sc.tl.leiden(work, **kw)
    if verbose:
        n = work.obs[key].nunique()
        print(f"  [viz] leiden(resolution={resolution}) on "
              f"{'Harmony-corrected' if batch_key else 'uncorrected'} graph -> {n} clusters")
    return work if return_adata else work.obs[key]


def score_modules(adata, modules, prefix="module"):
    """Add a per-cell expression score for each IDS module (scanpy score_genes).
    `modules` is the list returned by pipeline.extract_modules. Returns score names."""
    names = []
    for i, m in enumerate(modules):
        genes = [g for g in m["genes"] if g in adata.var_names]
        if len(genes) < 2:
            continue
        name = f"{prefix}{i+1}_{genes[0]}"          # label by top gene, recognizable
        sc.tl.score_genes(adata, gene_list=genes, score_name=name)
        names.append(name)
    return names


def plot_umap(adata, color, out_png, basis="X_umap", title=None, ncols=2, size=None):
    """Plot UMAP colored by one or more obs columns / genes / module scores."""
    color = [color] if isinstance(color, str) else list(color)
    color = [c for c in color if c in adata.obs.columns or c in adata.var_names]
    if not color:
        print("  [viz] nothing valid to color by; skipping", out_png)
        return None
    key = basis.replace("X_", "")
    if basis != "X_umap":
        adata.obsm["X_" + key] = adata.obsm[basis]
    sc.pl.embedding(adata, basis=key, color=color, ncols=ncols, show=False,
                    title=title, size=size, wspace=0.35)
    plt.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close("all")
    print(f"  [viz] wrote {out_png}")
    return out_png


def batch_mixing_score(adata, batch_key, use_rep, n_neighbors=30):
    """Quantify batch mixing: average fraction of each cell's neighbours that come
    from a DIFFERENT batch. Higher = better mixed. Compare before vs after Harmony
    to verify correction numerically, not just by eye."""
    from sklearn.neighbors import NearestNeighbors
    Z = adata.obsm[use_rep]
    b = adata.obs[batch_key].astype(str).values
    nn = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(Z)
    idx = nn.kneighbors(Z, return_distance=False)[:, 1:]
    return float((b[idx] != b[:, None]).mean())
