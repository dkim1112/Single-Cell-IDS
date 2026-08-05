"""
Single-cell IDS pipeline — a thin wrapper around the authors' library.

Design principle (per Dr. Patrick's ask): REUSE the repo's code. All dependence
math goes through the authors' `ids.compute_IDS` unchanged — including their
PyTorch backend (auto-selected when torch is installed). This file only adds
the single-cell-specific scaffolding the paper describes but the repo leaves
to the user: load scRNA-seq, preprocess, pick highly-variable genes, run
their function on the gene x gene matrix, and read modules out of the result.

Usage
-----
    # runnable now on built-in synthetic data:
    python sc_ids_pipeline.py

    # on the real GSE194315 data, which lives in this repo under data/
    # (kept local, gitignored). Run from the interdependence_scores/ folder:
    python sc_ids_pipeline.py \
        --data data/GSE194315_PBMC-01-07_processed_data_files \
        --prefix "PBMC-01-1." \
        --min-cells-frac 0.01 --max-mt 15 --min-genes 200 \
        --n-hvg 1000 --min-ids 0.5 --min-size 2 --cluster greedy

    # IMPORTANT: each 10x lane pools several subjects (demuxlet) + ~25% doublets.
    # Join the metadata to drop doublets and restrict to one subject before IDS:
    python sc_ids_pipeline.py \
        --data data/GSE194315_PBMC-01-07_processed_data_files --prefix "PBMC-01-1." \
        --metadata data/GSE194315_CellMetadata-AS_TotalCiteseq_20220711.tsv \
        --included-only --subset-col Subject --subset-value PSA26 \
        --min-cells-frac 0.02 --min-ids 0.5 --min-size 2 --cluster greedy
    # --included-only uses the SOURCE STUDY's own per-cell QC verdict
    # (IncludedInStudy==TRUE) — preferred over --singlets, which it subsumes.
    # NOTE: no --n-hvg -> HVG is OFF; all genes passing --min-cells-frac are run via
    # gene-blocking (auto-enabled). Runtime grows ~quadratically in gene count, so the
    # gene floor doubles as the compute knob: 0.15 (~4k genes) ~2 min; 0.02 (~11k) ~15 min
    # on CPU. metadata columns for --subset-col: Sample, Subject, Status, CellType, Cluster

The gene x gene IDS call is exactly the repo's documented API:
    ids.compute_IDS(X)                       -> (d, d) all-pairs
    ids.compute_IDS(X, p_val=True, ...)      -> (d, d), (d, d) p-values
"""
import argparse
import numpy as np
import scanpy as sc
import anndata as ad

import ids  # the authors' package — do NOT reimplement; call it directly.
from ids_batched import compute_IDS_blocked      # gene-blocked wrapper over ids.compute_IDS
from ids_streaming import compute_IDS_streaming  # cell-streaming + gene-blocking (out-of-core)

try:
    import scipy.sparse as sp
    HAVE_SP = True
except ImportError:
    HAVE_SP = False

try:
    import torch
    HAVE_TORCH = True
except ImportError:
    HAVE_TORCH = False


def minmax_scale_sparse(X, scale_hi=8.0):
    """Per-gene min-max to [0, scale_hi] that PRESERVES sparsity.

    For scRNA-seq, essentially every gene has at least one zero-count cell, so
    the per-gene minimum is 0 and min-max reduces to pure per-column scaling
    (x -> x * scale_hi / span), which keeps the matrix sparse. If some gene has
    a strictly positive minimum, subtracting it would turn stored zeros into
    -lo (destroying sparsity for that column), so in that case we densify and
    scale exactly, printing how many genes triggered it.
    Returns (X_scaled, was_sparse: bool)."""
    lo = np.asarray(X.min(axis=0).todense()).ravel()
    hi = np.asarray(X.max(axis=0).todense()).ravel()
    span = np.where((hi - lo) == 0, 1.0, hi - lo)
    if np.allclose(lo, 0.0):
        scale = (scale_hi / span).astype(np.float64)
        Xs = X.multiply(scale[None, :]).tocsr()      # stays sparse
        return Xs, True
    n_pos = int((lo > 0).sum())
    print(f"  [minmax] {n_pos} gene(s) have a nonzero minimum after "
          f"normalization; densifying to scale exactly (sparsity not preserved).")
    Xd = np.asarray(X.todense(), dtype=np.float64)
    Xd = (Xd - lo[None, :]) / span[None, :] * scale_hi
    return Xd, False


# ----------------------------------------------------------------------
# 1. Data
# ----------------------------------------------------------------------
def load_data(path, transpose=False, prefix=None):
    """Load a real scRNA-seq matrix as AnnData (rows = cells, columns = genes).

    Accepts the two formats you're most likely to be handed:
      * .h5ad                -> AnnData, the standard single-cell format (scanpy)
      * .csv / .tsv / .txt   -> a plain table: first column = cell IDs,
                                header row = gene names, values = counts.
    A 10x folder (matrix.mtx + barcodes + features) is also read if `path` is
    that folder. Use --transpose if your table is genes-as-rows.

    `prefix` is passed to scanpy's read_10x_mtx for GEO datasets whose 10x files
    are prefixed per-sample, e.g. GSE194315 files 'PBMC-01-1.matrix.mtx.gz',
    'PBMC-01-1.barcodes.tsv.gz', 'PBMC-01-1.features.tsv.gz' -> prefix='PBMC-01-1.'
    (keep the trailing dot; it must match the real filenames exactly).
    """
    import os
    if os.path.isdir(path):
        adata = sc.read_10x_mtx(path, prefix=prefix)       # 10x Cell Ranger output
    elif path.endswith(".h5ad"):
        adata = sc.read_h5ad(path)
    elif path.endswith((".csv", ".tsv", ".txt")):
        import pandas as pd
        sep = "," if path.endswith(".csv") else "\t"
        df = pd.read_csv(path, sep=sep, index_col=0)
        adata = ad.AnnData(df.astype("float32"))
        adata.obs_names = df.index.astype(str)
        adata.var_names = df.columns.astype(str)
    else:
        raise ValueError(f"Unsupported file type: {path} (use .h5ad, .csv/.tsv, or a 10x folder)")
    if transpose:
        adata = adata.T
    return adata


def synthetic_adata(n_cells=3000, n_genes=400, seed=0):
    """Synthetic scRNA-seq with two planted co-expression modules + a
    nonlinear gene pair, so the pipeline is end-to-end runnable and its
    output is checkable without real data.

    Design note: modules are kept small relative to n_genes so that a module
    switching on does not dominate a cell's total counts — otherwise standard
    library-size normalization would partly cancel the very co-expression we
    plant (a genuine scRNA-seq confounder, not just a synthetic artifact)."""
    rng = np.random.default_rng(seed)
    counts = rng.poisson(lam=rng.uniform(1.0, 5.0, n_genes), size=(n_cells, n_genes)).astype(np.float32)

    # Module A: genes 0-7 strongly co-expressed in a ~25% subpopulation (bimodal on/off).
    onA = rng.random(n_cells) < 0.25
    counts[np.ix_(onA, range(8))] += rng.poisson(30, size=(onA.sum(), 8))
    # Module B: genes 8-14 co-expressed in a different ~15% subpopulation.
    onB = rng.random(n_cells) < 0.15
    counts[np.ix_(onB, range(8, 15))] += rng.poisson(35, size=(onB.sum(), 7))
    # Nonlinear pair: gene 21 ~ (gene 20)^2, made high-variance so it survives HVG.
    g20 = rng.poisson(6, n_cells).astype(np.float32)
    counts[:, 20] = g20
    counts[:, 21] = (g20 - g20.mean()) ** 2 + rng.normal(0, 1, n_cells)

    A = ad.AnnData(np.clip(counts, 0, None))
    A.var_names = [f"gene{i}" for i in range(n_genes)]
    A.obs_names = [f"cell{i}" for i in range(n_cells)]
    A.uns["planted"] = {"moduleA": list(range(8)), "moduleB": list(range(8, 15)),
                        "nonlinear_pair": (20, 21)}
    return A


# ----------------------------------------------------------------------
# 2. Preprocessing (paper's recipe)
# ----------------------------------------------------------------------
def join_metadata(adata, path, barcode_col="CellName", sample=None,
                  singlets_only=False, singlet_col="DemuxletDropletType",
                  singlet_value="SNG", included_only=False,
                  included_col="IncludedInStudy", subset_col=None,
                  subset_value=None, verbose=True):
    """Join per-cell metadata (barcode -> Subject / Status / CellType / ...) onto
    adata, then optionally keep singlets only and/or subset to one group.

    WHY THIS MATTERS for GSE194315: each 10x lane (e.g. PBMC-01-1) is NOT one
    sample — it pools several genetically-distinct subjects that were demultiplexed
    with demuxlet, and ~25% of droplets are doublets. Running IDS on the raw pool
    therefore mixes donors and includes doublets, so apparent 'co-expression' can be
    driven by donor identity or doublet artifacts rather than biology. This join
    lets you (a) drop doublets (--singlets) and (b) restrict to one Subject / one
    CellType / etc. before computing IDS.

    BARCODE RECONCILIATION: scanpy's read_10x_mtx yields barcodes like
    'AAACCCAAGACCATAA-1'; this metadata's CellName is 'PBMC-01-1_AAACCCAAGACCATAA'
    (sample-prefixed, no gem suffix). We strip the '-<n>' suffix, prepend the
    sample name (taken from --prefix, e.g. 'PBMC-01-1.' -> 'PBMC-01-1'), and match
    that against CellName — which is globally unique, so the match is unambiguous."""
    import pandas as pd
    meta = pd.read_csv(path, sep="\t", low_memory=False)
    if barcode_col not in meta.columns:
        barcode_col = meta.columns[0]
        if verbose:
            print(f"  [metadata] '{barcode_col}' assumed as barcode column")

    bare = pd.Index(adata.obs_names).astype(str).str.replace(r"-\d+$", "", regex=True)
    if sample:
        keys = pd.Index(sample + "_" + bare)
    else:
        keys = pd.Index(bare)
        if verbose:
            print("  [metadata] no sample (from --prefix/--metadata-sample); "
                  "matching on bare barcode — may be ambiguous.")

    lut = meta.drop_duplicates(barcode_col).set_index(barcode_col)
    matched = keys.isin(lut.index)
    n = int(matched.sum())
    if verbose:
        print(f"  [metadata] matched {n}/{adata.n_obs} cells "
              f"({100 * n / max(adata.n_obs,1):.1f}%)")
    if n == 0:
        raise SystemExit("  [metadata] 0 cells matched — check --prefix / "
                         "--metadata-sample (CellName looks like 'PBMC-01-1_AAAC...').")

    aligned = lut.reindex(keys)
    for col in meta.columns:
        if col != barcode_col:
            adata.obs[col] = aligned[col].values
    adata = adata[matched].copy()                 # drop cells absent from metadata

    # The SOURCE STUDY'S OWN QC VERDICT. GSE194315's metadata carries a per-cell
    # IncludedInStudy flag: the cells the original authors kept after their full QC
    # pipeline. It is STRICTLY STRONGER than our singlet filter — on PBMC-01-1 it
    # keeps 17,858 cells vs 20,517 for singlets-only, i.e. it additionally drops
    # 2,659 singlets that failed their QC for other reasons (and CellType is
    # annotated only for included cells). Preferring this over hand-tuned cutoffs
    # means we inherit QC done by the people who generated the data.
    if included_only and included_col in adata.obs:
        b = adata.n_obs
        keep = adata.obs[included_col].astype(str).str.upper().isin(["TRUE", "T", "1"])
        adata = adata[keep.values].copy()
        if verbose:
            print(f"  [metadata] study QC only ({included_col}==TRUE): {b} -> {adata.n_obs}")

    if singlets_only and singlet_col in adata.obs:
        b = adata.n_obs
        adata = adata[adata.obs[singlet_col].astype(str) == singlet_value].copy()
        if verbose:
            print(f"  [metadata] singlets only ({singlet_col}=={singlet_value}): {b} -> {adata.n_obs}")

    if subset_col is not None and subset_value is not None:
        if subset_col not in adata.obs:
            raise SystemExit(f"  [metadata] --subset-col '{subset_col}' not found; "
                             f"available: {list(meta.columns)}")
        b = adata.n_obs
        adata = adata[adata.obs[subset_col].astype(str) == str(subset_value)].copy()
        if verbose:
            print(f"  [metadata] subset {subset_col}=={subset_value}: {b} -> {adata.n_obs}")
    return adata


def qc_filter(adata, min_cells_frac=0.0, max_mt=None, min_genes=None, verbose=True):
    """Standard single-cell quality control, applied to raw counts BEFORE
    normalization. All three are opt-in (defaults are no-ops).

      * min_genes    : drop cells expressing fewer than this many genes
                       (empty droplets / very low-quality cells).
      * max_mt       : drop cells whose mitochondrial-gene % exceeds this
                       (dying/stressed cells; MT-high). e.g. 10-15.
      * min_cells_frac: GENE FLOOR — drop genes detected in fewer than this
                       fraction of cells. Removes ultra-rare genes whose only
                       'co-expression' is co-detection noise (a gene seen in 10
                       of 27,000 cells can score a high IDS by chance). e.g. 0.01.

    Neither the authors' repo nor our earlier pipeline did any of this; on real
    PBMC data its absence let lowly-expressed lincRNAs and RBC/mito programs
    dominate the modules. Returns the filtered AnnData (a copy)."""
    adata = adata.copy()
    n_cells0, n_genes0 = adata.n_obs, adata.n_vars
    if min_genes is not None:
        sc.pp.filter_cells(adata, min_genes=min_genes)
    if max_mt is not None:
        adata.var["mt"] = adata.var_names.str.upper().str.startswith(("MT-", "MT."))
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                                   log1p=False, inplace=True)
        adata = adata[adata.obs["pct_counts_mt"] <= max_mt].copy()
    if min_cells_frac and min_cells_frac > 0:
        sc.pp.filter_genes(adata, min_cells=int(min_cells_frac * adata.n_obs))
    if verbose:
        print(f"  QC: cells {n_cells0} -> {adata.n_obs}, genes {n_genes0} -> {adata.n_vars}"
              f"  (min_genes={min_genes}, max_mt={max_mt}, gene_floor={min_cells_frac})")
    return adata


def preprocess(adata, n_hvg=500, target_sum=1e4, scale_hi=8.0, log1p=False,
               keep_sparse=False):
    """Library-size normalization + HVG selection, then the paper's per-gene
    min-max scaling to [0, scale_hi]. Returns (X, gene_names).

    HVG selection is what keeps IDS tractable: the authors' function forms a
    (k*d) x (k*d) correlation matrix, so d must be limited to hundreds/low
    thousands of genes, not all ~20k. If your data is already small (a few
    hundred genes), pass n_hvg >= n_vars to skip HVG entirely.

    CAVEAT worth knowing: scanpy's default HVG uses a mean-variance model, and a
    co-expression program that is OFF in most cells has low mean, so its (very
    real) bimodal signal may NOT rank as "highly variable" -> the genes you most
    want can get dropped before IDS ever sees them. On real data, sanity-check
    that known marker genes survive HVG, raise n_hvg, or select genes another way.

    IMPORTANT — log1p defaults to OFF, and that is deliberate. The paper's
    recipe is normalize -> min-max to [0,8]. Adding the usual log1p step BEFORE
    the [0,8] min-max collapses IDS: log spreads counts evenly across [0,8],
    pushing many cells into the range where IDS's Gaussian feature envelope
    exp(-B*x^2) saturates to ~0 (B=0.5 => exp(-0.5*8^2)~0). Empirically this
    drops a strong module's mean IDS from ~0.89 to ~0.09. Raw/normalized counts
    are right-skewed (most mass near 0, where the features are informative), so
    the [0,8] scaling is calibrated for them, not for log-space data. HVG
    selection is still done on log space (standard) when log1p=True is requested.
    """
    adata = adata.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)

    if n_hvg is not None and n_hvg < adata.n_vars:
        # HVG ranking wants log space; compute on a log copy, subset the original.
        logged = adata.copy(); sc.pp.log1p(logged)
        sc.pp.highly_variable_genes(logged, n_top_genes=n_hvg)
        adata = adata[:, logged.var["highly_variable"].values].copy()

    if log1p:
        sc.pp.log1p(adata)

    # Sparse path (for the blocked/large-scale route): scale WITHOUT densifying
    # the whole matrix, so peak memory stays at one gene block downstream.
    if keep_sparse and HAVE_SP and sp.issparse(adata.X) and not log1p:
        Xs, still_sparse = minmax_scale_sparse(adata.X, scale_hi=scale_hi)
        return Xs, list(adata.var_names)

    X = adata.X
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    X = X.astype(np.float32)

    # per-gene min-max -> [0, scale_hi]  (paper preprocessing)
    lo = X.min(0, keepdims=True)
    span = X.max(0, keepdims=True) - lo
    X = (X - lo) / np.where(span == 0, 1.0, span) * scale_hi
    return X, list(adata.var_names)


# ----------------------------------------------------------------------
# 3. IDS — calls the authors' function directly
# ----------------------------------------------------------------------
def gene_gene_ids(X, p_norm="max", p_val=False, num_tests=1000,
                  num_terms=6, bandwidth_term=0.5, use_torch=True):
    """Gene x gene IDS via the authors' `ids.compute_IDS`.

    Passing a torch tensor routes through their torch backend (GPU if avail);
    passing ndarray routes through their numpy backend. Either way the math is
    theirs, unchanged.
    """
    if use_torch and HAVE_TORCH:
        device = "cuda" if torch.cuda.is_available() else "cpu"
        X = torch.tensor(X, dtype=torch.float64, device=device)

    out = ids.compute_IDS(X, num_terms=num_terms, p_norm=p_norm,
                          p_val=p_val, num_tests=num_tests,
                          bandwidth_term=bandwidth_term)

    def to_np(v):
        return v.detach().cpu().numpy() if HAVE_TORCH and torch.is_tensor(v) else np.asarray(v)

    if p_val:
        C, P = out
        return to_np(C), to_np(P)
    return to_np(out)


# ----------------------------------------------------------------------
# 4. Read co-expression modules out of the IDS matrix
# ----------------------------------------------------------------------
def extract_modules(C, gene_names, P=None, alpha=0.01, min_ids=0.3, min_size=3,
                    method="components"):
    """Discover co-expression modules from a thresholded gene-gene IDS graph
    (bottom-up: no cell-type labels imposed). Most genes are background and
    belong to NO module, so we keep only strong+significant edges rather than
    forcing every gene into a cluster.

    method:
      * "components" — connected components of the thresholded graph. Simple, but
        chains of weak links can fuse distinct programs into one giant component.
      * "greedy" — greedy modularity community detection (networkx) on the same
        edges, weighted by IDS. Splits weakly-linked chains, so it recovers
        coherent modules instead of one blob. Recommended for real data.

    Edge kept between genes i,j if BOTH:
      - IDS[i,j] >= min_ids                 (effect size), and
      - permutation p-value P[i,j] < alpha  (significance, when computed).

    Why an effect-size floor and not just significance: with thousands of cells,
    library-size normalization induces weak but *statistically significant*
    dependencies across many genes, so a p-value filter alone lets background
    edges bridge real modules into one giant component. The IDS magnitude gate
    removes those. Choose min_ids from the data: histogram the off-diagonal IDS
    (real modules usually sit in a high mode separated by a gap from the
    background bulk) or set it at the max IDS seen under permutation. Here the
    background bulk tops out ~0.14 and modules sit ~0.9, so 0.3 sits in the gap.
    """
    d = C.shape[0]
    off = ~np.eye(d, dtype=bool)
    edge = (C >= min_ids) & off
    if P is not None:
        edge &= (P < alpha)
    edge = edge & edge.T  # keep mutual/symmetric edges

    if method == "greedy":
        try:
            import networkx as nx
            from networkx.algorithms.community import greedy_modularity_communities
        except ImportError:
            print("  [extract_modules] networkx not installed; falling back to "
                  "connected components. `pip install networkx` for method='greedy'.")
            method = "components"

    def build(members):
        members = np.asarray(sorted(members))
        sub = C[np.ix_(members, members)].copy(); np.fill_diagonal(sub, np.nan)
        return {"genes": [gene_names[i] for i in members],
                "idx": list(members), "mean_ids": float(np.nanmean(sub))}

    modules = []
    if method == "greedy":
        G = nx.Graph()
        G.add_nodes_from(range(d))
        ii, jj = np.where(np.triu(edge, k=1))
        for i, j in zip(ii.tolist(), jj.tolist()):
            G.add_edge(i, j, weight=float(C[i, j]))
        for comm in greedy_modularity_communities(G, weight="weight"):
            if len(comm) >= min_size:
                modules.append(build(comm))
    else:
        from scipy.sparse import csr_matrix
        from scipy.sparse.csgraph import connected_components
        n_comp, labels = connected_components(csr_matrix(edge), directed=False)
        for lab in range(n_comp):
            members = np.where(labels == lab)[0]
            if len(members) >= min_size:
                modules.append(build(members))
    return sorted(modules, key=lambda m: -m["mean_ids"])


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=None, help="path to .h5ad, .csv/.tsv, or 10x folder; omit for synthetic demo")
    ap.add_argument("--transpose", action="store_true", help="use if your table has genes as rows, cells as columns")
    ap.add_argument("--prefix", default=None, help="filename prefix for a 10x folder (e.g. 'PBMC-01-1.' for GSE194315); keep the trailing dot")
    ap.add_argument("--n-hvg", type=int, default=None,
                    help="highly-variable-gene count. DEFAULT None = NO HVG: keep all genes that "
                         "pass the --min-cells-frac floor and run them via gene-blocking (auto-enabled "
                         "when many genes remain). Set a number only to force HVG down-selection.")
    ap.add_argument("--p-norm", default="max")
    ap.add_argument("--p-val", action="store_true")
    ap.add_argument("--num-tests", type=int, default=200)
    ap.add_argument("--min-ids", type=float, default=0.3, help="effect-size threshold for calling a module edge; set from your data's IDS histogram (gap between background bulk and module mode)")
    ap.add_argument("--min-size", type=int, default=3, help="minimum genes for a reported module (use 2 to keep strong gene-pair programs)")
    ap.add_argument("--cluster", choices=["components", "greedy"], default="components",
                    help="module extraction: 'components' (connected components; can over-merge) "
                         "or 'greedy' (weighted modularity clustering; recommended for real data)")
    ap.add_argument("--log1p", action="store_true", help="apply log1p before min-max (usually hurts IDS; see preprocess docstring)")
    # --- single-cell QC (opt-in; applied to raw counts before normalization) ---
    ap.add_argument("--min-cells-frac", type=float, default=0.0,
                    help="GENE FLOOR: drop genes detected in fewer than this fraction of cells "
                         "(e.g. 0.01). Removes ultra-rare genes that produce co-detection-noise modules.")
    ap.add_argument("--max-mt", type=float, default=None,
                    help="CELL QC: drop cells with mitochondrial-gene %% above this (e.g. 15) — dying/stressed cells.")
    ap.add_argument("--min-genes", type=int, default=None,
                    help="CELL QC: drop cells expressing fewer than this many genes (e.g. 200) — empty droplets/low quality.")
    # --- metadata join / demultiplex / subset (opt-in via --metadata) ---
    ap.add_argument("--metadata", default=None,
                    help="path to CellMetadata TSV; joins barcode -> Subject/Status/CellType/... "
                         "and enables --singlets / --subset-col / --subset-value.")
    ap.add_argument("--metadata-barcode-col", default="CellName",
                    help="metadata column holding the cell barcode (default CellName).")
    ap.add_argument("--metadata-sample", default=None,
                    help="sample name used to match barcodes; default derived from --prefix "
                         "(e.g. 'PBMC-01-1.' -> 'PBMC-01-1').")
    ap.add_argument("--singlets", action="store_true",
                    help="with --metadata, keep only demultiplexed singlets "
                         "(DemuxletDropletType==SNG) — drops doublets/ambiguous droplets.")
    ap.add_argument("--included-only", action="store_true",
                    help="with --metadata, keep only cells the SOURCE STUDY kept after its own "
                         "QC (IncludedInStudy==TRUE). Stricter than --singlets (it also drops "
                         "singlets that failed their QC) and RECOMMENDED as the primary QC gate.")
    ap.add_argument("--subset-col", default=None,
                    help="with --metadata, metadata column to subset on (e.g. Subject, CellType).")
    ap.add_argument("--subset-value", default=None,
                    help="value of --subset-col to keep (e.g. PSA26, 'CD14 Mono').")
    ap.add_argument("--block-size", type=int, default=None,
                    help="genes per block for the memory-bounded batched IDS "
                         "(paper's d'). When set, IDS is assembled from gene "
                         "blocks and the matrix is densified one block at a "
                         "time (keeps sparse input sparse). Lets you skip HVG "
                         "and run all/many genes. p-values are not supported "
                         "with blocking.")
    ap.add_argument("--cell-batch", type=int, default=None,
                    help="cells per batch for the out-of-core STREAMING IDS. "
                         "When set, IDS is computed by streaming over cell "
                         "batches (feature matrix never fully materialized) and "
                         "blocking genes by --block-size (default 256). This is "
                         "the path for datasets with too many cells to hold in "
                         "RAM. p-values not supported.")
    ap.add_argument("--one-pass", action="store_true",
                    help="with --cell-batch, use the faster single-pass streaming "
                         "instead of the paper's numerically-robust two-pass "
                         "(mean pass, then centered correlate pass). Two-pass is "
                         "the default; use --one-pass only when speed matters more "
                         "than precision at very large cell counts.")
    ap.add_argument("--save-means", default=None,
                    help="with --cell-batch (two-pass), write the step-1 feature "
                         "means to this .npz (mirrors the paper's 'save running "
                         "mean to disk'; useful for future multi-file/resume runs).")
    args = ap.parse_args()

    # compute_IDS wants int 1/2 for IDS-1/IDS-2 ('max' stays a string);
    # argparse hands us strings, which would silently miss those branches.
    if args.p_norm in ("1", "2"):
        args.p_norm = int(args.p_norm)

    streaming = args.cell_batch is not None
    block_size = args.block_size
    blocked = block_size is not None
    if (blocked or streaming) and args.p_val:
        ap.error("--p-val is not supported with --block-size/--cell-batch "
                 "(consistent permutations across blocks are a separate "
                 "change). Drop one.")

    adata = load_data(args.data, transpose=args.transpose, prefix=args.prefix) if args.data else synthetic_adata()
    print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes"
          f"{' (synthetic)' if args.data is None else ''}")

    # metadata join / demultiplex / subset (before QC), if requested
    if args.metadata:
        sample = args.metadata_sample or (args.prefix.rstrip(".") if args.prefix else None)
        adata = join_metadata(adata, args.metadata,
                              barcode_col=args.metadata_barcode_col, sample=sample,
                              singlets_only=args.singlets,
                              included_only=args.included_only,
                              subset_col=args.subset_col, subset_value=args.subset_value)
        print(f"After metadata join/subset: {adata.n_obs} cells x {adata.n_vars} genes")

    # single-cell QC (no-op unless the flags are set)
    if args.min_cells_frac > 0 or args.max_mt is not None or args.min_genes is not None:
        adata = qc_filter(adata, min_cells_frac=args.min_cells_frac,
                          max_mt=args.max_mt, min_genes=args.min_genes)

    # Drop-HVG default: with no HVG and many genes, run ALL genes via gene-blocking
    # (the dense direct path can't build a (6d)x(6d) matrix for thousands of genes).
    # Pair with --min-cells-frac to trim rare genes first and keep this tractable.
    if not streaming and block_size is None and args.n_hvg is None and adata.n_vars > 2000:
        block_size = 1024
        blocked = True
        print(f"  [auto] no HVG and {adata.n_vars} genes -> computing IDS over ALL "
              f"genes via gene-blocking (block-size {block_size}).")

    X, genes = preprocess(adata, n_hvg=args.n_hvg, log1p=args.log1p,
                          keep_sparse=blocked or streaming)
    is_sparse = HAVE_SP and sp.issparse(X)
    print(f"After preprocessing: {X.shape[0]} cells x {X.shape[1]} genes "
          f"({'sparse' if is_sparse else 'dense'}; "
          f"backend: {'torch' if HAVE_TORCH else 'numpy'})")

    if streaming:
        bs = block_size or 256
        mode = "one-pass" if args.one_pass else "two-pass (paper: mean pass + centered correlate pass)"
        print(f"Computing IDS by streaming cells (batch {args.cell_batch}) "
              f"and blocking genes ({bs}); {mode}; peak memory ~ k^2*{bs}^2 + "
              f"{args.cell_batch}*k*{bs}")
        C = compute_IDS_streaming(X, block_size=bs, cell_batch=args.cell_batch,
                                  p_norm=args.p_norm, verbose=True,
                                  two_pass=not args.one_pass,
                                  save_means_to=args.save_means)
        P = None
    elif blocked:
        print(f"Computing IDS in gene blocks of {block_size} "
              f"(peak memory ~ k^2 * block_size^2)")
        C = compute_IDS_blocked(X, block_size=block_size,
                                p_norm=args.p_norm, verbose=True)
        P = None
    else:
        # unblocked path needs a dense matrix; densify if we kept it sparse
        if is_sparse:
            X = np.asarray(X.todense(), dtype=np.float32)
        res = gene_gene_ids(X, p_norm=args.p_norm, p_val=args.p_val,
                            num_tests=args.num_tests)
        C, P = res if args.p_val else (res, None)
    print(f"IDS matrix: {C.shape}, range [{C.min():.3f}, {np.nanmax(C[~np.eye(len(C),dtype=bool)]):.3f}]")

    modules = extract_modules(C, genes, P=P, min_ids=args.min_ids,
                              min_size=args.min_size, method=args.cluster)
    print(f"\nDiscovered {len(modules)} co-expression modules "
          f"(method={args.cluster}, min_ids={args.min_ids}, min_size={args.min_size}; top 5):")
    for m in modules[:5]:
        preview = ", ".join(m["genes"][:8]) + ("…" if len(m["genes"]) > 8 else "")
        print(f"  mean IDS {m['mean_ids']:.3f}  ({len(m['genes'])} genes): {preview}")

    # If synthetic, verify the planted structure was recovered.
    if "planted" in adata.uns:
        pl = adata.uns["planted"]
        name_to_idx = {g: i for i, g in enumerate(genes)}
        def block_mean(orig_idx):
            kept = [name_to_idx[f"gene{i}"] for i in orig_idx if f"gene{i}" in name_to_idx]
            if len(kept) < 2:
                return None, len(kept)
            sub = C[np.ix_(kept, kept)].copy(); np.fill_diagonal(sub, np.nan)
            return float(np.nanmean(sub)), len(kept)
        print("\nRecovery check vs planted ground truth:")
        for key in ("moduleA", "moduleB"):
            mean, k = block_mean(pl[key]);
            print(f"  {key}: mean IDS {mean:.3f} over {k} surviving genes" if mean else f"  {key}: too few HVGs survived")
        gi, gj = pl["nonlinear_pair"]
        gi, gj = f"gene{gi}", f"gene{gj}"
        if gi in name_to_idx and gj in name_to_idx:
            i, j = name_to_idx[gi], name_to_idx[gj]
            def col(k):  # dense 1-D column whether X is sparse or dense
                c = X[:, k]
                return np.asarray(c.todense()).ravel() if HAVE_SP and sp.issparse(X) else np.asarray(c).ravel()
            lin = abs(np.corrcoef(col(i), col(j))[0, 1])
            print(f"  nonlinear pair: IDS {C[i,j]:.3f} vs |Pearson| {lin:.3f}")
            if P is not None:
                print(f"                  permutation p-value {P[i,j]:.3f}")


if __name__ == "__main__":
    main()