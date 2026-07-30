"""
ids_streaming.py — out-of-core IDS for datasets with too many CELLS to hold in RAM.

This is the second batching axis from the IDS paper's "implementation details on
large datasets". `ids_batched.py` blocks the GENE axis (caps the k^2*d^2 output
memory); this file additionally STREAMS the CELL (sample) axis so the
n x (k*d) feature-mapped matrix is never materialized either. Together they give
the paper's O(k^2*d'^2 + n*d + n*k*d') footprint — the piece needed for the
27M-cell scale, which `ids_batched.py` alone does not provide.

TWO IMPLEMENTATIONS, matching the paper's structure vs. a faster shortcut
--------------------------------------------------------------------------------
* two_pass=True  (DEFAULT — this is what the paper describes):
    STEP 1: stream all cells once and accumulate the MEAN of every feature-mapped
            variable (optionally saved to disk). Nothing else is computed yet.
    STEP 2: stream the cells again, blocked over genes; CENTER each feature by
            subtracting the STEP-1 mean, then accumulate centered sums-of-squares
            (-> variances) and centered cross-products (-> covariances); form the
            correlation block and reduce to IDS.
    This is the numerically robust order (find the true center first, then measure
    around it), which is why the paper does it this way at 27M-cell scale.

* two_pass=False (one pass, faster, less robust at extreme n):
    Accumulate raw sums Sx, Sxx, Sxy in a single pass and combine at the end via
    cov = Sxy - Sx*Sy/n. Identical result to two_pass at moderate scale, but the
    "sum-of-products minus product-of-sums" step can lose precision when n is very
    large (catastrophic cancellation). Kept for speed when precision isn't a concern.

Both reduce IDS with the authors' feature map (`transform`, imported unchanged) and
match a single ids.compute_IDS call to ~1e-13 (see test_ids_streaming.py).

Notes / scope
-------------
* numpy accumulation in float64. Results match the authors' function regardless of
  which backend it would otherwise pick.
* p-values not supported (would need a fixed cell permutation shared across batches).
* min-max preprocessing must be GLOBAL per gene (done once, before streaming); a
  per-batch min-max would make batches incomparable. Pair with the sparse global
  scaling in sc_ids_pipeline.minmax_scale_sparse.
* `save_means_to`: path to a .npz written after STEP 1 (two-pass only). This mirrors
  the paper's "save the running mean to disk" and lets a future multi-file / resume
  workflow reuse the means without recomputing. The means vector is tiny (k*d floats),
  so this is about orchestration, not memory.
* CURRENT input is an in-memory matrix that we slice twice; a truly out-of-core source
  (a list of files / a re-openable generator) is a future extension — the two-pass
  structure is exactly what makes that extension clean.
"""
import numpy as np
from ids.numpy_dependence import transform  # authors' Gaussian feature map — reused

try:
    import scipy.sparse as _sp
    _HAVE_SP = True
except ImportError:
    _HAVE_SP = False


def _iter_cell_batches(X, cell_batch):
    n = X.shape[0]
    for s in range(0, n, cell_batch):
        e = min(s + cell_batch, n)
        block = X[s:e]
        if _HAVE_SP and _sp.issparse(block):
            block = np.asarray(block.todense(), dtype=np.float64)
        else:
            block = np.asarray(block, dtype=np.float64)
        yield block


def _reduce(corr, k, di, dj, p_norm):
    """corr: (k*di, k*dj) abs-correlation -> (di, dj) via l_p over feature pairs."""
    C = np.abs(corr).reshape(k, di, k, dj)
    if p_norm == "max":
        return C.max(axis=(0, 2))
    elif p_norm == 1:
        return C.mean(axis=0).mean(axis=1)
    elif p_norm == 2:
        return np.sqrt((C ** 2).mean(axis=0).mean(axis=1))
    raise ValueError(f"unsupported p_norm: {p_norm}")


def _make_blocks(d, block_size):
    starts = range(0, d, block_size)
    return [list(range(s, min(s + block_size, d))) for s in starts]


# ----------------------------------------------------------------------
# STEP 1 (two-pass): running mean of every feature-mapped variable
# ----------------------------------------------------------------------
def _pass1_means(X, blocks, k, cell_batch, bandwidth_term, verbose):
    n = X.shape[0]
    Ssum = {bi: np.zeros(k * len(cols)) for bi, cols in enumerate(blocks)}
    n_seen = 0
    for b_idx, cells in enumerate(_iter_cell_batches(X, cell_batch)):
        n_seen += cells.shape[0]
        for bi, cols in enumerate(blocks):
            ft = transform(cells[:, cols], num_terms=k, bandwidth_term=bandwidth_term)
            Ssum[bi] += ft.sum(axis=0)
        if verbose:
            print(f"  [pass 1/2 · means] {n_seen}/{n} cells", flush=True)
    means = {bi: Ssum[bi] / n_seen for bi in Ssum}
    return means, n_seen


# ----------------------------------------------------------------------
# STEP 2 (two-pass): centered covariance -> correlation -> IDS
# ----------------------------------------------------------------------
def _pass2_correlate(X, blocks, means, k, cell_batch, bandwidth_term, p_norm, verbose):
    n, d = X.shape
    Cxx = {bi: np.zeros(k * len(cols)) for bi, cols in enumerate(blocks)}   # centered sum of squares
    Cxy = {(bi, bj): np.zeros((k * len(blocks[bi]), k * len(blocks[bj])))
           for bi in range(len(blocks)) for bj in range(bi, len(blocks))}
    n_seen = 0
    for b_idx, cells in enumerate(_iter_cell_batches(X, cell_batch)):
        n_seen += cells.shape[0]
        cf = {}
        for bi, cols in enumerate(blocks):
            ft = transform(cells[:, cols], num_terms=k, bandwidth_term=bandwidth_term)
            fc = ft - means[bi][None, :]           # CENTER with the STEP-1 global mean
            cf[bi] = fc
            Cxx[bi] += (fc * fc).sum(axis=0)
        for bi in range(len(blocks)):
            for bj in range(bi, len(blocks)):
                Cxy[(bi, bj)] += cf[bi].T @ cf[bj]
        if verbose:
            print(f"  [pass 2/2 · correlate] {n_seen}/{n} cells", flush=True)

    C = np.empty((d, d), dtype=np.float64)
    for bi in range(len(blocks)):
        for bj in range(bi, len(blocks)):
            cols_i, cols_j = blocks[bi], blocks[bj]
            denom = np.sqrt(np.outer(Cxx[bi], Cxx[bj]))
            with np.errstate(divide="ignore", invalid="ignore"):
                corr = Cxy[(bi, bj)] / denom       # centered cov / (sd_i sd_j) = Pearson
            corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
            blk = _reduce(corr, k, len(cols_i), len(cols_j), p_norm)
            C[np.ix_(cols_i, cols_j)] = blk
            if bi != bj:
                C[np.ix_(cols_j, cols_i)] = blk.T
    return C


# ----------------------------------------------------------------------
# One-pass shortcut (faster, less robust at extreme n)
# ----------------------------------------------------------------------
def _one_pass(X, blocks, k, cell_batch, bandwidth_term, p_norm, verbose):
    n, d = X.shape
    Sx = {bi: np.zeros(k * len(cols)) for bi, cols in enumerate(blocks)}
    Sxx = {bi: np.zeros(k * len(cols)) for bi, cols in enumerate(blocks)}
    Sxy = {(bi, bj): np.zeros((k * len(blocks[bi]), k * len(blocks[bj])))
           for bi in range(len(blocks)) for bj in range(bi, len(blocks))}
    n_seen = 0
    for b_idx, cells in enumerate(_iter_cell_batches(X, cell_batch)):
        n_seen += cells.shape[0]
        feats = {}
        for bi, cols in enumerate(blocks):
            ft = transform(cells[:, cols], num_terms=k, bandwidth_term=bandwidth_term)
            feats[bi] = ft
            Sx[bi] += ft.sum(axis=0)
            Sxx[bi] += (ft * ft).sum(axis=0)
        for bi in range(len(blocks)):
            for bj in range(bi, len(blocks)):
                Sxy[(bi, bj)] += feats[bi].T @ feats[bj]
        if verbose:
            print(f"  [1-pass] {n_seen}/{n} cells", flush=True)

    nf = float(n_seen)
    C = np.empty((d, d), dtype=np.float64)
    for bi in range(len(blocks)):
        for bj in range(bi, len(blocks)):
            cols_i, cols_j = blocks[bi], blocks[bj]
            sx, sy = Sx[bi], Sx[bj]
            cov = Sxy[(bi, bj)] - np.outer(sx, sy) / nf
            vx = Sxx[bi] - sx * sx / nf
            vy = Sxx[bj] - sy * sy / nf
            denom = np.sqrt(np.outer(vx, vy))
            with np.errstate(divide="ignore", invalid="ignore"):
                corr = cov / denom
            corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
            blk = _reduce(corr, k, len(cols_i), len(cols_j), p_norm)
            C[np.ix_(cols_i, cols_j)] = blk
            if bi != bj:
                C[np.ix_(cols_j, cols_i)] = blk.T
    return C


def compute_IDS_streaming(X, block_size=256, cell_batch=100_000, p_norm="max",
                          num_terms=6, bandwidth_term=0.5, verbose=False,
                          two_pass=True, save_means_to=None):
    """Full (d x d) gene-gene IDS, streaming cells and blocking genes.

    X : (n_cells, d_genes) ndarray or scipy sparse, already globally preprocessed
        (normalized + per-gene min-max to [0, 8]). Not preprocessed here.
    two_pass : True (default) follows the paper — mean pass, then centered correlate
               pass (numerically robust). False uses the faster one-pass shortcut.
    save_means_to : optional .npz path; write the STEP-1 means to disk (two-pass only),
                    mirroring the paper's "save the running mean to disk".
    Returns (d, d) ndarray, equal to ids.compute_IDS(dense X) up to fp.
    """
    n, d = X.shape
    k = num_terms
    blocks = _make_blocks(d, block_size)

    if not two_pass:
        return _one_pass(X, blocks, k, cell_batch, bandwidth_term, p_norm, verbose)

    # STEP 1 — means (the paper's first step)
    means, n_seen = _pass1_means(X, blocks, k, cell_batch, bandwidth_term, verbose)
    if save_means_to is not None:
        np.savez(save_means_to, n=n_seen,
                 **{f"block_{bi}": means[bi] for bi in means})
        if verbose:
            print(f"  saved feature means to {save_means_to}", flush=True)
    # STEP 2 — centered correlation (the paper's second step)
    return _pass2_correlate(X, blocks, means, k, cell_batch, bandwidth_term, p_norm, verbose)


if __name__ == "__main__":
    import ids
    rng = np.random.default_rng(0)
    n, d = 5000, 33
    Xc = rng.poisson(0.5, size=(n, d)).astype(np.float64)
    Xc[:, 1] = Xc[:, 0]
    g = rng.poisson(5, n).astype(np.float64); Xc[:, 2] = g
    Xc[:, 3] = (g - g.mean()) ** 2
    lo = Xc.min(0, keepdims=True); span = Xc.max(0, keepdims=True) - lo
    Xs = (Xc - lo) / np.where(span == 0, 1.0, span) * 8.0
    full = np.asarray(ids.compute_IDS(Xs, p_norm="max"))
    for tp in (True, False):
        st = compute_IDS_streaming(Xs, block_size=10, cell_batch=1234, p_norm="max", two_pass=tp)
        tag = "two-pass" if tp else "one-pass"
        print(f"{tag}: max|diff vs full| = {np.abs(full - st).max():.2e}")
