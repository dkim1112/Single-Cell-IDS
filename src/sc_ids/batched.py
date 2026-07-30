"""
ids_batched.py — memory-bounded IDS over the GENE (variable) axis.

Implements the batching described in the IDS paper's supplementary section
"IDS implementation details on large datasets": batch the computation over the
number of variables d (block size d' << d) so you never materialize the full
n x (k*d) feature-mapped matrix or the (k*d) x (k*d) correlation matrix.

DESIGN PRINCIPLE (unchanged from the rest of this project): we NEVER reimplement
IDS. Every number is produced by the authors' `ids.compute_IDS`, called on
gene blocks. The full d x d IDS matrix is assembled from d'-sized blocks:

    diagonal block (i, i):     compute_IDS(X[:, block_i])
    off-diagonal   (i, j):     compute_IDS(X[:, block_i], Y=X[:, block_j])
    lower triangle (j, i):     transpose of (i, j)   [IDS is symmetric]

Why this is the right primitive: `compute_IDS(X, Y)` returns exactly the cross
block of the correlation tensor between two variable sets, reduced by the same
l_p norm. So assembling blocks reproduces a single full `compute_IDS(X)` call to
floating-point precision (validated in test_ids_batched.py), while peak feature
memory is O(k^2 * d'^2 + n*k*d') instead of O(k^2 * d^2 + n*k*d). That is what
makes "all genes" (not just a few hundred HVGs) tractable.

SPARSE HANDLING: real scRNA-seq is a scipy sparse matrix. The authors' code is
dense-only (its x**2 step reads a non-square sparse matrix as a matrix power and
raises "sparse matrix is not square"), and the Gaussian feature map fills in the
zeros anyway (exp(-B*0^2)*0^0 = 1), so IDS cannot stay sparse internally. We keep
X sparse in storage and densify ONE gene block at a time — the whole matrix is
never densified.

NOT handled here: permutation p-values (p_val=True). Consistent p-values under
blocking require permuting the SAME cell ordering across every block; that is a
separate change. For p-values, use the unblocked path, or ask for the
consistent-permutation version.
"""
import numpy as np
import ids
import ids.dependence as _dep

try:
    import scipy.sparse as _sp
    _HAVE_SP = True
except ImportError:
    _HAVE_SP = False

_BACKEND = getattr(_dep, "backend", "numpy")
if _BACKEND == "torch":
    import torch


def _to_backend(block):
    """Return `block` (dense n x d' float array) as the type the active
    compute_IDS backend expects: a torch tensor under the torch backend
    (the authors' torch path rejects ndarrays), else a float64 ndarray."""
    block = np.ascontiguousarray(block, dtype=np.float64)
    if _BACKEND == "torch":
        device = "cuda" if torch.cuda.is_available() else "cpu"
        return torch.tensor(block, dtype=torch.float64, device=device)
    return block


def _to_numpy(v):
    if _BACKEND == "torch" and torch.is_tensor(v):
        return v.detach().cpu().numpy()
    return np.asarray(v)


def _densify_cols(X, cols):
    """Columns `cols` of X (sparse or dense) as a dense float64 ndarray."""
    if _HAVE_SP and _sp.issparse(X):
        return np.asarray(X[:, cols].todense(), dtype=np.float64)
    return np.asarray(X[:, cols], dtype=np.float64)


def compute_IDS_blocked(X, block_size=256, p_norm="max", num_terms=6,
                        bandwidth_term=0.5, verbose=False):
    """Full (d x d) gene-gene IDS matrix, computed in gene blocks.

    Parameters
    ----------
    X : ndarray or scipy.sparse matrix, shape (n_cells, d_genes)
        Already preprocessed (normalized + per-gene min-max to [0, 8]); this
        function does NOT preprocess. Sparse is fine and preferred.
    block_size : int
        Number of genes per block (the paper's d'). Peak memory grows ~quadratically
        in this; correctness is independent of it. 256-1024 is a sane range.
    p_norm, num_terms, bandwidth_term : passed straight to ids.compute_IDS.

    Returns
    -------
    C : ndarray, shape (d, d)
        Identical to ids.compute_IDS(dense X) up to floating point.
    """
    n, d = X.shape
    if block_size >= d:
        # single block — just the authors' function, densifying if needed
        Xf = _to_backend(_densify_cols(X, list(range(d))))
        return _to_numpy(ids.compute_IDS(Xf, num_terms=num_terms, p_norm=p_norm,
                                         bandwidth_term=bandwidth_term))

    starts = range(0, d, block_size)
    blocks = [list(range(s, min(s + block_size, d))) for s in starts]
    C = np.empty((d, d), dtype=np.float64)

    for bi, cols_i in enumerate(blocks):
        Xi = _to_backend(_densify_cols(X, cols_i))
        Cii = _to_numpy(ids.compute_IDS(Xi, num_terms=num_terms, p_norm=p_norm,
                                        bandwidth_term=bandwidth_term))
        C[np.ix_(cols_i, cols_i)] = Cii
        for bj in range(bi + 1, len(blocks)):
            cols_j = blocks[bj]
            Xj = _to_backend(_densify_cols(X, cols_j))
            Cij = _to_numpy(ids.compute_IDS(Xi, Y=Xj, num_terms=num_terms,
                                            p_norm=p_norm,
                                            bandwidth_term=bandwidth_term))
            C[np.ix_(cols_i, cols_j)] = Cij
            C[np.ix_(cols_j, cols_i)] = Cij.T
        if verbose:
            print(f"  gene block {bi + 1}/{len(blocks)} "
                  f"({len(cols_i)} genes) done", flush=True)
    return C


if __name__ == "__main__":
    # quick self-check on synthetic data (numpy backend)
    rng = np.random.default_rng(0)
    n, d = 800, 37
    Xc = rng.poisson(3.0, size=(n, d)).astype(np.float64)
    Xc[:, 1] = Xc[:, 0]                      # linear dep
    Xc[:, 3] = (Xc[:, 2] - Xc[:, 2].mean()) ** 2  # nonlinear dep
    lo = Xc.min(0, keepdims=True); span = Xc.max(0, keepdims=True) - lo
    Xs = (Xc - lo) / np.where(span == 0, 1.0, span) * 8.0

    full = np.asarray(ids.compute_IDS(Xs, p_norm="max"))
    for bs in (8, 16, 100):
        blk = compute_IDS_blocked(Xs, block_size=bs, p_norm="max")
        print(f"block_size={bs:3d}  max|diff|={np.abs(full - blk).max():.2e}")
