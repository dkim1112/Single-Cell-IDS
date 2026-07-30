"""Exactness tests: blocked IDS must equal a single full ids.compute_IDS call,
for dense AND scipy-sparse input, across block sizes and p_norms."""
import numpy as np
import scipy.sparse as sp
import ids
from sc_ids.batched import compute_IDS_blocked


def make_data(n=1000, d=41, seed=0):
    rng = np.random.default_rng(seed)
    # sparse-ish counts: mostly zeros, like scRNA-seq
    X = rng.poisson(0.4, size=(n, d)).astype(np.float64)
    X[:, 1] = X[:, 0]                                  # linear dependence
    g = rng.poisson(5, n).astype(np.float64)
    X[:, 2] = g
    X[:, 3] = (g - g.mean()) ** 2 + rng.normal(0, 0.5, n)  # nonlinear dependence
    # per-gene min-max to [0,8] (the paper's preprocessing)
    lo = X.min(0, keepdims=True); span = X.max(0, keepdims=True) - lo
    Xs = (X - lo) / np.where(span == 0, 1.0, span) * 8.0
    return Xs


def main():
    X = make_data()
    d = X.shape[1]
    density = (X != 0).mean()
    print(f"data: {X.shape}, nonzero fraction={density:.2f}")

    worst = 0.0
    for p_norm in ("max", 1, 2):
        full = np.asarray(ids.compute_IDS(X, p_norm=p_norm))
        for bs in (7, 16, d):  # uneven blocks, even blocks, single block
            # dense input
            blk = compute_IDS_blocked(X, block_size=bs, p_norm=p_norm)
            dd = np.abs(full - blk).max()
            # sparse (CSR) input — same numbers, densified per block internally
            Xsp = sp.csr_matrix(X)
            blk_sp = compute_IDS_blocked(Xsp, block_size=bs, p_norm=p_norm)
            ds = np.abs(full - blk_sp).max()
            worst = max(worst, dd, ds)
            print(f"  p_norm={str(p_norm):>3}  block={bs:>3}  "
                  f"max|diff| dense={dd:.2e}  sparse={ds:.2e}")

    # sanity: known structure recovered (uses full matrix)
    C = np.asarray(ids.compute_IDS(X, p_norm="max"))
    lin = C[0, 1]
    nonlin = C[2, 3]
    pear = abs(np.corrcoef(X[:, 2], X[:, 3])[0, 1])
    print(f"\nrecovery: linear pair IDS={lin:.3f}, "
          f"nonlinear pair IDS={nonlin:.3f} vs |Pearson|={pear:.3f}")

    assert worst < 1e-6, f"blocked IDS diverged from full: {worst:.2e}"
    print(f"\nPASS — blocked == full to {worst:.1e} across dense+sparse, all p_norms")


if __name__ == "__main__":
    main()
