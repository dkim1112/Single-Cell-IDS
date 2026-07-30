"""Exactness tests for cell-streaming + gene-blocking IDS (both the paper's
two-pass form and the one-pass shortcut): each must equal a single full
ids.compute_IDS call for dense AND sparse input, across p_norms, cell-batch
sizes, and gene-block sizes. Plus a small illustration of WHY two-pass is the
numerically safer order at large n."""
import numpy as np
import scipy.sparse as sp
import ids
from sc_ids.streaming import compute_IDS_streaming


def make_data(n=4000, d=37, seed=3):
    rng = np.random.default_rng(seed)
    X = rng.poisson(0.4, size=(n, d)).astype(np.float64)   # sparse-like counts
    X[:, 1] = X[:, 0]
    g = rng.poisson(6, n).astype(np.float64); X[:, 2] = g
    X[:, 3] = (g - g.mean()) ** 2 + rng.normal(0, 0.5, n)
    lo = X.min(0, keepdims=True); span = X.max(0, keepdims=True) - lo
    return (X - lo) / np.where(span == 0, 1.0, span) * 8.0


def test_exactness():
    X = make_data()
    n, d = X.shape
    print(f"data {X.shape}, nonzero frac {(X != 0).mean():.2f}")
    worst = 0.0
    for two_pass in (True, False):
        for p_norm in ("max", 1, 2):
            full = np.asarray(ids.compute_IDS(X, p_norm=p_norm))
            for cb in (777, n):
                for bs in (9, d):
                    sd = compute_IDS_streaming(X, block_size=bs, cell_batch=cb,
                                               p_norm=p_norm, two_pass=two_pass)
                    ss = compute_IDS_streaming(sp.csr_matrix(X), block_size=bs, cell_batch=cb,
                                               p_norm=p_norm, two_pass=two_pass)
                    dd, ds = np.abs(full - sd).max(), np.abs(full - ss).max()
                    worst = max(worst, dd, ds)
        print(f"  two_pass={two_pass}: all combos within {worst:.1e}")
    assert worst < 1e-6, f"streaming diverged: {worst:.2e}"
    print(f"PASS exactness — both paths == full to {worst:.1e}\n")


def test_twopass_more_stable():
    """Illustrate the mechanism the paper's two-pass avoids. A near-constant
    feature (large mean, tiny variance) makes the one-pass variance
    'sum(x^2) - sum(x)^2/n' subtract two large, nearly-equal numbers -> lost
    digits. We reproduce it in float32 (whose ~7 digits stand in for float64 at
    very large n) and compare each variance formula to a float64 reference."""
    rng = np.random.default_rng(0)
    n = 200_000
    x = (1000.0 + rng.normal(0, 1e-3, n))          # huge mean, tiny spread
    ref = np.var(x)                                 # float64 reference

    xf = x.astype(np.float32)
    one_pass_var = np.float32((xf * xf).sum()) / n - (np.float32(xf.sum()) / n) ** 2
    mean = np.float32(xf.sum()) / n                 # two-pass: mean first...
    two_pass_var = np.float32(((xf - mean) ** 2).sum()) / n   # ...then center

    e1 = abs(one_pass_var - ref) / ref
    e2 = abs(two_pass_var - ref) / ref
    print(f"reference var          = {ref:.6e}")
    print(f"one-pass (float32) rel error = {e1:.2e}")
    print(f"two-pass (float32) rel error = {e2:.2e}")
    assert e2 <= e1, "two-pass should be at least as accurate"
    print(f"PASS stability — two-pass is ~{e1 / max(e2, 1e-30):.0e}x more accurate here\n")


if __name__ == "__main__":
    test_exactness()
    test_twopass_more_stable()
