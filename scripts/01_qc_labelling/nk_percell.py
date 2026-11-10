"""
nk_percell.py — is the NK sub-cluster one cell type, or a mixture?

The cluster-level marker table (rank_markers.py) gives one number per cluster per
tissue, so a cluster that is half NK and half gamma-delta T looks the same as a
cluster of cells that each co-express both. This script makes the call CELL BY
CELL, so the two cases can be told apart.

Per-cell call, using the lab's definition (NK = FCGR3A/NKG7/GZMB/PRF1 present,
CD3D absent) plus the T-lineage receptor genes:
  NK score   = how many of FCGR3A, NKG7, GZMB, PRF1 are detected (>= 2 counts as NK+)
  T score    = how many of CD3D, CD3E, TRAC are detected   (>= 1 counts as T+)
  NK only    NK+ and not T+        T only    T+ and not NK+
  both       NK+ and T+            undetermined  neither

A clean cell type is mostly one column. A real mixture splits across NK-only and
T-only. To know what "clean" looks like in THIS data, the same call is applied to
reference clusters -- sub-clusters the annotations agree are plain T cells -- and
their numbers are printed alongside.

Run on the UNFILTERED gene set: the gene filter deletes FCGR3A, PRF1 and TRDC.

Usage:
  python scripts/01_qc_labelling/nk_percell.py --lineage results/01_qc_labelling/clustered-theta4/step5_t_embedded.h5ad \
      --cluster 0 --reference 1,2 --out results/01_qc_labelling/markers-theta4
"""
import argparse, os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

NK_GENES = ["FCGR3A", "NKG7", "GZMB", "PRF1"]
T_GENES = ["CD3D", "CD3E", "TRAC"]
GD_GENES = ["TRDC", "TRGC1", "TRGC2"]
NK_MIN, T_MIN = 2, 1


def detected(a, genes):
    """cells x len(genes) boolean: gene has at least one read in that cell."""
    cols = [a.var_names.get_loc(g) for g in genes if g in a.var_names]
    if not cols:
        return np.zeros((a.n_obs, 0), dtype=bool)
    X = sp.csr_matrix(a.X)[:, cols].toarray()
    return X > 0


def calls(a):
    nk = detected(a, NK_GENES).sum(1)
    t = detected(a, T_GENES).sum(1)
    nk_pos, t_pos = nk >= NK_MIN, t >= T_MIN
    lab = np.full(a.n_obs, "undetermined", dtype=object)
    lab[nk_pos & ~t_pos] = "NK only"
    lab[t_pos & ~nk_pos] = "T only"
    lab[nk_pos & t_pos] = "both"
    return lab


def summarise(a, name):
    rows = []
    lab = calls(a)
    gd = detected(a, GD_GENES).any(1)
    for tis in ("PBMC", "skin"):
        m = (a.obs["tissue"].astype(str).values == tis)
        if m.sum() < 20:
            continue
        r = {"group": name, "tissue": tis, "n": int(m.sum())}
        for k in ("NK only", "T only", "both", "undetermined"):
            r[k] = int(((lab == k) & m).sum())
        r["gamma-delta receptor +"] = int((gd & m).sum())
        rows.append(r)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lineage", required=True)
    ap.add_argument("--cluster", required=True)
    ap.add_argument("--reference", default="", help="comma-separated clean-T clusters")
    ap.add_argument("--cluster-col", default="leiden_res_0.40")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    e = ad.read_h5ad(args.lineage)
    cl = e.obs[args.cluster_col].astype(str).values
    rows = summarise(e[cl == str(args.cluster)].copy(), f"cluster {args.cluster}")
    for c in [x for x in args.reference.split(",") if x.strip()]:
        rows += summarise(e[cl == c.strip()].copy(), f"cluster {c.strip()} (reference)")

    df = pd.DataFrame(rows)
    p = os.path.join(args.out, f"nk_percell_cluster{args.cluster}.csv")
    df.to_csv(p, index=False)
    print(df.to_string(index=False), flush=True)
    print(f"[out] {p}", flush=True)


if __name__ == "__main__":
    main()
