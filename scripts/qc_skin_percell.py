import os, sys, glob
import numpy as np, pandas as pd

REPO = os.path.expanduser("~/mnt/Single-Cell-IDS")
DATA = os.path.join(REPO, "data", "GSE173706")
FEAT = os.path.join(REPO, "data", "GSE194315_PBMC-01-07_processed_data_files",
                    "PBMC-01-1.features.tsv.gz")
PARTS = os.path.expanduser("~/qcskin/parts")

feat = pd.read_csv(FEAT, sep="\t", header=None, usecols=[0, 1], names=["e", "s"])
MT = set(feat.loc[feat["s"].astype(str).str.startswith("MT-"), "e"])

for path in sys.argv[1:]:
    base = os.path.basename(path)
    sample = base.split("_", 1)[1].replace(".csv.gz", "").replace("-", "_")
    out = os.path.join(PARTS, sample + ".csv")
    if os.path.exists(out):
        print("skip", sample, flush=True); continue
    header = pd.read_csv(path, nrows=0)
    cells = list(header.columns[1:])
    dt = {c: np.int32 for c in cells}
    n = len(cells)
    tot = np.zeros(n, np.int64); det = np.zeros(n, np.int64); mt = np.zeros(n, np.int64)
    for ch in pd.read_csv(path, index_col=0, dtype=dt, chunksize=3000):
        v = ch.to_numpy()
        tot += v.sum(0); det += (v > 0).sum(0)
        m = ch.index.isin(MT)
        if m.any(): mt += v[m].sum(0)
    pd.DataFrame({
        "sample": sample, "barcode": cells,
        "n_genes": det, "total_counts": tot,
        "pct_mt": np.where(tot > 0, 100.0 * mt / np.maximum(tot, 1), 0.0),
    }).to_csv(out, index=False)
    print("done", sample, n, flush=True)
