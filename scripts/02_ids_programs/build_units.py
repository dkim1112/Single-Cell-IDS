"""
build_units.py — one QC'd file per blood lane and per skin sample (see sc_ids/units.py).

FIRST STEP OF THE PIPELINE: every later script reads these files, one at a time, so memory is
set by the largest single lane (about 1.5 GB) and not by how many lanes are used. The cells are
exactly those of the merged objects this replaced (74,126 blood cells in the nine October lanes,
67,367 skin cells; checked cell for cell).

Filters, applied per unit:
  blood  source study's IncludedInStudy + demuxlet singlets, Status in Healthy/PSO/PSA,
         >= 200 genes, <= 10% mitochondrial, Platelet and Eryth cells removed
  skin   cells present in the source metadata, >= 200 genes, <= 10% mitochondrial
Checkpointed: a unit that already has a file is skipped, so the script can be re-run until done.

  python scripts/02_ids_programs/build_units.py --tissue blood [--units PBMC-01-1,PBMC-03-1] [--budget 150]
  python scripts/02_ids_programs/build_units.py --tissue skin
  python scripts/02_ids_programs/build_units.py --tissue blood --index      # write cells_/genes_ tables
"""
import argparse, glob, gzip, os, re, sys, time, resource
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from sc_ids import paths, units                                              # noqa: E402

KEEP_STATUS = {"Healthy", "PSO", "PSA"}
MIN_GENES, MAX_MT = 200, 10.0
EXCLUDE = {"Platelet", "Eryth"}
ALL_LANES = [f"PBMC-0{r}-{l}" for r in range(1, 8) for l in range(1, 5)]


def rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def qc(X, var_names, obs):
    mt = np.array([g.upper().startswith("MT-") for g in var_names])
    tot = np.asarray(X.sum(axis=1)).ravel()
    ng = np.asarray((X > 0).sum(axis=1)).ravel()
    pct = 100 * np.asarray(X[:, mt].sum(axis=1)).ravel() / np.maximum(tot, 1)
    keep = (ng >= MIN_GENES) & (pct <= MAX_MT)
    obs = obs.assign(total_counts=tot, n_genes=ng, pct_mt=pct)
    return X[keep], obs[keep]


def save(tissue, unit, X, var_names, obs):
    os.makedirs(os.path.dirname(units.unit_path(tissue, unit)), exist_ok=True)
    a = ad.AnnData(X=sp.csr_matrix(X, dtype=np.float32), obs=obs, var=pd.DataFrame(index=pd.Index(np.asarray(var_names, dtype=str), name="gene")))
    tmp = units.unit_path(tissue, unit) + ".tmp"
    a.write_h5ad(tmp, compression="lzf"); os.replace(tmp, units.unit_path(tissue, unit))


_META = {}
def blood_meta():
    if "b" not in _META:
        m = pd.read_csv(paths.PBMC_META, sep="\t", low_memory=False,
                        usecols=["CellName", "Sample", "Subject", "Status", "DemuxletDropletType",
                                 "IncludedInStudy", "CellType"])
        m = m[m.IncludedInStudy.astype(str).str.upper().isin(["TRUE", "T", "1"])
              & (m.DemuxletDropletType.astype(str) == "SNG") & m.Status.isin(KEEP_STATUS)]
        _META["b"] = m.drop_duplicates("CellName").set_index("CellName")
    return _META["b"]


def build_lane(lane):
    try:
        import fast_matrix_market as fmm          # seconds per lane
    except ImportError:
        import scipy.io as fmm                    # same result, about a minute per lane
    from anndata.utils import make_index_unique
    d = paths.PBMC_DIR
    with gzip.open(f"{d}/{lane}.matrix.mtx.gz", "rb") as f:
        X = sp.csr_matrix(fmm.mmread(f).T, dtype=np.float32)          # cells x genes
    bc = pd.read_csv(f"{d}/{lane}.barcodes.tsv.gz", header=None)[0].astype(str)
    ft = pd.read_csv(f"{d}/{lane}.features.tsv.gz", header=None, sep="\t")
    # CITE-seq lanes also carry antibody rows; keep RNA only (what scanpy's gex_only=True does).
    # Leaving them in inflates total counts, lowers the mitochondrial share and lets extra cells through QC.
    gex = (ft[2].astype(str) == "Gene Expression").values if ft.shape[1] > 2 else np.ones(len(ft), bool)
    X, ft = X[:, gex], ft[gex]
    var = make_index_unique(pd.Index(ft[1].astype(str).values))
    n_raw = X.shape[0]
    key = lane + "_" + bc.str.replace(r"-\d+$", "", regex=True)
    m = blood_meta()
    hit = key.isin(m.index).values
    X, key, bc = X[hit], key[hit], bc[hit]
    mm = m.loc[key.values]
    obs = pd.DataFrame({"donor": mm.Subject.astype(str).values, "subtype": mm.CellType.astype(str).values,
                        "disease_fine": mm.Status.astype(str).values, "unit": lane},
                       index=pd.Index(lane + "_" + bc.values, name="cell"))
    X, obs = qc(X, var, obs)
    ok = ~obs.subtype.isin(EXCLUDE).values
    X, obs = X[ok], obs[ok]
    save("blood", lane, X, var, obs)
    return n_raw, X.shape[0]


def skin_symbols():
    from anndata.utils import make_index_unique
    ft = pd.read_csv(f"{paths.PBMC_DIR}/PBMC-01-1.features.tsv.gz", header=None, sep="\t")
    return dict(zip(ft[0].astype(str), ft[1].astype(str))), make_index_unique


def build_skin_sample(f, meta):
    sample = re.sub(r"^GSM\d+_", "", os.path.basename(f)).replace(".csv.gz", "").replace("-", "_", 1)
    rows = meta[meta["orig.ident"] == sample]
    bare = pd.Index(rows.index).str.replace(r"-\d+$", "", regex=True)
    suffix = str(rows.index[0]).rsplit("-", 1)[-1]
    wanted = set(bare + "-1")
    with gzip.open(f, "rt") as fh:
        header = fh.readline().rstrip("\n").split(",")
    index_name = header[0] or "Unnamed: 0"
    take = [c for c in header[1:] if c.strip('"') in wanted]
    assert len(take) == len(wanted), f"{sample}: {len(wanted) - len(take)} metadata barcodes absent"
    chunks, idx, colnames = [], [], None
    for ch in pd.read_csv(f, index_col=0, usecols=[index_name] + take,
                          dtype={c: np.int32 for c in take}, chunksize=4000):
        colnames = ch.columns if colnames is None else colnames
        idx.append(ch.index.astype(str).to_numpy())
        chunks.append(sp.csr_matrix(ch.to_numpy(dtype=np.float32)))
    ens = np.concatenate(idx)
    X = sp.vstack(chunks, format="csr").T.tocsr()                         # cells x genes
    keys = pd.Index(colnames).str.replace(r"-\d+$", "", regex=True) + "-" + suffix
    sym, uniq = skin_symbols()
    var = uniq(pd.Index([sym.get(e) or e for e in ens]))
    mm = meta.loc[keys]
    obs = pd.DataFrame({"donor": mm["donor"].astype(str).values, "subtype": mm["fullsubtype"].astype(str).values,
                        "disease_fine": mm["disease"].astype(str).values, "unit": sample},
                       index=pd.Index(keys, name="cell"))
    n_raw = X.shape[0]
    X, obs = qc(X, var, obs)
    ok = ~obs.subtype.isin(EXCLUDE).values
    save("skin", sample, X[ok], var, obs[ok])
    return sample, n_raw, int(ok.sum())


def write_index(tissue, order):
    rows, g0 = [], None
    for u in order:
        p = units.unit_path(tissue, u)
        if not os.path.exists(p):
            continue
        a = ad.read_h5ad(p, backed="r")
        o = a.obs[["donor", "subtype", "disease_fine", "unit", "total_counts"]].copy()
        o.insert(0, "row", np.arange(len(o))); o.insert(0, "cell", a.obs_names.values)
        rows.append(o)
        g = list(a.var_names)
        assert g0 is None or g == g0, f"{u}: gene order differs"
        g0 = g
        a.file.close()
    c = pd.concat(rows, ignore_index=True)
    c.to_csv(os.path.join(units.WORK, f"cells_{tissue}.csv.gz"), index=False)
    open(os.path.join(units.WORK, f"genes_{tissue}.txt"), "w").write("\n".join(g0) + "\n")
    print(f"[index] {tissue}: {len(c):,} cells in {c.unit.nunique()} units, {len(g0):,} genes")
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", required=True, choices=["blood", "skin"])
    ap.add_argument("--units", default=None)
    ap.add_argument("--budget", type=float, default=150, help="stop starting new units after this many seconds")
    ap.add_argument("--index", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    if a.tissue == "blood":
        todo = a.units.split(",") if a.units else ALL_LANES
        if not a.index:
            for u in todo:
                if os.path.exists(units.unit_path("blood", u)):
                    continue
                if time.time() - t0 > a.budget:
                    print("[budget] stopping; re-run to continue"); break
                t = time.time(); n0, n1 = build_lane(u)
                print(f"[blood] {u}: {n0:,} droplets -> {n1:,} cells  {time.time()-t:.0f}s  peak {rss():.2f} GB", flush=True)
        else:
            write_index("blood", ALL_LANES)
    else:
        files = sorted(glob.glob(os.path.join(paths.SKIN_DIR, "*.csv.gz")))
        names = [re.sub(r"^GSM\d+_", "", os.path.basename(f)).replace(".csv.gz", "").replace("-", "_", 1) for f in files]
        if not a.index:
            meta = pd.read_csv(paths.SKIN_META, sep="\t", index_col=0, low_memory=False)
            meta.index = meta.index.astype(str)
            for f, nm in zip(files, names):
                if os.path.exists(units.unit_path("skin", nm)):
                    continue
                if time.time() - t0 > a.budget:
                    print("[budget] stopping; re-run to continue"); break
                t = time.time(); s, n0, n1 = build_skin_sample(f, meta)
                print(f"[skin] {s}: {n0:,} -> {n1:,} cells  {time.time()-t:.0f}s  peak {rss():.2f} GB", flush=True)
        else:
            write_index("skin", names)


if __name__ == "__main__":
    main()
