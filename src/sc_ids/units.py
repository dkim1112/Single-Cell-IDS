"""
units.py — the data as one QC'd file per UNIT (a blood lane or a skin sample), never merged.

Why: the earlier scripts loaded a whole tissue into memory, so memory grew with every lane
added and the 9-lane build was killed at 6 GB. Here every step reads one unit at a time, so
peak memory is set by the largest single lane, not by how many lanes are used.

Layout under WORK (env SC_IDS_WORK, default results/00_units):
  <tissue>/<unit>.h5ad      raw counts (CSR float32) + obs: donor, subtype, disease_fine
  cells_<tissue>.csv.gz     one row per cell: cell, unit, row, donor, subtype, disease_fine, total_counts
  genes_<tissue>.txt        gene order shared by every unit of the tissue

Scoring is two passes over the units and reproduces scanpy's `sc.tl.score_genes` on the
merged object: pass 1 gets each gene's mean log-expression over ALL cells, the control genes
are then picked once by scanpy's own binning code from those means, and pass 2 scores each
unit against the same controls.
"""
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
from . import paths

WORK = paths.UNITS
TISSUES = ("blood", "skin")
# disease label shared by the two tissues; the source label is kept as disease_fine
DISEASE = {"Healthy": "healthy", "PSO": "psoriasis", "PSA": "PsA", "NS": "healthy", "PN": "psoriasis", "PP": "psoriasis"}


def unit_path(tissue, unit):
    return os.path.join(WORK, tissue, f"{unit}.h5ad")


def cells(tissue):
    return pd.read_csv(os.path.join(WORK, f"cells_{tissue}.csv.gz"), dtype={"donor": str})


def genes(tissue):
    return [g.rstrip("\n") for g in open(os.path.join(WORK, f"genes_{tissue}.txt"))]


def read_unit(tissue, unit):
    import anndata as ad
    return ad.read_h5ad(unit_path(tissue, unit))


def lognorm(X, target=1e4):
    """normalize_total(target) over ALL genes, then log1p. CSR float32 in, CSR float32 out."""
    X = sp.csr_matrix(X, dtype=np.float32, copy=True)
    tot = np.asarray(X.sum(axis=1)).ravel()
    tot[tot == 0] = 1.0
    X.data *= np.repeat((target / tot).astype(np.float32), np.diff(X.indptr))
    np.log1p(X.data, out=X.data)
    return X


def read_cells(tissue, sel, gene_list=None, norm=None):
    """Counts for selected cells. `sel` is a slice of cells(tissue) (needs columns unit,row);
    rows come back in the order of `sel`. norm='total' divides by each cell's total over ALL
    genes (x1e4) before the gene subset is taken; None returns raw counts."""
    g_all = genes(tissue)
    cols = None if gene_list is None else pd.Index(g_all).get_indexer(gene_list)
    parts, order = [], []
    for unit, s in sel.groupby("unit", sort=False):
        a = read_unit(tissue, unit)
        X = sp.csr_matrix(a.X, dtype=np.float32)[s["row"].values]
        if norm == "total":
            tot = np.asarray(X.sum(axis=1)).ravel(); tot[tot == 0] = 1.0
            X = sp.diags((1e4 / tot).astype(np.float32)) @ X
        parts.append(sp.csr_matrix(X[:, cols]) if cols is not None else X)
        order.append(s.index.values)
        del a
    X = sp.vstack(parts, format="csr")
    pos = pd.Series(np.arange(len(sel)), index=np.concatenate(order)).loc[sel.index.values].values
    return X[pos]


# ---------------------------------------------------------------- scoring
def gene_means(tissue, units=None, cache=None):
    """Pass 1: mean log-normalized expression of every gene over all cells of the units."""
    if cache and os.path.exists(cache):
        z = np.load(cache); return z["mean"], int(z["n"])
    c = cells(tissue)
    units = list(units) if units is not None else list(c["unit"].unique())
    s, n = None, 0
    for u in units:
        X = lognorm(read_unit(tissue, u).X)
        v = np.asarray(X.sum(axis=0, dtype=np.float64)).ravel()
        s = v if s is None else s + v
        n += X.shape[0]
    m = s / n
    if cache:
        np.savez(cache, mean=m, n=n)
    return m, n


def control_sets(var_names, mean, programs, ctrl_size=50, n_bins=25, random_state=0):
    """Control genes per program, chosen by scanpy's own binning code from global means.
    `programs` is {id: [genes]}. Returns {id: (gene_idx, ctrl_idx)}; programs with <2 genes
    present are skipped, as in the earlier scripts."""
    from scanpy.tools._score_genes import _score_genes_bins
    vn = pd.Index(var_names)
    pool = vn.astype("string")
    row = np.asarray(mean, dtype=np.float64)[None, :]
    out = {}
    for pid, gl in programs.items():
        gl = pd.Index([g for g in gl if g in vn])
        if len(gl) < 2:
            continue
        gl = pd.Index(gl).intersection(vn)
        np.random.seed(random_state)
        ctrl = pd.Index([], dtype="string")
        for r in _score_genes_bins(gl, pool, ctrl_as_ref=True, ctrl_size=ctrl_size,
                                   n_bins=n_bins, get_subset=lambda genes: row):
            ctrl = ctrl.union(r)
        out[pid] = (vn.get_indexer(gl), vn.get_indexer(pd.Index(ctrl.astype(str))))
    return out


def score_matrix(X_lognorm, ctrl):
    """Pass 2 for one unit: cells x programs score = mean(program genes) - mean(controls)."""
    G = X_lognorm.shape[1]
    ids = list(ctrl)
    r, c, v = [], [], []
    for j, pid in enumerate(ids):
        gi, ci = ctrl[pid]
        r += list(gi) + list(ci); c += [j] * (len(gi) + len(ci))
        v += [1.0 / len(gi)] * len(gi) + [-1.0 / len(ci)] * len(ci)
    W = sp.csc_matrix((v, (r, c)), shape=(G, len(ids)), dtype=np.float64)
    return np.asarray((X_lognorm.astype(np.float64) @ W).todense()), ids


def stage_units(tissue, lanes=None):
    """Units a catalogue stage works on: the nine October blood lanes unless lanes='all'; every skin sample."""
    from .datasets.gse194315 import LANES_USED
    if tissue == "blood" and lanes != "all":
        return list(lanes or LANES_USED)
    return list(cells(tissue)["unit"].unique())


def iter_scores(tissue, programs, unit_list=None, cache=None):
    """Score `programs` ({id: [genes]}) in every cell, one unit at a time.
    Yields (unit, obs, S, ids): S is cells x programs for that unit, equal to what sc.tl.score_genes
    (ctrl_size 50, random_state 0) gives on the merged object of all `unit_list` units."""
    unit_list = list(unit_list) if unit_list is not None else stage_units(tissue)
    mean, _ = gene_means(tissue, unit_list, cache=cache)
    ctrl = control_sets(genes(tissue), mean, programs)
    for u in unit_list:
        a = read_unit(tissue, u)
        S, ids = score_matrix(lognorm(a.X), ctrl)
        yield u, a.obs, S, ids
