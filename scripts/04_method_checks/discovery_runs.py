"""
discovery_runs.py — rerun IDS program discovery under controlled changes.

Reproduces scripts/02_ids_programs/ids_programs_disease.py step for step (same saved gene
panel, same stratified 4,000-cell sample, same 2% detection floor, same [0,8] scaling, the
authors' compute_IDS through compute_IDS_blocked, greedy modularity, 3-50 genes), reading the
cells from the per-unit files, and lets two things vary:

  --seed   which random 4,000 cells are drawn         (October used seed 0 only)
  --norm   panel = library size over the panel genes that pass the floor (what October did)
           total = library size over ALL genes       (what the permutation test and the scoring use)
  --cap    how many cells are sampled per stratum    (October used 4,000)

Edge thresholds used in October, recovered by reproducing all 194 detections gene for gene:
  whole-tissue 0.50, T+NK 0.35, Myeloid 0.40, Keratinocyte 0.50   (MIN_IDS below)

Blood is restricted to the nine October lanes so that only the stated change differs.
The IDS matrix is cached, so thresholds can be tried without recomputing.

  python scripts/04_method_checks/discovery_runs.py --analysis whole-tissue --seeds 0 --norm panel --min-ids 0.5
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from sc_ids import paths, units, programs as P                                # noqa: E402
from sc_ids.pipeline import preprocess, extract_modules, minmax_scale_sparse  # noqa: E402
from sc_ids.batched import compute_IDS_blocked                                # noqa: E402
from sc_ids.lineage import LINEAGE                                            # noqa: E402
from sc_ids.datasets.gse194315 import LANES_USED                              # noqa: E402

RUNS = {"whole-tissue": paths.DISEASE_RUN, **paths.LINEAGE_RUNS}
STRATA = {"PBMC_Healthy": ("blood", "Healthy"), "PBMC_PSA": ("blood", "PSA"), "PBMC_PSO": ("blood", "PSO"),
          "skin_NS": ("skin", "NS"), "skin_PN": ("skin", "PN"), "skin_PP": ("skin", "PP")}
MIN_IDS = {"whole-tissue": 0.50, "T+NK": 0.35, "Myeloid": 0.40, "Keratinocyte": 0.50}
CACHE = os.path.join(units.WORK, "discovery")
_CELLS = {}


def stratum_cells(analysis, stratum):
    tissue, dis = STRATA[stratum]
    if tissue not in _CELLS:
        c = units.cells(tissue)
        _CELLS[tissue] = c[c.unit.isin(list(LANES_USED))].reset_index(drop=True) if tissue == "blood" else c
    c = _CELLS[tissue]
    m = (c.disease_fine == dis).values
    if analysis != "whole-tissue":
        m &= c.subtype.map(LINEAGE).fillna("Other").isin(set(analysis.split("+"))).values
    return tissue, c[m]


def subsample_idx(lab, cap, seed):
    """Same allocation and draws as ids_programs_disease.subsample; returns sorted positions."""
    n = len(lab)
    if n <= cap:
        return np.arange(n)
    rng = np.random.default_rng(seed)
    types, counts = np.unique(lab, return_counts=True)
    floor = max(1, cap // (len(types) * 8))
    alloc = np.minimum(counts, np.maximum(floor, (counts / counts.sum() * cap).astype(int)))
    while alloc.sum() > cap:
        alloc[np.argmax(alloc)] -= 1
    keep = []
    for t, k in zip(types, alloc):
        idx = np.where(lab == t)[0]
        keep.append(rng.choice(idx, size=min(k, len(idx)), replace=False))
    return np.sort(np.concatenate(keep))


def ids_matrix(analysis, stratum, seed, norm, cap=4000, min_cells=400, floor=0.02):
    tag = f"{analysis.replace('+','-')}__{stratum}__s{seed}__{norm}" + ("" if cap == 4000 else f"__cap{cap}")
    fC, fg = os.path.join(CACHE, tag + ".npy"), os.path.join(CACHE, tag + ".json")
    if os.path.exists(fC):
        return np.load(fC).astype(np.float64), json.load(open(fg))
    if os.path.exists(fg):                       # recorded as skipped
        return None, json.load(open(fg))
    os.makedirs(CACHE, exist_ok=True)
    tissue, c = stratum_cells(analysis, stratum)
    sel = c.iloc[subsample_idx(c.subtype.astype(str).values, cap, seed)]
    if len(sel) < min_cells:
        json.dump({"skipped": len(sel)}, open(fg, "w")); return None, {"skipped": len(sel)}
    panel = json.load(open(os.path.join(RUNS[analysis], "gene_panel.json")))["panel"]
    X = units.read_cells(tissue, sel, panel)                         # raw counts, panel genes
    keep = np.asarray((X > 0).sum(axis=0)).ravel() / X.shape[0] >= floor
    genes = [g for g, k in zip(panel, keep) if k]
    X = sp.csr_matrix(X[:, keep], dtype=np.float32)
    if norm == "panel":
        Xs, _ = preprocess(ad.AnnData(X=X, var=pd.DataFrame(index=genes)), n_hvg=None, keep_sparse=True)
    else:
        X = sp.diags((1e4 / np.maximum(sel.total_counts.values, 1)).astype(np.float32)) @ X
        Xs, _ = minmax_scale_sparse(sp.csr_matrix(X))
    C = compute_IDS_blocked(Xs, block_size=1024, p_norm="max", verbose=False)
    np.save(fC, C.astype(np.float32))
    json.dump({"genes": genes, "cells": int(len(sel))}, open(fg, "w"))
    return C, {"genes": genes, "cells": int(len(sel))}


def programs_from(C, genes, min_ids):
    mods = extract_modules(C, genes, min_ids=min_ids, min_size=2, method="greedy")
    mods = P.annotate_programs(mods, C, min_ids=min_ids)
    progs, _ = P.rank_programs(mods, verbose=False)
    return P.programs_table(progs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analysis", required=True, choices=list(RUNS))
    ap.add_argument("--seeds", default="0"); ap.add_argument("--norm", default="panel", choices=["panel", "total"])
    ap.add_argument("--min-ids", type=float, default=None)
    ap.add_argument("--cap", type=int, default=4000, help="cells sampled per stratum (October: 4000)")
    ap.add_argument("--budget", type=float, default=140)
    a = ap.parse_args()
    t0 = time.time()
    strata = list(pd.read_csv(os.path.join(RUNS[a.analysis], "ids_summary.csv")).stratum)
    for seed in [int(s) for s in a.seeds.split(",")]:
        for s in strata:
            if time.time() - t0 > a.budget:
                print("[budget] stopping; re-run to continue"); return
            t = time.time()
            C, info = ids_matrix(a.analysis, s, seed, a.norm, cap=a.cap)
            if C is None:
                print(f"{a.analysis} {s} seed {seed}: skipped ({info})"); continue
            msg = f"{a.analysis} {s} seed {seed} {a.norm}: {info['cells']} cells x {len(info['genes'])} genes"
            if a.min_ids is not None:
                msg += f" -> {len(programs_from(C, info['genes'], a.min_ids))} programs at {a.min_ids}"
            print(msg + f"  ({time.time()-t:.0f}s)", flush=True)
    print("done")


if __name__ == "__main__":
    main()
