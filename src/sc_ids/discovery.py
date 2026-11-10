"""
discovery.py — IDS program discovery on the per-lane files (sc_ids.units).

One implementation, used by scripts/02_ids_programs/ids_programs_disease.py (the catalogue run)
and scripts/04_method_checks/discovery_runs.py (the same run with one setting changed).

What a run does, per stratum (tissue x disease state, optionally inside one lineage):
  1. take the stratum's cells and draw at most `cap` of them, stratified by cell subtype
  2. keep the panel genes detected in at least `floor` of those cells
  3. library-size normalise, scale each gene to [0, 8] (no log1p), compute IDS with the
     authors' ids.compute_IDS through compute_IDS_blocked
  4. link genes whose IDS reaches `min_ids`, split the graph by greedy modularity, keep groups of 3-50 genes

IDS itself is never reimplemented here.

Settings behind the 194 catalogue detections (recovered by reproducing them gene for gene):
  cap 4000, floor 0.02, seed 0, norm "panel", nine blood lanes (LANES_USED), and MIN_IDS below.
"""
import json, os
import numpy as np, pandas as pd, scipy.sparse as sp
from . import units, programs as P
from .lineage import LINEAGE
from .datasets.gse194315 import LANES_USED

STRATA = {"PBMC_Healthy": ("blood", "Healthy"), "PBMC_PSA": ("blood", "PSA"), "PBMC_PSO": ("blood", "PSO"),
          "skin_NS": ("skin", "NS"), "skin_PN": ("skin", "PN"), "skin_PP": ("skin", "PP")}
# edge threshold per analysis; None = whole tissue. Read off each run's own IDS distribution in October.
MIN_IDS = {None: 0.50, "T+NK": 0.35, "Myeloid": 0.40, "Keratinocyte": 0.50}

# Added to the panel unconditionally: mean-variance HVG selection drops programs that are OFF in
# most cells, which are exactly the on/off programs IDS is good at finding.
CURATED = """
IL17A IL17F IL23A IL12B IL36G IL36A IL1B TNF IFNG IL22 IL19 IL20 IL6 IL10 TGFB1
S100A7 S100A8 S100A9 S100A12 KRT6A KRT6B KRT16 KRT17 SERPINB3 SERPINB4 DEFB4A
LCN2 PI3 CAMP SPRR2A SPRR2B SPRR1B
CXCL1 CXCL2 CXCL8 CXCL9 CXCL10 CXCL11 CXCL13 CCL2 CCL5 CCL17 CCL19 CCL20 CCL22
ISG15 IFI6 IFI27 IFI44L MX1 MX2 OAS1 OAS3 STAT1 IRF7 IFIT1 IFIT3 BST2 XAF1
CD3D CD3E CD3G TRAC TRBC2 IL7R FOXP3 CTLA4 PDCD1 RORC TBX21 GATA3 IL2RA CCR6
CCR7 SELL TCF7 LEF1 CD4 CD8A CD8B CD69 ITGAE CXCR4 KLRB1
GNLY NKG7 PRF1 GZMA GZMB GZMH GZMK KLRD1 FCGR3A NCAM1 KLRF1 TRDC TRGC1 TRGC2
KLRC1 KLRC2 TYROBP FCER1G CST7 CCL4
CD14 LYZ VCAN FCN1 CD1C CLEC10A CD207 CD1A FCER1A ITGAX CSF1R C1QA C1QB MRC1
CD68 CD163 AIF1 HLA-DRA HLA-DRB1 CIITA
MS4A1 CD79A CD79B IGHM IGKC JCHAIN TCL1A
HSPA1A HSPA1B HSPA8 DNAJB1 DNAJA1 HSPH1 HSP90AA1 JUN JUNB FOS FOSB EGR1 ZFP36
ATF3 NFKBIA TNFAIP3 SOCS3 DUSP1
MKI67 TOP2A PCNA TYMS STMN1 TUBA1B
KRT1 KRT2 KRT5 KRT10 KRT14 KRT15 FLG LOR IVL DSP CDSN
COL1A1 COL1A2 COL3A1 DCN LUM PDGFRB ACTA2 TAGLN MYL9
PECAM1 VWF CDH5 CLDN5 ACKR1 SELE
TPSAB1 TPSB2 CPA3 KIT MS4A2 CMA1
MLANA PMEL TYRP1 DCT TYR
""".split()
_CELLS = {}


def tissue_cells(tissue, lanes=None):
    """Cell table of a tissue. Blood defaults to the nine October lanes; pass lanes='all' for every lane."""
    key = (tissue, "all" if lanes == "all" else tuple(lanes or LANES_USED))
    if key not in _CELLS:
        c = units.cells(tissue)
        if tissue == "blood" and lanes != "all":
            c = c[c.unit.isin(list(lanes or LANES_USED))].reset_index(drop=True)
        _CELLS[key] = c
    return _CELLS[key]


def stratum_cells(lineage, stratum, lanes=None):
    """(tissue, cells) of one stratum, optionally restricted to one coarse lineage ('T+NK' keeps both)."""
    tissue, dis = STRATA[stratum]
    c = tissue_cells(tissue, lanes)
    m = (c.disease_fine == dis).values
    if lineage:
        m &= c.subtype.map(LINEAGE).fillna("Other").isin(set(lineage.split("+"))).values
    return tissue, c[m]


def subsample_idx(lab, cap, seed=0):
    """Positions of at most `cap` cells, stratified by label so rare types survive (floor of cap/(8*types))."""
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


def build_panel(lineage, n_hvg, out, lanes=None, verbose=True):
    """ONE gene panel for every stratum: union of each stratum's highly variable genes (chosen among
    genes present in both tissues) plus the curated list. Comparing programs across strata is only
    fair if every stratum was offered the same genes."""
    import scanpy as sc, anndata as ad
    shared = sorted(set(units.genes("blood")) & set(units.genes("skin")))
    chosen, per = set(), {}
    for s in sorted(STRATA):
        tissue, c = stratum_cells(lineage, s, lanes)
        if len(c) < 200:
            if verbose:
                print(f"  [panel] {s:<12} only {len(c)} cells in this lineage - skipped", flush=True)
            continue
        sub = ad.AnnData(X=units.read_cells(tissue, c, shared), var=pd.DataFrame(index=shared))
        sc.pp.normalize_total(sub, target_sum=1e4)
        sc.pp.log1p(sub)
        sc.pp.highly_variable_genes(sub, n_top_genes=min(n_hvg, sub.n_vars - 1))
        per[s] = list(sub.var_names[sub.var["highly_variable"].values])
        chosen |= set(per[s])
        if verbose:
            print(f"  [panel] {s:<12} {sub.n_obs:>7,} cells -> {len(per[s]):,} HVGs", flush=True)
        del sub
    hvg = sorted(chosen)
    cur = sorted((set(CURATED) & set(shared)) - chosen)
    panel = sorted(chosen | set(cur))
    if out:
        pd.DataFrame({"gene": panel, "source": ["HVG" if g in chosen else "curated" for g in panel]}
                     ).to_csv(os.path.join(out, "gene_panel.csv"), index=False)
        json.dump({"panel": panel, "per_stratum_hvg": per, "curated_added": cur,
                   "curated_absent_from_data": sorted(set(CURATED) - set(shared))},
                  open(os.path.join(out, "gene_panel.json"), "w"))
    if verbose:
        print(f"  [panel] HVG union {len(hvg):,} + curated not already in it {len(cur):,} = {len(panel):,} genes", flush=True)
    return panel


def compute_ids(lineage, stratum, panel, seed=0, norm="panel", cap=4000, floor=0.02,
                min_cells=400, block_size=1024, lanes=None):
    """IDS matrix of one stratum. Returns dict(C, genes, dropped, cells), or dict(skipped=n) when the
    stratum has fewer than `min_cells` cells (IDS estimates were unstable below roughly 500-1000).

    norm='panel'  library size over the panel genes that pass the floor (the catalogue run)
    norm='total'  library size over all genes
    """
    import anndata as ad
    from .pipeline import preprocess, minmax_scale_sparse
    from .batched import compute_IDS_blocked
    tissue, c = stratum_cells(lineage, stratum, lanes)
    sel = c.iloc[subsample_idx(c.subtype.astype(str).values, cap, seed)]
    if len(sel) < min_cells:
        return {"skipped": int(len(sel))}
    X = units.read_cells(tissue, sel, panel)                         # raw counts, panel genes
    keep = np.asarray((X > 0).sum(axis=0)).ravel() / X.shape[0] >= floor
    genes = [g for g, k in zip(panel, keep) if k]
    dropped = [g for g, k in zip(panel, keep) if not k]
    X = sp.csr_matrix(X[:, keep], dtype=np.float32)
    if norm == "panel":
        Xs, _ = preprocess(ad.AnnData(X=X, var=pd.DataFrame(index=genes)), n_hvg=None, keep_sparse=True)
    else:
        X = sp.diags((1e4 / np.maximum(sel.total_counts.values, 1)).astype(np.float32)) @ X
        Xs, _ = minmax_scale_sparse(sp.csr_matrix(X))
    C = compute_IDS_blocked(Xs, block_size=block_size, p_norm="max", verbose=False)
    return {"C": C, "genes": genes, "dropped": dropped, "cells": int(len(sel))}


def programs_from(C, genes, min_ids, min_genes=P.PROGRAM_SIZE_MIN, max_genes=P.PROGRAM_SIZE_MAX, with_dropped=False):
    """Ranked program table from an IDS matrix (greedy modularity on edges >= min_ids, size filter)."""
    from .pipeline import extract_modules
    mods = extract_modules(C, genes, min_ids=min_ids, min_size=2, method="greedy")
    mods = P.annotate_programs(mods, C, min_ids=min_ids)
    progs, dropped = P.rank_programs(mods, min_genes=min_genes, max_genes=max_genes, verbose=False)
    t = P.programs_table(progs)
    return (t, dropped) if with_dropped else t


def support(lineage, stratum, tbl, lanes=None):
    """How many cells of the WHOLE stratum (not the sample) detect each program's genes."""
    tissue, c = stratum_cells(lineage, stratum, lanes)
    need = sorted({g for gs in tbl["genes"] for g in gs.split(", ")})
    D = (units.read_cells(tissue, c, need) > 0).tocsc()
    pos = {g: i for i, g in enumerate(need)}
    n_all, n_half = [], []
    for gs in tbl["genes"]:
        cols = [pos[g] for g in gs.split(", ")]
        k = np.asarray(D[:, cols].sum(axis=1)).ravel()
        n_all.append(int((k == len(cols)).sum())); n_half.append(int((k >= max(2, len(cols) / 2)).sum()))
    tbl = tbl.copy()
    tbl["cells_all_genes"], tbl["cells_half_genes"] = n_all, n_half
    tbl["frac_cells_half_genes"] = [round(v / len(c), 4) for v in n_half]
    tbl["stratum_cells"] = len(c)
    return tbl
