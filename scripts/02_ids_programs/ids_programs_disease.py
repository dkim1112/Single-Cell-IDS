"""
ids_programs_disease.py — IDS gene programs per (tissue x disease) stratum.

Six strata: PBMC x {Healthy, PSO, PSA} and skin x {NS, PN, PP}. Each gets its
own gene-gene IDS matrix, so a program can later be called tissue-shared,
tissue-specific, disease-shared or disease-specific by comparing strata.

TWO DESIGN CHOICES THAT MATTER, both about comparability:

  1. ONE FIXED GENE PANEL for every stratum. Comparing programs across strata is
     comparing gene SETS, which is only meaningful if every stratum had the
     chance to use the same genes. Per-stratum HVG selection would make a
     program look tissue-specific merely because its genes were never offered to
     the other tissue. The panel is the union of each stratum's highly variable
     genes plus a curated marker list added unconditionally -- mean-variance HVG
     selection systematically misses programs that are OFF in most cells, which
     are exactly the bimodal programs IDS is good at finding.

  2. A COMMON CELL CAP per stratum. The strata differ in size by more than an
     order of magnitude. Without a cap the largest stratum would dominate the
     panel and the IDS estimates would not be comparable across strata at equal
     precision. Subsampling is stratified by cell type so rare types survive.

Representation is the paper's: normalize_total -> per-gene min-max to [0, 8],
NO log1p (log1p before the [0,8] scaling pushes cells into the saturated region
of the Gaussian feature envelope and collapses IDS). IDS itself is the authors'
`ids.compute_IDS`, called block-wise over the gene axis, never reimplemented.

Pass 1 computes and caches the matrices and prints each stratum's off-diagonal
IDS percentiles. min_ids is then read off those percentiles -- the value where
the background bulk ends and the program mode begins -- and pass 2 is run with
--min-ids and --reuse to extract programs without recomputing IDS.

Usage:
  python scripts/02_ids_programs/ids_programs_disease.py \
      --blood results/02_ids_programs/ids-disease/blood_qc.h5ad \
      --skin  results/02_ids_programs/ids-disease/skin_qc.h5ad \
      --out results/02_ids_programs/ids-disease --n-hvg 1200 --cell-cap 4000
  # one coarse lineage only (writes ids-lin-TNK / -MY / -KC):  add  --lineage T+NK
  python scripts/02_ids_programs/ids_programs_disease.py ... --reuse --min-ids 0.5
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from sc_ids.pipeline import preprocess, extract_modules      # noqa: E402
from sc_ids.batched import compute_IDS_blocked               # noqa: E402
from sc_ids import programs as P                             # noqa: E402

from sc_ids.lineage import LINEAGE                                           # noqa: E402

# Added to the panel unconditionally. Rationale in the module docstring: HVG
# selection drops low-mean bimodal programs, which is the failure mode that would
# quietly remove the psoriasis biology we are looking for.
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


def lineage_mask(obs, lineage):
    """Cells of one coarse lineage. 'T+NK' keeps both, since they are the pair the
    cross-tissue cytotoxic programme spans and separating them would split it."""
    if not lineage:
        return np.ones(len(obs), dtype=bool)
    want = set(lineage.split("+"))
    lin = obs["subtype"].astype(str).map(LINEAGE).fillna("Other")
    return lin.isin(want).values


def build_panel(src, shared, n_hvg, out, lineage=None):
    """Union of per-stratum HVGs, plus the curated list. One panel for all strata.

    HVGs are selected on the SHARED gene set, so a gene cannot enter the panel
    for one tissue while being absent from the other -- that alone would
    manufacture a tissue-specific program.
    """
    chosen, per, all_names = set(), {}, set(shared)
    pairs = [(s, t) for t, a in src.items()
             for s in sorted(a.obs["stratum"].astype(str).unique())]
    for s, t in pairs:
        a = src[t]
        m = np.where((a.obs["stratum"].astype(str) == s).values
                     & lineage_mask(a.obs, lineage))[0]
        if len(m) < 200:
            print(f"  [panel] {s:<12} only {len(m)} cells in this lineage \u2014 skipped",
                  flush=True)
            continue
        cols = [a.var_names.get_loc(g) for g in shared]
        sub = ad.AnnData(X=a.X[m][:, cols], obs=a.obs.iloc[m].copy(),
                         var=pd.DataFrame(index=shared))
        sc.pp.normalize_total(sub, target_sum=1e4)
        sc.pp.log1p(sub)
        sc.pp.highly_variable_genes(sub, n_top_genes=min(n_hvg, sub.n_vars - 1))
        g = list(sub.var_names[sub.var["highly_variable"].values])
        per[s] = g
        chosen |= set(g)
        print(f"  [panel] {s:<12} {sub.n_obs:>7,} cells -> {len(g):,} HVGs", flush=True)
        del sub
    hvg_only = sorted(chosen)
    cur = sorted((set(CURATED) & all_names) - set(hvg_only))
    panel = sorted(set(hvg_only) | set(cur))
    print(f"  [panel] HVG union {len(hvg_only):,}  + curated not already in it "
          f"{len(cur):,}  = {len(panel):,} genes", flush=True)
    missing = sorted(set(CURATED) - all_names)
    pd.DataFrame({"gene": panel,
                  "source": ["HVG" if g in set(hvg_only) else "curated" for g in panel]}
                 ).to_csv(os.path.join(out, "gene_panel.csv"), index=False)
    json.dump({"panel": panel, "per_stratum_hvg": per,
               "curated_added": cur, "curated_absent_from_data": missing},
              open(os.path.join(out, "gene_panel.json"), "w"))
    if missing:
        print(f"  [panel] curated genes absent from the data ({len(missing)}): "
              f"{', '.join(missing[:12])}{' ...' if len(missing) > 12 else ''}", flush=True)
    return panel


def subsample(sub, cap, seed=0):
    """Cap cells, stratified by cell type so rare types are not lost."""
    if sub.n_obs <= cap:
        return sub
    rng = np.random.default_rng(seed)
    lab = sub.obs["subtype"].astype(str).values
    keep = []
    types, counts = np.unique(lab, return_counts=True)
    # proportional allocation with a floor, so a 30-cell type keeps all 30
    floor = max(1, cap // (len(types) * 8))
    alloc = np.minimum(counts, np.maximum(floor, (counts / counts.sum() * cap).astype(int)))
    while alloc.sum() > cap:                      # trim the largest first
        alloc[np.argmax(alloc)] -= 1
    for t, k in zip(types, alloc):
        idx = np.where(lab == t)[0]
        keep.append(rng.choice(idx, size=min(k, len(idx)), replace=False))
    return sub[np.sort(np.concatenate(keep))].copy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blood", required=True)
    ap.add_argument("--skin", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-hvg", type=int, default=1200)
    ap.add_argument("--cell-cap", type=int, default=4000)
    ap.add_argument("--block-size", type=int, default=1024)
    ap.add_argument("--min-cells", type=int, default=400,
                    help="skip a stratum with fewer cells than this; IDS estimates were "
                         "found unstable below roughly 500-1000 cells")
    ap.add_argument("--min-cells-frac", type=float, default=0.02,
                    help="detection floor applied within each stratum (registered default 0.02)")
    ap.add_argument("--lineage", default=None,
                    help="restrict every stratum to one coarse lineage (e.g. T+NK, "
                         "Keratinocyte). Panel, detection floor and IDS are then all "
                         "computed WITHIN the lineage, which is the point: a T-cell "
                         "cytokine diluted across whole tissue falls below the floor, "
                         "while inside T cells it does not.")
    ap.add_argument("--reuse", action="store_true")
    ap.add_argument("--min-ids", type=float, default=None)
    ap.add_argument("--min-genes", type=int, default=P.PROGRAM_SIZE_MIN)
    ap.add_argument("--max-genes", type=int, default=P.PROGRAM_SIZE_MAX)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    O = lambda *p: os.path.join(args.out, *p)

    # One file per tissue (see build_disease_object.py: concatenating them is
    # not survivable at this cell count). Loaded BACKED, so only the cells and
    # genes a stratum actually needs are ever read into memory.
    src = {"PBMC": ad.read_h5ad(args.blood, backed="r"),
           "skin": ad.read_h5ad(args.skin, backed="r")}
    shared = sorted(set(src["PBMC"].var_names) & set(src["skin"].var_names))
    strata = {}
    for t, a in src.items():
        for s in sorted(a.obs["stratum"].astype(str).unique()):
            strata[s] = t
    print("[in] " + "  ".join(f"{t}={a.n_obs:,}" for t, a in src.items())
          + f"  shared genes {len(shared):,}  strata {sorted(strata)}", flush=True)

    def load_stratum(s, genes):
        """Cells of one stratum, on `genes`, as an in-memory AnnData."""
        a = src[strata[s]]
        m = np.where((a.obs["stratum"].astype(str) == s).values
                     & lineage_mask(a.obs, args.lineage))[0]
        cols = [a.var_names.get_loc(g) for g in genes]
        return ad.AnnData(X=a.X[m][:, cols], obs=a.obs.iloc[m].copy(),
                          var=pd.DataFrame(index=genes))

    pj = O("gene_panel.json")
    if args.reuse and os.path.exists(pj):
        panel = json.load(open(pj))["panel"]
        print(f"  [panel] reusing {len(panel):,} genes from {pj}", flush=True)
    else:
        panel = build_panel(src, shared, args.n_hvg, args.out, lineage=args.lineage)

    rows = []
    for s in sorted(strata):
        cC, cg = O(f"ids_{s}.npy"), O(f"ids_{s}_genes.json")
        if args.reuse and os.path.exists(cC):
            C = np.load(cC).astype(np.float64)
            genes = json.load(open(cg))
            print(f"[{s}] reusing cached IDS {C.shape}", flush=True)
            n_cells = json.load(open(O("ids_runlog.json"))).get(s, {}).get("cells", -1)
        else:
            sub = subsample(load_stratum(s, panel), args.cell_cap)
            n_cells = int(sub.n_obs)
            if n_cells < args.min_cells:
                print(f"[{s}] only {n_cells:,} cells \u2014 below the --min-cells floor "
                      f"of {args.min_cells:,}; skipped", flush=True)
                continue
            # DETECTION FLOOR, APPLIED WITHIN THE STRATUM. Without it the ranking
            # is topped by 3-gene modules of lncRNAs and pseudogenes at median IDS
            # 1.00: a gene detected in a handful of cells coincides with another
            # such gene in those same few cells, and IDS saturates. The floor is
            # the dataset's registered default (min_cells_frac, 2%).
            #
            # Applied PER STRATUM, not to the shared panel, and that distinction
            # is deliberate. The panel stays identical across strata so every
            # stratum is OFFERED the same genes; the floor then removes, in each
            # stratum, the genes that stratum cannot support an estimate for. A
            # program can therefore be tissue-specific because the gene is not
            # expressed in the other tissue -- a biological reason -- rather than
            # because panel construction never showed it to the other tissue.
            det = np.asarray((sub.X > 0).sum(axis=0)).ravel() / max(sub.n_obs, 1)
            keep = det >= args.min_cells_frac
            dropped_genes = [g for g, k in zip(sub.var_names, keep) if not k]
            print(f"[{s}] detection floor {args.min_cells_frac:.0%}: "
                  f"{len(panel):,} -> {int(keep.sum()):,} genes", flush=True)
            json.dump(dropped_genes, open(O(f"genes_below_floor_{s}.json"), "w"))
            sub = sub[:, keep].copy()
            X, genes = preprocess(sub, n_hvg=None, keep_sparse=True)
            print(f"[{s}] IDS on {X.shape[0]:,} cells x {X.shape[1]:,} genes", flush=True)
            C = compute_IDS_blocked(X, block_size=args.block_size, p_norm="max",
                                    verbose=False)
            np.save(cC, C.astype(np.float32))
            json.dump(genes, open(cg, "w"))
            del sub, X

        off = C[~np.eye(len(C), dtype=bool)]
        q = {f"p{v}": float(np.percentile(off, v)) for v in (50, 90, 99, 99.9)}
        r = {"stratum": s, "cells": n_cells, "genes": len(genes),
             **{k: round(v, 4) for k, v in q.items()}, "max": round(float(off.max()), 4)}
        rows.append(r)
        print(f"  [{s}] IDS off-diagonal " +
              "  ".join(f"{k}={v:.3f}" for k, v in q.items()) +
              f"  max={off.max():.3f}", flush=True)

        if args.min_ids is not None:
            mods = extract_modules(C, genes, min_ids=args.min_ids, min_size=2,
                                   method="greedy")
            mods = P.annotate_programs(mods, C, min_ids=args.min_ids)
            progs, dropped = P.rank_programs(mods, min_genes=args.min_genes,
                                             max_genes=args.max_genes, verbose=False)
            tbl = P.programs_table(progs)
            tbl.insert(0, "stratum", s)
            # HOW MANY CELLS ACTUALLY SUPPORT EACH PROGRAM. A high median IDS says
            # the genes move together; it does not say how many cells that claim
            # rests on. Reported so a 3-gene module resting on 20 cells cannot be
            # mistaken for one resting on 2,000.
            need = sorted({g for gs in tbl["genes"] for g in gs.split(", ")})
            cells = load_stratum(s, need)
            D = np.asarray((cells.X > 0).todense()) if hasattr(cells.X, "todense") \
                else np.asarray(cells.X > 0)
            pos = {g: i for i, g in enumerate(need)}
            n_all, n_half, frac_half = [], [], []
            for gs in tbl["genes"]:
                cols = [pos[g] for g in gs.split(", ")]
                k = D[:, cols].sum(1)
                n_all.append(int((k == len(cols)).sum()))
                n_half.append(int((k >= max(2, len(cols) / 2)).sum()))
                frac_half.append(round(float(n_half[-1]) / cells.n_obs, 4))
            tbl["cells_all_genes"] = n_all
            tbl["cells_half_genes"] = n_half
            tbl["frac_cells_half_genes"] = frac_half
            tbl["stratum_cells"] = cells.n_obs
            del cells, D
            tbl.to_csv(O(f"programs_{s}.csv"), index=False)
            pd.DataFrame([{"stratum": s, "reason": k, "n_genes": n,
                           "first_genes": ", ".join(g)}
                          for k, items in dropped.items() for n, g in items]
                         ).to_csv(O(f"programs_dropped_{s}.csv"), index=False)
            r["n_programs"] = len(progs)
            r["n_dropped_large"] = len(dropped["too_large"])
            print(f"  [{s}] {len(progs)} programs kept "
                  f"({len(dropped['too_small'])} too small, "
                  f"{len(dropped['too_large'])} too large)", flush=True)
            for _, x in tbl.head(8).iterrows():
                g = x["genes"] if len(x["genes"]) < 95 else x["genes"][:92] + "..."
                print(f"      {x['rank']:>2}. n={x['n_genes']:<3} "
                      f"medIDS={x['median_ids']:.3f} dens={x['density']:.2f}  {g}",
                      flush=True)
        del C

    df = pd.DataFrame(rows)
    df.to_csv(O("ids_summary.csv"), index=False)
    json.dump({r["stratum"]: r for r in rows}, open(O("ids_runlog.json"), "w"))
    print()
    print(df.to_string(index=False), flush=True)
    if args.min_ids is None:
        print("\n--min-ids not given. Read it off the percentiles above: pick a value\n"
              "above the background bulk (p99) and below the program mode, then rerun\n"
              "with --reuse --min-ids <value> to extract programs without recomputing IDS.")


if __name__ == "__main__":
    main()
