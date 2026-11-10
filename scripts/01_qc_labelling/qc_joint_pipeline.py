"""
qc_joint_pipeline.py — Agenda 1 / Step 1: lab QC on both tissues, gene QC per
tissue, intersect, then Harmony + UMAP.

The hypothesis being tested: skin and PBMC may be separating on the UMAP simply
because each tissue contains genes the other does not. Filtering genes inside
each tissue and keeping only the intersection removes that explanation.

Cell QC = the lab pipeline (CutaneousBioinf/single-cell_pipeline, run_QC.py):
    n_genes >= 200, pct_mt <= 10, doublet_score <= 0.15.
  * PBMC has no doublet score; demuxlet singlet is used instead, and the source
    study already required it, so it removes nothing here.
  * Skin has NO doublet information at all — not multiplexed, run_Doublet.py
    never run. Two of the three lab filters are applied there. This is the one
    real gap.

Gene QC is NOT in the lab pipeline — there is no sc.pp.filter_genes or min_cells
anywhere in it. The threshold here is our own choice, applied per tissue and
then intersected. Sizes at several thresholds are reported so the choice is
visible rather than buried.

Everything downstream (normalize 1e6 -> log1p -> HVG 2000 -> PCA 20 ->
Harmony -> neighbors -> UMAP) follows the lab pipeline's own steps.

Usage:
  python scripts/01_qc_labelling/qc_joint_pipeline.py \
      --pbmc results/01_qc_labelling/PSA26/embedded.h5ad --skin <GSE173706 skin h5ad> \
      --gene-filter 0.01 --out results/01_qc_labelling/qc-joint
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from sc_ids.viz import compute_embedding  # noqa: E402

MIN_GENES, MAX_MT = 200, 10.0
THRESHOLDS = [("3 cells", 3), ("10 cells", 10), ("0.1%", 0.001), ("1%", 0.01), ("2%", 0.02)]
TISSUE_COLORS = {"PBMC": "#1f77b4", "skin": "#ff7f0e"}
PBMC_MAIN = {"T Cell": "T cell", "Myeloid Cell": "Myeloid cell", "B Cell": "B cell",
             "NK Cell": "NK cell", "Other": "Other"}


def qc_metrics(X, var_names):
    """n_genes and pct_mt straight off the counts, the same way run_QC.py does."""
    X = sp.csr_matrix(X)
    mt = np.array([str(g).startswith("MT-") for g in var_names])
    total = np.asarray(X.sum(1)).ravel()
    n_genes = np.asarray((X > 0).sum(1)).ravel()
    mt_counts = np.asarray(X[:, mt].sum(1)).ravel() if mt.any() else np.zeros_like(total)
    pct_mt = np.where(total > 0, 100.0 * mt_counts / np.maximum(total, 1), 0.0)
    return n_genes, pct_mt, int(mt.sum())


def load_pbmc(path):
    p = ad.read_h5ad(path)
    X = sp.csr_matrix(p.layers["counts"] if "counts" in p.layers else p.X)
    q = ad.AnnData(X=X, obs=pd.DataFrame(index=[f"{b}-0" for b in p.obs_names]),
                   var=pd.DataFrame(index=p.var_names))
    main = p.obs["CellTypeMain"].astype(str).map(PBMC_MAIN).fillna("Other")
    main[p.obs["CellType"].astype(str) == "ILC"] = "ILC"
    q.obs["tissue"] = "PBMC"
    q.obs["main"] = main.values
    q.obs["subtype"] = p.obs["CellType"].astype(str).values
    q.obs["donor"] = p.obs["Subject"].astype(str).values
    q.obs["singlet"] = (p.obs["DemuxletDropletType"].astype(str) == "SNG").values
    return q


def load_skin(path):
    s = ad.read_h5ad(path)
    q = ad.AnnData(X=sp.csr_matrix(s.X), obs=pd.DataFrame(index=[f"{b}-1" for b in s.obs_names]),
                   var=pd.DataFrame(index=s.var_names))
    q.obs["tissue"] = "skin"
    q.obs["main"] = s.obs["celltype"].astype(str).values
    q.obs["subtype"] = s.obs["fullsubtype"].astype(str).values
    q.obs["donor"] = s.obs["donor"].astype(str).values
    q.obs["singlet"] = True          # no doublet call exists for this dataset
    return q


def cell_qc(a, name, log):
    n_genes, pct_mt, n_mt = qc_metrics(a.X, a.var_names)
    a.obs["n_genes"], a.obs["pct_mt"] = n_genes, pct_mt
    keep = (n_genes >= MIN_GENES) & (pct_mt <= MAX_MT) & a.obs["singlet"].values
    log[name] = dict(
        cells_in=int(a.n_obs), mt_genes_found=n_mt,
        fail_min_genes=int((n_genes < MIN_GENES).sum()),
        fail_pct_mt=int((pct_mt > MAX_MT).sum()),
        fail_singlet=int((~a.obs["singlet"].values).sum()),
        cells_kept=int(keep.sum()), pct_kept=round(100 * keep.mean(), 1))
    print(f"[qc] {name}: {a.n_obs:,} -> {keep.sum():,} ({100*keep.mean():.1f}%)", flush=True)
    return a[keep].copy()


def genes_passing(a, thr):
    n = np.asarray((sp.csr_matrix(a.X) > 0).sum(0)).ravel()
    need = thr if thr >= 1 else max(1, int(np.ceil(thr * a.n_obs)))
    return set(a.var_names[n >= need]), need


def mixing(e, keys, n=30):
    idx = {r: NearestNeighbors(n_neighbors=n + 1).fit(e.obsm[r]).kneighbors(
        e.obsm[r], return_distance=False)[:, 1:] for r in ("X_pca", "X_pca_harmony")}
    rows = []
    for k in keys:
        lab = e.obs[k].astype(str).values
        p = pd.Series(lab).value_counts(normalize=True).values
        perfect = 1 - float((p ** 2).sum())
        row = {"key": k, "perfect": round(perfect, 3)}
        for r, nm in (("X_pca", "before"), ("X_pca_harmony", "after")):
            raw = float((lab[idx[r]] != lab[:, None]).mean())
            row[f"{nm}_raw"] = round(raw, 3)
            row[f"{nm}_of_perfect"] = round(raw / perfect, 3) if perfect else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def plot(e, key, colors, out_png, title):
    rng = np.random.default_rng(0)
    order = rng.permutation(e.n_obs)
    lab = e.obs[key].astype(str).values[order]
    if colors is None:
        cats = sorted(pd.unique(e.obs[key].astype(str)))
        cmap = plt.get_cmap("tab20")
        colors = {c: matplotlib.colors.to_hex(cmap(i % 20)) for i, c in enumerate(cats)}
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    for ax, basis, t in zip(axes, ("X_umap_before", "X_umap"), ("before Harmony", "after Harmony")):
        U = e.obsm[basis][order]
        ax.scatter(U[:, 0], U[:, 1], s=1.2, c=[colors[x] for x in lab], linewidths=0, rasterized=True)
        ax.set_title(f"{title} — {t}", fontsize=11)
        ax.set_xticks([]); ax.set_yticks([])
        for s_ in ax.spines.values():
            s_.set_color("#BBBBBB")
    counts = e.obs[key].astype(str).value_counts()
    fig.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=colors[g], ms=6,
                                   label=f"{g} ({counts[g]:,})") for g in counts.index],
               loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=9)
    fig.tight_layout(); fig.savefig(out_png, dpi=140, bbox_inches="tight"); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pbmc", required=True)
    ap.add_argument("--skin", required=True)
    ap.add_argument("--gene-filter", type=float, default=0.01,
                    help="fraction of cells (<1) or absolute cell count (>=1), per tissue")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    O = lambda *p: os.path.join(args.out, *p)
    log = {"gene_filter": args.gene_filter}

    pbmc, skin = load_pbmc(args.pbmc), load_skin(args.skin)
    log["shared_genes_before_any_filter"] = len(pbmc.var_names.intersection(skin.var_names))
    pbmc = cell_qc(pbmc, "PBMC", log)
    skin = cell_qc(skin, "skin", log)

    # per-cell QC table, and what QC did to each cell type
    rows = []
    for a in (pbmc, skin):
        d = a.obs[["tissue", "main", "subtype", "n_genes", "pct_mt"]].copy()
        rows.append(d)
    pd.concat(rows).to_csv(O("cells_passing_qc.csv"))

    # ---- gene QC per tissue, several thresholds, then intersect
    tab = []
    for name, thr in THRESHOLDS + ([("chosen", args.gene_filter)] if
                                   args.gene_filter not in [t for _, t in THRESHOLDS] else []):
        gp, needp = genes_passing(pbmc, thr)
        gs, needs = genes_passing(skin, thr)
        tab.append(dict(threshold=name, pbmc_min_cells=needp, skin_min_cells=needs,
                        genes_pbmc=len(gp), genes_skin=len(gs), intersection=len(gp & gs),
                        pbmc_only=len(gp - gs), skin_only=len(gs - gp)))
    tab = pd.DataFrame(tab)
    tab.to_csv(O("gene_filter_thresholds.csv"), index=False)
    print(tab.to_string(index=False), flush=True)

    gp, _ = genes_passing(pbmc, args.gene_filter)
    gs, _ = genes_passing(skin, args.gene_filter)
    keep = sorted(gp & gs)
    json.dump({"pbmc_only": sorted(gp - gs)[:2000], "skin_only": sorted(gs - gp)[:2000],
               "kept": keep}, open(O("gene_lists.json"), "w"))
    log["genes_kept"] = len(keep)
    print(f"[genes] keeping {len(keep):,}", flush=True)

    j = ad.concat([pbmc[:, keep], skin[:, keep]], join="inner")
    for c in ("tissue", "main", "subtype", "donor"):
        j.obs[c] = j.obs[c].astype("category")
    log["joint_cells"], log["joint_genes"] = int(j.n_obs), int(j.n_vars)

    e = compute_embedding(j, batch_key="tissue", n_pcs=20, n_hvg=2000, target_sum=1e6,
                          store_raw=False, copy=False, verbose=True)
    e.obsm["X_umap_before"] = e.obsm.pop("X_umap_before_harmony")
    e.write_h5ad(O("qc_joint_embedded.h5ad"))

    mixing(e, ["tissue", "main", "subtype", "donor"]).to_csv(O("mixing_scores.csv"), index=False)
    plot(e, "tissue", TISSUE_COLORS, O("umap_tissue.png"), "after QC + gene intersection")
    plot(e, "main", None, O("umap_main.png"), "after QC + gene intersection")
    (e.obs.groupby(["tissue", "main"], observed=True).size().rename("n").reset_index()
       .to_csv(O("celltype_counts_after_qc.csv"), index=False))
    json.dump(log, open(O("run_log.json"), "w"), indent=2)
    print(json.dumps(log, indent=2), flush=True)


if __name__ == "__main__":
    main()
