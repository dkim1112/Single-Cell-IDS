"""
joint_umaps.py — skin + PBMC on one UMAP: stranded cells, expanded cell sets,
per-lineage Harmony, and per-tissue subtype overlays.

Four embeddings, each computed before and after Harmony on `tissue`:

  original  T + Myeloid, both tissues. NOT recomputed: read from the existing
            joint object so the plots are exactly the ones already shown.
  expanded  PBMC T, Myeloid, B, NK, ILC  +  skin T, Myeloid, Mast
  myeloid   Myeloid only, both tissues, re-embedded from scratch
  tcell     T only, both tissues, re-embedded from scratch

Re-embedding a subset (new HVGs, new PCA, new Harmony) is the point of
`myeloid` / `tcell`: the question is whether the two lineages were competing
for the same 20 PCs and degrading each other's correction.

Recipe = sc_ids.viz.compute_embedding, the same one the original joint used:
normalize_total(1e6) -> log1p -> HVG(2000, seurat) -> PCA -> 20 PCs
-> Harmony(tissue) -> neighbors(15) -> UMAP(seed 0). Harmony corrects the
embedding only; nothing here touches expression or IDS.

"Stranded" = a cell whose 30 nearest neighbours on the UMAP are mostly (>50%)
a different main cell type. On the original embedding this agrees exactly with
the simpler "which of the two islands is it on" rule (15 = 15 cells).

Usage:
  python scripts/joint_umaps.py \
      --pbmc    results/PSA26/embedded.h5ad \
      --skin    GSE173706_skin.h5ad \
      --original results/joint/joint_embedded.h5ad \
      --pbmc-qc results/qc-compare/pbmc_percell_qc.csv.gz \
      --skin-qc results/qc-compare/skin_percell_qc.csv.gz \
      --out     results/joint-v2
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from sc_ids.viz import compute_embedding  # noqa: E402

GREY = "#D4D4D4"
MAIN_COLORS = {"T cell": "#ff7f0e", "Myeloid cell": "#1f77b4", "B cell": "#9467bd",
               "NK cell": "#2ca02c", "ILC": "#8c564b", "Mast cell": "#e377c2"}
TISSUE_COLORS = {"PBMC": "#1f77b4", "skin": "#ff7f0e"}
PBMC_MAIN = {"T Cell": "T cell", "Myeloid Cell": "Myeloid cell", "B Cell": "B cell",
             "NK Cell": "NK cell"}
T_MARK = ["CD3E", "CD3D", "CD2", "TRAC"]
MY_MARK = ["LYZ", "CST3", "AIF1", "CD68", "FCER1G", "TYROBP"]


# ------------------------------------------------------------------ building
def load_pbmc(path):
    p = ad.read_h5ad(path)
    X = p.layers["counts"] if "counts" in p.layers else p.X
    q = ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame(index=p.obs_names), var=pd.DataFrame(index=p.var_names))
    main = p.obs["CellTypeMain"].astype(str).map(PBMC_MAIN)
    main[p.obs["CellType"].astype(str) == "ILC"] = "ILC"     # ILC sits in "Other" in the map
    q.obs["tissue"] = "PBMC"
    q.obs["main"] = main.values
    q.obs["subtype"] = p.obs["CellType"].astype(str).values
    q.obs["donor"] = p.obs["Subject"].astype(str).values
    q.obs_names = [f"{b}-0" for b in p.obs_names]
    return q


def load_skin(path, keep):
    s = ad.read_h5ad(path, backed="r")
    mask = s.obs["celltype"].astype(str).isin(keep).values
    s = s[mask].to_memory()
    q = ad.AnnData(X=sp.csr_matrix(s.X), obs=pd.DataFrame(index=s.obs_names), var=pd.DataFrame(index=s.var_names))
    q.obs["tissue"] = "skin"
    q.obs["main"] = s.obs["celltype"].astype(str).values
    q.obs["subtype"] = s.obs["fullsubtype"].astype(str).values
    q.obs["donor"] = s.obs["donor"].astype(str).values
    q.obs_names = [f"{b}-1" for b in s.obs_names]
    return q


def build(pbmc, skin, pbmc_keep, skin_keep):
    a = pbmc[pbmc.obs["main"].isin(pbmc_keep).values]
    b = skin[skin.obs["main"].isin(skin_keep).values]
    genes = a.var_names.intersection(b.var_names)
    j = ad.concat([a[:, genes], b[:, genes]], join="inner")
    j.obs_names_make_unique()
    j.layers["counts"] = j.X.copy()
    for c in ("tissue", "main", "subtype", "donor"):
        j.obs[c] = j.obs[c].astype("category")
    return j


def embed(j):
    e = compute_embedding(j, batch_key="tissue", n_pcs=20, n_hvg=2000, target_sum=1e6,
                          store_raw=False, copy=True, verbose=True)
    e.obsm["X_umap_before"] = e.obsm.pop("X_umap_before_harmony")
    return e


# ------------------------------------------------------------------ measuring
def mixing(e, keys, n=30):
    """Fraction of 30 nearest neighbours from a different group. 'perfect' is the
    value expected if groups were fully intermingled: 1 - sum(p_i^2)."""
    rows = []
    idx = {r: NearestNeighbors(n_neighbors=n + 1).fit(e.obsm[r]).kneighbors(
        e.obsm[r], return_distance=False)[:, 1:] for r in ("X_pca", "X_pca_harmony")}
    for k in keys:
        lab = e.obs[k].astype(str).values
        p = pd.Series(lab).value_counts(normalize=True).values
        row = {"key": k, "perfect": round(1 - float((p ** 2).sum()), 3)}
        for r, name in (("X_pca", "before"), ("X_pca_harmony", "after")):
            raw = float((lab[idx[r]] != lab[:, None]).mean())
            row[f"{name}_raw"] = round(raw, 3)
            # share of perfect mixing: 0 = fully separated, 1 = fully intermingled.
            # This is the number the Sep-14 table reported (0.039 -> 0.322).
            row[f"{name}_of_perfect"] = round(raw / row["perfect"], 3) if row["perfect"] else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def other_tissue_neighbours(e, n=30):
    idx = NearestNeighbors(n_neighbors=n + 1).fit(e.obsm["X_pca_harmony"]).kneighbors(
        e.obsm["X_pca_harmony"], return_distance=False)[:, 1:]
    t = e.obs["tissue"].astype(str).values
    frac = (t[idx] != t[:, None]).mean(1)
    df = pd.DataFrame({"tissue": t, "main": e.obs["main"].astype(str).values, "f": frac})
    out = df.groupby(["tissue", "main"]).agg(n=("f", "size"), pct_other_tissue=("f", "mean")).reset_index()
    share = df["tissue"].value_counts(normalize=True)
    out["expected_if_mixed"] = out["tissue"].map(lambda x: 1 - share[x])
    out[["pct_other_tissue", "expected_if_mixed"]] = (100 * out[["pct_other_tissue", "expected_if_mixed"]]).round(1)
    return out


def partner_subtypes(e, n=30, min_cells=20):
    """For each subtype: what share of its neighbours (Harmony space) come from the
    other tissue, and which other-tissue subtype they mostly are. This is the check
    that correction pairs like with like (skin DC -> blood DC, not blood monocyte)."""
    idx = NearestNeighbors(n_neighbors=n + 1).fit(e.obsm["X_pca_harmony"]).kneighbors(
        e.obsm["X_pca_harmony"], return_distance=False)[:, 1:]
    t, s, m = (e.obs[c].astype(str).values for c in ("tissue", "subtype", "main"))
    rows = []
    for tis in ("skin", "PBMC"):
        for sub in pd.unique(s[t == tis]):
            ii = np.where((t == tis) & (s == sub))[0]
            if len(ii) < min_cells:
                continue
            nb = idx[ii].ravel()
            other = nb[t[nb] != tis]
            top = pd.Series(s[other]).value_counts(normalize=True)
            rows.append(dict(tissue=tis, main=m[ii[0]], subtype=sub, n=len(ii),
                             pct_neighbours_other_tissue=round(100 * len(other) / len(nb), 1),
                             top_partner_subtype=top.index[0] if len(top) else "",
                             top_partner_share=round(100 * top.iloc[0], 1) if len(top) else np.nan))
    return pd.DataFrame(rows)


def stranded(e, basis="X_umap", n=30):
    lab = e.obs["main"].astype(str).values
    idx = NearestNeighbors(n_neighbors=n + 1).fit(e.obsm[basis]).kneighbors(
        e.obsm[basis], return_distance=False)[:, 1:]
    nb = lab[idx]
    frac_other = (nb != lab[:, None]).mean(1)
    majority = pd.DataFrame(nb).mode(axis=1)[0].values
    return pd.DataFrame({"frac_other": frac_other, "neighbour_majority": majority,
                         "stranded": frac_other > 0.5}, index=e.obs_names)


def marker_verdict(e, cells):
    X = sp.csr_matrix(e.layers["counts"][e.obs_names.get_indexer(cells)])
    var = e.var_names

    def det(genes):
        g = [x for x in genes if x in var]
        return (X[:, var.get_indexer(g)] > 0).toarray().mean(1)

    t, m = det(T_MARK), det(MY_MARK)
    main = e.obs.loc[cells, "main"].astype(str).values
    out = []
    for ti, mi, lab in zip(t, m, main):
        own, oth = (ti >= .5, mi >= .5) if lab == "T cell" else (mi >= .5, ti >= .5)
        if lab not in ("T cell", "Myeloid cell"):
            out.append("n/a")
        elif own and oth:
            out.append("both lineages (doublet-like)")
        elif oth:
            out.append("other lineage (likely mislabeled)")
        elif own:
            out.append("own lineage (label fits)")
        else:
            out.append("neither (low quality)")
    return pd.DataFrame({"T_markers_detected": t.round(2), "myeloid_markers_detected": m.round(2),
                         "n_genes": np.asarray((X > 0).sum(1)).ravel(),
                         "UMIs": np.asarray(X.sum(1)).ravel().astype(int), "verdict": out}, index=cells)


# ------------------------------------------------------------------ plotting
def _ax_clean(ax, title):
    ax.set_title(title, fontsize=11)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color("#BBBBBB")


def plot_groups(e, key, colors, out_png, title):
    """Before | after, colored by `key`, drawn in RANDOM order so no group hides another."""
    rng = np.random.default_rng(0)
    order = rng.permutation(e.n_obs)
    lab = e.obs[key].astype(str).values[order]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    for ax, basis, t in zip(axes, ("X_umap_before", "X_umap"), ("before Harmony", "after Harmony")):
        U = e.obsm[basis][order]
        ax.scatter(U[:, 0], U[:, 1], s=1.5, c=[colors.get(x, "#777777") for x in lab],
                   linewidths=0, rasterized=True)
        _ax_clean(ax, f"{title} — {t}")
    counts = e.obs[key].astype(str).value_counts()
    handles = [plt.Line2D([], [], marker="o", ls="", color=colors.get(g, "#777777"), ms=6,
                          label=f"{g} ({counts[g]:,})") for g in counts.index]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)


# 26 distinct, saturated hues (no greys — grey is reserved for the other tissue).
SUBTYPE_HUES = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
    "#17becf", "#bcbd22", "#393b79", "#e6550d", "#31a354", "#843c39", "#7b4173",
    "#3182bd", "#fd8d3c", "#74c476", "#ad494a", "#a55194", "#6baed6", "#b5cf6b",
    "#e7ba52", "#ce6dbd", "#637939", "#9c9ede", "#d6616b",
]


def subtype_palette(obs, tissue):
    """One fixed colour per subtype, the same in every embedding (colour follows
    the subtype, never its position). Ordered by lineage so the legend groups."""
    d = obs[obs.tissue.astype(str) == tissue][["main", "subtype"]].astype(str).drop_duplicates()
    subs = d.sort_values(["main", "subtype"])["subtype"].drop_duplicates().tolist()
    if len(subs) > len(SUBTYPE_HUES):
        raise ValueError(f"{len(subs)} subtypes > {len(SUBTYPE_HUES)} colours")
    return dict(zip(subs, SUBTYPE_HUES))


def plot_subtypes(e, tissue, pal, out_png, title, min_label=40):
    """Other tissue in light grey underneath; this tissue's subtypes on top."""
    rng = np.random.default_rng(0)
    is_t = (e.obs["tissue"].astype(str) == tissue).values
    sub = e.obs["subtype"].astype(str).values
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6))
    for ax, basis, t in zip(axes, ("X_umap_before", "X_umap"), ("before Harmony", "after Harmony")):
        U = e.obsm[basis]
        bg = np.where(~is_t)[0]
        ax.scatter(U[bg, 0], U[bg, 1], s=1.2, c=GREY, linewidths=0, rasterized=True, zorder=1)
        fg = rng.permutation(np.where(is_t)[0])
        ax.scatter(U[fg, 0], U[fg, 1], s=2.2, c=[pal[x] for x in sub[fg]], linewidths=0,
                   rasterized=True, zorder=2)
        # direct labels at each subtype's median position
        for s_ in pd.unique(sub[is_t]):
            m = is_t & (sub == s_)
            if m.sum() < min_label:
                continue
            x, y = np.median(U[m, 0]), np.median(U[m, 1])
            ax.text(x, y, s_, fontsize=7, ha="center", va="center", zorder=3, color="#1A1A1A",
                    path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])
        _ax_clean(ax, f"{title} — {tissue} subtypes, {t}")
    counts = pd.Series(sub[is_t]).value_counts()
    main_of = e.obs.loc[is_t, ["subtype", "main"]].astype(str).drop_duplicates().set_index("subtype")["main"]
    order = sorted(counts.index, key=lambda s: (main_of.get(s, ""), s))
    other = "skin" if tissue == "PBMC" else "PBMC"
    handles = [plt.Line2D([], [], marker="o", ls="", color=GREY, ms=6,
                          label=f"{other} (all, {int((~is_t).sum()):,})")]
    handles += [plt.Line2D([], [], marker="o", ls="", color=pal[s], ms=6,
                           label=f"{s} ({counts[s]:,})") for s in order]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False,
               fontsize=8, title=f"{tissue} subtype (n)", title_fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_stranded(e, st, out_png):
    rng = np.random.default_rng(0)
    order = rng.permutation(e.n_obs)
    U = e.obsm["X_umap"]
    lab = e.obs["main"].astype(str).values
    fig, ax = plt.subplots(figsize=(7.5, 5.6))
    ax.scatter(U[order, 0], U[order, 1], s=1.5, c=[MAIN_COLORS[x] for x in lab[order]],
               linewidths=0, rasterized=True, zorder=1)
    s = st[st.stranded]
    for tis, mk in (("PBMC", "o"), ("skin", "s")):
        c = s.index[s.tissue == tis]
        ii = e.obs_names.get_indexer(c)
        ax.scatter(U[ii, 0], U[ii, 1], s=55, marker=mk, facecolors="none", edgecolors="black",
                   linewidths=1.0, zorder=3)
    _ax_clean(ax, f"Stranded cells, after Harmony ({len(s)} of {e.n_obs:,})")
    handles = [plt.Line2D([], [], marker="o", ls="", color=MAIN_COLORS[m], ms=6, label=m)
               for m in ("T cell", "Myeloid cell")]
    handles += [plt.Line2D([], [], marker="o", ls="", mfc="none", mec="black", ms=8, label="stranded — PBMC"),
                plt.Line2D([], [], marker="s", ls="", mfc="none", mec="black", ms=8, label="stranded — skin")]
    ax.legend(handles=handles, loc="best", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pbmc", required=True)
    ap.add_argument("--skin", required=True)
    ap.add_argument("--original", required=True)
    ap.add_argument("--pbmc-qc", default=None)
    ap.add_argument("--skin-qc", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--replot", action="store_true",
                    help="reuse <out>/<name>_embedded.h5ad instead of re-embedding")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    O = lambda *p: os.path.join(args.out, *p)
    log = {}

    # ---- original (read, not recomputed)
    print("[original] reading", args.original, flush=True)
    orig = ad.read_h5ad(args.original)
    if "X_umap_before" not in orig.obsm:
        orig.obsm["X_umap_before"] = orig.obsm["X_umap_before_harmony"]

    # ---- Agenda 1a: stranded cells on the original after-Harmony UMAP
    st = stranded(orig)
    st = st.join(orig.obs[["tissue", "main", "subtype"]].astype(str))
    s = st[st.stranded].copy()
    s = s.join(marker_verdict(orig, s.index))
    if args.pbmc_qc and args.skin_qc:
        pq = pd.read_csv(args.pbmc_qc, low_memory=False)
        pq = pq[pq["sample"] == "PBMC-01-1"].set_index("barcode")["pct_mt"]
        sq = pd.read_csv(args.skin_qc).set_index("key")["pct_mt"]
        s["pct_mt"] = [(pq if t == "PBMC" else sq).get(i.rsplit("-", 1)[0], np.nan)
                       for i, t in zip(s.index, s.tissue)]
    s.sort_values(["main", "tissue", "subtype"]).to_csv(O("stranded_cells.csv"))
    summ = (s.groupby(["main", "tissue", "subtype"]).size().rename("n_stranded").reset_index()
            .merge(orig.obs.astype(str).groupby(["main", "tissue", "subtype"]).size()
                   .rename("n_total").reset_index(), on=["main", "tissue", "subtype"]))
    summ["pct_of_subtype"] = (100 * summ.n_stranded / summ.n_total).round(2)
    summ.sort_values("n_stranded", ascending=False).to_csv(O("stranded_by_subtype.csv"), index=False)
    plot_stranded(orig, st, O("umap_after_stranded.png"))
    log["stranded_original"] = int(st.stranded.sum())
    print(f"[stranded] {st.stranded.sum()} cells", flush=True)

    # ---- build the other three
    pbmc = load_pbmc(args.pbmc)
    skin = load_skin(args.skin, keep={"T cell", "Myeloid cell", "Mast cell"})
    sets = {
        "expanded": (["T cell", "Myeloid cell", "B cell", "NK cell", "ILC"], ["T cell", "Myeloid cell", "Mast cell"]),
        "myeloid": (["Myeloid cell"], ["Myeloid cell"]),
        "tcell": (["T cell"], ["T cell"]),
    }
    emb = {"original": orig}
    for name, (pk, skk) in sets.items():
        cached = O(f"{name}_embedded.h5ad")
        if args.replot and os.path.exists(cached):
            emb[name] = ad.read_h5ad(cached)
            continue
        print(f"[{name}] building", flush=True)
        j = build(pbmc, skin, pk, skk)
        print(f"[{name}] {j.n_obs:,} cells x {j.n_vars:,} genes", flush=True)
        emb[name] = embed(j)
        emb[name].write_h5ad(cached)

    # ---- measurements + standard plots for every embedding
    all_obs = pd.concat([x.obs[["tissue", "main", "subtype"]] for x in emb.values()])
    pals = {tis: subtype_palette(all_obs, tis) for tis in ("PBMC", "skin")}
    mix_all, nb_all, strand_all, partner_all = [], [], [], []
    for name, e in emb.items():
        keys = ["tissue", "main", "subtype", "donor"] if e.obs["main"].nunique() > 1 else ["tissue", "subtype", "donor"]
        m = mixing(e, keys); m.insert(0, "embedding", name); mix_all.append(m)
        nb = other_tissue_neighbours(e); nb.insert(0, "embedding", name); nb_all.append(nb)
        pt = partner_subtypes(e); pt.insert(0, "embedding", name); partner_all.append(pt)
        if e.obs["main"].nunique() > 1:
            sti = stranded(e)
            if name != "original":   # original's breakdown is written above, with marker verdicts
                os.makedirs(O(name), exist_ok=True)
                x = sti[sti.stranded].join(e.obs[["tissue", "main", "subtype"]].astype(str))
                (x.groupby(["main", "tissue", "subtype", "neighbour_majority"]).size()
                   .rename("n").reset_index().sort_values("n", ascending=False)
                   .to_csv(O(name, "stranded_by_subtype.csv"), index=False))
            strand_all.append({"embedding": name, "n_cells": e.n_obs, "n_stranded": int(sti.stranded.sum()),
                               "pct": round(100 * sti.stranded.mean(), 3)})
        d = O(name); os.makedirs(d, exist_ok=True)
        if name != "original":
            plot_groups(e, "tissue", TISSUE_COLORS, os.path.join(d, "umap_tissue.png"), name)
            if e.obs["main"].nunique() > 1:
                plot_groups(e, "main", MAIN_COLORS, os.path.join(d, "umap_main.png"), name)
        # Agenda 3: subtype overlays, one file per tissue
        for tis in ("PBMC", "skin"):
            plot_subtypes(e, tis, pals[tis], os.path.join(d, f"umap_subtypes_{tis}.png"), name)
        print(f"[{name}] plotted", flush=True)
    pd.concat(mix_all).to_csv(O("mixing_scores.csv"), index=False)
    pd.concat(nb_all).to_csv(O("other_tissue_neighbours.csv"), index=False)
    pd.DataFrame(strand_all).to_csv(O("stranded_counts_by_embedding.csv"), index=False)
    pd.concat(partner_all).to_csv(O("partner_subtypes.csv"), index=False)
    log["cells"] = {k: int(v.n_obs) for k, v in emb.items()}
    json.dump(log, open(O("run_log.json"), "w"), indent=2)
    print("done", flush=True)


if __name__ == "__main__":
    main()
