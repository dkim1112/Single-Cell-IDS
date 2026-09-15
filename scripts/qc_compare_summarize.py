import os, glob, numpy as np, pandas as pd
REPO = os.path.expanduser("~/mnt/Single-Cell-IDS")
OUT = os.path.join(REPO, "results", "qc-compare"); os.makedirs(OUT, exist_ok=True)
MIN_GENES, MAX_MT = 200, 10.0

# ---------------- SKIN ----------------
meta = pd.read_csv(os.path.join(REPO, "data", "GSE173706_psorSkinMeta.tsv"),
                   sep="\t", index_col=0, low_memory=False)
suf = pd.Series([str(b).rsplit("-", 1)[-1] for b in meta.index], index=meta.index)
tab = pd.crosstab(suf, meta["orig.ident"])
assert (tab > 0).sum(1).max() == 1 and (tab > 0).sum(0).max() == 1
s2s = {c: tab.index[tab[c] > 0][0] for c in tab.columns}
annotated = set(meta.index)

parts = []
for f in sorted(glob.glob(os.path.expanduser("~/qcskin/parts/*.csv"))):
    d = pd.read_csv(f)
    s = d["sample"].iloc[0]
    bare = d["barcode"].astype(str).str.replace(r"-\d+$", "", regex=True)
    d["key"] = bare + "-" + str(s2s[s])
    d["annotated"] = d["key"].isin(annotated)
    parts.append(d)
sk = pd.concat(parts, ignore_index=True)
sk["pass_min_genes"] = sk.n_genes >= MIN_GENES
sk["pass_pct_mt"] = sk.pct_mt <= MAX_MT
sk["pass_lab"] = sk.pass_min_genes & sk.pass_pct_mt
sk.to_csv(os.path.join(OUT, "skin_percell_qc.csv.gz"), index=False)

print("=" * 62); print("SKIN — GSE173706"); print("=" * 62)
print(f"cells in GEO files              {len(sk):,}")
print(f"annotated (study kept)          {sk.annotated.sum():,}")
print(f"unannotated (study dropped)     {(~sk.annotated).sum():,}")
print(f"metadata rows matched           {sk.annotated.sum():,} / {len(meta):,}")
a = sk[sk.annotated]
print(f"\nOF THE {len(a):,} STUDY-KEPT CELLS:")
for lab, col in [("min_gene >= 200", "pass_min_genes"),
                 ("pct_mt <= 10", "pass_pct_mt"),
                 ("BOTH (lab QC)", "pass_lab")]:
    k = a[col].sum()
    print(f"  {lab:<20} {k:>7,}  {100*k/len(a):5.1f}%")
f = a[~a.pass_lab]
print(f"\n  of the {len(f):,} that fail:")
print(f"    fail mito only      {(f.pass_min_genes & ~f.pass_pct_mt).sum():>7,}")
print(f"    fail genes only     {(~f.pass_min_genes & f.pass_pct_mt).sum():>7,}")
print(f"    fail both           {(~f.pass_min_genes & ~f.pass_pct_mt).sum():>7,}")
u = sk[~sk.annotated]
print(f"\nOF THE {len(u):,} STUDY-DROPPED CELLS: lab would keep "
      f"{u.pass_lab.sum():,} ({100*u.pass_lab.mean():.1f}%)")
print(f"\nmito percentiles (study-kept): "
      + "  ".join(f"p{p}={np.percentile(a.pct_mt, p):.1f}" for p in (50, 75, 90, 95, 99)))

# per-sample
ps = a.groupby("sample").agg(n=("pass_lab", "size"), pass_lab=("pass_lab", "sum"),
                             med_mt=("pct_mt", "median")).reset_index()
ps["pct"] = 100 * ps.pass_lab / ps.n
ps.sort_values("pct").to_csv(os.path.join(OUT, "skin_qc_by_sample.csv"), index=False)
print("\nworst 5 samples:"); print(ps.sort_values("pct").head(5).to_string(index=False))
print("best 5 samples:");  print(ps.sort_values("pct").tail(5).to_string(index=False))

# ---------------- PBMC ----------------
print()
print("=" * 62); print("PBMC — GSE194315 (all 28 lanes)"); print("=" * 62)
pb = pd.read_csv(os.path.expanduser("~/qcpbmc/pbmc_percell_qc.tsv"), sep="\t",
                 header=None, names=["sample", "barcode", "n_genes", "total_counts", "pct_mt"])
pb["CellName"] = pb["sample"] + "_" + pb["barcode"].str.replace(r"-\d+$", "", regex=True)
pm = pd.read_csv(os.path.join(REPO, "data",
                 "GSE194315_CellMetadata-AS_TotalCiteseq_20220711.tsv"),
                 sep="\t", low_memory=False)
pb = pb.merge(pm[["CellName", "IncludedInStudy", "DemuxletDropletType", "CellType"]],
              on="CellName", how="left")
pb["pass_min_genes"] = pb.n_genes >= MIN_GENES
pb["pass_pct_mt"] = pb.pct_mt <= MAX_MT
pb["singlet"] = pb.DemuxletDropletType.eq("SNG")
pb["pass_lab"] = pb.pass_min_genes & pb.pass_pct_mt & pb.singlet
pb["pass_lab_nodbl"] = pb.pass_min_genes & pb.pass_pct_mt
pb.to_csv(os.path.join(OUT, "pbmc_percell_qc.csv.gz"), index=False)

print(f"cells in GEO barcode files      {len(pb):,}")
print(f"matched to metadata             {pb.IncludedInStudy.notna().sum():,}")
inc = pb[pb.IncludedInStudy == True]
print(f"IncludedInStudy == True         {len(inc):,}")
print(f"\nOF THE {len(inc):,} STUDY-KEPT CELLS:")
for lab, col in [("min_gene >= 200", "pass_min_genes"), ("pct_mt <= 10", "pass_pct_mt"),
                 ("demuxlet singlet", "singlet"),
                 ("genes+mito only", "pass_lab_nodbl"), ("ALL THREE (lab QC)", "pass_lab")]:
    k = inc[col].sum()
    print(f"  {lab:<20} {k:>7,}  {100*k/len(inc):5.1f}%")
fi = inc[~inc.pass_lab_nodbl]
print(f"\n  of the {len(fi):,} failing genes+mito:")
print(f"    fail mito only      {(fi.pass_min_genes & ~fi.pass_pct_mt).sum():>7,}")
print(f"    fail genes only     {(~fi.pass_min_genes & fi.pass_pct_mt).sum():>7,}")
print(f"    fail both           {(~fi.pass_min_genes & ~fi.pass_pct_mt).sum():>7,}")
print(f"\nmito percentiles (study-kept): "
      + "  ".join(f"p{p}={np.percentile(inc.pct_mt, p):.1f}" for p in (50, 75, 90, 95, 99)))

pls = inc.groupby("sample").agg(n=("pass_lab_nodbl", "size"),
                                keep=("pass_lab_nodbl", "sum"),
                                med_mt=("pct_mt", "median")).reset_index()
pls["pct"] = 100 * pls.keep / pls.n
pls.sort_values("pct").to_csv(os.path.join(OUT, "pbmc_qc_by_lane.csv"), index=False)
print("\nworst 5 lanes (genes+mito):"); print(pls.sort_values("pct").head(5).to_string(index=False))
print("best 5 lanes:"); print(pls.sort_values("pct").tail(5).to_string(index=False))

# side by side
sb = pd.DataFrame([
    dict(tissue="PBMC (GSE194315)", study_kept=len(inc),
         pass_genes=int(inc.pass_min_genes.sum()), pass_mito=int(inc.pass_pct_mt.sum()),
         pass_genes_mito=int(inc.pass_lab_nodbl.sum()),
         pct_genes_mito=round(100*inc.pass_lab_nodbl.mean(), 1),
         pass_all_lab=int(inc.pass_lab.sum()), pct_all_lab=round(100*inc.pass_lab.mean(), 1)),
    dict(tissue="skin (GSE173706)", study_kept=len(a),
         pass_genes=int(a.pass_min_genes.sum()), pass_mito=int(a.pass_pct_mt.sum()),
         pass_genes_mito=int(a.pass_lab.sum()),
         pct_genes_mito=round(100*a.pass_lab.mean(), 1),
         pass_all_lab=None, pct_all_lab=None),
])
sb.to_csv(os.path.join(OUT, "labqc_vs_study_both_tissues.csv"), index=False)
print("\n" + "=" * 62); print(sb.to_string(index=False))
