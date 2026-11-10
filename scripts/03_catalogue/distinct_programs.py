"""
distinct_programs.py — from the 194 detections to the distinct gene programs.

A "detection" is one gene program found in one search (analysis x tissue x disease state). The
same gene set is detected again in other conditions, so the 194 detections contain far fewer
distinct programs. This script applies three gates, in order, to the detections, then merges repeats.

  Gate 0  technical   mitochondrial / ribosomal / haemoglobin / platelet / sex-linked / housekeeping genes
                      (flag from program_catalogue.py, plus a gene-content check: ribosomal + mitochondrial
                      genes >= 25% of the set, or XIST in a set of <= 6 genes). Applied once, here.
  Gate 1  permutation fraction of within-program gene pairs that beat all 60 shuffles (p = 0 of 60; the catalogue column is named p < 0.01) must be >= PERM_CUT
                      (authors' test; the cutoff is set below, default 0.9; PERM_CUT=1 means every pair)
  Gate 2  activity    mean score >= 0.2 in at least one cell subtype (>= 50 cells) of either tissue
  Merge   detections whose gene sets overlap (Jaccard >= 0.5) are one distinct program.
  Specificity is not a gate: tau (0 = equal everywhere, 1 = one group) is reported per axis.
"""
import os, sys, json
import numpy as np, pandas as pd
import marker_rules as cw
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from sc_ids import paths                                                      # noqa: E402
from sc_ids.lineage import LINEAGE                                           # noqa: E402

OUT = paths.CATALOGUE
ACTIVE = 0.2
JACC = 0.5
PERM_CUT = float(os.environ.get("PERM_CUT", 0.9))

IDENTITY = {"T-cell identity / receptor signalling", "B-cell identity", "conventional dendritic cell (cDC2)",
            "plasmacytoid dendritic cell", "monocyte / macrophage core", "melanocyte identity / pigmentation",
            "mast cell identity", "lymphatic endothelium", "vascular endothelium", "fibroblast / extracellular matrix",
            "smooth muscle / pericyte contractility", "sweat gland / secretory", "nerve / Schwann cell",
            "basal keratinocyte / adhesion to basement membrane"}
STRESS = {"AP-1 / immediate-early activation", "heat-shock / protein folding (also a dissociation artefact)"}


def tau(x):
    x = np.clip(np.asarray(x, float), 0, None)
    if len(x) < 2 or x.max() <= 0:
        return np.nan
    return float((1 - x / x.max()).sum() / (len(x) - 1))


def content_technical(genes):
    gs = [g.strip() for g in genes.split(",")]
    n = sum(g.startswith(("RPL", "RPS", "MT-", "MTRNR")) for g in gs)
    return n / len(gs) >= 0.25 or ("XIST" in gs and len(gs) <= 6)


def run(cut, write=True):
    c = pd.read_csv(f"{OUT}/program_catalogue.csv")
    act = pd.read_csv(f"{OUT}/program_activity_long.csv.gz")
    con = pd.read_csv(f"{OUT}/disease_contrasts.csv")
    sub = act[(act.axis == "subtype") & (act.n_cells >= 50)]

    c["technical"] = c.technical_flag.notna() | c.genes.map(content_technical)
    c["perm_ok"] = c["frac_pairs_beat_all_shuffles"] >= cut
    mx = sub.groupby("program_id").mean_score.max()
    c["active"] = c.program_id.map(mx).fillna(-9) >= ACTIVE
    c["status"] = np.where(c.technical, "removed: technical",
                  np.where(~c.perm_ok, "removed: failed permutation test",
                  np.where(~c.active, "removed: not active (score below 0.2)", "kept")))
    keep = c[c.status == "kept"].reset_index(drop=True)

    G = [set(x.strip() for x in g.split(",")) for g in keep.genes]
    n = len(keep); parent = list(range(n))
    def f(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    for i in range(n):
        for j in range(i + 1, n):
            if len(G[i] & G[j]) / len(G[i] | G[j]) >= JACC:
                parent[f(i)] = f(j)
    keep["grp"] = [f(i) for i in range(n)]

    DM = pd.read_csv(f"{OUT}/donor_means.csv.gz"); DM["lane"] = DM["lane"].astype(str)
    DM = DM[DM.lineage != "Other"]; DM["w"] = DM.score * DM.n_cells
    DD = DM.groupby(["tissue", "program_id", "lineage", "disease", "donor"], as_index=False)[["w", "n_cells"]].sum()
    DD["score"] = DD.w / DD.n_cells
    ORDER = {"blood": ["Healthy", "PSO", "PSA"], "skin": ["NS", "PN", "PP"]}

    def dz(ids, tissue, force=None):
        x = DD[(DD.tissue == tissue) & DD.program_id.isin(ids)]
        nd = x.groupby("lineage").donor.nunique()
        x = x[x.lineage.isin(nd[nd >= 5].index)]
        if x.empty:
            return "", np.nan
        top = force if force in set(x.lineage) else x.groupby("lineage").score.mean().idxmax()
        y = x[x.lineage == top]; parts, means = [], []
        for d in ORDER[tissue]:
            v = y[y.disease == d].groupby("donor").score.mean()
            means.append(v.mean() if len(v) else np.nan)
            parts.append(f"{d} {v.mean():.2f} (n={len(v)})" if len(v) else f"{d} -")
        return f"{top}: " + " | ".join(parts), tau([m for m in means if not np.isnan(m)])

    rows = []
    for _, m in keep.groupby("grp"):
        ids = set(m.program_id)
        gl = [set(x.strip() for x in g.split(",")) for g in m.genes]
        genes = sorted(set().union(*gl))
        labs = [cw.annotate(g)[0] for g in m.genes]
        labs = [l for l in labs if l] or ([cw.annotate(", ".join(genes))[0]] if cw.annotate(", ".join(genes))[0] else [])
        label = pd.Series(labs).mode().iloc[0] if labs else ""
        a = sub[sub.program_id.isin(ids)]
        best = a.groupby(["tissue", "group"]).mean_score.max().reset_index()
        best["lineage"] = best.group.map(LINEAGE).fillna("Other")
        on = best[best.mean_score >= ACTIVE]
        tis = sorted(on.tissue.unique())
        tissue = "both" if len(tis) == 2 else tis[0] + " only"
        cts = lambda t: "; ".join(f"{r.group} ({r.mean_score:.2f})" for r in on[on.tissue == t].sort_values("mean_score", ascending=False).itertuples())
        lins = sorted(on.lineage.unique())
        d = con[con.program_id.isin(ids) & (con.tissue == "skin") & (con.p_adj < 0.05) & (con.delta.abs() >= ACTIVE)]
        d = d.sort_values("delta", key=abs, ascending=False).drop_duplicates(["lineage", "contrast"])
        ms, ts = dz(ids, "skin", d.lineage.iloc[0] if len(d) else None)
        mb, tb = dz(ids, "blood")
        tt = tau([best[best.tissue == t].mean_score.max() if (best.tissue == t).any() else 0 for t in ("blood", "skin")])
        tc = tau(best.groupby("group").mean_score.max().values)
        rows.append(dict(
            n_detections=len(m), detected_in=", ".join(sorted(set(m.stratum))), analyses=", ".join(sorted(set(m.analysis))),
            n_genes=len(genes), genes=", ".join(genes), label_rule=label, tissue_call=tissue,
            celltypes_blood=cts("blood"), celltypes_skin=cts("skin"), n_subtypes=on.group.nunique(), n_lineages=len(lins), lineages=", ".join(lins),
            tau_tissue=round(tt, 2), tau_celltype=round(tc, 2), tau_disease_skin=round(ts, 2) if ts == ts else np.nan,
            tau_disease_blood=round(tb, 2) if tb == tb else np.nan,
            skin_disease_change="; ".join(f"{r.lineage}: {r.contrast} {r.delta:+.2f}" for r in d.itertuples()),
            disease_means_skin=ms, disease_means_blood=mb,
            median_IDS=round(float(m.median_IDS.median()), 2), member_ids=" | ".join(sorted(ids))))
    F = pd.DataFrame(rows)
    F["celltype_call"] = np.where(F.n_lineages >= 3, "broad (3+ lineages)", np.where(F.n_lineages == 2, "two lineages",
                          np.where(F.n_subtypes >= 2, "one lineage, several subtypes", "one subtype")))
    F["category"] = np.where(F.label_rule.isin(IDENTITY), "cell-identity",
                     np.where(F.label_rule.isin(STRESS), "stress / possible handling artefact",
                     np.where(F.label_rule.fillna("") == "", "unlabelled", "functional")))
    F = F.sort_values(["tissue_call", "n_lineages"], ascending=[True, False]).reset_index(drop=True)
    F.insert(0, "program_no", [f"GP{i+1:02d}" for i in range(len(F))])

    chain = dict(perm_cut=cut, detections=len(c), technical=int(c.technical.sum()),
                 after_technical=int((~c.technical).sum()),
                 failed_permutation=int((c.status == "removed: failed permutation test").sum()),
                 after_permutation=int((~c.technical & c.perm_ok).sum()),
                 not_active=int((c.status == "removed: not active (score below 0.2)").sum()),
                 kept_detections=len(keep), distinct_programs=len(F))
    if write:
        print(json.dumps(chain))
    if write:
        print(F.tissue_call.value_counts().to_dict(), F.category.value_counts().to_dict())
    if write:
        pid2gp = {p: r.program_no for r in F.itertuples() for p in r.member_ids.split(" | ")}
        c["program_no"] = c.program_id.map(pid2gp).fillna("")
        c.to_csv(f"{OUT}/detections_status.csv", index=False)
        F.to_csv(f"{OUT}/distinct_programs.csv", index=False)
        json.dump(chain, open(f"{OUT}/selection_chain.json", "w"))
    return chain, F


def main():
    chain, F = run(PERM_CUT)
    alt = 1.0 if PERM_CUT < 1 else 0.9
    ach, _ = run(alt, write=False)
    chain["alt_cut"], chain["alt_distinct_programs"] = alt, ach["distinct_programs"]
    json.dump(chain, open(f"{OUT}/selection_chain.json", "w"))


if __name__ == "__main__":
    main()
