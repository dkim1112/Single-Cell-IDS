"""
program_catalogue.py — the systematic catalogue of IDS gene programs.

One row per program, with everything needed to decide whether it is worth
reporting: its genes, where it is active (cell type, tissue, disease state), how
concentrated that activity is, whether it is shared or distinct, whether its
within-program dependence survives permutation testing, and whether it is
biology or a technical artefact.

FOUR ANALYSES ARE POOLED, and the catalogue records which one each program came
from, because they answer different questions:
  whole-tissue   programs visible across a whole tissue, dominated by the
                 abundant populations
  T+NK, Myeloid, Keratinocyte
                 programs found INSIDE one lineage, where a gene expressed only
                 by that lineage is no longer diluted below the detection floor

THREE SPECIFICITY SCORES, one per axis, all the same index so they are
comparable. Tau = sum(1 - x_i / x_max) / (n - 1) over the group means, which is
0 when a program is equally active everywhere and 1 when it is active in exactly
one group. Reported alongside the top group, because a high tau is only
interesting once you know what it is specific TO.

PERMUTATION P-VALUES use the authors' own `ids.compute_IDS(p_val=True)`, which is
the filter their paper provides and which our blocked path could not run. They
are computed on the union of each stratum's program genes rather than the whole
matrix. That is exact, not an approximation: the IDS between two genes depends
only on those two genes' feature maps, so removing other genes from the matrix
cannot change their pairwise value or its null distribution.

TECHNICAL FLAGS. A high-IDS program is not automatically biology. Mitochondrial,
ribosomal, haemoglobin, platelet and sex-linked gene sets recur across strata and
reflect ambient RNA, library composition or donor sex rather than disease. Note
in particular that removing platelet-ANNOTATED cells does not remove platelet
ambient transcripts, which still form a coherent programme inside T cells.

Usage:
  python scripts/03_catalogue/program_catalogue.py --out results/03_catalogue
"""
import argparse, glob, os, sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
import ids                                                          # noqa: E402
from sc_ids import paths                                            # noqa: E402
from sc_ids.lineage import LINEAGE                                  # noqa: E402

ANALYSES = {"whole-tissue": paths.DISEASE_RUN, **paths.LINEAGE_RUNS}

# --- technical signatures. Membership rule, then the share of the programme's
# genes that match; a programme is flagged when the share is at least half.
TECHNICAL = {
    "mitochondrial": lambda g: g.startswith("MT-") or g.startswith("MTRNR"),
    "ribosomal": lambda g: g.startswith("RPL") or g.startswith("RPS"),
    "haemoglobin": lambda g: g in {"HBA1", "HBA2", "HBB", "HBD", "HBM", "ALAS2"},
    "platelet (ambient)": lambda g: g in {"PF4", "PPBP", "NRGN", "SDPR", "TUBB1",
                                          "GNG11", "RGS18", "CLU", "PRKAR2B",
                                          "ACRBP", "CAVIN2", "TMEM40", "SPARC"},
    "sex-linked (donor)": lambda g: g in {"XIST", "RPS4Y1", "DDX3Y", "UTY", "USP9Y",
                                          "EIF1AY", "KDM5D", "NLGN4Y", "TXLNGY"},
    "cytoskeleton/housekeeping": lambda g: g in {"ACTB", "ACTG1", "PFN1", "TMSB4X",
                                                 "TMSB10", "B2M", "GAPDH", "EEF1A1"},
    "dissociation stress": lambda g: g in {"HSPA1A", "HSPA1B", "HSPA6", "HSPA8",
                                           "HSPH1", "HSP90AA1", "DNAJB1", "DNAJA1",
                                           "HSPB1", "HSPE1"},
}

EXPRESSING = 0.0
MIN_CELLS = 20


def tau(x):
    """Tissue-specificity index on group means. 0 = everywhere, 1 = one group."""
    x = np.asarray([v for v in x if np.isfinite(v)], dtype=float)
    x = np.clip(x, 0, None)                 # score <= 0 means "not above background"
    if len(x) < 2 or x.max() <= 0:
        return np.nan
    return float((1 - x / x.max()).sum() / (len(x) - 1))


def classify(genes):
    hits = {}
    for name, rule in TECHNICAL.items():
        share = np.mean([bool(rule(g)) for g in genes])
        if share >= 0.5:
            hits[name] = round(float(share), 2)
    if not hits:
        return "", ""
    best = max(hits, key=hits.get)
    return best, "; ".join(f"{k} {v:.0%}" for k, v in sorted(hits.items(),
                                                             key=lambda kv: -kv[1]))


def load_all():
    rows = []
    for analysis, d in ANALYSES.items():
        for f in sorted(glob.glob(os.path.join(d, "programs_*.csv"))):
            if "dropped" in f:
                continue
            t = pd.read_csv(f)
            if t.empty:
                continue
            for _, r in t.iterrows():
                genes = [g.strip() for g in r["genes"].split(",")]
                tech, tech_detail = classify(genes)
                rows.append({
                    "program_id": f"{analysis}|{r['stratum']}|P{int(r['rank']):02d}",
                    "analysis": analysis,
                    "found_in_tissue": "blood" if r["stratum"].startswith("PBMC") else "skin",
                    "found_in_disease": {"Healthy": "control", "PSO": "psoriasis",
                                         "PSA": "PsA", "NS": "control",
                                         "PN": "psoriasis (non-lesional)",
                                         "PP": "psoriasis (lesional)"}[
                                             r["stratum"].split("_", 1)[1]],
                    "stratum": r["stratum"],
                    "n_genes": int(r["n_genes"]),
                    "median_IDS": round(float(r["median_ids"]), 3),
                    "edge_density": round(float(r["density"]), 2),
                    "cells_supporting": int(r.get("cells_half_genes", -1)),
                    "stratum_cells": int(r.get("stratum_cells", -1)),
                    "technical_flag": tech,
                    "technical_detail": tech_detail,
                    "genes": ", ".join(genes),
                    "_genes": genes})
    return rows


def permutation_pvalues(rows, out, n_tests=100):
    """Authors' permutation test, per stratum, on the union of programme genes.

    CHECKPOINTED. This environment recycles the container roughly every 20 minutes,
    so any step that cannot finish inside that window never finishes at all. Each
    stratum writes its own file and is skipped if that file already exists, so a
    restart resumes instead of starting over.
    """
    pdir = os.path.join(out, "perm")
    os.makedirs(pdir, exist_ok=True)
    by = {}
    for r in rows:
        by.setdefault((r["analysis"], r["stratum"]), []).append(r)

    def cpath(analysis, stratum):
        return os.path.join(pdir, f"{analysis.replace('+','-')}__{stratum}.csv")

    todo = [k for k in sorted(by) if not os.path.exists(cpath(*k))]
    for k in sorted(by):
        if os.path.exists(cpath(*k)):
            done = (pd.read_csv(cpath(*k))
                    .rename(columns={"frac_pairs_p_lt_0.01": "frac_pairs_beat_all_shuffles"})  # pre-rename checkpoints
                    .set_index("program_id"))
            for r in by[k]:
                if r["program_id"] in done.index:
                    for c in ("perm_p_median", "perm_p_max", "frac_pairs_beat_all_shuffles"):
                        r[c] = done.loc[r["program_id"], c]
    if not todo:
        print(f"  [perm] all {len(by)} strata already done", flush=True)
        return
    print(f"  [perm] {len(todo)} strata to do, {len(by) - len(todo)} cached", flush=True)
    src = {}
    for (analysis, stratum), rs in [(k, by[k]) for k in todo]:
        t = "PBMC" if stratum.startswith("PBMC") else "skin"
        if t not in src:
            src[t] = ad.read_h5ad(paths.DISEASE_RUN + "/"
                                  + ("blood_qc.h5ad" if t == "PBMC" else "skin_qc.h5ad"))
            if src[t].X.dtype != np.float32:
                src[t].X = src[t].X.astype(np.float32)
            sc.pp.normalize_total(src[t], target_sum=1e4)
        a = src[t]
        m = (a.obs["stratum"].astype(str) == stratum).values
        if analysis != "whole-tissue":
            want = set(analysis.split("+"))
            lin = a.obs["subtype"].astype(str).map(LINEAGE).fillna("Other")
            m &= lin.isin(want).values
        genes = sorted({g for r in rs for g in r["_genes"]})
        genes = [g for g in genes if g in a.var_names]
        idx = np.where(m)[0]
        if len(idx) > 4000:
            idx = np.random.default_rng(0).choice(idx, 4000, replace=False)
        X = a[idx][:, genes].X
        X = np.asarray(X.todense() if hasattr(X, "todense") else X, dtype=np.float64)
        lo, span = X.min(0, keepdims=True), X.max(0, keepdims=True) - X.min(0, keepdims=True)
        X = (X - lo) / np.where(span == 0, 1, span) * 8.0      # same [0,8] as discovery
        C, Pv = ids.compute_IDS(X, num_terms=6, p_norm="max", p_val=True,
                                num_tests=n_tests, bandwidth_term=0.5)
        pos = {g: i for i, g in enumerate(genes)}
        for r in rs:
            gi = [pos[g] for g in r["_genes"] if g in pos]
            if len(gi) < 2:
                continue
            sub = Pv[np.ix_(gi, gi)]
            iu = np.triu_indices_from(sub, k=1)
            vals = sub[iu]
            r["perm_p_median"] = round(float(np.median(vals)), 4)
            r["perm_p_max"] = round(float(vals.max()), 4)
            r["frac_pairs_beat_all_shuffles"] = round(float((vals < 0.01).mean()), 3)
        pd.DataFrame([{"program_id": r["program_id"],
                       "perm_p_median": r.get("perm_p_median"),
                       "perm_p_max": r.get("perm_p_max"),
                       "frac_pairs_beat_all_shuffles": r.get("frac_pairs_beat_all_shuffles")}
                      for r in rs]).to_csv(cpath(analysis, stratum), index=False)
        print(f"  [perm] {analysis} {stratum}: {len(genes)} genes, "
              f"{len(rs)} programmes \u2192 cached", flush=True)
    del src


def score_and_summarise(rows, out, block=40):
    """Also checkpointed: each tissue writes its summary once computed, and within
    a tissue the programmes are scored in blocks so a restart loses one block."""
    long = []
    for tag, f in (("blood", os.path.join(paths.DISEASE_RUN, "blood_qc.h5ad")),
                   ("skin", os.path.join(paths.DISEASE_RUN, "skin_qc.h5ad"))):
        cp = os.path.join(out, f"activity_{tag}.csv.gz")
        if os.path.exists(cp):
            print(f"  [score] {tag}: cached", flush=True)
            long.append(pd.read_csv(cp))
            continue
        a = ad.read_h5ad(f)
        if a.X.dtype != np.float32:
            a.X = a.X.astype(np.float32)
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        names, scores = [], {}
        for i, r in enumerate(rows):
            g = [x for x in r["_genes"] if x in a.var_names]
            if len(g) < 2:
                continue
            sc.tl.score_genes(a, g, score_name="_s", ctrl_size=50, random_state=0)
            scores[r["program_id"]] = a.obs.pop("_s").values
            names.append(r["program_id"])
            if (i + 1) % block == 0:
                print(f"  [score] {tag}: {len(names)} programmes so far", flush=True)
        for p, v in scores.items():
            a.obs[p] = v
        del scores
        print(f"  [score] {tag}: {len(names)} programmes", flush=True)
        o = a.obs
        lin = o["subtype"].astype(str).map(LINEAGE).fillna("Other")
        key = pd.DataFrame({"lineage": lin.values,
                            "subtype": o["subtype"].astype(str).values,
                            "disease": o["disease"].astype(str).values,
                            "donor": o["donor"].astype(str).values})
        this = []
        for p in names:
            key["s"] = a.obs[p].values
            for axis, cols in (("lineage", ["lineage"]), ("subtype", ["subtype"]),
                               ("disease", ["disease"])):
                g = key.groupby(cols, observed=True)["s"].agg(
                    ["mean", "size", lambda x: (x > EXPRESSING).mean()])
                g.columns = ["mean_score", "n_cells", "prop_expressing"]
                g = g[g["n_cells"] >= MIN_CELLS].reset_index()
                g["program_id"], g["tissue"], g["axis"] = p, tag, axis
                g = g.rename(columns={cols[0]: "group"})
                this.append(g)
        del a
        T = pd.concat(this, ignore_index=True)
        T.to_csv(cp, index=False, compression="gzip")      # checkpoint this tissue
        print(f"  [score] {tag}: cached to {cp}", flush=True)
        long.append(T)
    L = pd.concat(long, ignore_index=True)
    L.to_csv(os.path.join(out, "program_activity_long.csv.gz"), index=False,
             compression="gzip")
    return L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=paths.CATALOGUE)
    ap.add_argument("--perm-tests", type=int, default=100)
    ap.add_argument("--skip-perm", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    rows = load_all()
    print(f"[in] {len(rows)} programmes across {len(ANALYSES)} analyses")
    print(pd.Series([r["analysis"] for r in rows]).value_counts().to_string())
    print(f"  technical-flagged: {sum(1 for r in rows if r['technical_flag'])}")

    if not args.skip_perm:
        permutation_pvalues(rows, args.out, n_tests=args.perm_tests)

    L = score_and_summarise(rows, args.out)

    # ---- fold the activity summary back into one row per programme
    for r in rows:
        p = r["program_id"]
        sub = L[L["program_id"] == p]
        s = sub[sub["axis"] == "lineage"]
        if not s.empty:
            grp = s.groupby("group", observed=True)["mean_score"].mean()
            r["tau_celltype"] = round(tau(grp.values), 3)
            r["top_celltype"] = grp.idxmax()
            r["top_celltype_score"] = round(float(grp.max()), 2)
            act = grp[grp > 0.2]
            r["celltypes_active"] = "; ".join(sorted(act.index.astype(str)))
            r["n_celltypes_active"] = int(len(act))

        # DISEASE AXIS IS COMPUTED WITHIN TISSUE, NEVER POOLED. Skin has no PsA
        # samples, so pooling lets a skin programme score its maximum in blood PsA
        # cells and be reported as "PsA-specific" -- which would be false.
        sd = sub[sub["axis"] == "disease"]
        n_active_total = 0
        for t in ("blood", "skin"):
            st = sd[sd["tissue"] == t]
            if st.empty:
                continue
            grp = st.groupby("group", observed=True)["mean_score"].mean()
            r[f"{t}_tau_disease"] = round(tau(grp.values), 3)
            r[f"{t}_top_disease"] = grp.idxmax()
            r[f"{t}_top_disease_score"] = round(float(grp.max()), 2)
            act = grp[grp > 0.2]
            r[f"{t}_diseases_active"] = "; ".join(sorted(act.index.astype(str)))
            n_active_total += len(act)
        r["n_diseases_active"] = n_active_total
        s = sub[sub["axis"] == "lineage"]
        if not s.empty:
            per_t = s.groupby("tissue", observed=True)["mean_score"].max()
            r["tau_tissue"] = round(tau(per_t.values), 3)
            r["top_tissue"] = per_t.idxmax()
            r["blood_max_score"] = round(float(per_t.get("blood", np.nan)), 2)
            r["skin_max_score"] = round(float(per_t.get("skin", np.nan)), 2)
            r["tissue_pattern"] = ("both" if (per_t > 0.2).sum() > 1 else
                                   f"{per_t.idxmax()} only" if (per_t > 0.2).any()
                                   else "neither")
        sb = sub[(sub["axis"] == "subtype")]
        if not sb.empty:
            top = sb.sort_values("mean_score", ascending=False).head(3)
            r["top_subtypes"] = "; ".join(
                f"{x.group} ({x.mean_score:.1f}, {100*x.prop_expressing:.0f}%)"
                for x in top.itertuples())

    df = pd.DataFrame([{k: v for k, v in r.items() if not k.startswith("_")}
                       for r in rows])
    order = ["program_id", "analysis", "found_in_tissue", "found_in_disease",
             "n_genes", "median_IDS", "edge_density", "cells_supporting",
             "stratum_cells", "perm_p_median", "perm_p_max", "frac_pairs_beat_all_shuffles",
             "technical_flag", "tissue_pattern", "top_tissue", "blood_max_score",
             "skin_max_score", "tau_tissue", "top_celltype", "top_celltype_score",
             "tau_celltype", "n_celltypes_active", "celltypes_active",
             "top_subtypes", "top_disease", "top_disease_score", "tau_disease",
             "n_diseases_active", "diseases_active", "technical_detail", "genes"]
    df = df[[c for c in order if c in df.columns]
            + [c for c in df.columns if c not in order]]
    df = df.sort_values(["technical_flag", "median_IDS"], ascending=[True, False])
    df.to_csv(os.path.join(args.out, "program_catalogue.csv"), index=False)
    print(f"[out] {os.path.join(args.out, 'program_catalogue.csv')}: "
          f"{len(df)} rows x {len(df.columns)} columns")
    return df


if __name__ == "__main__":
    main()
