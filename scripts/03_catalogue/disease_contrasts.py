"""
disease_contrasts.py — donor-level disease contrasts on the catalogue programmes.

The unit of replication is the DONOR, one value per donor. Blood donors were
sequenced in several lanes (44 of 53 appear in more than one), so scores are first
averaged over a donor's cells ACROSS lanes (cell-weighted). Counting each lane as a
separate donor, which an earlier version of this script did, inflated every blood
result: PsA vs psoriasis went from 0 to 52 "significant" programmes that way.

Designs, and what each can support
  blood  Healthy vs PSO   lanes 03 and 05 only. These are the only lanes that hold
                          both groups. The PBMC-06 lanes contain psoriasis and PsA
                          donors but no healthy ones; including them makes the
                          psoriasis group come from different lanes than the healthy
                          group, and produced 154 false positives.
  blood  PSO vs PsA       (a) PBMC-06 lanes only: 12 PsA and 6 PSO donors in the SAME
                          lanes. (b) all lanes. A real effect appears in both.
  blood  Healthy vs PsA   cross-lane by construction (no lane has both groups in
                          usable numbers). Reported, marked, not interpretable.
  skin   NS vs PN, NS vs PP, PN vs PP (paired within donor). Each skin donor is
         one sample, so no collapsing is needed.

Scoring (slow, needs the full matrices) is cached in donor_means.csv.gz; delete that
file to rescore.
"""
import os, sys
import numpy as np
import pandas as pd
from scipy import stats
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from sc_ids import paths                                                      # noqa: E402
from sc_ids.lineage import LINEAGE                                           # noqa: E402

OUT = paths.CATALOGUE
CACHE = f"{OUT}/donor_means.csv.gz"
MIN_CELLS = 20
MIN_EFFECT = 0.2   # score units; same cutoff the catalogue uses for "active"
BOTH_GROUP_LANES = ("PBMC-03", "PBMC-05")
LANE06 = ("PBMC-06",)


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.empty_like(p)
    adj[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)


def score_donors(cat):
    """Score every biological, significant programme per cell, then average per
    donor x lineage x disease x lane. Cached because it needs both full matrices."""
    import scanpy as sc, anndata as ad
    keep = cat[cat.technical_flag.isna() & (cat["frac_pairs_beat_all_shuffles"] >= 0.9)]
    progs = {r.program_id: [g.strip() for g in r.genes.split(",")] for r in keep.itertuples()}
    rows = []
    for tag, f in (("blood", os.path.join(paths.DISEASE_RUN, "blood_qc.h5ad")),
                   ("skin", os.path.join(paths.DISEASE_RUN, "skin_qc.h5ad"))):
        a = ad.read_h5ad(f)
        sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
        names = []
        for p, g in progs.items():
            g = [x for x in g if x in a.var_names]
            if len(g) < 2:
                continue
            sc.tl.score_genes(a, g, score_name=p, ctrl_size=50, random_state=0)
            names.append(p)
        o = a.obs
        key = pd.DataFrame({"donor": o["donor"].astype(str).values,
                            "lineage": o["subtype"].astype(str).map(LINEAGE).fillna("Other").values,
                            "disease": o["disease_fine"].astype(str).values,
                            "lane": o["lane"].astype(str).values if "lane" in o else tag})
        for p in names:
            key["s"] = a.obs[p].values
            g = (key.groupby(["donor", "lineage", "disease", "lane"], observed=True)["s"]
                    .agg(["mean", "size"]).reset_index())
            g = g[g["size"] >= MIN_CELLS]
            g["program_id"], g["tissue"] = p, tag
            rows.append(g.rename(columns={"mean": "score", "size": "n_cells"}))
        del a
    D = pd.concat(rows, ignore_index=True)
    D.to_csv(CACHE, index=False, compression="gzip")
    return D


def collapse_lanes(d):
    """One value per donor: cell-weighted mean over the lanes kept."""
    d = d.assign(w=d["score"] * d["n_cells"])
    g = d.groupby(["tissue", "program_id", "lineage", "disease", "donor"], observed=True)
    out = g[["w", "n_cells"]].sum().reset_index()
    out["score"] = out["w"] / out["n_cells"]
    return out.drop(columns="w")


def contrast(d, a_, b_, paired):
    out = []
    for (p, lin), s in d.groupby(["program_id", "lineage"], observed=True):
        x = s[s.disease == a_].set_index("donor")["score"]
        y = s[s.disease == b_].set_index("donor")["score"]
        if paired:
            both = x.index.intersection(y.index)
            if len(both) < 4:
                continue
            pv = stats.wilcoxon(x[both], y[both]).pvalue
            na = nb = len(both)
            dl = float((y[both] - x[both]).mean())
            ma, mb = float(x[both].mean()), float(y[both].mean())
        else:
            if len(x) < 3 or len(y) < 3:
                continue
            pv = stats.mannwhitneyu(x, y, alternative="two-sided").pvalue
            na, nb = len(x), len(y)
            dl = float(y.mean() - x.mean())
            ma, mb = float(x.mean()), float(y.mean())
        out.append(dict(program_id=p, lineage=lin, n_a=na, n_b=nb,
                        mean_a=round(ma, 3), mean_b=round(mb, 3),
                        delta=round(dl, 3), p=pv))
    return pd.DataFrame(out)


def main():
    cat = pd.read_csv(f"{OUT}/program_catalogue.csv")
    D = pd.read_csv(CACHE) if os.path.exists(CACHE) else score_donors(cat)
    D["lane"] = D["lane"].astype(str)

    blood = D[D.tissue == "blood"]
    skin = collapse_lanes(D[D.tissue == "skin"])
    in_lanes = lambda pre: blood[blood.lane.str.startswith(pre)]

    tests = [
        ("blood", "Healthy vs PSO", "lanes 03+05 (both groups present)",
         collapse_lanes(in_lanes(BOTH_GROUP_LANES)), "Healthy", "PSO", False),
        ("blood", "PSO vs PSA", "PBMC-06 lanes (same lanes)",
         collapse_lanes(in_lanes(LANE06)), "PSO", "PSA", False),
        ("blood", "PSO vs PSA", "all lanes",
         collapse_lanes(blood), "PSO", "PSA", False),
        ("blood", "Healthy vs PSA", "CROSS-LANE, not interpretable",
         collapse_lanes(blood), "Healthy", "PSA", False),
        ("skin", "NS vs PN", "all skin", skin, "NS", "PN", False),
        ("skin", "NS vs PP", "all skin", skin, "NS", "PP", False),
        ("skin", "PN vs PP", "paired within donor", skin, "PN", "PP", True),
    ]
    frames = []
    for tissue, name, design, d, a_, b_, paired in tests:
        r = contrast(d, a_, b_, paired)
        r["p_adj"] = bh(r["p"].values) if len(r) else []
        r.insert(0, "tissue", tissue); r.insert(1, "contrast", name)
        r.insert(2, "design", design); r.insert(3, "paired", paired)
        frames.append(r)
    R = pd.concat(frames, ignore_index=True)
    R = R.merge(cat[["program_id", "genes", "n_genes", "analysis"]], on="program_id", how="left")
    R = R.sort_values(["tissue", "contrast", "design", "p"])
    R.to_csv(f"{OUT}/disease_contrasts.csv", index=False)

    sig = R[R.p_adj < 0.05]
    summ = (R.groupby(["tissue", "contrast", "design"])
              .agg(donors_a=("n_a", "max"), donors_b=("n_b", "max"), tests=("p", "size"),
                   raw_p_below_0_05=("p", lambda x: int((x < 0.05).sum())))
              .join(sig.groupby(["tissue", "contrast", "design"]).size().rename("significant"))
              .fillna({"significant": 0}).astype({"significant": int}))
    summ["expected_by_chance"] = (summ["tests"] * 0.05).round(0).astype(int)
    summ.to_csv(f"{OUT}/disease_contrasts_summary.csv")
    print(summ.to_string())

    # one line per programme for the workbook
    def ev(sub):
        s = sub[(sub.p_adj < 0.05) & (sub.delta.abs() >= MIN_EFFECT)]
        return {p: "; ".join(f"{r.lineage}: {r.contrast.split(' vs ')[0]}->{r.contrast.split(' vs ')[1]} "
                             f"({r.delta:+.2f})" for r in g.itertuples())
                for p, g in s.groupby("program_id")}
    sk = ev(R[R.tissue == "skin"])
    bl = ev(R[(R.tissue == "blood") & ~R.design.str.startswith("CROSS")])
    pd.DataFrame({"program_id": cat.program_id,
                  "skin_disease_evidence": cat.program_id.map(sk).fillna(""),
                  "blood_disease_evidence": cat.program_id.map(bl).fillna("")}
                 ).to_csv(f"{OUT}/disease_evidence.csv", index=False)

    # Does what changes in lesional skin also change in blood? Same program in the
    # same cell type, skin normal-vs-lesional against blood healthy-vs-psoriasis
    # (lanes 03+05 only, the valid blood design).
    sk_ = R[(R.tissue == "skin") & (R.contrast == "NS vs PP")]
    bl_ = R[(R.tissue == "blood") & (R.contrast == "Healthy vs PSO")]
    m = sk_.merge(bl_, on=["program_id", "lineage"], suffixes=("_skin", "_blood"))
    m[["program_id", "lineage", "delta_skin", "p_adj_skin", "delta_blood", "p_adj_blood"]
      ].to_csv(f"{OUT}/skin_to_blood_check.csv", index=False)
    s_sig = m[m.p_adj_skin < 0.05]
    print(f"skin-to-blood: {len(m)} matched program x cell-type pairs; "
          f"{len(s_sig)} significant in skin; "
          f"{int((s_sig.p_adj_blood < 0.05).sum())} of those significant in blood")


if __name__ == "__main__":
    main()
