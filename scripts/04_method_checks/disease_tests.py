"""
disease_tests.py — donor-level disease tests, on any set of lanes and at either cell-type level.

Same statistics as scripts/03_catalogue/disease_contrasts.py (one value per donor, cell-weighted
over the lanes kept; Mann-Whitney, paired Wilcoxon for PN vs PP; Benjamini-Hochberg within each
comparison; a donor x group x lane needs >= 20 cells), generalised so that two choices can be varied:
  which lanes go into a comparison   (more lanes = more donors)
  the cell-type level                (13 broad lineages, or the source subtypes)

Reads results/04_method_checks/scores_*.csv.gz (from score_units.py).
Writes tests_summary.csv (one row per comparison) and tests_detail.csv.gz (every test).
"""
import os, sys
from math import comb
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from sc_ids import paths                                                      # noqa: E402
from sc_ids.lineage import LINEAGE                                           # noqa: E402

OUT = os.path.join(paths.RESULTS, "04_method_checks")
MIN_CELLS, MIN_EFFECT = 20, 0.2


def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); adj = np.empty_like(p)
    adj[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)


def tested_programs():
    cat = pd.read_csv(os.path.join(paths.CATALOGUE, "program_catalogue.csv"))
    keep = cat[cat.technical_flag.isna() & (cat["frac_pairs_beat_all_shuffles"] >= 0.9)]
    return list(keep.program_id)                     # the same 93 detections tested in October


def donor_table(S, ids, level, units_keep=None):
    """One score per donor x group: >= 20 cells per donor x group x unit, then cell-weighted over units."""
    S = S if units_keep is None else S[S.unit.isin(units_keep)]
    grp = S.subtype.map(LINEAGE).fillna("Other") if level == "lineage" else S.subtype
    W = S[ids].multiply(S.n_cells, axis=0)
    W["donor"], W["group"], W["disease"], W["unit"], W["n"] = S.donor.values, grp.values, S.disease.values, S.unit.values, S.n_cells.values
    a = W.groupby(["donor", "group", "disease", "unit"], observed=True).sum()
    a = a[a.n >= MIN_CELLS]
    b = a.groupby(["donor", "group", "disease"], observed=True).sum()
    D = b[ids].div(b.n, axis=0); D["n_cells"] = b.n
    return D.reset_index()


def contrast(D, ids, a_, b_, paired=False):
    rows = []
    for grp, s in D.groupby("group", observed=True):
        xa, xb = s[s.disease == a_].set_index("donor"), s[s.disease == b_].set_index("donor")
        if paired:
            both = xa.index.intersection(xb.index)
            if len(both) < 4:
                continue
            xa, xb = xa.loc[both], xb.loc[both]
        elif len(xa) < 3 or len(xb) < 3:
            continue
        A, B = xa[ids].to_numpy(float), xb[ids].to_numpy(float)
        for j, p in enumerate(ids):
            x, y = A[:, j], B[:, j]
            if paired:
                pv = 1.0 if np.allclose(x, y) else stats.wilcoxon(x, y).pvalue
            else:
                pv = stats.mannwhitneyu(x, y, alternative="two-sided").pvalue
            sd = np.sqrt(((len(x) - 1) * x.var(ddof=1) + (len(y) - 1) * y.var(ddof=1)) / (len(x) + len(y) - 2))
            rows.append((p, grp, len(x), len(y), x.mean(), y.mean(), y.mean() - x.mean(), sd, pv))
    R = pd.DataFrame(rows, columns=["program_id", "group", "n_a", "n_b", "mean_a", "mean_b", "delta", "pooled_sd", "p"])
    if len(R):
        R["p_adj"] = bh(R.p.values)
    return R


def summarise(R, paired):
    if not len(R):
        return dict(tests=0)
    na, nb, m = int(R.n_a.max()), int(R.n_b.max()), len(R)
    z = stats.norm.ppf(1 - 0.05 / (2 * m)) + stats.norm.ppf(0.8)          # Bonferroni-level, 80% power
    mdd = (z * R.pooled_sd * np.sqrt(1 / R.n_a + 1 / R.n_b)).median()
    min_p = 2 / 2 ** na if paired else 2 / comb(na + nb, na)
    sig = R[R.p_adj < 0.05]
    return dict(donors_a=na, donors_b=nb, groups=R.group.nunique(), tests=m,
                raw_p_below_0_05=int((R.p < 0.05).sum()), expected_by_chance=int(round(m * 0.05)),
                significant=len(sig), significant_and_change_0_2=int((sig.delta.abs() >= MIN_EFFECT).sum()),
                programs_significant=sig.program_id.nunique(),
                smallest_possible_p=min_p, largest_abs_change=round(float(R.delta.abs().max()), 3),
                median_detectable_change=round(float(mdd), 3))


def lanes(*runs):
    return [f"PBMC-0{r}-{l}" for r in runs for l in range(1, 5)]


def main():
    ids = tested_programs()
    nine = pd.read_csv(f"{OUT}/scores_blood_nine.csv.gz", dtype={"donor": str})
    all28 = pd.read_csv(f"{OUT}/scores_blood_all28.csv.gz", dtype={"donor": str})
    skin = pd.read_csv(f"{OUT}/scores_skin_all.csv.gz", dtype={"donor": str})
    L0305 = ["PBMC-03-1", "PBMC-03-3", "PBMC-05-1", "PBMC-05-2"]
    # (tissue, lanes label, score table, lanes kept, a, b, paired, note)
    plan = [
        ("blood", "October: lanes 03-1, 03-3, 05-1, 05-2", nine, L0305, "Healthy", "PSO", False, "valid: both groups in the same lanes"),
        ("blood", "October: PBMC-06 lanes", nine, lanes(6), "PSO", "PSA", False, "valid: same lanes"),
        ("blood", "October: all nine lanes", nine, None, "PSO", "PSA", False, "valid"),
        ("blood", "October: all nine lanes", nine, None, "Healthy", "PSA", False, "NOT interpretable: groups in different lanes"),
        ("blood", "more lanes: all 16 lanes of runs 02-05", all28, lanes(2, 3, 4, 5), "Healthy", "PSO", False, "valid: both groups in every lane"),
        ("blood", "more lanes: runs 06 + 07 (8 lanes)", all28, lanes(6, 7), "PSO", "PSA", False, "valid: both groups in every lane (run 07 has 1 PSO donor)"),
        ("blood", "more lanes: all 28 lanes", all28, None, "PSO", "PSA", False, "valid, but most PsA and PSO cells sit in different runs"),
        ("blood", "more lanes: run 01 only (4 lanes)", all28, lanes(1), "Healthy", "PSA", False, "same run, DIFFERENT lanes: lane and disease cannot be separated"),
        ("blood", "more lanes: all 28 lanes", all28, None, "Healthy", "PSA", False, "NOT interpretable: groups in different lanes"),
        ("skin", "all skin", skin, None, "NS", "PN", False, "valid"),
        ("skin", "all skin", skin, None, "NS", "PP", False, "valid"),
        ("skin", "all skin, same donors", skin, None, "PN", "PP", True, "paired"),
    ]
    summ, detail = [], []
    for level in ("lineage", "subtype"):
        for tissue, label, S, keep, a_, b_, paired, note in plan:
            use = [p for p in ids if p in S.columns]
            R = contrast(donor_table(S, use, level, keep), use, a_, b_, paired)
            row = dict(tissue=tissue, cell_type_level=level, comparison=f"{a_} vs {b_}", lanes=label, **summarise(R, paired), note=note)
            summ.append(row)
            R.insert(0, "lanes", label); R.insert(0, "comparison", f"{a_} vs {b_}"); R.insert(0, "cell_type_level", level); R.insert(0, "tissue", tissue)
            detail.append(R)
            print(f"{level:8s} {tissue:5s} {a_:>7s} vs {b_:<4s} {label:42s} donors {row.get('donors_a')}/{row.get('donors_b')} "
                  f"tests {row['tests']:5d} sig {row.get('significant')}", flush=True)
    pd.DataFrame(summ).to_csv(f"{OUT}/tests_summary.csv", index=False)
    pd.concat(detail, ignore_index=True).round(6).to_csv(f"{OUT}/tests_detail.csv.gz", index=False, compression="gzip")


if __name__ == "__main__":
    main()
