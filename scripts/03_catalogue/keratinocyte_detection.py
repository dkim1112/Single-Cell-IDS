"""
keratinocyte_detection.py — how many keratinocytes detect the antimicrobial and barrier
genes, by skin state, pooled over cells and per donor. Descriptive companion to the
program-score contrasts: it shows the antimicrobial genes switching on in a fraction of
non-lesional cells. Per-gene donor tests are unadjusted and exploratory.
Writes results/03_catalogue/keratinocyte_detection.csv
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from sc_ids import paths, units                                               # noqa: E402
from sc_ids.lineage import LINEAGE                                           # noqa: E402
LIN = LINEAGE
genes = ["S100A7", "S100A8", "S100A9", "DSC1", "DSG1"]
o = units.cells("skin")                                    # every skin cell, read sample by sample
X = units.read_cells("skin", o, genes).toarray()
k = (o.subtype.astype(str).map(LIN) == "Keratinocyte").values
D = pd.DataFrame(X[k] > 0, columns=genes)
D["any_S100A7/8/9"] = D[["S100A7", "S100A8", "S100A9"]].any(axis=1)
D["disease"] = o.disease_fine.values[k]; D["donor"] = o.donor.values[k]
cols = ["S100A8", "any_S100A7/8/9", "DSC1", "DSG1"]
pooled = D.groupby("disease", observed=True)[cols].mean() * 100
n = D.groupby("disease", observed=True).size().rename("keratinocytes")
don = D.groupby(["disease", "donor"], observed=True)[cols].mean().reset_index()
don = don.merge(D.groupby(["disease", "donor"], observed=True).size().rename("n").reset_index())
don = don[don.n >= 20]
med = don.groupby("disease", observed=True)[cols].median() * 100
rows = []
for g in cols:
    for b in ("PN", "PP"):
        x, y = don[don.disease == "NS"][g], don[don.disease == b][g]
        rows.append(dict(gene=g, comparison=f"NS vs {b}", donors_a=len(x), donors_b=len(y),
                         pct_cells_NS=pooled.loc["NS", g], pct_cells_b=pooled.loc[b, g],
                         median_donor_pct_NS=med.loc["NS", g], median_donor_pct_b=med.loc[b, g],
                         p_unadjusted=stats.mannwhitneyu(x, y).pvalue))
out = pd.DataFrame(rows).round(4)
out.to_csv(os.environ.get("KC_OUT") or os.path.join(paths.CATALOGUE, "keratinocyte_detection.csv"), index=False)
print(n.to_dict()); print(out.to_string(index=False))
