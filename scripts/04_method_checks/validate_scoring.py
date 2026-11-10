"""
validate_scoring.py — is lane-by-lane scoring the same as scanpy on the merged object?

Merges a few units in memory, scores programs with sc.tl.score_genes (the old route), then
scores the same cells unit by unit through sc_ids.units (the new route) and reports the
largest difference. Expect rounding error only.

  python scripts/04_method_checks/validate_scoring.py --tissue blood --units PBMC-03-1,PBMC-03-3 --every 4
"""
import argparse, os, sys, time
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from sc_ids import paths, units                                              # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--tissue", default="blood"); ap.add_argument("--units", required=True)
ap.add_argument("--every", type=int, default=4)
a = ap.parse_args()
us = a.units.split(",")
cat = pd.read_csv(os.path.join(paths.CATALOGUE, "program_catalogue.csv"))
progs = {r.program_id: [g.strip() for g in r.genes.split(",")] for r in cat.iloc[::a.every].itertuples()}

t = time.time()
m = ad.concat([units.read_unit(a.tissue, u) for u in us], join="inner")
sc.pp.normalize_total(m, target_sum=1e4); sc.pp.log1p(m)
old = {}
for pid, g in progs.items():
    g = [x for x in g if x in m.var_names]
    if len(g) < 2:
        continue
    sc.tl.score_genes(m, g, score_name="_s", ctrl_size=50, random_state=0)
    old[pid] = m.obs.pop("_s").values
print(f"scanpy on merged object: {m.n_obs:,} cells, {len(old)} programs, {time.time()-t:.0f}s")
del m

t = time.time()
mean, n = units.gene_means(a.tissue, us)
ctrl = units.control_sets(units.genes(a.tissue), mean, progs)
new = []
for u in us:
    S, ids = units.score_matrix(units.lognorm(units.read_unit(a.tissue, u).X), ctrl)
    new.append(S)
new = np.vstack(new)
print(f"unit by unit: {new.shape[0]:,} cells, {len(ids)} programs, {time.time()-t:.0f}s")
d = np.array([np.abs(new[:, j] - old[p]).max() for j, p in enumerate(ids)])
print(f"largest absolute difference over all cells and programs: {d.max():.2e}  (scores span "
      f"{min(v.min() for v in old.values()):.2f} to {max(v.max() for v in old.values()):.2f})")
print("programs differing by more than 1e-4:", int((d > 1e-4).sum()), "of", len(d))
