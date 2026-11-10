"""
score_units.py — score every detection in every cell, one unit at a time, and keep the
donor x cell-subtype means (never the per-cell matrix).

  python scripts/04_method_checks/score_units.py --tissue blood --set nine     # the 9 lanes used in October
  python scripts/04_method_checks/score_units.py --tissue blood --set all28    # every lane
  python scripts/04_method_checks/score_units.py --tissue skin  --set all

Output: results/04_method_checks/scores_<tissue>_<set>.csv.gz
  one row per unit x donor x subtype: n_cells, then one mean-score column per detection.
Checkpointed per unit; re-run until it prints "done".
"""
import argparse, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from sc_ids import paths, units                                              # noqa: E402
from sc_ids.datasets.gse194315 import LANES_USED                             # noqa: E402

OUT = os.path.join(paths.RESULTS, "04_method_checks")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", required=True); ap.add_argument("--set", required=True)
    ap.add_argument("--budget", type=float, default=140)
    ap.add_argument("--programs", default=None, help="CSV with program_id, genes (default: the 194 October detections)")
    ap.add_argument("--tag", default="", help="suffix for the output name when --programs is given")
    a = ap.parse_args()
    t0 = time.time()
    c = units.cells(a.tissue)
    us = list(LANES_USED) if (a.tissue == "blood" and a.set == "nine") else list(c.unit.unique())
    name = f"{a.tissue}_{a.set}" + (f"_{a.tag}" if a.tag else "")
    part = os.path.join(units.WORK, "scores", name); os.makedirs(part, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)
    cat = pd.read_csv(a.programs or os.path.join(paths.CATALOGUE, "program_catalogue.csv"))
    progs = {r.program_id: [g.strip() for g in r.genes.split(",")] for r in cat.itertuples()}

    mean, n = units.gene_means(a.tissue, us, cache=os.path.join(part, "_gene_means.npz"))
    ctrl = units.control_sets(units.genes(a.tissue), mean, progs)
    print(f"[pass 1] gene means over {n:,} cells in {len(us)} units; {len(ctrl)} of {len(progs)} detections scorable "
          f"({time.time()-t0:.0f}s)", flush=True)
    for u in us:
        f = os.path.join(part, f"{u}.csv.gz")
        if os.path.exists(f):
            continue
        if time.time() - t0 > a.budget:
            print("[budget] stopping; re-run to continue"); return
        ad_ = units.read_unit(a.tissue, u)
        S, ids = units.score_matrix(units.lognorm(ad_.X), ctrl)
        df = pd.DataFrame(S, columns=ids)
        df.insert(0, "disease", ad_.obs["disease_fine"].astype(str).values)
        df.insert(0, "subtype", ad_.obs["subtype"].astype(str).values)
        df.insert(0, "donor", ad_.obs["donor"].astype(str).values)
        g = df.groupby(["donor", "subtype", "disease"], observed=True)
        o = g[ids].mean(); o.insert(0, "n_cells", g.size())
        o = o.reset_index(); o.insert(0, "unit", u)
        o.to_csv(f + ".tmp", index=False, compression="gzip"); os.replace(f + ".tmp", f)
    allp = pd.concat([pd.read_csv(os.path.join(part, f"{u}.csv.gz"), dtype={"donor": str}) for u in us], ignore_index=True)
    dst = os.path.join(OUT, f"scores_{name}.csv.gz")
    allp.to_csv(dst, index=False, compression="gzip")
    print(f"done: {dst}  {len(allp):,} rows, {int(allp.n_cells.sum()):,} cells ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
