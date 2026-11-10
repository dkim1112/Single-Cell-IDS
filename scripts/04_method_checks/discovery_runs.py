"""
discovery_runs.py — rerun IDS program discovery with ONE setting changed.

Same code as the catalogue run (src/sc_ids/discovery.py, driven in stage 02 by
ids_programs_disease.py) and the same saved gene panels, so only the stated change differs:

  --seeds  which random 4,000 cells are drawn         (catalogue: seed 0)
  --norm   panel = library size over the panel genes that pass the floor (catalogue)
           total = library size over ALL genes       (what the permutation test and the scoring use)
  --cap    how many cells are sampled per stratum    (catalogue: 4,000)

Matrices are cached under the work folder, so thresholds can be tried without recomputing
(discovery_compare.py moves the edge threshold by +/- 0.05 that way).

  python scripts/04_method_checks/discovery_runs.py --analysis whole-tissue --seeds 0,1,2,3,4
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from sc_ids import paths, units, discovery                                    # noqa: E402

RUNS = {"whole-tissue": paths.DISEASE_RUN, **paths.LINEAGE_RUNS}
MIN_IDS = {a: discovery.MIN_IDS[None if a == "whole-tissue" else a] for a in RUNS}
CACHE = os.path.join(units.WORK, "discovery")
programs_from = discovery.programs_from


def lineage_of(analysis):
    return None if analysis == "whole-tissue" else analysis


def stratum_cells(analysis, stratum):
    return discovery.stratum_cells(lineage_of(analysis), stratum)


def strata_of(analysis):
    return list(pd.read_csv(os.path.join(RUNS[analysis], "ids_summary.csv")).stratum)


def ids_matrix(analysis, stratum, seed, norm, cap=4000):
    """(C, info) for one search, cached. C is None when the stratum is below the 400-cell floor."""
    tag = f"{analysis.replace('+','-')}__{stratum}__s{seed}__{norm}" + ("" if cap == 4000 else f"__cap{cap}")
    fC, fg = os.path.join(CACHE, tag + ".npy"), os.path.join(CACHE, tag + ".json")
    if os.path.exists(fC):
        return np.load(fC).astype(np.float64), json.load(open(fg))
    if os.path.exists(fg):
        return None, json.load(open(fg))
    os.makedirs(CACHE, exist_ok=True)
    panel = json.load(open(os.path.join(RUNS[analysis], "gene_panel.json")))["panel"]
    r = discovery.compute_ids(lineage_of(analysis), stratum, panel, seed=seed, norm=norm, cap=cap)
    if "skipped" in r:
        json.dump(r, open(fg, "w")); return None, r
    np.save(fC, r["C"].astype(np.float32))
    info = {"genes": r["genes"], "cells": r["cells"]}
    json.dump(info, open(fg, "w"))
    return r["C"], info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analysis", required=True, choices=list(RUNS))
    ap.add_argument("--seeds", default="0"); ap.add_argument("--norm", default="panel", choices=["panel", "total"])
    ap.add_argument("--cap", type=int, default=4000)
    ap.add_argument("--min-ids", type=float, default=None)
    ap.add_argument("--budget", type=float, default=140)
    a = ap.parse_args()
    t0 = time.time()
    for seed in [int(s) for s in a.seeds.split(",")]:
        for s in strata_of(a.analysis):
            if time.time() - t0 > a.budget:
                print("[budget] stopping; re-run to continue"); return
            t = time.time()
            C, info = ids_matrix(a.analysis, s, seed, a.norm, cap=a.cap)
            if C is None:
                print(f"{a.analysis} {s} seed {seed}: skipped ({info})"); continue
            msg = f"{a.analysis} {s} seed {seed} {a.norm}: {info['cells']} cells x {len(info['genes'])} genes"
            if a.min_ids is not None:
                msg += f" -> {len(programs_from(C, info['genes'], a.min_ids))} programs at {a.min_ids}"
            print(msg + f"  ({time.time()-t:.0f}s)", flush=True)
    print("done")


if __name__ == "__main__":
    main()
