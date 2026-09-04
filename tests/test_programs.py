"""
Tests for Agenda 1 (cell-type mapping + pre-registered program ranking) and
Agenda 3 (cell-type enrichment).

Run:  python tests/test_programs.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
import pandas as pd
import anndata as ad

from sc_ids.celltypes import (load_celltype_map, add_celltype_main,
                              coverage_report, MAIN_ORDER)
from sc_ids.datasets import get as get_dataset

PBMC_MAP = get_dataset("gse194315").celltype_map_path()
from sc_ids import programs as P


def test_mapping_file():
    m = load_celltype_map(PBMC_MAP, verbose=False)
    assert len(m) == 30, f"expected 30 granular labels, got {len(m)}"
    # the supplied file's 'K Cell' typo must be repaired, so NK subtypes do not
    # split off into a phantom sixth lineage
    assert "K Cell" not in set(m.values()), "typo 'K Cell' survived load"
    assert m["NK_CD56bright"] == "NK Cell" and m["NK Proliferating"] == "NK Cell"
    assert m["NK"] == "NK Cell"
    assert set(m.values()) <= set(MAIN_ORDER), set(m.values())
    assert len(set(m.values())) == 5, sorted(set(m.values()))
    # opting out of the repair must give back exactly what the file says
    raw = load_celltype_map(PBMC_MAP, fix_typos=False, verbose=False)
    assert raw["NK_CD56bright"] == "K Cell"
    print(f"  mapping: 30 -> {sorted(set(m.values()))}  [ok]")


def test_add_celltype_main():
    labels = ["CD4 TCM", "CD14 Mono", "B naive", "NK", "NK Proliferating",
              "Platelet", "SomethingNew"]
    A = ad.AnnData(np.zeros((len(labels), 3), dtype=np.float32))
    A.obs["CellType"] = labels
    add_celltype_main(A, path=PBMC_MAP, verbose=False)
    got = list(A.obs["CellTypeMain"].astype(str))
    assert got == ["T Cell", "Myeloid Cell", "B Cell", "NK Cell", "NK Cell",
                   "Other", "Unmapped"], got
    # unmapped stays visible as its own category rather than becoming NaN
    assert A.obs["CellTypeMain"].isna().sum() == 0
    # categories follow MAIN_ORDER so colours are stable across figures
    cats = list(A.obs["CellTypeMain"].cat.categories)
    assert cats == [c for c in MAIN_ORDER if c in cats], cats
    missing, unused = coverage_report(load_celltype_map(PBMC_MAP, verbose=False), labels)
    assert missing == ["SomethingNew"] and len(unused) == 24, (missing, len(unused))
    print("  add_celltype_main + coverage_report  [ok]")


def _toy_C(blocks, d, fill=0.05, seed=0):
    rng = np.random.default_rng(seed)
    C = rng.uniform(0, fill, size=(d, d)); C = (C + C.T) / 2
    for idx, vals in blocks:
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                C[idx[a], idx[b]] = C[idx[b], idx[a]] = vals[(a, b)]
    np.fill_diagonal(C, 1.0)
    return C


def test_rank_uses_median_not_mean():
    """A program carried by one near-duplicate pair must NOT outrank a program
    where every pair is strong. This is the concrete failure of ranking on mean."""
    d = 8
    names = [f"g{i}" for i in range(d)]
    # A: genes 0-3, all pairs 0.70            -> mean 0.700, median 0.70
    a_idx = [0, 1, 2, 3]
    a_vals = {(i, j): 0.70 for i in range(4) for j in range(i + 1, 4)}
    # B: genes 4-7, two near-duplicate pairs at 1.00, the other four at 0.60
    #    -> mean 0.733 (BEATS A), median 0.600 (LOSES to A)
    b_idx = [4, 5, 6, 7]
    b_vals = {(i, j): 0.60 for i in range(4) for j in range(i + 1, 4)}
    b_vals[(0, 1)] = 1.00
    b_vals[(2, 3)] = 1.00
    C = _toy_C([(a_idx, a_vals), (b_idx, b_vals)], d)
    mods = [{"genes": [names[i] for i in a_idx], "idx": a_idx},
            {"genes": [names[i] for i in b_idx], "idx": b_idx}]
    mods = P.annotate_programs(mods, C, min_ids=0.5)
    byname = {m["genes"][0]: m for m in mods}
    assert byname["g4"]["mean_ids"] > byname["g0"]["mean_ids"], "test setup wrong"
    assert byname["g4"]["median_ids"] < byname["g0"]["median_ids"], "test setup wrong"
    ranked, _ = P.rank_programs(mods, verbose=False)
    assert ranked[0]["genes"][0] == "g0", [m["genes"][0] for m in ranked]
    assert ranked[0]["rank"] == 1 and ranked[0]["name"].startswith("P01_")
    print(f"  median-ranked A(med .70) over B(mean {byname['g4']['mean_ids']:.3f})  [ok]")


def test_size_filter_and_reporting():
    d = 60
    names = [f"g{i}" for i in range(d)]
    small = {"genes": names[:2], "idx": [0, 1]}                 # a pair, not a program
    ok = {"genes": names[2:8], "idx": list(range(2, 8))}
    huge = {"genes": names[8:60], "idx": list(range(8, 60))}    # 52 genes: a compartment
    C = _toy_C([], d, fill=0.6)
    mods = P.annotate_programs([small, ok, huge], C, min_ids=0.5)
    kept, dropped = P.rank_programs(mods, min_genes=3, max_genes=50, verbose=False)
    assert [m["genes"][0] for m in kept] == ["g2"], [m["genes"][0] for m in kept]
    assert len(dropped["too_small"]) == 1 and len(dropped["too_large"]) == 1
    # the filter is a REPORTING rule, so nothing is silently deleted: both
    # exclusions are returned for the record
    assert dropped["too_large"][0][0] == 52
    print("  size filter 3-50 with an audit trail  [ok]")


def _toy_adata_with_types(seed=0):
    """600 cells, 3 cell types, 800 genes. Genes 0-4 are ON only in type A;
    genes 5-9 are ON in every cell (a housekeeping decoy).

    The background genes are given a SPREAD of expression levels on purpose.
    sc.tl.score_genes builds its control set from genes in the same expression
    bin, so a toy matrix with a flat background and two equally-huge gene sets
    makes each program the other's control and both scores collapse to ~0. That
    is a property of score_genes, not of this code, but it is worth knowing
    before trusting a score computed on a heavily gene-filtered object.
    """
    rng = np.random.default_rng(seed)
    n, g = 600, 800
    X = rng.poisson(lam=rng.uniform(0.5, 25.0, g), size=(n, g)).astype(np.float32)
    types = np.array(["A"] * 200 + ["B"] * 200 + ["C"] * 200)
    X[types == "A", :5] += rng.poisson(80, size=((types == "A").sum(), 5))
    X[:, 5:10] += rng.poisson(80, size=(n, 5))
    A = ad.AnnData(np.clip(X, 0, None))
    A.var_names = [f"g{i}" for i in range(g)]
    A.obs["CellType"] = types
    import scanpy as sc
    sc.pp.normalize_total(A, target_sum=1e6); sc.pp.log1p(A)
    return A


def test_enrichment_finds_the_right_type_and_flags_the_decoy():
    A = _toy_adata_with_types()
    progs = [{"genes": [f"g{i}" for i in range(5)], "idx": list(range(5)), "name": "P01_specific"},
             {"genes": [f"g{i}" for i in range(5, 10)], "idx": list(range(5, 10)), "name": "P02_ubiquitous"}]
    names = P.score_programs(A, progs, verbose=False)
    assert names == ["P01_specific", "P02_ubiquitous"]
    enr = P.celltype_enrichment(A, names, celltype_col="CellType")

    top = enr[enr["program"] == "P01_specific"].iloc[0]
    assert top["cell_type"] == "A" and top["rank_in_program"] == 1, top.to_dict()
    assert top["prop_expressing"] > 0.95, top["prop_expressing"]
    assert top["enrichment"] > 1.5, top["enrichment"]
    # Wilson interval must bracket the point estimate and stay in [0,1]
    assert top["ci_low"] <= top["prop_expressing"] <= top["ci_high"] <= 1.0

    spec = P.specificity_summary(enr)
    s = spec.set_index("program")
    assert bool(s.loc["P01_specific", "specific"]) is True
    # the housekeeping decoy is expressed everywhere: high proportion, no specificity.
    # This is exactly the case where "sort cell types by proportion" is misleading
    # on its own, which is why enrichment is reported beside it.
    assert bool(s.loc["P02_ubiquitous", "specific"]) is False
    assert s.loc["P02_ubiquitous", "top_prop"] > 0.9
    assert abs(s.loc["P02_ubiquitous", "top_enrichment"] - 1.0) < 0.15
    print(f"  enrichment: specific->A (prop {top['prop_expressing']:.2f}, "
          f"enrich {top['enrichment']:.2f}); decoy flagged diffuse  [ok]")


def test_saturation_tiebreak_and_disagreement():
    """The two failure modes of "sort cell types by proportion", both seen on the
    real PBMC run. Scores are injected directly so the arithmetic is exact."""
    n = 900
    A = ad.AnnData(np.zeros((n, 4), dtype=np.float32))
    A.obs["CellType"] = np.array(["A"] * 300 + ["B"] * 300 + ["C"] * 300)
    s = np.empty(n)
    s[:300] = 1.0            # type A: every cell expresses, weakly   -> prop 1.00
    s[300:600] = 5.0         # type B: every cell expresses, strongly -> prop 1.00
    s[600:870] = -1.0        # type C: bimodal - 10% enormous, 90% off
    s[870:] = 100.0          #                                        -> prop 0.10
    A.obs["P01"] = s
    enr = P.celltype_enrichment(A, ["P01"], celltype_col="CellType")
    order = list(enr.sort_values("rank_in_program")["cell_type"])
    # 1. SATURATION: A and B both sit at prop 1.000, so proportion alone cannot
    #    order them; the mean-score tie-break must put B first, not insertion order.
    assert order == ["B", "A", "C"], order
    # 2. DISAGREEMENT: by strength the answer is C, which proportion ranks LAST.
    by_mean = list(enr.sort_values("rank_by_mean_score")["cell_type"])
    assert by_mean[0] == "C", by_mean
    spec = P.specificity_summary(enr).iloc[0]
    assert bool(spec["saturated"]) and int(spec["n_tied_at_top_prop"]) == 2
    assert spec["top_cell_type"] == "B" and spec["top_by_mean_score"] == "C"
    assert not bool(spec["agree"])
    print("  saturation tie-break + prop/mean disagreement both flagged  [ok]")


def test_wilson_widths_scale_with_n():
    """The whole reason the CI is in the table: 55 cells and 55,000 cells give
    the same proportion very different credibility."""
    lo_s, hi_s = P._wilson(np.array([44.]), np.array([55.]))
    lo_b, hi_b = P._wilson(np.array([44000.]), np.array([55000.]))
    assert (hi_s - lo_s) > 8 * (hi_b - lo_b)
    assert 0 <= lo_s and hi_s <= 1
    print(f"  Wilson width: n=55 -> {(hi_s-lo_s)[0]:.3f}, "
          f"n=55000 -> {(hi_b-lo_b)[0]:.4f}  [ok]")


if __name__ == "__main__":
    print("Agenda 1a — cell-type mapping")
    test_mapping_file(); test_add_celltype_main()
    print("Agenda 1b — pre-registered program ranking")
    test_rank_uses_median_not_mean(); test_size_filter_and_reporting()
    print("Agenda 3 — cell-type enrichment")
    test_enrichment_finds_the_right_type_and_flags_the_decoy()
    test_saturation_tiebreak_and_disagreement()
    test_wilson_widths_scale_with_n()
    print("\nAll tests passed.")
