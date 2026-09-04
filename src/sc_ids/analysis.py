"""
analysis.py — the end-to-end run that turns one IDS stratum into reportable
gene programs plus their cell-type assignment (Agendas 1 and 3).

    python -m sc_ids.analysis \
        --data data/GSE194315_PBMC-01-07_processed_data_files --prefix "PBMC-01-1." \
        --metadata data/GSE194315_CellMetadata-AS_TotalCiteseq_20220711.tsv \
        --discovery-col Subject --discovery-value PSA26 \
        --min-cells-frac 0.02 --min-ids 0.5 --out results/PSA26

DISCOVERY SET vs EVALUATION SET — the important structural choice here.
IDS runs on the DISCOVERY set (one subject inside one file: batch and donor both
held constant, all cell types kept so identity programs stay visible). The
resulting programs are then SCORED on the EVALUATION set — every included cell in
the lane, all nine subjects. Two reasons this is not a detail:

  * Cell-type labels never touch program discovery. The programs are found
    bottom-up from gene-gene dependence alone, so the enrichment table in step 6
    is a genuine out-of-sample check, not a circular restatement of the input.
  * A single subject has too few cells of the rare types (ASDC, HSPC, ILC) to say
    anything about them. Scoring on the whole lane restores those counts, which is
    exactly what Agenda 3 needs to resolve granular subtypes.

Everything is written to --out: two CSVs (the ranked program table, the full
enrichment table), a specificity summary, and the figures.
"""
from __future__ import annotations

import argparse
import gc
import json
import resource
import os

import numpy as np
import scanpy as sc

from .pipeline import load_data, join_metadata, qc_filter, preprocess, extract_modules
from .batched import compute_IDS_blocked
from .celltypes import (add_celltype_main, load_celltype_map, coverage_report,
                         map_from_obs)
from .datasets import get as get_dataset, available as available_datasets
from . import programs as P


def _rss(tag):
    """Peak RSS so far. A 67k-cell cohort sits close enough to the memory ceiling
    that a silent OOM kill is a realistic failure mode; printing the peak makes it
    obvious which step is responsible instead of leaving a truncated log."""
    print(f"  [mem] {tag}: peak RSS "
          f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:.2f} GB",
          flush=True)


def build(args):
    os.makedirs(args.out, exist_ok=True)
    ds = get_dataset(args.dataset)
    print(f"Dataset: {ds.name} — {ds.description}")
    log = {"dataset": ds.name}

    # ---------------------------------------------------------------- 1. load
    full = load_data(args.data, prefix=args.prefix)
    # 10x features files carry duplicate gene symbols; without this, var_names
    # lookups in score_genes are ambiguous and scanpy raises or silently picks one.
    full.var_names_make_unique()
    print(f"Loaded {full.n_obs} cells x {full.n_vars} genes")
    if args.metadata:
        sample = args.metadata_sample or (args.prefix.rstrip(".") if args.prefix else None)
        full = join_metadata(full, args.metadata, barcode_col=ds.barcode_col,
                             sample=sample, included_only=ds.qc_gate is not None,
                             included_col=ds.qc_gate[0] if ds.qc_gate else None)
    else:
        # A pre-built .h5ad already carries the metadata and the study's filtering
        # (see datasets/psoriasis_skin.build_h5ad). Joining again would be a no-op
        # at best, so require the columns instead of silently continuing without.
        need = [c for c in (ds.celltype_col, ds.celltype_main_col, *ds.stratum_cols)
                if c and c not in full.obs.columns]
        if need:
            raise SystemExit(f"--metadata not given and the object is missing "
                             f"{need}; pass --metadata or build the .h5ad first")
        print(f"  [metadata] already joined on the input object "
              f"({full.n_obs} cells)")
    log["evaluation_cells"] = int(full.n_obs)

    # ------------------------------------------------- 2. Agenda 1a: main types
    if args.celltype_map:
        mapping = load_celltype_map(args.celltype_map)
    elif ds.celltype_main_col and ds.celltype_main_col in full.obs.columns:
        mapping = map_from_obs(full, ds.celltype_col, ds.celltype_main_col)
    else:
        mapping = load_celltype_map(ds.celltype_map_path())
    missing, unused = coverage_report(mapping, full.obs[ds.celltype_col].unique())
    print(f"  [celltypes] labels in data but not in map: {missing or 'none'}")
    print(f"  [celltypes] labels in map but not in data: {unused or 'none'}")
    add_celltype_main(full, mapping=mapping, source_col=ds.celltype_col)
    log["celltype_map_missing"], log["celltype_map_unused"] = missing, unused
    full.obs[[ds.celltype_col, "CellTypeMain"]].value_counts().to_csv(
        os.path.join(args.out, "celltype_counts.csv"))

    # ---------------------------------------------- 3. discovery set -> IDS
    disc = full
    if args.discovery_col:
        disc = full[full.obs[args.discovery_col].astype(str) == str(args.discovery_value)].copy()
        print(f"  [discovery] {args.discovery_col}=={args.discovery_value}: "
              f"{full.n_obs} -> {disc.n_obs} cells")
    disc = qc_filter(disc, min_cells_frac=args.min_cells_frac)
    log["discovery_cells"], log["discovery_genes"] = int(disc.n_obs), int(disc.n_vars)
    if disc.n_obs < 500:
        print(f"  [discovery] WARNING: {disc.n_obs} cells — below the ~500-1000 floor "
              f"at which IDS estimates were found to be stable; treat as provisional.")

    cache_C = os.path.join(args.out, "ids_matrix.npy")
    cache_g = os.path.join(args.out, "ids_genes.json")
    if args.reuse_ids and os.path.exists(cache_C) and os.path.exists(cache_g):
        # min_ids is chosen by looking at the IDS histogram, which means a second
        # pass. Recomputing IDS to change a downstream threshold would be pure
        # waste, so the matrix is cached and reused verbatim.
        C = np.load(cache_C).astype(np.float64)
        genes = json.load(open(cache_g))
        print(f"  [ids] reusing cached matrix {C.shape} from {cache_C}")
    else:
        X, genes = preprocess(disc, n_hvg=None, keep_sparse=True)
        print(f"IDS input: {X.shape[0]} cells x {X.shape[1]} genes")
        C = compute_IDS_blocked(X, block_size=args.block_size, p_norm="max", verbose=True)
        np.save(cache_C, C.astype(np.float32))
        with open(cache_g, "w") as fh:
            json.dump(genes, fh)

    off = C[~np.eye(len(C), dtype=bool)]
    qs = {f"p{q}": float(np.percentile(off, q))
          for q in (50, 90, 99, 99.9, 99.99)}
    log.update({f"ids_{k}": v for k, v in qs.items()})
    log["ids_max"] = float(off.max())
    print("IDS off-diagonal: " + "  ".join(f"{k}={v:.3f}" for k, v in qs.items())
          + f"  max={log['ids_max']:.3f}")

    if args.min_ids is None:
        # DELIBERATELY NOT GUESSED. min_ids is where the background bulk ends and
        # real programs begin, and that gap sits in a different place in every
        # dataset. The matrix is now cached, so choosing it costs one cheap rerun.
        raise SystemExit(
            "\n--min-ids was not given and this dataset registers no default.\n"
            "Read it off the percentiles above: pick a value above the background\n"
            "bulk and below the program mode, then rerun with --reuse-ids to skip\n"
            "recomputing IDS:\n"
            f"    python -m sc_ids.analysis ... --out {args.out} "
            f"--reuse-ids --min-ids <value>\n"
            f"(for reference, GSE194315's background reached p99={0.195:.3f} and "
            f"0.5 was chosen; this run reached p99={qs['p99']:.3f})")

    # ------------------------------------- 4. Agenda 1b: rank, do not cherry-pick
    mods = extract_modules(C, genes, min_ids=args.min_ids, min_size=2, method="greedy")
    mods = P.annotate_programs(mods, C, min_ids=args.min_ids)
    progs, dropped = P.rank_programs(mods, min_genes=args.min_genes,
                                     max_genes=args.max_genes, top_n=args.top_n)
    tbl = P.programs_table(progs)
    tbl.to_csv(os.path.join(args.out, "programs_ranked.csv"), index=False)
    print(f"\nRanked programs (size {args.min_genes}-{args.max_genes}, all reported):")
    for _, r in tbl.iterrows():
        g = r["genes"] if len(r["genes"]) < 90 else r["genes"][:87] + "..."
        print(f"  {r['rank']:>2}. {r['program']:<22} n={r['n_genes']:<3} "
              f"medIDS={r['median_ids']:.3f} dens={r['density']:.2f}  {g}")
    log["n_programs_reported"] = len(progs)
    log["n_dropped_too_small"] = len(dropped["too_small"])
    log["n_dropped_too_large"] = len(dropped["too_large"])
    # AUDIT TRAIL. The size filter is a reporting rule, so what it excluded is part
    # of the result and is written out every run. On real data the boundary does
    # real work: a 57-gene platelet community sits just above the 50-gene cap. The
    # correct response to seeing that is NOT to move the cap afterwards — that is
    # the cherry-picking this rule exists to prevent — but to change it in advance
    # for all future runs, or to raise --min-ids so greedy modularity stops
    # agglomerating, and to say which was done.
    import pandas as _pd
    _pd.DataFrame(
        [{"reason": r, "n_genes": n, "first_genes": ", ".join(gs)}
         for r, items in dropped.items() for n, gs in items]
    ).to_csv(os.path.join(args.out, "programs_dropped_by_size_filter.csv"), index=False)
    if dropped["too_large"] or dropped["too_small"]:
        print("  [programs] excluded by the size filter (recorded, not hidden):")
        for r, items in dropped.items():
            for n, gs in items:
                print(f"      {r:<10} n={n:<4} {', '.join(gs)}...")
    if not progs:
        raise SystemExit("No programs survived the size filter — lower --min-ids "
                         "or widen --min-genes/--max-genes, and say which you changed.")

    # ------------------------------- 5. evaluation representation (log space)
    # The IDS matrix and the discovery subset are finished with. C alone is
    # d x d float64 (~1 GB at 11,400 genes) and is held for nothing while the
    # evaluation object is built beside it.
    del C
    for _n in ("disc", "X"):
        if _n in locals():
            del locals()[_n]
    gc.collect()
    _rss("before evaluation copy")

    # No counts layer here: nothing downstream reads it, and on a 67k-cell cohort
    # it is another full copy of the matrix for no benefit.
    ev = full.copy()
    sc.pp.normalize_total(ev, target_sum=1e6)
    sc.pp.log1p(ev)
    _rss("evaluation object normalized")
    names = P.score_programs(ev, progs)
    _rss("programs scored")

    # ------------------------------------------ 6. Agenda 3: enrichment table
    # Granular labels only. The whole reason this exists is that granular subtypes
    # are the ones a UMAP cannot separate; a 5-main-type version answers the
    # question the UMAP already answers, so it is not produced.
    enr = P.celltype_enrichment(ev, names, celltype_col=ds.celltype_col)
    enr.to_csv(os.path.join(args.out, "enrichment_granular.csv"), index=False)
    # Heatmap columns grouped by lineage, not by average value: a value-sorted axis
    # scatters the B subtypes across the figure and undoes Agenda 1's whole point.
    order = P.lineage_column_order(enr["cell_type"].unique(), mapping)
    P.plot_enrichment_heatmap(
        enr, os.path.join(args.out, "enrichment_granular.png"), column_order=order,
        title="Proportion of each CellType expressing each program (score > 0)")

    spec = P.specificity_summary(enr)
    spec.to_csv(os.path.join(args.out, "specificity_summary.csv"), index=False)
    print("\nAgenda 3 — top cell type per program (granular labels, "
          "types with >=100 cells):")
    for _, r in spec.iterrows():
        flags = []
        if r["saturated"]:
            flags.append(f"SATURATED({int(r['n_tied_at_top_prop'])} types at "
                         f"prop {r['top_prop']:.2f})")
        if not r["agree"]:
            flags.append(f"DISAGREES: by mean score -> {r['top_by_mean_score']}")
        flags.append("specific" if r["specific"] else "diffuse")
        print(f"  {r['program']:<22} {r['top_cell_type']:<18} "
              f"prop={r['top_prop']:.2f} (n={r['top_n_cells']}) "
              f"mean={r['top_mean_score']:.2f} "
              f"enrich={r['top_enrichment']:.2f}  {'; '.join(flags)}")

    # ------------------------------------------------------- 7. UMAP figures
    if not args.no_umap:
        from .viz import compute_embedding, plot_umap
        # Keep only the scores and drop the evaluation object: compute_embedding
        # copies the input AND stores a counts layer, so holding a second full
        # matrix here is what pushes a large cohort into swap.
        scores = ev.obs[names].copy()
        del ev
        gc.collect()
        _rss("before embedding")
        # `full` is not read after this point, so the embedding is built in place
        # rather than on a copy.
        emb = compute_embedding(full, batch_key=None, n_pcs=args.n_pcs,
                                with_umap=True, store_raw=False, copy=False)
        full = None
        _rss("embedding done")
        # Carry the per-cell program scores into the saved object so embedded.h5ad
        # is self-contained and can be re-plotted without re-running IDS. Nothing
        # in this script plots them — Agenda 3's enrichment table is where programs
        # are read against cell types.
        for n in names:
            emb.obs[n] = scores[n].values
        plot_umap(emb, "CellTypeMain", os.path.join(args.out, "umap_celltype_main.png"),
                  title="Main cell type (Agenda 1)")
        plot_umap(emb, ds.celltype_col, os.path.join(args.out, "umap_celltype_granular.png"),
                  title="Granular cell type (30 labels — the readability problem)")
        # The metadata join brings in object columns with mixed types (Cluster
        # holds ints and NaN), which h5py cannot serialise. Stringify them rather
        # than dropping them, so the saved object still carries the full metadata.
        for c in emb.obs.columns:
            if emb.obs[c].dtype == object:
                emb.obs[c] = emb.obs[c].astype(str)
        emb.write(os.path.join(args.out, "embedded.h5ad"))

    with open(os.path.join(args.out, "run_log.json"), "w") as fh:
        json.dump({**log, **{k: v for k, v in vars(args).items()}}, fh, indent=2, default=str)
    print(f"\nWrote everything to {args.out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", default="gse194315",
                    help="which dataset's particulars to use (barcode key, QC gate, "
                         "cell-type column and mapping, stratum columns, starting "
                         "parameters). See sc_ids/datasets/.")
    ap.add_argument("--data", required=True)
    ap.add_argument("--prefix", default=None)
    ap.add_argument("--metadata", default=None,
                    help="per-cell metadata TSV to join. Omit when --data is an "
                         ".h5ad that already carries it.")
    ap.add_argument("--reuse-ids", action="store_true",
                    help="load ids_matrix.npy from --out instead of recomputing. "
                         "Use for the second pass once --min-ids is chosen.")
    ap.add_argument("--metadata-sample", default=None)
    ap.add_argument("--celltype-map", default=None,
                    help="granular<TAB>main TSV. Defaults to the mapping registered "
                         "for --dataset; pass a path to override.")
    ap.add_argument("--discovery-col", default="Subject",
                    help="obs column defining the IDS stratum (default Subject). "
                         "Programs are DISCOVERED here and SCORED on the whole file.")
    ap.add_argument("--discovery-value", default=None)
    # These two are NOT universal constants. The defaults come from the chosen
    # dataset, where they were derived from that dataset's own IDS histogram.
    ap.add_argument("--min-cells-frac", type=float, default=None)
    ap.add_argument("--min-ids", type=float, default=None)
    ap.add_argument("--min-genes", type=int, default=P.PROGRAM_SIZE_MIN)
    ap.add_argument("--max-genes", type=int, default=P.PROGRAM_SIZE_MAX)
    ap.add_argument("--top-n", type=int, default=None,
                    help="report only the top N ranked programs (default: all that pass)")
    ap.add_argument("--block-size", type=int, default=None)
    ap.add_argument("--n-pcs", type=int, default=None)
    ap.add_argument("--no-umap", action="store_true")
    ap.add_argument("--out", required=True,
                    help="output directory for the ranked-program table, the two "
                         "enrichment tables, the specificity summary and the figures")
    args = ap.parse_args()

    # Fill unset parameters from the dataset, and say which came from where, so a
    # run log never leaves it ambiguous whether a threshold was chosen or inherited.
    try:
        _ds = get_dataset(args.dataset)
    except KeyError as e:
        ap.error(f"{e}. Registered datasets: {available_datasets()}")
    for k, v in _ds.defaults.items():
        if getattr(args, k, None) is None:
            setattr(args, k, v)
            print(f"  [defaults] {k}={v} (from dataset '{_ds.name}')")

    if args.discovery_col and args.discovery_value is None:
        ap.error("--discovery-col given without --discovery-value "
                 "(pass --discovery-col '' to discover on the whole file)")
    build(args)


if __name__ == "__main__":
    main()
