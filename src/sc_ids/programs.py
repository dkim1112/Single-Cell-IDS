"""
programs.py — turn a raw IDS module list into a RANKED, REPORTABLE set of gene
programs (Agenda 1b) and quantify which cell types express each one (Agenda 3).

Two things live here because they are two halves of one claim. Agenda 1b decides
*which* programs we are allowed to talk about; Agenda 3 decides *what each one
means*. Separating them lets the selection rule be fixed before any biology is
looked at, which is the entire point of the exercise.

----------------------------------------------------------------------
THE PRE-REGISTERED SELECTION RULE (Agenda 1b)
----------------------------------------------------------------------
The previous session picked figures with `[m for m in mods if len(m) <= 60][:6]`
— an undisclosed 60-gene cap and a top-6 cut on MEAN IDS. Two problems: the cap
was chosen after seeing the modules, and mean IDS rewards a large module carried
by a few very strong pairs.

Replaced by a rule stated in advance and applied without exception:

  1. SIZE FILTER: keep programs with PROGRAM_SIZE_MIN <= n_genes <= PROGRAM_SIZE_MAX
     (default 3-50). Lower bound: a 2-gene "program" is a gene pair, reportable as
     a pair but not as a program. Upper bound: above ~50 genes a greedy-modularity
     community is a compartment (ribosome, mitochondria) rather than a program, and
     its members no longer share a single interpretation.
  2. RANK by (recurrence, median IDS) descending, in that order — recurrence first
     because a program seen in 15/18 strata at IDS 0.6 is a better finding than one
     seen once at 0.95. Ties broken by edge density.
  3. REPORT ALL of the top N (default: all that pass), never a hand-picked subset.

MEDIAN, not mean: the median within-program IDS is the level at which HALF the
gene pairs sit, so a program cannot be inflated by a couple of near-duplicate
genes (e.g. two haemoglobin chains at IDS 0.99 inside an otherwise weak set).

RECURRENCE IS UNAVAILABLE FROM A SINGLE RUN. It is defined across strata (Agenda
2). With one run, `rank_programs` ranks on median IDS alone and says so out loud
rather than pretending the primary key was applied.

----------------------------------------------------------------------
WHICH REPRESENTATION IS SCORED (Agenda 3)
----------------------------------------------------------------------
A THIRD representation of the same counts, on purpose (the repo already keeps two
— see viz.py). IDS runs on normalized counts min-max scaled to [0,8], no log1p.
Per-cell program scores use `sc.tl.score_genes`, which subtracts the mean of a
matched control gene set and therefore expects the standard log1p-normalized
representation. Scoring in the [0,8] IDS space would make the control subtraction
meaningless. Programs are DISCOVERED in IDS space and SCORED in log space.

"Expressing" is defined as score > 0 and that cutoff is fixed in advance. It is
not arbitrary: score_genes is already background-subtracted against a control set
matched on expression bin, so 0 is the value expected when a program's genes are
no higher than equivalently-expressed genes elsewhere. It needs no tuning and is
comparable across programs of different size and mean expression.
"""
from __future__ import annotations

import numpy as np

PROGRAM_SIZE_MIN = 3
PROGRAM_SIZE_MAX = 50
EXPRESSING_CUTOFF = 0.0     # sc.tl.score_genes is background-subtracted; 0 is the null


# ----------------------------------------------------------------------
# Agenda 1b — annotate and rank
# ----------------------------------------------------------------------
def annotate_programs(modules, C, min_ids=None):
    """Attach distribution statistics of the within-program IDS values.

    `mean_ids` alone (what extract_modules returns) hides the shape of the
    distribution. We add the median (the ranking key), the quartiles, the weakest
    pair, and the edge density — the fraction of within-program pairs that would
    themselves have passed `min_ids`. A program with density 0.2 is a chain held
    together by a few links; density near 1 is a genuine clique.
    """
    out = []
    for m in modules:
        idx = np.asarray(m["idx"])
        sub = C[np.ix_(idx, idx)].astype(float)
        vals = sub[np.triu_indices_from(sub, k=1)]
        vals = vals[np.isfinite(vals)]
        d = dict(m)
        d.update(
            n_genes=len(idx),
            median_ids=float(np.median(vals)) if vals.size else float("nan"),
            mean_ids=float(np.mean(vals)) if vals.size else float("nan"),
            q25_ids=float(np.percentile(vals, 25)) if vals.size else float("nan"),
            min_pair_ids=float(vals.min()) if vals.size else float("nan"),
            max_pair_ids=float(vals.max()) if vals.size else float("nan"),
            density=float((vals >= min_ids).mean()) if (vals.size and min_ids is not None) else float("nan"),
        )
        out.append(d)
    return out


def rank_programs(modules, min_genes=PROGRAM_SIZE_MIN, max_genes=PROGRAM_SIZE_MAX,
                  top_n=None, verbose=True):
    """Apply the pre-registered rule above. Returns (kept, dropped_summary).

    `modules` must already carry median_ids (run annotate_programs first).
    Recurrence is used as the primary key only if EVERY module carries a finite
    'recurrence'; otherwise ranking falls back to median IDS and says so.
    """
    have_rec = bool(modules) and all(np.isfinite(m.get("recurrence", np.nan)) for m in modules)

    kept, dropped = [], {"too_small": [], "too_large": []}
    for m in modules:
        n = m.get("n_genes", len(m["genes"]))
        if n < min_genes:
            dropped["too_small"].append((n, m["genes"][:4]))
        elif n > max_genes:
            dropped["too_large"].append((n, m["genes"][:4]))
        else:
            kept.append(m)

    if have_rec:
        kept.sort(key=lambda m: (-m["recurrence"], -m["median_ids"], -m.get("density", 0.0)))
    else:
        kept.sort(key=lambda m: (-m["median_ids"], -m.get("density", 0.0)))

    for rank, m in enumerate(kept, 1):
        m["rank"] = rank
        m["name"] = f"P{rank:02d}_{m['genes'][0]}"

    if top_n is not None:
        kept = kept[:top_n]

    if verbose:
        key = ("recurrence, then median IDS" if have_rec
               else "median IDS only (recurrence needs >1 stratum — see Agenda 2)")
        print(f"  [programs] size filter {min_genes}-{max_genes} genes: kept {len(kept)}, "
              f"dropped {len(dropped['too_small'])} too small / "
              f"{len(dropped['too_large'])} too large")
        print(f"  [programs] ranked by {key}"
              + (f"; reporting top {top_n}" if top_n else "; reporting all"))
    return kept, dropped


def programs_table(programs):
    """Ranked programs as a tidy DataFrame (one row per program)."""
    import pandas as pd
    rows = []
    for m in programs:
        rows.append({
            "rank": m.get("rank"),
            "program": m.get("name"),
            "n_genes": m.get("n_genes", len(m["genes"])),
            "median_ids": m.get("median_ids"),
            "q25_ids": m.get("q25_ids"),
            "mean_ids": m.get("mean_ids"),
            "min_pair_ids": m.get("min_pair_ids"),
            "density": m.get("density"),
            "recurrence": m.get("recurrence", np.nan),
            "genes": ", ".join(m["genes"]),
        })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# Agenda 3 — per-cell scoring and cell-type enrichment
# ----------------------------------------------------------------------
def score_programs(adata, programs, ctrl_size=50, verbose=True):
    """Per-cell score for each program via sc.tl.score_genes, on the CURRENT
    adata.X — which must be the log1p-normalized representation (see module
    docstring). Returns the list of obs column names created, in rank order.

    Genes absent from adata.var_names are dropped with a count; a program that
    loses more than half its genes is skipped rather than silently reported on a
    remnant.
    """
    import scanpy as sc
    if adata.n_vars < 1000:
        print(f"  [programs] WARNING: scoring on only {adata.n_vars} genes. "
              f"score_genes matches controls by expression bin; on a heavily "
              f"gene-filtered object the control pool can contain the very "
              f"programs being scored, which drives scores toward 0. Score on "
              f"the unfiltered object where possible.")
    names = []
    for m in programs:
        present = [g for g in m["genes"] if g in adata.var_names]
        if len(present) < 2 or len(present) < 0.5 * len(m["genes"]):
            if verbose:
                print(f"  [programs] skip {m.get('name')}: only {len(present)}/"
                      f"{len(m['genes'])} genes present in the scored object")
            continue
        name = m.get("name") or f"P_{m['genes'][0]}"
        # ctrl_size must be bounded by the genes actually AVAILABLE, not by the
        # program's size: score_genes draws its control set from expression bins,
        # and asking for 50 controls out of a few hundred genes empties the bins
        # and raises. Real scRNA-seq (~20k genes) always takes the 50 default.
        sc.tl.score_genes(adata, gene_list=present, score_name=name,
                          ctrl_size=max(5, min(ctrl_size, adata.n_vars // 10)))
        m["scored_genes"] = present
        names.append(name)
    if verbose:
        print(f"  [programs] scored {len(names)} programs on {adata.n_obs} cells")
    return names


def _wilson(k, n, z=1.959963984540054):
    """Wilson 95% CI for a proportion. Used because cell-type sizes here span
    three orders of magnitude (ASDC ~55 cells vs CD14 Mono ~65,000): a raw
    proportion from 55 cells is not comparable to one from 65,000, and sorting a
    table by proportion without showing that is how a rare subtype lands at the
    top of the list on nothing but sampling noise."""
    n = np.asarray(n, dtype=float)
    k = np.asarray(k, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = k / n
        denom = 1 + z ** 2 / n
        centre = (p + z ** 2 / (2 * n)) / denom
        half = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / denom
    return np.clip(centre - half, 0, 1), np.clip(centre + half, 0, 1)


def celltype_enrichment(adata, score_names, celltype_col="CellType",
                        cutoff=EXPRESSING_CUTOFF, min_cells=0):
    """AGENDA 3. For each (program, cell type): what proportion of that cell
    type's cells express the program, how strongly, and how that compares to the
    dataset as a whole.

    Returns a tidy DataFrame sorted by program, then descending proportion — so
    the head of each program's block is the cell type most associated with it.

    Columns
    -------
    prop_expressing : fraction of that cell type's cells with score > cutoff.
                      THE SORT KEY, per the agenda: "how widespread".
    ci_low/ci_high  : Wilson 95% interval on that proportion. Read it before
                      believing an ordering — overlapping intervals mean the rank
                      difference is not supported by the number of cells.
    mean_score      : "how strong", among ALL cells of the type (not only the
                      expressing ones), so it is not inflated by a small tail.
    prop_overall    : the same proportion across every cell in the object.
    enrichment      : prop_expressing / prop_overall. A program expressed by 85%
                      of Platelets and 80% of everything else has a high
                      proportion but an enrichment of ~1 and is NOT platelet
                      biology. Proportion answers "widespread"; enrichment is what
                      makes the ranking mean something. Report both.

    `min_cells` drops cell types with too few cells to interpret; the default of 0
    keeps everything and leaves the judgement to the CI column.
    """
    import pandas as pd

    if celltype_col not in adata.obs.columns:
        raise KeyError(f"obs['{celltype_col}'] not found")
    ct = adata.obs[celltype_col].astype(str).values

    rows = []
    for name in score_names:
        s = np.asarray(adata.obs[name].values, dtype=float)
        expressing = s > cutoff
        prop_overall = float(expressing.mean())
        for t in pd.unique(ct):
            m = ct == t
            n = int(m.sum())
            if n < min_cells:
                continue
            k = int(expressing[m].sum())
            rows.append({
                "program": name,
                "cell_type": t,
                "n_cells": n,
                "n_expressing": k,
                "prop_expressing": k / n if n else np.nan,
                "mean_score": float(s[m].mean()) if n else np.nan,
                "median_score": float(np.median(s[m])) if n else np.nan,
                "prop_overall": prop_overall,
                "enrichment": (k / n) / prop_overall if (n and prop_overall > 0) else np.nan,
            })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    lo, hi = _wilson(df["n_expressing"].values, df["n_cells"].values)
    df["ci_low"], df["ci_high"] = lo, hi
    # Sort by proportion (the agenda's key), TIE-BROKEN BY MEAN SCORE. The
    # tie-break is not cosmetic: with score > 0 as the cutoff, a strong program
    # saturates at prop = 1.000 in every type that expresses it at all, and on real
    # data that was five cell types at once. Without the tie-break the head of the
    # list is then whichever type pandas happened to encounter first.
    df = df.sort_values(["program", "prop_expressing", "mean_score"],
                        ascending=[True, False, False])
    df["rank_in_program"] = df.groupby("program").cumcount() + 1
    # A SECOND ranking by strength alone. Report both: proportion answers "in how
    # many cells of this type is the program on", mean score answers "how strongly".
    # They can disagree sharply for a bimodal program (measured: the haemoglobin
    # program has its HIGHEST mean score in Eryth while ranking near the BOTTOM of
    # Eryth by proportion, because it is huge in a few Eryth cells and absent in the
    # rest). Whichever single column you sort on, the other one is the check.
    df["rank_by_mean_score"] = (df.groupby("program")["mean_score"]
                                  .rank(ascending=False, method="first").astype(int))
    cols = ["program", "rank_in_program", "rank_by_mean_score", "cell_type",
            "n_cells", "n_expressing", "prop_expressing", "ci_low", "ci_high",
            "mean_score", "median_score", "prop_overall", "enrichment"]
    return df[cols].reset_index(drop=True)


def specificity_summary(enrichment_df, min_cells=100):
    """One row per program: is this ranking actually telling us anything?

    A program whose top cell type is expressed at 0.92 and whose median cell type
    is at 0.88 is ubiquitous — the sort order is noise. `gap_to_second` and
    `top_enrichment` are what separate a cell-identity program from a
    housekeeping one. Cell types below `min_cells` are excluded from the
    'top' calculation so a 55-cell subtype cannot win on sampling noise, but they
    remain in the full table.
    """
    import pandas as pd
    df = enrichment_df[enrichment_df["n_cells"] >= min_cells]
    rows = []
    for prog, g in df.groupby("program", sort=False):
        g = g.sort_values(["prop_expressing", "mean_score"], ascending=[False, False])
        top = g.iloc[0]
        second = g.iloc[1] if len(g) > 1 else None
        by_mean = g.sort_values("mean_score", ascending=False).iloc[0]
        n_tied = int((g["prop_expressing"] >= top["prop_expressing"] - 1e-12).sum())
        rows.append({
            "program": prog,
            "top_cell_type": top["cell_type"],
            "top_prop": top["prop_expressing"],
            "top_n_cells": int(top["n_cells"]),
            "gap_to_second": (top["prop_expressing"] - second["prop_expressing"]) if second is not None else np.nan,
            "second_cell_type": second["cell_type"] if second is not None else None,
            "median_prop_across_types": float(g["prop_expressing"].median()),
            "top_enrichment": top["enrichment"],
            # How many cell types share the top proportion. >1 means the agenda's
            # sort key cannot separate them and the ordering came from the
            # mean-score tie-break, not from the proportion.
            "n_tied_at_top_prop": n_tied,
            "saturated": bool(n_tied > 1),
            # The independent answer from strength rather than breadth. When this
            # disagrees with top_cell_type, neither column is safe to quote alone.
            "top_by_mean_score": by_mean["cell_type"],
            "top_mean_score": float(by_mean["mean_score"]),
            "agree": bool(by_mean["cell_type"] == top["cell_type"]),
            "specific": bool(top["enrichment"] >= 1.5 and
                             (top["prop_expressing"] - float(g["prop_expressing"].median())) >= 0.2),
        })
    return pd.DataFrame(rows).sort_values("top_enrichment", ascending=False).reset_index(drop=True)


def lineage_column_order(cell_types, mapping):
    """Order granular cell types so that lineages sit together on the heatmap.

    Ordering columns by their average value (the obvious default) interleaves B
    subtypes with monocytes with T subtypes, which is precisely the readability
    problem Agenda 1 set out to fix — reintroduced one figure later. Grouping by
    the CellTypeMain mapping instead lets you read a lineage block at a glance and
    still see which subtype inside it breaks rank.
    """
    from .celltypes import MAIN_ORDER
    return sorted(cell_types,
                  key=lambda t: (MAIN_ORDER.index(mapping.get(t, "Unmapped"))
                                 if mapping.get(t, "Unmapped") in MAIN_ORDER
                                 else len(MAIN_ORDER), str(t)))


def plot_enrichment_heatmap(enrichment_df, out_png, value="prop_expressing",
                            programs=None, title=None, figsize=None,
                            column_order=None):
    """Programs (rows) x cell types (columns) heatmap of `value`.

    This is the figure that replaces squinting at a 30-colour UMAP: every
    granular subtype gets its own column whether or not it forms a visible island
    in two dimensions.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    piv = enrichment_df.pivot(index="program", columns="cell_type", values=value)
    if programs is not None:
        piv = piv.reindex([p for p in programs if p in piv.index])
    if column_order is not None:
        piv = piv[[c for c in column_order if c in piv.columns]]
    else:
        piv = piv[piv.mean(axis=0).sort_values(ascending=False).index]

    h = max(2.2, 0.55 * len(piv) + 1.0) ; w = max(5.0, 0.34 * piv.shape[1] + 3)
    fig, ax = plt.subplots(figsize=figsize or (w, h))
    im = ax.imshow(piv.values, aspect="auto", cmap="viridis",
                   vmin=0, vmax=1 if value == "prop_expressing" else None)
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels(piv.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels(piv.index, fontsize=7)
    ax.set_title(title or f"{value} by program x cell type", fontsize=9)
    ax.set_xticks(np.arange(-0.5, piv.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, piv.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.6)
    ax.tick_params(which="minor", length=0)
    fig.colorbar(im, ax=ax, shrink=0.7, label=value)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [programs] wrote {out_png}")
    return out_png
