"""
celltypes.py — collapse granular cell-subtype labels onto a handful of MAIN cell
types, so a UMAP is legible (Agenda 1).

WHY: an Azimuth-style annotation has dozens of levels (for PBMC: CD4 TCM, CD8 TEM,
B intermediate, ... 30 in all). Colouring a UMAP by 30 categories is unreadable —
the palette recycles, the legend dominates the figure, and a lineage continuum
reads as noise. Mapping to ~5 main types makes the global structure visible. The
granular labels are NOT discarded: the enrichment table in programs.py is the tool
for resolving subtypes, and it does that better than a UMAP can.

The mapping is DATA, not code. It lives in resources/<dataset>_celltype_mapping.txt
as a two-column TSV (granular <TAB> main) supplied by the lab, and the path must be
passed explicitly — THERE IS NO DEFAULT. A PBMC mapping applied to skin would map
nothing and silently label every cell "Unmapped"; a mapping is only ever valid for
the annotation it was written against. Use `sc_ids.datasets.get(<name>).celltype_map_path()`.
"""
from __future__ import annotations

import warnings

# Verbatim typos in the supplied mapping file that we repair on load rather than
# editing the file (the file is the lab's artefact; we keep it byte-identical so
# a future re-send diffs cleanly). Every repair is announced, never silent.
_TYPO_FIXES = {"K Cell": "NK Cell"}

# Plot order for the collapsed labels: lineages first, catch-all last, so the
# legend and the colour cycle are stable across figures and across runs.
MAIN_ORDER = ["T Cell", "NK Cell", "B Cell", "Myeloid Cell", "Other", "Unmapped"]


def load_celltype_map(path, fix_typos=True, verbose=True):
    """Read the granular -> main cell-type mapping.

    Parameters
    ----------
    path : str
        Two-column TSV, `granular <TAB> main`. Required — see the module
        docstring for why there is no default.
    fix_typos : bool
        Apply _TYPO_FIXES to the MAIN column. The file as supplied assigns
        NK_CD56bright and NK Proliferating to "K Cell", which would otherwise
        appear as a sixth category holding two NK subtypes while "NK" itself sat
        elsewhere — a silent split of the NK lineage in every figure.

    Returns
    -------
    dict  granular label -> main label
    """
    mapping, fixed = {}, {}
    with open(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = [p.strip() for p in line.split("\t")]
            if len(parts) < 2:
                parts = [p.strip() for p in line.rsplit("  ", 1)]
            if len(parts) < 2:
                raise ValueError(f"{path}:{lineno}: expected 'granular<TAB>main', got {line!r}")
            granular, main = parts[0], parts[1]
            if fix_typos and main in _TYPO_FIXES:
                fixed.setdefault(_TYPO_FIXES[main], []).append(granular)
                main = _TYPO_FIXES[main]
            mapping[granular] = main
    if verbose:
        print(f"  [celltypes] loaded {len(mapping)} granular labels -> "
              f"{len(set(mapping.values()))} main types from "
              f"{path.rsplit('/', 1)[-1]}")
        for main, subs in fixed.items():
            print(f"  [celltypes] WARNING: repaired typo in mapping file -> "
                  f"'{main}' for {subs} (file left unchanged)")
    return mapping


def add_celltype_main(adata, mapping=None, path=None, source_col="CellType",
                      out_col="CellTypeMain", unmapped="Unmapped", verbose=True):
    """Add obs[out_col] = main cell type, derived from obs[source_col].

    Unmapped labels are kept as an explicit 'Unmapped' category rather than NaN,
    so they stay visible on the UMAP instead of silently vanishing. The number of
    unmapped cells is reported — a nonzero count means the mapping file is stale
    relative to the annotation and should be reconciled before the figure is used.

    Returns the same adata (mutated in place) for chaining.
    """
    import pandas as pd

    if source_col not in adata.obs.columns:
        raise KeyError(f"obs['{source_col}'] not found — run the metadata join "
                       f"first (--metadata); CellType is annotated only for "
                       f"IncludedInStudy cells.")
    if mapping is None:
        if path is None:
            raise ValueError("pass either `mapping` or `path` — there is no default "
                             "cell-type mapping (see module docstring)")
        mapping = load_celltype_map(path, verbose=verbose)

    src = adata.obs[source_col].astype(str)
    main = src.map(mapping).fillna(unmapped)

    present = [c for c in MAIN_ORDER if c in set(main)]
    present += sorted(set(main) - set(present))
    adata.obs[out_col] = pd.Categorical(main.values, categories=present, ordered=False)

    if verbose:
        missing = sorted(set(src[main == unmapped]))
        n_un = int((main == unmapped).sum())
        print(f"  [celltypes] obs['{out_col}']: {src.nunique()} granular -> "
              f"{len(present)} main ({', '.join(present)})")
        if n_un:
            warnings.warn(f"{n_un} cells ({n_un / len(main):.1%}) have a CellType "
                          f"absent from the mapping: {missing}", RuntimeWarning)
            print(f"  [celltypes] WARNING: {n_un} cells unmapped, labels: {missing}")
    return adata


def coverage_report(mapping, observed_labels):
    """Two-sided reconciliation between the mapping file and the data.

    Returns (missing_from_map, unused_in_map). A nonempty second list is not an
    error but is worth knowing: it means the file was written for a superset of
    this dataset, so a label you expect to see may simply be absent here.
    """
    obs = set(map(str, observed_labels))
    return sorted(obs - set(mapping)), sorted(set(mapping) - obs)


def map_from_obs(adata, granular_col, main_col, verbose=True):
    """Granular -> main mapping read off the data, for studies that annotate both.

    Majority vote per granular label. A handful of cells usually carry a
    different main label than their siblings (annotation noise — in GSE173706
    seven granular labels have 1-4 such cells), and taking the mode rather than
    the first value keeps those from silently redirecting a whole subtype.
    Disagreements are reported, never swallowed.
    """
    import pandas as pd
    df = adata.obs[[granular_col, main_col]].astype(str)
    counts = df.groupby([granular_col, main_col]).size().rename("n").reset_index()
    top = counts.sort_values("n", ascending=False).drop_duplicates(granular_col)
    mapping = dict(zip(top[granular_col], top[main_col]))
    if verbose:
        noisy = counts.groupby(granular_col).size()
        noisy = noisy[noisy > 1]
        print(f"  [celltypes] mapping read from obs['{granular_col}'] -> "
              f"obs['{main_col}']: {len(mapping)} granular -> "
              f"{len(set(mapping.values()))} main")
        if len(noisy):
            n_minor = int(counts.merge(top[[granular_col, main_col]], how="left",
                                       on=[granular_col, main_col], indicator=True)
                          .query("_merge == 'left_only'")["n"].sum())
            print(f"  [celltypes] {len(noisy)} granular label(s) span more than one "
                  f"main type ({n_minor} cells total); majority vote applied")
    return mapping
