"""
build_disease_object.py — one joint blood + skin object carrying DISEASE STATE.

OLDER DATA LAYOUT. This builds one merged object per tissue, which needs more memory with
every lane added (the nine-lane blood build was killed at 6 GB). For new work use
scripts/04_method_checks/build_units.py, which writes one file per lane with the same
filters and gives exactly the same cells. Kept because stages 02 and 03 still read the
merged objects (blood_qc.h5ad, skin_qc.h5ad).

Why this exists: every earlier object in this project had tissue and cell type
but no disease label, and the single PBMC lane we had been using (PBMC-01-1) is
100% psoriatic arthritis, so no disease contrast was possible in blood at all.

What goes in:
  BLOOD — GSE194315, nine 10x lanes (PBMC-01-1, 03-1, 03-3, 05-1, 05-2, 06-1 to 06-4;
    roles in sc_ids.datasets.gse194315.LANES_USED and docs/data_inventory.md). Each
    lane pools nine genetically distinct subjects, demultiplexed with demuxlet,
    so a lane is NOT a sample. We keep the source study's own verdict
    (IncludedInStudy) and demuxlet singlets, then keep only subjects whose
    Status is Healthy, PSO or PSA. AXI and PSX are dropped — different diseases.
    Lanes 03 and 05 each contain Healthy, PSO and PSA subjects TOGETHER, which
    is the point of including them: disease is not confounded with lane.
  SKIN — GSE173706, 22 donors, disease PP (psoriasis lesional), PN (psoriasis
    non-lesional) and NS (normal skin).

Harmonised disease column, so the two tissues can be compared:
  healthy   <- blood Healthy,  skin NS
  psoriasis <- blood PSO,      skin PP and PN (kept apart in `disease_fine`)
  PsA       <- blood PSA.      NO SKIN COUNTERPART — the PsA-vs-psoriasis
                               contrast is blood-only, by construction.

Cell QC is the lab pipeline's (min_gene_per_cell 200, max_mt_percentage 10).
No doublet-score filter: it does not exist for the skin data, and for blood the
demuxlet singlet gate already removes cross-donor doublets.

Cell-type labels are the SOURCE annotations, used as given. They are the
"reckless" labels: good enough to run IDS against, and the programs are
discovered without them, so they can be revised afterwards without rerunning
anything upstream.

Usage:
  python scripts/02_ids_programs/build_disease_object.py --out results/02_ids_programs/ids-disease
"""
import argparse, os, sys
import scanpy as sc
import anndata as ad
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from sc_ids.pipeline import load_data, join_metadata                      # noqa: E402
from sc_ids.datasets.psoriasis_skin import build_h5ad                     # noqa: E402
from sc_ids.datasets.gse194315 import LANES_USED                          # noqa: E402
from sc_ids import paths                                                      # noqa: E402

PBMC_DIR, PBMC_META = paths.PBMC_DIR, paths.PBMC_META
SKIN_DIR, SKIN_META = paths.SKIN_DIR, paths.SKIN_META
FEATURES = f"{PBMC_DIR}/PBMC-01-1.features.tsv.gz"
LANES = list(LANES_USED)
KEEP_STATUS = {"Healthy", "PSO", "PSA"}
DISEASE_MAP = {"Healthy": "healthy", "PSO": "psoriasis", "PSA": "PsA",
               "NS": "healthy", "PP": "psoriasis", "PN": "psoriasis"}
MIN_GENES, MAX_MT = 200, 10.0
# Non-nucleated cells the lab thresholds cannot catch: they filter on gene count and
# mitochondrial share, not on cell identity, and a platelet can clear 200 genes on
# residual transcripts alone. Removing the CELLS does not remove platelet or
# haemoglobin AMBIENT RNA, which still forms programmes inside other cells -- those
# are flagged as technical in the catalogue rather than deleted.
EXCLUDE_SUBTYPES = {"Platelet", "Eryth"}


def lab_qc(a, tag):
    """The lab pipeline's cell QC: min_gene_per_cell 200, max_mt_percentage 10."""
    n0 = a.n_obs
    a.var["mt"] = a.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    keep = (a.obs["n_genes_by_counts"] >= MIN_GENES) & (a.obs["pct_counts_mt"] <= MAX_MT)
    a = a[keep.values].copy()
    print(f"  [qc] {tag}: {n0:,} -> {a.n_obs:,} cells "
          f"({n0 - a.n_obs:,} removed)", flush=True)
    return a


def build_blood(scratch=None):
    """Lane by lane, each written to disk and freed before the next is opened.

    Holding all nine lanes in memory to concatenate them was OOM-killed at 6 GB:
    ad.concat needs every input plus the result resident at once. QC and the
    non-nucleated exclusion are per-cell, so doing them per lane is equivalent and
    keeps peak memory at one lane. The lanes are then merged ON DISK.
    """
    import tempfile
    from anndata.experimental import concat_on_disk
    scratch = tempfile.mkdtemp(prefix="blood_", dir=scratch or os.environ.get("SCRATCH"))
    parts = []
    for lane in LANES:
        a = load_data(PBMC_DIR, prefix=lane + ".")
        a.var_names_make_unique()
        a = join_metadata(a, PBMC_META, barcode_col="CellName", sample=lane,
                          singlets_only=True, included_only=True, verbose=False)
        a = a[a.obs["Status"].isin(KEEP_STATUS)].copy()
        a.obs["lane"] = lane
        a.obs = a.obs.rename(columns={"Subject": "donor", "CellType": "subtype",
                                      "Status": "disease_source"})
        a.obs["tissue"] = "PBMC"
        a.obs_names = [f"{lane}_{i}" for i in a.obs_names]
        n_in = a.n_obs
        a = lab_qc(a, lane)
        a = a[~a.obs["subtype"].astype(str).isin(EXCLUDE_SUBTYPES)].copy()
        # Keep only the columns we use, as strings. The source metadata carries
        # mixed-type columns (Cluster holds both labels and NaN) that h5py cannot
        # write, and concat_on_disk needs consistent dtypes across lanes anyway.
        keep = ["donor", "subtype", "disease_source", "lane", "tissue",
                "n_genes_by_counts", "total_counts", "pct_counts_mt"]
        a.obs = a.obs[[c for c in keep if c in a.obs.columns]].copy()
        for c in ("donor", "subtype", "disease_source", "lane", "tissue"):
            a.obs[c] = a.obs[c].astype(str)
        # concat_on_disk refuses CSC; scanpy's mtx reader gives CSC here.
        if sp.issparse(a.X) and a.X.format == "csc":
            a.X = a.X.tocsr()
        f = os.path.join(scratch, f"{lane}.h5ad")
        a.write_h5ad(f)
        print(f"  [blood] {lane}: {n_in:,} -> {a.n_obs:,} cells after QC  "
              + str(a.obs["disease_source"].value_counts().to_dict()), flush=True)
        parts.append(f)
        del a
    merged = os.path.join(scratch, "merged.h5ad")
    concat_on_disk(parts, merged, join="inner", label=None, index_unique=None)
    b = ad.read_h5ad(merged)
    print(f"  [blood] merged on disk: {b.n_obs:,} cells x {b.n_vars:,} genes", flush=True)
    return b


def build_skin(cache):
    if not os.path.exists(cache):
        print("  [skin] building from the GEO CSVs (one-off, slow)", flush=True)
        build_h5ad(SKIN_DIR, SKIN_META, FEATURES, cache,
                   scratch_dir=os.environ.get("SCRATCH", "/tmp"))
    s = ad.read_h5ad(cache)
    s.obs = s.obs.rename(columns={"disease": "disease_source"})
    s.obs["tissue"] = "skin"
    keep = [c for c in ("donor", "disease_source", "celltype", "subtype",
                        "fullsubtype", "orig.ident") if c in s.obs.columns]
    s.obs = s.obs[keep]
    if "celltype" in s.obs.columns:
        s.obs["subtype"] = s.obs["fullsubtype"] if "fullsubtype" in s.obs.columns else s.obs["subtype"]
    s.obs["tissue"] = "skin"
    print(f"  [skin] {s.n_obs:,} cells  "
          + str(s.obs["disease_source"].value_counts().to_dict()), flush=True)
    return lab_qc(s, "skin")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--skin-cache", default=paths.SKIN_CACHE)
    ap.add_argument("--tissue", required=True, choices=["blood", "skin"])
    ap.add_argument("--scratch", default=None)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    # ONE TISSUE PER PROCESS, AND NO GENE SUBSET. Two OOM kills at 6 GB taught
    # this layout. Concatenating the tissues needs both inputs plus the result
    # resident; even holding both tissues while column-subsetting one of them is
    # too much at 67k cells x 33k genes. Nothing downstream needs them combined:
    # every IDS stratum sits inside one tissue, and the cross-tissue step compares
    # GENE SETS, not cells. The gene intersection is applied later, when the
    # ~3k-gene panel is built, where it is free.
    if args.tissue == "blood":
        a = build_blood(scratch=args.scratch)        # QC and exclusion done per lane
    else:
        a = build_skin(args.skin_cache)
        n0 = a.n_obs
        a = a[~a.obs["subtype"].astype(str).isin(EXCLUDE_SUBTYPES)].copy()
        if a.n_obs != n0:
            print(f"  [exclude] {n0 - a.n_obs} non-nucleated cells removed", flush=True)
    a.obs["disease"] = a.obs["disease_source"].astype(str).map(DISEASE_MAP).astype(str)
    a.obs["disease_fine"] = a.obs["disease_source"].astype(str)
    a.obs["stratum"] = (a.obs["tissue"].astype(str) + "_"
                        + a.obs["disease_fine"].astype(str))
    for c in ("tissue", "disease", "disease_fine", "donor", "subtype", "stratum"):
        a.obs[c] = a.obs[c].astype(str)

    p = os.path.join(args.out, f"{args.tissue}_qc.h5ad")
    a.write_h5ad(p)
    print(f"[out] {p}: {a.n_obs:,} cells x {a.n_vars:,} genes", flush=True)

    obs = a.obs[["tissue", "disease", "disease_fine", "donor", "subtype", "stratum"]]
    obs.to_csv(os.path.join(args.out, f"cells_{args.tissue}.csv"))
    tab = (obs.groupby(["tissue", "disease", "disease_fine"], observed=True)
             .agg(cells=("donor", "size"), donors=("donor", "nunique")).reset_index())
    tab.to_csv(os.path.join(args.out, f"design_{args.tissue}.csv"), index=False)
    print()
    print(tab.to_string(index=False))


if __name__ == "__main__":
    main()
