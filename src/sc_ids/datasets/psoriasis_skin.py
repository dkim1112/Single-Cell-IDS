"""
GSE173706 — psoriasis skin, single-cell RNA-seq.
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173706

33 samples, three phenotypes:
    NS  normal skin, healthy donors
    PN  non-lesional skin, psoriasis patients
    PP  lesional skin, psoriasis patients

WHAT IS DIFFERENT FROM GSE194315, and why each difference needed code
---------------------------------------------------------------------
1. BARCODE KEY. GSE194315 keyed metadata on '<Sample>_<barcode>'. Here the
   metadata is a single merged Seurat object, so every cell is
   '<barcode>-<n>' where n is that sample's index in the merge (1..33) — NOT
   the 10x '-1' that each per-sample file carries. Verified from the metadata:
   the suffix<->sample map is exactly 1:1 (33 suffixes, 33 samples, no
   collisions), so `suffix_key` below rebuilds the key deterministically.
   Getting this wrong does not error, it silently matches ~1/33 of the cells.

2. NO QC FLAG COLUMN. GSE194315 published `IncludedInStudy`. Here the gate is
   membership: only cells that passed the source study's filtering appear in the
   metadata at all, so the join itself is the QC step (dropping ~30% — 96,057
   raw cells to 67,378). `qc_gate=None` is therefore correct, not an omission.

3. THE CELL-TYPE MAPPING IS ALREADY IN THE METADATA. No file needed from the
   lab: `fullsubtype` holds 26 granular labels and `celltype` holds 10 main ones.
   Seven fullsubtypes appear under more than one celltype, but only for 1-4
   stray cells each (e.g. one 'Spinous Keratinocyte' labelled 'Eccrine gland'),
   so the parent is taken by majority vote — see `celltype_map_from_metadata`.

4. GENE IDs ARE ENSEMBL. GSE194315 gave symbols (HBB, MS4A1); these CSVs give
   ENSG00000244734. Programs reported as ENSG IDs are unreadable, and — more
   importantly — cannot be compared with the blood programs at all. The PBMC
   features.tsv.gz is itself an ENSG->symbol table over the same CellRanger
   reference family (33,694 vs 33,538 genes), so it is used as the map. That
   also guarantees both tissues land in ONE symbol space, which is what any
   skin<->blood comparison requires.

5. GENES ARE ROWS. The CSVs are genes x cells; load with transpose=True.

6. CHEMISTRY IS CONFOUNDED WITH PHENOTYPE. NS is 100% v3_lib; PN and PP are
   ~2/3 v2_lib. Any NS-vs-PP difference is therefore partly a v2-vs-v3
   difference. Donor 30696 appears as both '30696' (v2) and '30696V3' (v3) and
   is the control for measuring how large that effect is. Stratifying on
   `orig.ident` holds chemistry constant WITHIN a run; it is only the
   ACROSS-phenotype comparisons that are exposed.

7. PN/PP ARE PAIRED WITHIN DONOR. Eleven donors contribute both non-lesional
   and lesional samples (929, 8940, 31170, 31277, 369PC, 5851, 7802ED, 8659ED,
   9709PC, 30696, 30696V3). GSE194315 had no analogue. This is the strongest
   design feature in the dataset and the closest thing here to a directional
   handle.
"""
from . import Dataset, register

#: Sample sizes after the source study's filtering, from the metadata. 33 samples;
#: 25 clear 1,000 cells, 32 clear 500, only NS_AR008 (42 cells) is unusable.
SMALLEST_USABLE_STRATUM = 500

#: Donors contributing both PN and PP — the paired lesional/non-lesional design.
PAIRED_DONORS = ("30696", "30696V3", "31170", "31277", "369PC", "5851",
                 "7802ED", "8659ED", "8940", "929", "9709PC")


def suffix_key(obs_names, sample=None):
    """Rebuild the merged-object barcode key: '<bare barcode>-<sample suffix>'.

    `sample` here is the integer suffix this file's cells carry in the metadata,
    passed through as a string. The per-file 10x suffix ('-1') is stripped first.
    """
    import pandas as pd
    bare = pd.Index(obs_names).astype(str).str.replace(r"-\d+$", "", regex=True)
    if sample is None:
        return bare
    return pd.Index(bare + "-" + str(sample))


def suffix_for_sample(metadata_path, sample):
    """Which merged-object suffix belongs to `sample` (e.g. 'PP_929' -> '2').

    Read from the metadata rather than assumed, because the merge order is not
    recoverable from anything else. Raises if the mapping is not 1:1, which would
    mean the assumption behind `suffix_key` no longer holds for this file.
    """
    import pandas as pd
    m = pd.read_csv(metadata_path, sep="\t", index_col=0, usecols=[0, 1],
                    low_memory=False)
    suf = pd.Series([str(b).rsplit("-", 1)[-1] for b in m.index], index=m.index)
    tab = pd.crosstab(suf, m["orig.ident"])
    if (tab > 0).sum(axis=1).max() > 1 or (tab > 0).sum(axis=0).max() > 1:
        raise ValueError("suffix <-> orig.ident is no longer 1:1 in this metadata; "
                         "the barcode key rule in psoriasis_skin.py must be revisited")
    if sample not in tab.columns:
        raise KeyError(f"sample '{sample}' not in metadata; have {sorted(tab.columns)[:6]}...")
    return tab.index[tab[sample] > 0][0]


def celltype_map_from_metadata(metadata_path, granular="fullsubtype", main="celltype"):
    """Granular -> main mapping, taken from the metadata instead of a hand file.

    Majority vote, because seven granular labels carry a handful of cells under a
    different main label (annotation noise, 1-4 cells each). Returns a dict in the
    same shape `sc_ids.celltypes` expects.
    """
    import pandas as pd
    m = pd.read_csv(metadata_path, sep="\t", index_col=0,
                    usecols=[0, granular, main], low_memory=False)
    counts = m.groupby([granular, main]).size().rename("n").reset_index()
    top = counts.sort_values("n", ascending=False).drop_duplicates(granular)
    return dict(zip(top[granular].astype(str), top[main].astype(str)))


def ensembl_to_symbol(features_tsv_gz):
    """ENSG -> gene symbol from a 10x features.tsv.gz (col 1 = id, col 2 = symbol).

    GSE194315's features file is used for this. Same CellRanger reference family,
    and reusing it means skin and blood programs are named in one symbol space.
    """
    import gzip
    mapping = {}
    with gzip.open(features_tsv_gz, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] and parts[1]:
                mapping[parts[0]] = parts[1]
    return mapping


PSORIASIS_SKIN = register(Dataset(
    name="psoriasis_skin",
    description="GSE173706 — psoriasis skin, 33 samples, NS/PN/PP, genes as rows, Ensembl IDs",

    barcode_col=None,          # the metadata's index IS the barcode column
    make_key=suffix_key,
    qc_gate=None,              # membership in the metadata is the gate — see note 2
    singlet_gate=None,         # not multiplexed; one sample per file

    celltype_col="fullsubtype",          # 26 granular labels
    celltype_map_file=None,              # no hand file needed — see note 3
    celltype_main_col="celltype",        # 10 main labels, already in the metadata

    # One sample fixes donor, phenotype AND chemistry simultaneously, so it is a
    # cleaner stratum than GSE194315's (lane x subject) — nothing else to hold.
    stratum_cols=("orig.ident",),

    defaults=dict(
        # DELIBERATELY NO min_ids. PBMC's 0.5 came from PBMC's own IDS histogram
        # (background topping out ~0.20, programs ~0.5-0.9). Skin has different
        # cell counts, a different gene set and far more keratinocytes; the gap
        # sits somewhere else. Run once, histogram the off-diagonal, then set it.
        min_cells_frac=0.02,
        block_size=1024,
        n_pcs=20,
    ),
))


# ----------------------------------------------------------------------
# One-time conversion: 33 gene-x-cell CSVs + metadata  ->  one sparse .h5ad
# ----------------------------------------------------------------------
def build_h5ad(raw_dir, metadata_path, features_tsv_gz, out_h5ad,
               scratch_dir=None, verbose=True):
    """Turn the GEO archive into a single AnnData the pipeline can read.

    Done once, because the raw form is hostile to repeated use: 33 DENSE
    gene-x-cell CSVs, 33,538 rows each, Ensembl-ID rows, and a barcode suffix
    that differs from the metadata's. The result is cells x genes, sparse, gene
    symbols, with the metadata joined and the source study's filtering applied.

    Two memory notes, both load-bearing on a laptop:
      * columns are subset AT PARSE TIME (`usecols`) to the cells the metadata
        keeps, so the ~30% of droplets the source study dropped are never read
        into memory at all;
      * each file is converted to sparse before the next is opened, so peak
        memory is one sample, not the cohort.
    """
    import glob, gzip, os, re, tempfile, shutil, resource
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import anndata as ad

    meta = pd.read_csv(metadata_path, sep="\t", index_col=0, low_memory=False)
    meta.index = meta.index.astype(str)
    files = sorted(glob.glob(os.path.join(raw_dir, "*.csv.gz")))
    if not files:
        raise FileNotFoundError(f"no *.csv.gz under {raw_dir}")

    # TWO PASSES, because one is not survivable on a laptop. Holding all 33 sparse
    # blocks AND the stacked result needs roughly twice the final matrix (~2.2 GB
    # here) and gets the process OOM-killed with no traceback. Pass 1 spills each
    # sample to scratch; pass 2 allocates the final arrays exactly once, from the
    # now-known total non-zero count.
    # Spill to LOCAL scratch. The spill is ~2 GB; writing it into a network- or
    # bridge-mounted folder buffers in RAM and will get the process OOM-killed on a
    # small machine. `scratch_dir` should be a local disk, never the output folder.
    scratch = tempfile.mkdtemp(prefix="gse173706_", dir=scratch_dir)
    obs_keys, genes, parts = [], None, []
    try:
        for i, f in enumerate(files, 1):
            sample = re.sub(r"^GSM\d+_", "", os.path.basename(f)).replace(".csv.gz", "")
            sample = sample.replace("-", "_", 1)            # 'NS-AR001' -> 'NS_AR001'
            rows = meta[meta["orig.ident"] == sample]
            if rows.empty:
                raise KeyError(f"{sample} has no cells in the metadata")

            # metadata key '<bare>-<n>'  ->  this file's column name '<bare>-1'
            bare = pd.Index(rows.index).str.replace(r"-\d+$", "", regex=True)
            wanted = set(bare + "-1")
            # dtype must be given per COLUMN, not globally: a global dtype is applied
            # to the index too, and the index holds Ensembl IDs. Reading the header
            # first also lets us fail loudly if a metadata barcode is missing.
            with gzip.open(f, "rt") as fh:
                header = fh.readline().rstrip("\n").split(",")
            # The gene column has an EMPTY header in these files; pandas renames it
            # to 'Unnamed: 0', so usecols must ask for that name, not for "".
            index_name = header[0] or "Unnamed: 0"
            cols = header[1:]
            take = [c for c in cols if c in wanted]
            if len(take) != len(wanted):
                raise ValueError(f"{sample}: {len(wanted) - len(take)} metadata "
                                 f"barcodes absent from the matrix")
            # Read in GENE-ROW CHUNKS. Materializing a whole sample dense costs
            # genes x cells x 4 bytes (~500 MB for the largest here) and then again
            # for the float32 copy scipy needs, which is what killed the process on
            # PP_7802ED. A chunk is ~60 MB and the sparse accumulation is ~80 MB.
            chunks, index_parts, colnames = [], [], None
            for chunk in pd.read_csv(f, index_col=0, usecols=[index_name] + take,
                                     dtype={c: np.int32 for c in take},
                                     chunksize=4000):
                if colnames is None:
                    colnames = chunk.columns
                index_parts.append(chunk.index.astype(str))
                chunks.append(sp.csr_matrix(chunk.to_numpy(dtype=np.float32)))
                del chunk
            g_index = pd.Index(np.concatenate([p.to_numpy() for p in index_parts]))

            if genes is None:
                genes = g_index
            elif not g_index.equals(genes):
                raise ValueError(f"{sample}: gene order differs from the first file")

            # (genes x cells) sparse -> transpose is a free CSC->CSR relabel
            X = sp.vstack(chunks, format="csr").T.tocsr()        # cells x genes
            del chunks
            keys = pd.Index(colnames).str.replace(r"-\d+$", "", regex=True) + \
                "-" + str(rows.index[0]).rsplit("-", 1)[-1]
            part = os.path.join(scratch, f"{i:02d}.npz")
            np.savez(part, data=X.data, indices=X.indices, indptr=X.indptr,
                     shape=np.array(X.shape))
            parts.append(part)
            obs_keys.append(keys)
            if verbose:
                rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
                print(f"  [{i:>2}/{len(files)}] {sample:<12} {X.shape[0]:>5} cells, "
                      f"{X.nnz / 1e6:.1f}M non-zero, peak rss {rss:.2f} GB", flush=True)
            del X

        # pass 2 — one exact allocation
        nnz = n_cells = 0
        for part in parts:
            with np.load(part) as z:
                nnz += len(z["data"]); n_cells += int(z["shape"][0])
        n_genes = len(genes)
        if verbose:
            print(f"  [assemble] {n_cells} cells x {n_genes} genes, "
                  f"{nnz / 1e6:.0f}M non-zero -> "
                  f"{nnz * 8 / 1e9:.2f} GB", flush=True)

        data = np.empty(nnz, dtype=np.float32)
        indices = np.empty(nnz, dtype=np.int32)
        indptr = np.zeros(n_cells + 1, dtype=np.int64)
        at_nnz = at_row = 0
        for part in parts:
            with np.load(part) as z:
                d, ix, ip = z["data"], z["indices"], z["indptr"]
            k, r = len(d), len(ip) - 1
            data[at_nnz:at_nnz + k] = d
            indices[at_nnz:at_nnz + k] = ix
            indptr[at_row + 1:at_row + 1 + r] = ip[1:] + at_nnz
            at_nnz += k; at_row += r
            del d, ix, ip
        X = sp.csr_matrix((data, indices, indptr), shape=(n_cells, n_genes))
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    adata = ad.AnnData(X)
    adata.obs_names = pd.Index(np.concatenate([k.to_numpy() for k in obs_keys]))
    adata.var_names = genes

    missing = adata.obs_names.difference(meta.index)
    if len(missing):
        raise ValueError(f"{len(missing)} cells not found in the metadata — the "
                         f"barcode key rule is wrong; first: {list(missing[:3])}")
    for col in meta.columns:
        adata.obs[col] = meta.loc[adata.obs_names, col].values

    # Ensembl -> symbol, so skin and blood programs share one namespace.
    sym = ensembl_to_symbol(features_tsv_gz)
    hit = adata.var_names.map(lambda g: sym.get(g))
    n_unmapped = int(pd.isna(hit).sum())
    adata.var["ensembl_id"] = adata.var_names
    adata.var_names = pd.Index([s if isinstance(s, str) and s else e
                                for s, e in zip(hit, adata.var["ensembl_id"])])
    adata.var_names_make_unique()
    if verbose:
        print(f"  [symbols] {adata.n_vars - n_unmapped}/{adata.n_vars} Ensembl IDs "
              f"mapped to symbols; {n_unmapped} kept as ENSG")
        print(f"  [build] {adata.n_obs} cells x {adata.n_vars} genes, "
              f"{adata.X.nnz / (adata.n_obs * adata.n_vars):.1%} non-zero")

    adata.write(out_h5ad)
    if verbose:
        print(f"  [build] wrote {out_h5ad} "
              f"({os.path.getsize(out_h5ad) / 1e6:.0f} MB)")
    return adata


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Build the GSE173706 h5ad (run once).")
    ap.add_argument("--raw-dir", default="data/GSE173706")
    ap.add_argument("--metadata", default="data/GSE173706_psorSkinMeta.tsv")
    ap.add_argument("--features",
                    default="data/GSE194315_PBMC-01-07_processed_data_files/"
                            "PBMC-01-1.features.tsv.gz")
    ap.add_argument("--out", default="data/GSE173706_skin.h5ad")
    ap.add_argument("--scratch", default=None,
                    help="local disk for the ~2 GB spill; must NOT be a mounted "
                         "or network folder. Defaults to the system temp dir.")
    a = ap.parse_args()
    build_h5ad(a.raw_dir, a.metadata, a.features, a.out, scratch_dir=a.scratch)
