"""
GSE194315 — PBMC CITE-seq, psoriatic arthritis and healthy controls.
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194315

Everything here was established empirically on this dataset; see
docs/agenda1_and_3_results.md and the handoff for the measurements behind it.

Structure worth remembering when reading the field values below:
  * Each 10x lane (PBMC-01-1, PBMC-01-2, ...) pools NINE genetically distinct
    subjects, demultiplexed with demuxlet, and ~25% of droplets are doublets.
    A lane is therefore NOT a sample.
  * PBMC-01-1 and PBMC-01-2 carry the SAME nine subjects on two lanes — paired
    technical replicates, which is why recurrence across the 18 (lane x subject)
    strata must be scored hierarchically rather than flat (Agenda 2, section 2).
  * `CellName` is '<Sample>_<barcode>' with no '-1' suffix; scanpy produces
    '<barcode>-1'. Reconstructing the key matched 27,299/27,299 cells.
  * `IncludedInStudy` is the source study's own QC verdict and is strictly
    stronger than a singlet filter: 17,858 kept vs 20,517 for singlets alone on
    PBMC-01-1. `CellType` is annotated ONLY for included cells.
"""
from . import Dataset, register, strip_gem_suffix

GSE194315 = register(Dataset(
    name="gse194315",
    description="PBMC CITE-seq, 28 samples, demuxlet-multiplexed 9 subjects per lane",

    barcode_col="CellName",
    make_key=strip_gem_suffix,
    qc_gate=("IncludedInStudy", {"TRUE", "T", "1"}),
    singlet_gate=("DemuxletDropletType", "SNG"),

    celltype_col="CellType",
    celltype_map_file="gse194315_celltype_mapping.txt",

    # One lane removes batch (one file = one machine run); one subject removes
    # donor, which was measured to be the LARGER effect (same donor across lanes
    # r=0.43 vs different donors within a lane r=0.23). All cell types are kept,
    # so cell-identity programs stay visible.
    stratum_cols=("Sample", "Subject"),

    defaults=dict(
        # Starting points, not transferable settings. min_ids came from this
        # dataset's IDS histogram: background bulk tops out ~0.14-0.20, real
        # programs sit ~0.5-0.9, so 0.5 lands in the gap. Re-derive per dataset.
        min_cells_frac=0.02,
        min_ids=0.5,
        block_size=1024,
        n_pcs=20,
    ),
))
