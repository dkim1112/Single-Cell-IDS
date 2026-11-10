"""Where things live. One place, so no script hardcodes a path.

  data/         INPUT. Raw GEO downloads. Read-only, gitignored.
  results/      OUTPUT. One folder per pipeline stage, one subfolder per run.
  deliverables/ OUTPUT meant for other people (workbooks sent to Dr. Patrick, etc.).
  docs/         Notes, reports, the abstract.

Override the locations with environment variables when the data live elsewhere:
  SC_IDS_DATA       folder holding GSE194315_* and GSE173706_*  (default: <repo>/data)
  SC_IDS_PBMC_DIR   folder holding the PBMC-0X-Y.{barcodes,features,matrix} files
"""
import os

ROOT = os.environ.get("SC_IDS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", ".."))

# ---- inputs
DATA = os.environ.get("SC_IDS_DATA") or os.path.join(ROOT, "data")
PBMC_DIR = os.environ.get("SC_IDS_PBMC_DIR") or os.path.join(
    DATA, "GSE194315_PBMC-01-07_processed_data_files")
PBMC_META = os.path.join(DATA, "GSE194315_CellMetadata-AS_TotalCiteseq_20220711.tsv")
SKIN_DIR = os.path.join(DATA, "GSE173706")
SKIN_META = os.path.join(DATA, "GSE173706_psorSkinMeta.tsv")

# ---- outputs, one folder per stage
RESULTS = os.path.join(ROOT, "results")
QC = os.path.join(RESULTS, "01_qc_labelling")       # QC, Harmony/UMAP, clustering, labels, NK check
IDS = os.path.join(RESULTS, "02_ids_programs")      # disease object, IDS matrices, programs per stratum
CATALOGUE = os.path.join(RESULTS, "03_catalogue")   # filtered, merged, annotated programs
DELIVERABLES = os.path.join(ROOT, "deliverables")

# ---- named runs inside the stages
DISEASE_RUN = os.path.join(IDS, "ids-disease")      # whole-tissue run: blood_qc.h5ad, skin_qc.h5ad, programs_*.csv
LINEAGE_RUNS = {"T+NK": os.path.join(IDS, "ids-lin-TNK"),
                "Myeloid": os.path.join(IDS, "ids-lin-MY"),
                "Keratinocyte": os.path.join(IDS, "ids-lin-KC")}
ATLAS_RUN = os.path.join(IDS, "ids-atlas")
SKIN_CACHE = os.path.join(IDS, "skin_raw.h5ad")
