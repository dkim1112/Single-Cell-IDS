"""
datasets/ — everything that is true of ONE dataset, kept out of the generic code.

WHY THIS EXISTS. The pipeline was written against GSE194315 (PBMC) and absorbed
its particulars: barcodes are matched by rebuilding '<Sample>_<barcode>', the QC
gate is `IncludedInStudy == TRUE`, doublets are flagged by `DemuxletDropletType`,
cell types live in `CellType`, and a batch-free stratum is (Sample x Subject)
because each 10x lane pools nine demuxlet-multiplexed donors. Not one of those
sentences will be true of a skin dataset.

So a `Dataset` here is a small description of those particulars. The generic code
(pipeline, programs, celltypes, viz) takes them as arguments and knows nothing
about any specific study. Adding a second dataset means adding ONE file here, not
editing the core.

    from sc_ids.datasets import get
    ds = get("gse194315")

What a Dataset does NOT contain: thresholds you have to justify. `defaults` holds
starting values that worked on that dataset, not settings that transfer. The gene
floor and `min_ids` must be re-derived from the new data's own IDS histogram —
the background level depends on cell count, gene count and tissue.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Callable

RESOURCES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                         "resources")


def resource_path(filename):
    """Absolute path to a file shipped under sc_ids/resources/."""
    return os.path.join(RESOURCES, filename)


def strip_gem_suffix(obs_names, sample=None):
    """Default barcode key: drop 10x's '-1' suffix, optionally prefix the sample.

    Cell Ranger writes 'AAACCCAAGACCATAA-1'. Most per-cell metadata files key on
    the bare barcode, some on '<sample>_<barcode>'. Override on the Dataset when a
    study does something else.
    """
    import pandas as pd
    bare = pd.Index(obs_names).astype(str).str.replace(r"-\d+$", "", regex=True)
    return pd.Index(sample + "_" + bare) if sample else bare


@dataclass
class Dataset:
    """A dataset's particulars. Every field is consumed by generic code."""

    name: str
    description: str = ""

    # --- metadata join -------------------------------------------------
    barcode_col: str = "CellName"
    make_key: Callable = strip_gem_suffix
    # (column, {values meaning KEEP}). The source study's own per-cell QC verdict
    # when it publishes one — preferable to hand-tuned cutoffs. None if it doesn't.
    qc_gate: tuple | None = None
    # (column, value) marking a droplet as a genuine singlet, for multiplexed lanes.
    singlet_gate: tuple | None = None

    # --- cell types ----------------------------------------------------
    celltype_col: str = "CellType"
    celltype_map_file: str | None = None      # filename under resources/
    # Some studies publish BOTH levels of annotation, in which case the mapping is
    # read off the data instead of a hand-written file. Set this to the column
    # holding the main type; celltype_col holds the granular one.
    celltype_main_col: str | None = None

    # --- stratification (Agenda 2) -------------------------------------
    # Columns whose combination defines one batch-free, donor-free stratum.
    # Order matters: outermost first (file/lane, then donor).
    stratum_cols: tuple = ()

    # --- starting parameters, NOT transferable settings ----------------
    defaults: dict = field(default_factory=dict)

    def celltype_map_path(self):
        if self.celltype_map_file is None:
            raise ValueError(
                f"dataset '{self.name}' has no cell-type mapping file. Obtain one "
                f"from the lab (granular<TAB>main, one line per label) and put it "
                f"in sc_ids/resources/, then set celltype_map_file."
            )
        return resource_path(self.celltype_map_file)


_REGISTRY = {}


def register(ds):
    _REGISTRY[ds.name] = ds
    return ds


def get(name):
    if name not in _REGISTRY:
        raise KeyError(f"unknown dataset '{name}'; available: {sorted(_REGISTRY)}")
    return _REGISTRY[name]


def available():
    return sorted(_REGISTRY)


from . import gse194315        # noqa: E402,F401  (registers on import)
from . import psoriasis_skin   # noqa: E402,F401
