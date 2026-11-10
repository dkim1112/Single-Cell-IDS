"""Coarse lineage axis shared by blood and skin: source cell-type label -> lineage.

Single source of truth. Scripts import LINEAGE from here; do not copy it."""

# Coarse lineage axis shared by the two tissues. Blood labels come from
# GSE194315's CellType, skin from GSE173706's fullsubtype; they use different
# vocabularies for the same lineages, which is precisely why the program
# framing is being used, but a common axis is still needed to SUMMARISE on.
LINEAGE = {
    # --- blood
    "B intermediate": "B", "B memory": "B", "B naive": "B", "Plasmablast": "B",
    "CD14 Mono": "Myeloid", "CD16 Mono": "Myeloid", "ASDC": "Myeloid",
    "cDC1": "Myeloid", "cDC2": "Myeloid", "pDC": "pDC",
    "CD4 CTL": "T", "CD4 Naive": "T", "CD4 Proliferating": "T", "CD4 TCM": "T",
    "CD4 TEM": "T", "CD8 Naive": "T", "CD8 Proliferating": "T", "CD8 TCM": "T",
    "CD8 TEM": "T", "MAIT": "T", "Treg": "T", "dnT": "T", "gdT": "T",
    "NK": "NK", "NK Proliferating": "NK", "NK_CD56bright": "NK", "ILC": "NK",
    "Platelet": "Platelet", "Eryth": "Erythroid", "HSPC": "Other",
    # --- skin
    "Basal Keratinocyte": "Keratinocyte", "Spinous Keratinocyte": "Keratinocyte",
    "Supraspinous Keratinocyte": "Keratinocyte",
    "COL11A1+ Fibroblast": "Fibroblast", "RAMP1+ Fibroblast": "Fibroblast",
    "SFRP2+ Fibroblast": "Fibroblast", "SFRP4+ Fibroblast": "Fibroblast",
    "Endothelial": "Endothelial", "Smooth muscle": "Smooth muscle",
    "Melanocyte": "Melanocyte", "Mast cell": "Mast", "Nerve cell": "Nerve",
    "Eccrine gland": "Eccrine gland",
    "Cytotoxic T cell": "T", "IFN T cell": "T", "Naive T cell": "T",
    "Stressed T cell": "T", "Tc17 cell": "T", "yd T cell": "T",
    "LC": "Myeloid", "PVM": "Myeloid", "CD16+ DC": "Myeloid",
    "cDC2A": "Myeloid", "cDC2B": "Myeloid",
}
