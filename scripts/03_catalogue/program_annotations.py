"""
program_annotations.py — hand-curated annotation for the programs.

Each entry matches a distinct program when ALL listed genes are in its gene set (first match
wins, most specific first). Nothing here is auto-generated. 'checked' = the statement in
`implication` was read in the cited source; 'gene function only' = standard gene function, no
psoriasis source was read; 'not checked' = no annotation attempted.
"""
SRC = {
 "amp":   ("Antimicrobial peptides in psoriasis review (PMC3527011)", "https://pmc.ncbi.nlm.nih.gov/articles/PMC3527011"),
 "krt":   ("Zhang X et al., Cells 2019;8:807", "https://www.mdpi.com/2073-4409/8/8/807"),
 "il36":  ("IL36G in psoriatic skin (PMC9510771)", "https://pmc.ncbi.nlm.nih.gov/articles/PMC9510771"),
 "chemo": ("The Role of Chemokines in Psoriasis, Medicina 2021;57:754", "https://www.mdpi.com/1648-9144/57/8/754"),
 "apc":   ("Antigen-Presenting Cells in Psoriasis, Life 2022;12:234", "https://pmc.ncbi.nlm.nih.gov/articles/PMC8880330"),
 "endo":  ("ICAM-1, E-selectin, VCAM-1 on endothelium in psoriasis, Acta Derm Venereol", "https://medicaljournalssweden.se/actadv/article/view/5399"),
 "ifn_psa": ("Type-I interferon pathway in PsA blood (PMC10731026)", "https://pmc.ncbi.nlm.nih.gov/articles/PMC10731026"),
 "cd8_psa": ("Durham LE et al., Arthritis Rheumatol 2026;78:152 (via ReachMD summary)", "https://reachmd.com/programs/living-rheum/shared-cd8-t-cell-clones-link-skin-and-joint-inflammation-in-psoriatic-arthritis/56345/"),
 "diss":  ("O'Flanagan 2019 Genome Biol; van den Brink 2017 Nat Methods", "https://pubmed.ncbi.nlm.nih.gov/31623682/"),
}

# (required genes, short name, implication, source key, evidence)
ANN = [
 ({"S100A7", "S100A8", "S100A9"}, "Keratinocyte antimicrobial / alarmin response",
  "S100A7/8/9 and elafin (PI3) are abundant in psoriatic lesions and have antimicrobial activity. Marks the inflamed keratinocyte state. "
  "It also reads as 'shared with blood' because S100A8/A9 are normal monocyte genes, so the blood signal is monocytes, not the same biology. Rises in non-keratinocyte cell types of lesional skin may be keratinocyte RNA leaking into other cells; ambient RNA was not ruled out.",
  "amp", "checked"),
 ({"KRT6A", "KRT16"}, "Hyperproliferation keratins (KRT6/16)",
  "KRT6/16/17 are high in psoriatic suprabasal keratinocytes and act as early barrier alarmins; T-cell cytokines (IFN-g, IL-17A, IL-22) help keep them on. "
  "The proliferation effect is shown for KRT17, which is not in this set.",
  "krt", "checked"),
 ({"IL36G"}, "IL-36 / keratinocyte stress-response set",
  "IL36G is up in psoriatic skin, in suprabasal keratinocytes, and tracks PASI. SPRR2A/D, APOBEC3A and HBEGF: no source read.",
  "il36", "checked (IL36G only)"),
 ({"CXCL9", "CXCL10", "CXCL11"}, "CXCR3 chemokines (T-cell recruiting)",
  "CXCL10 is in keratinocytes and infiltrate of active plaques and falls with treatment; CXCL11 comes mainly from basal keratinocytes and recruits T cells. "
  "CXCL9 is only named as inflammatory.",
  "chemo", "checked"),
 ({"IFIT1", "IFIT3", "ISG15"}, "Interferon-stimulated genes (skin)",
  "IL36G co-occurs with ISG15, MX1, IFI44L in lesional suprabasal keratinocytes but not in normal skin. "
  "In PsA blood a published IFN score (IFIT1, MX1, IFI44) was lower than in controls, so ISG here is not expected to be raised in blood.",
  "il36", "checked"),
 ({"IFIT1", "IFIT2", "IFIT3"}, "Interferon-stimulated genes (core)",
  "Type-I interferon response genes. In PsA blood a published IFN score (IFIT1, MX1, IFI44) was lower than in controls.",
  "ifn_psa", "checked"),
 ({"IFI6", "IFITM3", "LY6E"}, "Interferon-stimulated genes (small)",
  "Interferon response genes; same context as the other ISG sets.", "ifn_psa", "gene function only"),
 ({"CCL3", "CCL3L3", "CCL4", "CCL4L2", "CCL5", "GZMA"}, "Cytotoxic lymphocyte program with MIP-1 chemokines",
  "CCL3/CCL4 blood levels are raised in psoriasis; CCL4 pulls Th1 cells, DCs and monocytes into lesions, CCL3 favours CD8 T cells. "
  "In PsA, CD8 clones shared between skin and joint show GZMB, GNLY, PRF1.",
  "chemo", "checked"),
 ({"CCL3", "CCL3L3", "CCL4", "CCL4L2"}, "MIP-1 chemokines (CCL3/CCL4)",
  "CCL3/CCL4 blood levels are raised in psoriasis; CCL4 pulls Th1 cells, DCs and monocytes into lesions, CCL3 favours CD8 T cells.",
  "chemo", "checked"),
 ({"CXCL8", "IL1B", "CXCL2"}, "Neutrophil / monocyte inflammatory chemokines",
  "CXCL8 with CXCL1 attracts neutrophils to the epidermis. IL1B, CXCL2/3: no source read.",
  "chemo", "checked (CXCL8 only)"),
 ({"CXCL8", "IL1B", "IL17RA"}, "Inflammatory monocyte program (with IL-17 receptor A)",
  "CXCL8 attracts neutrophils to the epidermis (source). IL17RA is the receptor for IL-17A; no source read on this set as a whole.",
  "chemo", "checked (CXCL8 only)"),
 ({"GNLY", "GZMB", "NKG7", "CCL5"}, "Cytotoxic effector core (GNLY, GZMB, NKG7, CCL5)",
  "NK / CD8 killing machinery. CCL5 is also made by lesional keratinocytes. In PsA, CD8 clones shared between skin and joint show GZMB, GNLY, PRF1.",
  "cd8_psa", "checked"),
 ({"GNLY", "GZMA", "GZMH", "CD3D"}, "Cytotoxic T-cell program (with T-cell receptor genes)",
  "NK / CD8 killing machinery. In PsA, CD8 clones shared between skin and joint show GZMB, GNLY, PRF1.", "cd8_psa", "checked"),
 ({"GZMA", "GZMB", "GNLY", "FGFBP2"}, "Cytotoxic NK / CD8 effector program (large)",
  "NK / CD8 killing machinery. In PsA, CD8 clones shared between skin and joint show GZMB, GNLY, PRF1.", "cd8_psa", "checked"),
 ({"GNLY", "KLRD1", "TYROBP"}, "NK-cell cytotoxic set",
  "NK / CD8 killing machinery. In PsA, CD8 clones shared between skin and joint show GZMB, GNLY, PRF1.", "cd8_psa", "checked"),
 ({"CCL5", "GZMA", "NKG7"}, "Cytotoxic core (small)",
  "NK / CD8 killing machinery. In PsA, CD8 clones shared between skin and joint show GZMB, GNLY, PRF1.", "cd8_psa", "checked"),
 ({"GZMK", "NKG7", "HAVCR2"}, "GZMK+ cytotoxic / memory T-cell set",
  "Granzyme K / NKG7 T-cell state. No psoriasis source read.", "", "gene function only"),
 ({"HLA-DRA", "CD74", "HLA-DPA1"}, "Antigen-presenting cell program (MHC class II)",
  "Marks dendritic cells, monocytes, Langerhans cells and B cells. In psoriatic lesions the review links DCs, Langerhans cells and macrophages to IL-23 and "
  "to presenting self-antigens (LL37, ADAMTSL5); it does not describe MHC-II levels, so this set's change is not covered by a source.",
  "apc", "checked (cell role only)"),
 ({"CD79A", "CD79B", "MS4A1"}, "B-cell program (with MHC class II)",
  "B-cell identity and antigen presentation. No psoriasis source read.", "", "gene function only"),
 ({"SELE", "ICAM1"}, "Endothelial activation (E-selectin, ICAM-1)",
  "ICAM-1 and E-selectin are up on endothelium in lesional AND non-lesional psoriatic skin.", "endo", "checked"),
 ({"FOS", "JUN"}, "Immediate-early / AP-1 stress response",
  "Dissociation of tissue can itself switch these genes on (FOS, JUN and heat-shock genes). Treat any change here as a possible handling effect until "
  "checked against the handling of each sample.", "diss", "checked (artefact risk)"),
 ({"CD83", "JUNB", "TNFAIP3"}, "NF-kB / immediate-early activation",
  "Activation-response genes (CD83, JUNB, NFKBIA, TNFAIP3). Same handling-effect risk as the AP-1 sets.", "diss", "gene function only"),
 ({"MT1E", "MT1X", "MT2A"}, "Metallothioneins (metal / oxidative stress)",
  "Metal-binding stress-response genes. The one skin source read covers atopic dermatitis, not psoriasis.", "", "no psoriasis source found"),
 ({"HIST1H4C", "NUSAP1", "STMN1"}, "Cell-cycle / proliferation",
  "S/G2-M phase genes (histones, NUSAP1, STMN1, tubulin). Gene function only.", "", "gene function only"),
 ({"DSC1", "DSG1"}, "Desmosome genes of differentiating keratinocytes",
  "Desmosomal adhesion genes of upper epidermis. In our data DSC1 is detected in fewer lesional keratinocytes (61% to 32%). No psoriasis source read.",
  "", "our data only"),
 ({"FLG", "CDSN"}, "Cornified envelope / terminal differentiation",
  "Skin-barrier genes of the upper epidermis (filaggrin, involucrin, corneodesmosin). Lower in lesional keratinocytes in our data. No source read.",
  "", "our data only"),
 ({"IVL", "CALML5"}, "Cornified envelope / terminal differentiation (variant)",
  "Upper-epidermis differentiation genes (involucrin, calmodulin-like 5, DSC1, KRT2). No source read.", "", "our data only"),
 ({"CLDN1", "FOSL1", "SFN"}, "Keratinocyte stress / differentiation set (small)", "No source read.", "", "not checked"),
 ({"PLAU", "PLAUR"}, "Plasminogen-activator system (PLAU / PLAUR)", "Tissue-remodelling genes. No source read.", "", "not checked"),
 ({"CD8A", "CD8B"}, "CD8 T-cell identity", "CD8 co-receptor genes.", "", "gene function only"),
 ({"ACTB", "TMSB10", "IL32"}, "Cytoskeletal / activation genes (small)", "ACTB, TMSB10, S100A4, LGALS1 are broadly expressed; no specific function claimed.", "", "not checked"),
]

def annotate_program(genes):
    g = set(x.strip() for x in genes.split(","))
    for need, name, impl, src, ev in ANN:
        if need <= g:
            s = SRC.get(src, ("", ""))
            return name, impl, s[0], s[1], ev
    return "", "", "", "", "not annotated"
