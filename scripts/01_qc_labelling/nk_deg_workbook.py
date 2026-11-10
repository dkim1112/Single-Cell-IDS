"""
nk_deg_workbook.py — package the NK-cluster DEG tables as one shareable workbook.

Reads the CSVs written by nk_deg_export.py / rank_markers.py / nk_percell.py and
writes NK_cluster_DEG_for_review.xlsx. No new statistics are computed here.
"""
import os, sys
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from sc_ids import paths                                                      # noqa: E402

R = os.path.join(paths.QC, "markers-theta4") + os.sep
OUT = os.path.join(paths.DELIVERABLES, "NK_cluster_DEG_for_review.xlsx")
ARIAL = "Arial"
HEAD = PatternFill("solid", fgColor="E8E8E8")
WARN = PatternFill("solid", fgColor="FDE2E2")
OKF = PatternFill("solid", fgColor="EAF3EA")
THIN = Side(style="thin", color="C8C8C8")
BOX = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)


def sheet(wb, name, df, widths, notes=None, flag_col=None):
    ws = wb.create_sheet(name)
    r0 = 1
    if notes:
        for i, n in enumerate(notes):
            ws.cell(row=1 + i, column=1, value=n).font = Font(name=ARIAL, size=10, italic=True, color="666666")
        r0 = 1 + len(notes) + 1
    for j, col in enumerate(df.columns, start=1):
        c = ws.cell(row=r0, column=j, value=str(col))
        c.font = Font(name=ARIAL, size=10, bold=True)
        c.fill = HEAD
        c.border = BOX
        c.alignment = Alignment(wrap_text=True, vertical="center")
    for i, row in enumerate(df.itertuples(index=False), start=r0 + 1):
        for j, v in enumerate(row, start=1):
            c = ws.cell(row=i, column=j, value=(None if pd.isna(v) else v))
            c.font = Font(name=ARIAL, size=10)
            c.border = BOX
            col = df.columns[j - 1]
            if isinstance(v, float):
                c.number_format = "0.00E+00" if ("p_value" in col or "p_adj" in col) else "0.00"
            if flag_col and col == flag_col:
                c.fill = WARN if str(v) == "NO" else OKF
                c.font = Font(name=ARIAL, size=10, bold=(str(v) == "NO"))
    for j, w in enumerate(widths, start=1):
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.freeze_panes = ws.cell(row=r0 + 1, column=1)
    return ws


def main():
    up = pd.read_csv(R + "NK_cluster0_upregulated.csv")
    sv = pd.read_csv(R + "NK_cluster0_skin_vs_blood.csv")
    panel = pd.read_csv(R + "panel_t_theta4_by_cluster_tissue.csv")
    perc = pd.read_csv(R + "nk_percell_cluster0.csv")

    wb = Workbook()
    wb.remove(wb.active)

    # ---------------- Read me
    ws = wb.create_sheet("Read me")
    lines = [
        ("NK / gamma-delta T cluster — differentially expressed genes for review", 14, True, None),
        ("", 10, False, None),
        ("What this is", 11, True, None),
        ("Sub-cluster 0 of the T/NK lineage in a joint blood + skin single-cell object. We believe it is NK cells.", 10, False, None),
        ("The open question is whether the skin cells inside it are NK or gamma-delta T, which we cannot settle on our own.", 10, False, None),
        ("", 10, False, None),
        ("Cells in the cluster", 11, True, None),
        ("blood (PBMC): 1,295 cells     skin: 159 cells     total: 1,454", 10, False, None),
        ("Source datasets: PBMC lane PBMC-01-1 (lab) and skin GSE173706. Both passed the lab QC pipeline.", 10, False, None),
        ("", 10, False, None),
        ("How the genes were ranked", 11, True, None),
        ("scanpy rank_genes_groups, method='wilcoxon' — the same call as the lab pipeline's run_Ranking.py.", 10, False, None),
        ("Cluster 0 against all other T/NK sub-clusters (1, 2, 3) pooled, both tissues together.", 10, False, None),
        ("Kept at the lab config thresholds: adjusted p < 0.01 and log2 fold change > 1.", 10, False, None),
        ("Clustering ran on the Harmony-corrected embedding at theta = 4, Leiden resolution 0.40.", 10, False, None),
        ("Expression is CPM-normalised and log1p-transformed, so log2FC is on that scale.", 10, False, None),
        ("", 10, False, None),
        ("Sheets", 11, True, None),
        ("Upregulated genes — the main table: 547 genes, with log2FC, p-value, adjusted p-value, and the", 10, False, None),
        ("    share of cells each gene is detected in, split into the blood half and the skin half of the cluster.", 10, False, None),
        ("Removed by gene filter — the 72 of those 547 that a 1%-of-cells gene filter would delete.", 10, False, None),
        ("Marker panel — detection rates for the panel we were given (NK: FCGR3A, NKG7, GZMB, PRF1 present,", 10, False, None),
        ("    CD3D absent), for every sub-cluster and tissue, so cluster 0 can be compared against the plain T clusters.", 10, False, None),
        ("Per-cell calls — the same test applied one cell at a time, with two plain-T clusters as a reference.", 10, False, None),
        ("Skin vs blood in cluster — differential expression between the two tissue halves of cluster 0 itself.", 10, False, None),
        ("", 10, False, None),
        ("On the gene filter", 11, True, None),
        ("The filter keeps a gene only if it is detected in at least 1% of cells, applied per tissue then intersected:", 10, False, None),
        ("10,126 genes kept of 18,250. Of the 547 genes up in this cluster, 475 survive and 72 are deleted.", 10, False, None),
        ("The 72 include KLRD1, PRF1, KLRF1, FCGR3A, NCAM1, KLRC1, KLRC2, NCR3, TBX21, TRDC and TRGC1 —", 10, False, None),
        ("that is, the cluster would still be recognisable as cytotoxic (GNLY, NKG7, GZMB, GZMA, TYROBP all survive),", 10, False, None),
        ("but the genes that separate NK from gamma-delta T specifically are the ones removed.", 10, False, None),
        ("Filter status is marked per gene in the main table, so the effect can be judged gene by gene.", 10, False, None),
        ("", 10, False, None),
        ("What we are asking", 11, True, None),
        ("Given these genes, is cluster 0 NK, gamma-delta T, or a mixture — and does the answer differ between", 10, False, None),
        ("the blood half and the skin half?", 10, False, None),
    ]
    for i, (t, sz, bold, _) in enumerate(lines, start=1):
        c = ws.cell(row=i, column=1, value=t)
        c.font = Font(name=ARIAL, size=sz, bold=bold)
    ws.column_dimensions["A"].width = 118

    # ---------------- main DEG table
    sheet(wb, "Upregulated genes", up,
          [7, 14, 10, 13, 13, 15, 16, 16, 15],
          notes=["T/NK sub-cluster 0 vs sub-clusters 1-3 pooled. Wilcoxon (scanpy rank_genes_groups), "
                 "adjusted p < 0.01 and log2FC > 1. Sorted by test statistic.",
                 "Red in the last column = the gene would be deleted by a 1%-of-cells gene filter."],
          flag_col="survives 1% gene filter")

    lost = up[up["survives 1% gene filter"] == "NO"].drop(columns=["survives 1% gene filter"])
    sheet(wb, "Removed by gene filter", lost, [7, 14, 10, 13, 13, 15, 16, 16],
          notes=["The 72 up-regulated genes a 1%-of-cells filter (applied per tissue, then intersected) would delete.",
                 "These are the genes whose loss matters for telling NK from gamma-delta T."])

    sheet(wb, "Marker panel", panel, [11, 9, 8] + [11] * (len(panel.columns) - 3),
          notes=["Share of cells in which each gene is detected (at least one read), per sub-cluster and tissue.",
                 "Clusters 1, 2 and 3 are T-cell clusters and serve as the comparison for cluster 0."])

    sheet(wb, "Per-cell calls", perc, [26, 9, 11, 11, 10, 11, 20],
          notes=["Counts of cells, not percentages. A cell is NK+ if at least 2 of FCGR3A, NKG7, GZMB, PRF1 are "
                 "detected, and T+ if any of CD3D, CD3E, TRAC are.",
                 "The two reference rows show what a single cell type looks like under this test."])

    sv2 = sv.rename(columns={"group": "up in tissue", "names": "gene", "logfoldchanges": "log2FC",
                             "pvals": "p_value", "pvals_adj": "p_adj", "scores": "wilcoxon_score",
                             "survives_1pct_filter": "survives 1% gene filter"})
    sv2 = sv2[["up in tissue", "gene", "log2FC", "p_value", "p_adj", "wilcoxon_score",
               "survives 1% gene filter"]]
    sheet(wb, "Skin vs blood in cluster", sv2, [13, 14, 10, 13, 13, 15, 15],
          notes=["Differential expression between the skin cells and the blood cells INSIDE cluster 0.",
                 "Same thresholds. This is what makes the two halves look like different populations."],
          flag_col="survives 1% gene filter")

    # a recalculating count so the workbook is not purely static
    ws2 = wb["Upregulated genes"]
    n = len(up) + 3
    ws2.cell(row=n + 2, column=1, value="genes listed:").font = Font(name=ARIAL, size=10, bold=True)
    ws2.cell(row=n + 2, column=2, value=f"=COUNTA(B4:B{n - 1})").font = Font(name=ARIAL, size=10)
    ws2.cell(row=n + 3, column=1, value="removed by filter:").font = Font(name=ARIAL, size=10, bold=True)
    ws2.cell(row=n + 3, column=2, value=f'=COUNTIF(I4:I{n - 1},"NO")').font = Font(name=ARIAL, size=10)

    os.makedirs(paths.DELIVERABLES, exist_ok=True)
    wb.save(OUT)
    print(f"[out] {OUT}  sheets: {wb.sheetnames}")
    print(f"      {len(up)} up genes, {len(lost)} removed by filter, {len(sv2)} skin-vs-blood genes")


if __name__ == "__main__":
    main()
