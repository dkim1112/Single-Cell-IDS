"""
catalogue_workbook.py — workbook for Matthew's question: which gene programs are shared or
specific across tissue, cell type and disease state, and what do they mean.
One row per distinct program (not per detection).

Inputs: distinct_programs.csv, detections_status.csv, selection_chain.json (from distinct_programs.py),
disease_contrasts*.csv, program_annotations.py
"""
import os, sys, json
import numpy as np, pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from sc_ids import paths                                                      # noqa: E402
import program_annotations as PA                                              # noqa: E402

OUT = paths.CATALOGUE
ARIAL = "Arial"
HEAD = PatternFill("solid", fgColor="E8E8E8")
T1 = PatternFill("solid", fgColor="E3F1E3")
T2 = PatternFill("solid", fgColor="FFF4CC")
FLAG = PatternFill("solid", fgColor="FDE2E2")
THIN = Side(style="thin", color="C8C8C8")
BOX = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)


def shorten(txt, n=5):
    if not isinstance(txt, str) or not txt:
        return ""
    parts = txt.split("; ")
    return "; ".join(parts[:n]) + (f"; (+{len(parts) - n} more)" if len(parts) > n else "")


def write_table(ws, df, r0, widths, wrap=(), fills=None):
    for j, col in enumerate(df.columns, 1):
        c = ws.cell(row=r0, column=j, value=str(col))
        c.font = Font(name=ARIAL, size=10, bold=True); c.fill = HEAD; c.border = BOX
        c.alignment = Alignment(wrap_text=True, vertical="center")
    for i, row in enumerate(df.itertuples(index=False), r0 + 1):
        for j, v in enumerate(row, 1):
            v = None if (isinstance(v, float) and np.isnan(v)) else v
            c = ws.cell(row=i, column=j, value=v)
            c.font = Font(name=ARIAL, size=10); c.border = BOX
            c.alignment = Alignment(vertical="top", wrap_text=(df.columns[j - 1] in wrap))
            if fills and df.columns[j - 1] in fills and v in fills[df.columns[j - 1]]:
                c.fill = fills[df.columns[j - 1]][v]
    for j, col in enumerate(df.columns, 1):
        ws.column_dimensions[get_column_letter(j)].width = widths.get(col, 14)
    ws.freeze_panes = ws.cell(row=r0 + 1, column=3)
    ws.auto_filter.ref = f"A{r0}:{get_column_letter(len(df.columns))}{r0 + len(df)}"


def main():
    det = pd.read_csv(f"{OUT}/detections_status.csv")
    P0 = pd.read_csv(f"{OUT}/distinct_programs.csv").fillna("")
    con = pd.read_csv(f"{OUT}/disease_contrasts.csv")
    summ = pd.read_csv(f"{OUT}/disease_contrasts_summary.csv")
    ch = json.load(open(f"{OUT}/selection_chain.json"))

    ann = P0.genes.map(PA.annotate_program)
    P0["name"] = [a[0] for a in ann]; P0["implication"] = [a[1] for a in ann]
    P0["source"] = [a[2] for a in ann]; P0["source_url"] = [a[3] for a in ann]; P0["evidence"] = [a[4] for a in ann]
    P0["name"] = np.where(P0.name == "", "Unnamed (" + P0.genes.str.split(", ").str[:3].str.join(", ") + "…)", P0.name)
    P0.loc[P0.name.str.contains("identity"), "category"] = "cell-identity"
    P0["changes_skin"] = P0.skin_disease_change != ""
    note = " Our donor-level test: no lesional-vs-non-lesional change above 0.2 in skin."
    P0["implication"] = np.where(P0.evidence.str.startswith("checked") & ~P0.changes_skin & ~P0.category.str.startswith("stress"),
                                 P0.implication + note, P0.implication)

    def rel(r):
        if r.category.startswith("stress"): return "check handling artefact first"
        lit = r.evidence in ("checked", "checked (IL36G only)", "checked (CXCL8 only)")
        if r.changes_skin and lit: return "1"
        if r.changes_skin or (lit and r.tissue_call == "both"): return "2"
        return "3"
    P0["relevance"] = P0.apply(rel, axis=1)
    def skin_pattern(txt):
        """Direction of the donor-level skin changes (adjusted p < 0.05, |change| >= 0.2), from skin_disease_change."""
        import re
        if not txt: return "no change"
        found = re.findall(r"([^;:]+): (NS vs PP|NS vs PN|PN vs PP) ([+-][0-9.]+)", txt)
        by = {}
        for lin, con_, d in found:
            by.setdefault(con_, []).append((lin.strip(), float(d.rstrip("."))))
        def word(con_, name):
            v = by.get(con_)
            if not v: return None
            ups = sorted({l for l, d in v if d > 0}); dns = sorted({l for l, d in v if d < 0})
            if ups and not dns: return f"up in {name}"
            if dns and not ups: return f"down in {name}"
            return f"mixed in {name} (up: {', '.join(ups)}; down: {', '.join(dns)})"
        parts = [w for w in (word("NS vs PP", "lesional"), word("NS vs PN", "non-lesional"), word("PN vs PP", "lesional vs non-lesional")) if w]
        return "; ".join(parts) if parts else "no change"
    P0["pattern_skin"] = P0.skin_disease_change.map(skin_pattern)
    P0["disease_change"] = P0.apply(lambda r: f"Skin: {shorten(r.skin_disease_change, 4) if r.skin_disease_change else 'none above 0.2 score units'}. Blood: none detected.", axis=1)
    P0["spread"] = P0.apply(lambda r: f"{r.celltype_call}: {r.lineages}", axis=1)
    for k in ("tau_tissue", "tau_celltype", "tau_disease_skin", "tau_disease_blood"):
        P0[k] = pd.to_numeric(P0[k], errors="coerce")
    def spec_txt(r):
        f2 = lambda v: "n/a" if v != v else f"{v:.2f}"
        return (f"tissue {f2(r.tau_tissue)} | cell type {f2(r.tau_celltype)} | skin disease {f2(r.tau_disease_skin)} "
                f"| blood disease {f2(r.tau_disease_blood)}")
    P0["specificity"] = P0.apply(spec_txt, axis=1)

    # top activity score: highest subtype mean score (>= 50 cells) over the program's kept detections
    act = pd.read_csv(f"{OUT}/program_activity_long.csv.gz")
    act = act[(act.axis == "subtype") & (act.n_cells >= 50)]
    def top_act(r):
        ids = [x.strip() for x in r.member_ids.replace(" | ", "|||").split("|||")]
        x = act[act.program_id.isin(ids)]
        if x.empty: return pd.Series({"top_score": np.nan, "top_where": ""})
        b = x.loc[x.mean_score.idxmax()]
        return pd.Series({"top_score": round(float(b.mean_score), 2), "top_where": f"{b.group} ({b.tissue})"})
    P0[["top_score", "top_where"]] = P0.apply(top_act, axis=1)

    ident = P0[P0.category == "cell-identity"]
    prog = P0[P0.category != "cell-identity"].copy()
    order = {"1": 0, "2": 1, "3": 2, "check handling artefact first": 3}
    prog["_o"] = prog.relevance.map(order)
    prog["score_rank"] = prog.top_score.rank(ascending=False, method="min").astype(int)
    prog = prog.sort_values(["_o", "top_score"], ascending=[True, False])

    P = pd.DataFrame({
        "ID": prog.program_no, "Disease relevance (1 = highest)": prog.relevance,
        "Score rank (1 = highest of 40)": prog.score_rank, "Top activity score": prog.top_score, "Top activity: cell type (tissue)": prog.top_where,
        "Program": prog.name, "Genes": prog.genes,
        "Tissue": prog.tissue_call, "Cell types active, blood": prog.celltypes_blood.map(shorten),
        "Cell types active, skin": prog.celltypes_skin.map(shorten), "Cell-type spread": prog.spread,
        "Specificity (0 = equal everywhere, 1 = one group only)": prog.specificity,
        "Pattern in skin": prog.pattern_skin,
        "Skin: mean score by state (strongest cell type)": prog.disease_means_skin,
        "Blood: mean score by state (strongest cell type)": prog.disease_means_blood,
        "Changes with disease state?": prog.disease_change,
        "What it does / implication": prog.implication, "Source": prog.source, "Link": prog.source_url, "Evidence": prog.evidence})
    I = pd.DataFrame({"ID": ident.program_no, "Cell-type signature": np.where(ident.label_rule == "", ident.name, ident.label_rule),
                      "Genes": ident.genes, "Tissue": ident.tissue_call, "Cell types active, blood": ident.celltypes_blood.map(shorten),
                      "Cell types active, skin": ident.celltypes_skin.map(shorten), "Skin change with disease": ident.disease_change})

    det["program_no"] = det.program_no.fillna("").replace("", "—")
    A = det[["program_id", "status", "program_no", "analysis", "stratum", "n_genes", "genes", "median_IDS",
             "frac_pairs_beat_all_shuffles", "cells_supporting", "technical_flag", "top_celltype", "tissue_pattern"]].rename(columns={
        "program_id": "detection", "status": "outcome of the gates", "program_no": "distinct program",
        "stratum": "detected in (tissue_state)", "frac_pairs_beat_all_shuffles": "share of gene pairs that beat all 60 shuffles",
        "top_celltype": "highest-scoring cell type", "tissue_pattern": "tissue pattern (score ≥0.2)"})
    pid2gp = dict(zip(det.program_id, det.program_no))

    sk = con[(con.tissue == "skin") & (con.p_adj < 0.05) & (con.delta.abs() >= 0.2)].copy()
    sk["program_no"] = sk.program_id.map(pid2gp)
    sk = sk[sk.program_no.notna() & (sk.program_no != "—")]
    sk = sk.sort_values("delta", key=abs, ascending=False).drop_duplicates(["program_no", "lineage", "contrast"])
    D = sk[["program_no", "lineage", "contrast", "n_a", "n_b", "mean_a", "mean_b", "delta", "p_adj"]].rename(columns={
        "program_no": "distinct program", "lineage": "cell type", "contrast": "comparison (first vs second)", "n_a": "donors, first",
        "n_b": "donors, second", "mean_a": "mean score, first", "mean_b": "mean score, second",
        "delta": "change (second − first)", "p_adj": "adjusted p"}).sort_values(["distinct program", "adjusted p"])
    S = summ.copy()
    S["what it can tell us"] = S.design.map(lambda d: "not interpretable (donors sit in different lanes)" if d.startswith("CROSS") else
                                             ("limited: 9 pairs cannot reach adjusted p<0.05" if "paired" in d else "valid"))
    S = S.rename(columns={"donors_a": "donors, first group", "donors_b": "donors, second group", "tests": "programs × cell types tested",
                          "raw_p_below_0_05": "p<0.05 before correction", "significant": "significant after correction",
                          "expected_by_chance": "expected by chance at 0.05"})

    nprog = len(P0); nfun = len(prog)
    t1 = int((prog.relevance == "1").sum()); t2 = int((prog.relevance == "2").sum()); t3 = int((prog.relevance == "3").sum())
    nart = int((prog.relevance == "check handling artefact first").sum())
    nboth = int((P0.tissue_call == "both").sum()); nskin = int((P0.tissue_call == "skin only").sum())
    nchg = int(P0.changes_skin.sum()); ncheck = int(prog.evidence.str.startswith("checked").sum())
    cut = int(round(ch["perm_cut"] * 100))

    wb = Workbook(); wb.remove(wb.active)
    ws = wb.create_sheet("Read me")
    L = [
     ("IDS gene programs: shared or specific across tissue, cell type and disease state", 14, True),
     ("Psoriasis blood (healthy, psoriasis, psoriatic arthritis) and skin (non-lesional-healthy, non-lesional psoriasis, lesional).", 10, False),
     ("", 10, False),
     ("The question this sheet answers", 12, True),
     ("For each gene program: which genes, in which cell types, in which tissue, in which disease states, and what it does.", 10, False),
     ("Then a plain call on each axis: shared (many) or specific (few).", 10, False),
     ("", 10, False),
     ("How the programs were chosen (one chain, counts in order)", 12, True),
     (f"{ch['detections']} detections. A detection is one gene program found in one search (one analysis, one tissue, one disease state). The same program is often detected several times.", 10, False),
     (f"Step 0, remove technical: {ch['technical']} removed (mitochondrial, ribosomal, haemoglobin, platelet, sex-linked, housekeeping genes), leaving {ch['after_technical']}.", 10, False),
     (f"Filter 1, permutation test (the authors' own): for each gene pair the data are shuffled 60 times and the pair's score is compared with the shuffled scores (p = 0 means no shuffle beat it). A detection stays if at least {cut}% of its gene pairs beat all 60 shuffles. {ch['failed_permutation']} removed, leaving {ch['after_permutation']}.", 10, False),
     ("   With 60 shuffles the smallest non-zero p is 0.017, so 'beat all 60' supports p below roughly 0.05, not 0.01. 100 or more shuffles would be firmer.", 10, False),
     (f"   The {cut}% cutoff is our choice, not a tested value. With a cutoff of {int(round(ch['alt_cut'] * 100))}% of pairs instead, there would be {ch['alt_distinct_programs']} distinct programs rather than {nprog}.", 10, False),
     (f"Filter 2, activity: the program must score at least 0.2 in at least one cell subtype (50+ cells). {ch['not_active']} removed, leaving {ch['kept_detections']} detections.", 10, False),
     (f"Merge repeats: detections whose gene sets overlap by at least half are one distinct program. The {ch['kept_detections']} detections are {nprog} distinct programs.", 10, False),
     ("Specificity is not a filter. It is reported for every program (see below) so you can read shared vs specific.", 10, False),
     ("", 10, False),
     ("Permutation, activity and specificity in plain words", 12, True),
     ("Permutation test: is the dependence between the program's genes real, or would shuffled data give it too?", 10, False),
     ("Activity: does the program actually switch on (mean score 0.2 or more) in some cell type, or is it only a statistical link?", 10, False),
     ("Specificity: how concentrated the program is. 0 = equally active in every group, 1 = active in exactly one group. Shown for tissue, cell type and disease state.", 10, False),
     ("", 10, False),
     ("How to read a row on the 'Programs' sheet", 12, True),
     ("1. Tissue: 'both' = the program scores 0.2 or more in at least one cell subtype of blood AND of skin. 'skin only' = skin cell types only.", 10, False),
     ("2. Cell-type spread: how many broad lineages (T, NK, Myeloid, B, Keratinocyte…) carry it. 'broad' (3 or more) = shared; 'one subtype' = specific.", 10, False),
     ("3. Disease: 'Pattern in skin' says up, down, mixed or no change (normal vs lesional first, then non-lesional); filter on it to find programs specific to one disease state. Then compare the three mean scores for the strongest cell type. 'Changes with disease state?' says whether the donor-level test found a real difference.", 10, False),
     ("4. 'What it does' and 'Evidence' tell you how far to trust the meaning. Only 'checked' rows had a source read for the statement.", 10, False),
     ("5. Sort order: first 'Disease relevance' (1 = changes in lesional skin AND a psoriasis source was read for these genes; 2 = one of the two; 3 = neither), then 'Top activity score' (highest first). 'Score rank' ranks all 40 programs by that score alone, so you can also sort by it directly.", 10, False),
     ("", 10, False),
     ("Sheets", 12, True),
     (f"Programs: {nfun} distinct programs ({t1} relevance 1, {t2} relevance 2, {t3} relevance 3, {nart} flagged as a possible tissue-handling effect).", 10, False),
     (f"Cell-identity programs: {len(ident)} programs that mostly mark a cell type (T cell, pDC, melanocyte, sweat gland…). They show a cell type exists, so they are not ranked.", 10, False),
     (f"All {ch['detections']} detections: every detection with the outcome of the gates and the distinct program it belongs to. Use this to trace any number back.", 10, False),
     (f"Technical: the {ch['technical']} detections removed in step 0.", 10, False),
     ("Disease tests: the donor-level results behind 'Changes with disease state?', and what each comparison can and cannot show.", 10, False),
     ("", 10, False),
     ("Counts", 12, True),
     (f"{nprog} distinct programs: {nboth} active in both tissues, {nskin} in skin only. {nchg} change in lesional skin. {ncheck} of the {nfun} on the Programs sheet have a source read for the statement.", 10, False),
     ("", 10, False),
     ("What is solid and what is not", 12, True),
     ("Skin: lesional vs non-lesional differences are tested donor by donor (7 non-lesional-healthy, 10 non-lesional psoriasis, 13 lesional donors) and hold for the programs listed.", 10, False),
     ("Blood: in the valid comparisons no program differs (healthy 16 vs psoriasis 8 donors, lanes 03+05; psoriasis 6 vs PsA 12, PBMC-06 lanes; psoriasis 14 vs PsA 23, all lanes).", 10, False),
     ("   That means 'not detectable with these donor numbers', not 'no difference'. Healthy vs PsA shows 86 hits but donors sit in different lanes, so it is not interpretable.", 10, False),
     ("Skin-to-blood: of the 72 program × cell-type pairs that change in lesional skin and can be matched in blood, 0 change in blood (psoriasis only; PsA not covered).", 10, False),
     ("Skin PN vs PP (paired, 9 donors) cannot reach significance after correction, so it is not a finding either way.", 10, False),
     ("'Shared' is loose: a program shared with blood can be a monocyte gene that skin keratinocytes also use (S100A8/A9). Read the cell types before calling it the same biology.", 10, False),
     ("A program that rises in many unrelated skin cell types at once (e.g. antimicrobial S100 genes in T cells) may reflect RNA leaking from the dominant cell type; this was not tested.", 10, False),
     ("Immediate-early / AP-1 programs can be produced by handling the tissue. Their drop in lesional skin is not yet interpretable.", 10, False),
     ("The merge rule (gene sets overlapping by at least half) is a judgment call, not a standard method.", 10, False),
    ]
    for i, (t, sz, b) in enumerate(L, 1):
        ws.cell(row=i, column=1, value=t).font = Font(name=ARIAL, size=sz, bold=b)
    ws.column_dimensions["A"].width = 150

    W = {"ID": 7, "Disease relevance (1 = highest)": 12, "Score rank (1 = highest of 40)": 11, "Top activity score": 11, "Top activity: cell type (tissue)": 28, "Pattern in skin": 30, "Program": 34, "Genes": 55, "Tissue": 11, "Cell types active, blood": 34,
         "Cell types active, skin": 34, "Cell-type spread": 34, "Specificity (0 = equal everywhere, 1 = one group only)": 30,
         "Skin: mean score by state (strongest cell type)": 36, "Blood: mean score by state (strongest cell type)": 36,
         "Changes with disease state?": 44, "What it does / implication": 70, "Source": 34, "Link": 30, "Evidence": 18}
    ws = wb.create_sheet("Programs")
    ws.cell(row=1, column=1, value=f"{len(P)} distinct programs").font = Font(name=ARIAL, size=12, bold=True)
    ws.cell(row=2, column=1, value="Scores are mean log-expression gene-set scores (score_genes); 0.2 or more counts as active. Rows are sorted by Disease relevance, then by Top activity score (highest first). Filter on Tissue, Disease relevance or Evidence.").font = Font(name=ARIAL, size=10, italic=True, color="666666")
    write_table(ws, P, 4, W, wrap=tuple(c for c in P.columns if c not in ("ID", "Tissue", "Link", "Evidence", "Disease relevance (1 = highest)", "Score rank (1 = highest of 40)", "Top activity score")),
        fills={"Disease relevance (1 = highest)": {"1": T1, "2": T2, "check handling artefact first": FLAG}})
    ws = wb.create_sheet("Cell-identity programs")
    ws.cell(row=1, column=1, value=f"{len(I)} programs that mostly mark a cell type").font = Font(name=ARIAL, size=12, bold=True)
    write_table(ws, I, 3, {"ID": 7, "Cell-type signature": 40, "Genes": 55, "Tissue": 11, "Cell types active, blood": 34,
                           "Cell types active, skin": 34, "Skin change with disease": 40}, wrap=tuple(I.columns))
    ws = wb.create_sheet(f"All {ch['detections']} detections")
    ws.cell(row=1, column=1, value=f"All {ch['detections']} detections and what the gates did to each").font = Font(name=ARIAL, size=12, bold=True)
    ws.cell(row=2, column=1, value="'Detected in' is where the program was found, not where it is active. Filter 'outcome of the gates' to see each step.").font = Font(name=ARIAL, size=10, italic=True, color="666666")
    write_table(ws, A, 4, {"detection": 32, "outcome of the gates": 34, "distinct program": 10, "analysis": 13, "detected in (tissue_state)": 14,
                           "n_genes": 8, "genes": 60, "median_IDS": 10, "share of gene pairs that beat all 60 shuffles": 14, "cells_supporting": 12,
                           "technical_flag": 14, "highest-scoring cell type": 22, "tissue pattern (score ≥0.2)": 14}, wrap=("genes",))
    ws = wb.create_sheet("Technical")
    T = det[det.status == "removed: technical"][["program_id", "technical_flag", "technical_detail", "n_genes", "genes"]].fillna("")
    T["technical_flag"] = np.where(T.technical_flag == "", "gene content (ribosomal / mitochondrial / sex-linked)", T.technical_flag)
    ws.cell(row=1, column=1, value=f"{len(T)} detections removed as technical").font = Font(name=ARIAL, size=12, bold=True)
    ws.cell(row=2, column=1, value="Mitochondrial, ribosomal, haemoglobin, platelet, sex-linked, housekeeping. Platelet ambient RNA still forms a program inside T cells after platelet cells are removed.").font = Font(name=ARIAL, size=10, italic=True, color="666666")
    write_table(ws, T, 4, {"program_id": 32, "technical_flag": 30, "technical_detail": 26, "n_genes": 9, "genes": 80}, wrap=("genes",))
    ws = wb.create_sheet("Disease tests")
    ws.cell(row=1, column=1, value="Donor-level disease tests").font = Font(name=ARIAL, size=12, bold=True)
    ws.cell(row=2, column=1, value="One value per donor. Mann-Whitney (paired Wilcoxon for PN vs PP), Benjamini-Hochberg correction. Detail below lists skin results with adjusted p<0.05 and a change of at least 0.2 score units.").font = Font(name=ARIAL, size=10, italic=True, color="666666")
    write_table(ws, S, 4, {"tissue": 8, "contrast": 18, "design": 34, "donors, first group": 12, "donors, second group": 12,
                           "programs × cell types tested": 14, "p<0.05 before correction": 14, "significant after correction": 14,
                           "expected by chance at 0.05": 14, "what it can tell us": 50}, wrap=("design", "what it can tell us"))
    ws.freeze_panes = None
    r = 4 + len(S) + 3
    ws.cell(row=r - 1, column=1, value=f"Skin results by distinct program ({len(D)} rows)").font = Font(name=ARIAL, size=10, bold=True)
    for j, col in enumerate(D.columns, 1):
        c = ws.cell(row=r, column=j, value=col); c.font = Font(name=ARIAL, size=10, bold=True); c.fill = HEAD; c.border = BOX
        c.alignment = Alignment(wrap_text=True)
    for i, row in enumerate(D.itertuples(index=False), r + 1):
        for j, v in enumerate(row, 1):
            c = ws.cell(row=i, column=j, value=(float(v) if isinstance(v, (np.floating,)) else v)); c.font = Font(name=ARIAL, size=10); c.border = BOX
            if isinstance(v, (float, np.floating)):
                c.number_format = "0.000" if abs(v) < 1 else "0.00"

    os.makedirs(paths.DELIVERABLES, exist_ok=True)
    out = os.path.join(paths.DELIVERABLES, "IDS_program_catalogue_v2.xlsx")
    wb.save(out)
    print("[out]", out)
    print(json.dumps(ch))
    print(f"programs sheet {len(P)} (rel1 {t1}, rel2 {t2}, rel3 {t3}, artefact {nart}) | identity {len(I)} | both {nboth} | skin only {nskin} | skin change {nchg} | checked {ncheck} | disease-test rows {len(D)}")


if __name__ == "__main__":
    main()
