# Where things moved (2026-10-26 reorganisation)

Older reports and notes refer to the old paths on the left. Nothing was renamed inside a run folder.

| Old path | New path |
|---|---|
| `scripts/qc_joint_pipeline.py`, `cluster_pipeline.py`, `theta_tuning.py`, `rank_markers.py`, `nk_deg_export.py`, `nk_percell.py`, `nk_deg_workbook.py` | `scripts/01_qc_labelling/` |
| `scripts/build_disease_object.py`, `ids_programs_disease.py`, `ids_program_atlas.py` | `scripts/02_ids_programs/` |
| `scripts/program_catalogue.py`, `disease_contrasts.py`, `distinct_programs.py`, `keratinocyte_detection.py`, `catalogue_workbook.py`, `program_annotations.py`, `marker_rules.py` | `scripts/03_catalogue/` |
| `results/PSA26`, `qc-joint`, `qc-joint-minfilter`, `clustered-nofilter`, `clustered-theta4`, `theta`, `markers-theta4` | `results/01_qc_labelling/<same name>` |
| `results/ids-disease`, `ids-lin-TNK`, `ids-lin-MY`, `ids-lin-KC`, `ids-atlas`, `skin_raw.h5ad` | `results/02_ids_programs/<same name>` |
| `results/catalogue` | `results/03_catalogue` |
| `docs/IDS_program_catalogue_v2.xlsx`, `docs/NK_cluster_DEG_for_review.xlsx` | `deliverables/` |
| `docs/Agenda1_QC_report.docx`, `Agenda2_clustering_report.docx`, `Theta_and_labelling_final.docx`, `Research_Agenda_4_skin.docx` | `docs/reports/` |
| `docs/abstract_notes.md` | `docs/abstract/abstract_notes.md` |
| `docs/agenda1_and_3_results.md`, `agenda2_combining_strata.md`, `scope_skin_extension.md`, `qc_labpipeline_both_tissues.md`, `architecture.html` | `docs/notes/` |
| Everything else listed in `_archive/MANIFEST.md` | `_archive/<same relative path>` |

Code no longer hardcodes any of these: it reads them from `src/sc_ids/paths.py`.

## November 2026: one data layout

| Before | Now |
|---|---|
| `scripts/02_ids_programs/build_disease_object.py` (one merged object per tissue: `blood_qc.h5ad`, `skin_qc.h5ad`) | `scripts/02_ids_programs/build_units.py` (one file per lane / sample under `results/00_units/`, read through `sc_ids.units`) |
| discovery code inside `ids_programs_disease.py` | `src/sc_ids/discovery.py`, shared by stage 02 and stage 04 |
| `scripts/02_ids_programs/ids_program_atlas.py`, `results/02_ids_programs/ids-atlas` | archived: superseded by stage 03, and its disease tests counted each lane as a separate donor |
