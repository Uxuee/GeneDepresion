"""Copy original project-generated manuscript figures with provenance tracking.

The script avoids assistant draft figures from 14_Preliminary_Manuscript_Package/02_figures_main
and copies original pipeline PNG/PDF outputs into 02_figures_main_original and
03_figures_supplementary_original. Figure 1 is regenerated only by the deterministic
project script code/make_workflow_schematic.py.
"""
from __future__ import annotations

import hashlib
import shutil
import subprocess
from pathlib import Path
from typing import Iterable

import pandas as pd

ROOT = Path(r"D:/genes/GeneDepresion")
CLEAN = ROOT / "outputs/length_gt_100_clean_pipeline"
PKG = CLEAN / "14_Preliminary_Manuscript_Package"
MAIN = PKG / "02_figures_main_original"
SUPP = PKG / "03_figures_supplementary_original"
DRAFT_MAIN = PKG / "02_figures_main"
MANIFEST = PKG / "manuscript_figure_provenance_manifest.csv"
REPORT = PKG / "Original_Figure_Copy_Report.md"

for folder in (MAIN, SUPP):
    folder.mkdir(parents=True, exist_ok=True)

DRAFT_EXCLUDE = {
    "Figure_1_workflow_schematic.png",
    "Figure_2_Tier2_ML1_candidate_summary.png",
    "Figure_3_AHBA_context.png",
    "Figure_3_AHBA_context_summary.png",
    "Figure_4_broad_region_WGCNA_networks.png",
    "Figure_4_broad_region_WGCNA_prioritized_networks.png",
    "Figure_6_Module_Level_ML2_performance.png",
    "Figure_6_module_level_ML2_performance.png",
}

rows: list[dict] = []


def md5(path: Path) -> str:
    if not path.exists():
        return ""
    h = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def is_draft(path: Path) -> bool:
    try:
        path.relative_to(DRAFT_MAIN)
        return True
    except ValueError:
        return path.name in DRAFT_EXCLUDE


def find_files(search_roots: Iterable[Path], patterns: Iterable[str]) -> list[Path]:
    found: list[Path] = []
    for root in search_roots:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if path.suffix.lower() not in {".png", ".pdf"}:
                continue
            if is_draft(path):
                continue
            lower = path.name.lower()
            if any(p.lower().replace("*", "") in lower for p in patterns):
                found.append(path)
    return sorted(set(found), key=lambda p: (p.suffix.lower() != ".png", len(str(p)), str(p)))


def matching_pdf(png: Path) -> Path | None:
    pdf = png.with_suffix(".pdf")
    return pdf if pdf.exists() and not is_draft(pdf) else None


def copy_record(fig: str, panel: str, src: Path | None, dst_name: str, use: str, source_type: str, notes: str) -> None:
    if src is None or not src.exists():
        rows.append({
            "figure_id": fig,
            "panel_id": panel,
            "copied_filename": dst_name,
            "source_filepath": "",
            "source_type": "missing",
            "png_found": False,
            "pdf_found": False,
            "md5_hash": "",
            "recommended_use": use,
            "notes": notes,
        })
        return
    dst = MAIN / dst_name if use == "main" else SUPP / dst_name
    if src.resolve() != dst.resolve():
        shutil.copy2(src, dst)
    rows.append({
        "figure_id": fig,
        "panel_id": panel,
        "copied_filename": dst.name,
        "source_filepath": str(src),
        "source_type": source_type,
        "png_found": src.suffix.lower() == ".png",
        "pdf_found": src.suffix.lower() == ".pdf",
        "md5_hash": md5(dst),
        "recommended_use": use,
        "notes": notes,
    })


def copy_png_pdf(fig: str, panel: str, png: Path | None, out_stem: str, notes: str, source_type: str = "original_pipeline_output") -> None:
    if png is None or not png.exists():
        copy_record(fig, panel, None, out_stem + ".png", "main", "missing", notes)
        return
    copy_record(fig, panel, png, out_stem + ".png", "main", source_type, notes)
    pdf = matching_pdf(png)
    if pdf:
        copy_record(fig, panel, pdf, out_stem + ".pdf", "main", source_type, notes + " Matching PDF found.")
    else:
        rows.append({
            "figure_id": fig,
            "panel_id": panel,
            "copied_filename": out_stem + ".pdf",
            "source_filepath": "",
            "source_type": "missing",
            "png_found": True,
            "pdf_found": False,
            "md5_hash": "",
            "recommended_use": "main",
            "notes": notes + " Matching PDF missing.",
        })


def choose_first(paths: list[Path], preferred: Iterable[str]) -> Path | None:
    if not paths:
        return None
    for needle in preferred:
        for path in paths:
            if needle.lower() in path.name.lower() and path.suffix.lower() == ".png":
                return path
    for path in paths:
        if path.suffix.lower() == ".png":
            return path
    return paths[0]


# Figure 1: deterministic project-script workflow.
subprocess.run(["python", str(ROOT / "code/make_workflow_schematic.py")], cwd=ROOT, check=True)
copy_record(
    "Figure 1",
    "workflow",
    MAIN / "Figure_1_workflow_schematic_REGENERATED_BY_SCRIPT.png",
    "Figure_1_workflow_schematic_REGENERATED_BY_SCRIPT.png",
    "main",
    "regenerated_by_project_script",
    "Deterministic text-box/arrow schematic generated by code/make_workflow_schematic.py; no generative AI imagery.",
)
copy_record(
    "Figure 1",
    "workflow",
    MAIN / "Figure_1_workflow_schematic_REGENERATED_BY_SCRIPT.pdf",
    "Figure_1_workflow_schematic_REGENERATED_BY_SCRIPT.pdf",
    "main",
    "regenerated_by_project_script",
    "Deterministic text-box/arrow schematic generated by code/make_workflow_schematic.py; no generative AI imagery.",
)

# Figure 2: original Tier2 + ML #1 candidate summary.
fig2_candidates = find_files(
    [CLEAN / "10_ML_prioritized_Tier2_corrected", CLEAN / "07_Tier2_candidate_gene_sets"],
    ["Tier2", "ML", "count", "high_moderate", "candidate", "MLpriority", "counts_by_region", "summary"],
)
fig2 = choose_first(fig2_candidates, ["barplot_high_moderate", "mlprioritized_count", "counts_by_region", "summary"])
copy_png_pdf(
    "Figure 2",
    "Tier2_ML1_summary",
    fig2,
    "Figure_2_original_Tier2_ML1_candidate_summary",
    "Original pipeline figure selected for Tier2 + ML #1 candidate summary. If missing, create a deterministic project-script figure from CSV tables later.",
)

# Figure 3: original AHBA/context figure.
fig3_candidates = find_files(
    [CLEAN / "11_AHBA_modules_enrichment", CLEAN / "12_MDD_WGCNA_modules"],
    ["ahba_M26_evidence_across_gene_sets", "ahba_enrichment_evidence_ladder_by_gene_set", "published_ahba_module_enrichment_dotplot", "AHBA_overlap_annotation_summary_plot", "ahba", "dotplot", "overlap"],
)
fig3 = choose_first(fig3_candidates, ["ahba_m26_evidence", "ahba_enrichment_evidence_ladder", "published_ahba_module_enrichment_dotplot", "AHBA_overlap_annotation_summary_plot"])
copy_png_pdf(
    "Figure 3",
    "AHBA_context",
    fig3,
    "Figure_3_original_AHBA_context",
    "Original AHBA/context figure. AHBA is healthy-brain coexpression context, not disease validation.",
)

# Figure 4: broad-region WGCNA networks, prefer original multi-panel summary.
summary = CLEAN / "12_MDD_WGCNA_modules/plots_cluster_visualization/full_module_network_summary_top_modules.png"
if summary.exists():
    copy_png_pdf("Figure 4", "broad_WGCNA_summary", summary, "Figure_4_original_broad_region_WGCNA_network_summary", "Original multi-panel broad-region full-module network summary.")
else:
    targets = [
        "full_module_network_Ventral_Striatum_MElavenderblush3",
        "full_module_network_Ventral_Striatum_MElightcyan1",
        "full_module_network_Ventral_Striatum_MEdarkred",
        "full_module_network_Hippocampal_Formation_MEgreen",
        "full_module_network_Insula_MEorange",
        "full_module_network_labeled_overlaps_Ventral_Striatum_MElavenderblush3",
        "full_module_network_labeled_overlaps_Ventral_Striatum_MElightcyan1",
        "full_module_network_labeled_overlaps_Ventral_Striatum_MEdarkred",
        "full_module_network_labeled_overlaps_Hippocampal_Formation_MEgreen",
        "full_module_network_labeled_overlaps_Insula_MEorange",
    ]
    for target in targets:
        png = CLEAN / f"12_MDD_WGCNA_modules/plots_cluster_visualization/{target}.png"
        copy_png_pdf("Figure 4", target, png if png.exists() else None, f"Figure_4_original_{target}", "Original individual broad-region WGCNA module network.")

# Figure 5: exact original orange-overlap frontal networks.
front_dir = CLEAN / "12_MDD_WGCNA_modules/frontal_subregion_wgcna/plots/cluster_visualizations_orange_overlap"
fig5_files = [
    "dlPFC_BA8_9_MEblue_Tier2_overlap_orange_labeled.png",
    "dlPFC_BA8_9_MEorange_Tier2_overlap_orange_labeled.png",
    "CG25_MEwhite_Tier2_overlap_orange_labeled.png",
    "OFC_BA11_MEtan_highModerateML_overlap_orange_labeled.png",
]
for i, name in enumerate(fig5_files, 1):
    copy_png_pdf("Figure 5", f"panel_{i}", front_dir / name, f"Figure_5_original_panel_{i}_{Path(name).stem}", "Exact original frontal subregion orange-overlap network requested for manuscript package.")

# Supplementary original copies: keep useful originals separate from main.
for src_root in [
    CLEAN / "12_MDD_WGCNA_modules/plots_cluster_visualization",
    CLEAN / "12_MDD_WGCNA_modules/frontal_subregion_wgcna/plots/cluster_visualizations_orange_overlap",
    CLEAN / "11_AHBA_modules_enrichment/published_modules/figures",
    CLEAN / "12_MDD_WGCNA_modules/biological_annotation/plots",
]:
    if not src_root.exists():
        continue
    for path in src_root.glob("*"):
        if path.suffix.lower() not in {".png", ".pdf"} or is_draft(path):
            continue
        dst_name = path.name
        dst = SUPP / dst_name
        if dst.exists():
            dst_name = f"{src_root.name}_{path.name}"
        copy_record("Supplementary", path.stem, path, dst_name, "supplementary", "original_pipeline_output", "Original pipeline output copied to supplementary original figure folder.")

pd.DataFrame(rows).to_csv(MANIFEST, index=False)

missing = [r for r in rows if r["source_type"] == "missing"]
main_rows = [r for r in rows if r["recommended_use"] == "main"]
main_present_rows = [r for r in main_rows if r["source_type"] != "missing"]
pairs = pd.DataFrame(main_present_rows).groupby(["figure_id", "panel_id"], dropna=False).agg(
    png_found=("png_found", "max"),
    pdf_found=("pdf_found", "max"),
    files=("copied_filename", lambda x: "; ".join([str(v) for v in x if str(v).lower() != "nan"])),
).reset_index()
missing_pairs = pd.DataFrame(main_rows).groupby(["figure_id", "panel_id"], dropna=False).agg(
    any_missing=("source_type", lambda x: any(v == "missing" for v in x)),
    missing_files=("copied_filename", lambda x: "; ".join([str(v) for v in x if str(v).lower().endswith(".pdf")]))
).reset_index()
if not pairs.empty and not missing_pairs.empty:
    pairs = pairs.merge(missing_pairs, on=["figure_id", "panel_id"], how="left")
    pairs["notes"] = pairs.apply(lambda r: "Missing expected file(s): " + str(r["missing_files"]) if r.get("any_missing", False) else "Complete PNG/PDF pair copied.", axis=1)
    pairs = pairs.drop(columns=["any_missing", "missing_files"])

report_lines = [
    "# Original Figure Copy Report",
    "",
    "Assistant draft figures in `02_figures_main/` were excluded from the manuscript-original copy set.",
    "Figures in `02_figures_main_original/` were copied from original project outputs when available.",
    "Figure 1 was regenerated by deterministic repository script `code/make_workflow_schematic.py` using only text boxes and arrows.",
    "No generative AI images were used.",
    "No original pipeline files were deleted or overwritten.",
    "",
    "## Main Figure PNG/PDF Availability",
    "",
    pairs.to_markdown(index=False),
    "",
    "## Missing Items",
    "",
]
if missing:
    report_lines.append(pd.DataFrame(missing).to_markdown(index=False))
else:
    report_lines.append("No required main figure files were missing. Some optional matching PDFs may be marked missing only if the original pipeline did not create one.")
report_lines += [
    "",
    "## Provenance Manifest",
    "",
    f"Full copy provenance is recorded in `{MANIFEST}`.",
]
REPORT.write_text("\n".join(report_lines), encoding="utf-8")
print(f"Wrote {MANIFEST}")
print(f"Wrote {REPORT}")
print(f"Main rows: {len(main_rows)}; supplementary rows: {len(rows)-len(main_rows)}; missing rows: {len(missing)}")



