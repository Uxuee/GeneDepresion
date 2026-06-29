"""Annotate published AHBA modules using Hawrylycz et al. 2015 supplementary enrichments."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
MANUAL_DIR = PROJECT_ROOT / "data" / "ahba" / "published_modules"
PUB_DIR = PROJECT_ROOT / "outputs" / "length_gt_100_clean_pipeline" / "ahba_modules" / "published_modules"
ENR_DIR = PUB_DIR / "enrichment"
OUT_DIR = PUB_DIR / "annotations"
SUPP8_URL = "https://static-content.springer.com/esm/art%3A10.1038%2Fnn.4171/MediaObjects/41593_2015_BFnn4171_MOESM99_ESM.csv"
SUPP8_LOCAL = MANUAL_DIR / "Hawrylycz_2015_likely_Supplementary_Table_8_MOESM99.csv"
SOURCE_REF = "Hawrylycz et al. 2015, Canonical genetic signatures of the adult human brain, Supplementary Table 8 / MOESM99"

CATEGORY_PRIORITY = {
    "GO: Biological Process": 1,
    "GO: Molecular Function": 2,
    "GO: Cellular Component": 3,
    "Brain Gene Expression": 4,
    "Human Disease": 5,
    "Mouse Phenotype": 6,
}


def download_if_needed() -> tuple[Path | None, str]:
    MANUAL_DIR.mkdir(parents=True, exist_ok=True)
    if SUPP8_LOCAL.exists():
        return SUPP8_LOCAL, "local"
    try:
        import urllib.request

        print(f"Attempting official Supplementary Table 8 download: {SUPP8_URL}")
        urllib.request.urlretrieve(SUPP8_URL, SUPP8_LOCAL)
        return SUPP8_LOCAL, "downloaded"
    except Exception as exc:  # noqa: BLE001 - leave manual fallback.
        msg = (
            "Supplementary Table 8 / MOESM99 could not be downloaded automatically.\n"
            f"Official URL tried: {SUPP8_URL}\n"
            f"Place the CSV at: {SUPP8_LOCAL}\n"
            f"Error: {exc}\n"
        )
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        (OUT_DIR / "AHBA_MODULE_ANNOTATIONS_NOT_CREATED.txt").write_text(msg, encoding="utf-8")
        print(msg)
        return None, "missing"


def relevant_modules() -> set[str]:
    corrected = pd.read_csv(ENR_DIR / "published_ahba_enrichment_corrected_results_table.csv")
    nominal = pd.read_csv(ENR_DIR / "published_ahba_enrichment_nominal_results_table.csv")
    modules = set(corrected["AHBA_module_ID"].astype(str))
    modules.update(nominal.head(20)["AHBA_module_ID"].astype(str))
    return modules


def module_color_table() -> pd.DataFrame:
    assign_path = PUB_DIR / "ahba_published_gene_to_module_assignments.csv"
    if not assign_path.exists():
        return pd.DataFrame(columns=["AHBA_module_ID", "module_label_or_color"])
    assign = pd.read_csv(assign_path)
    color_col = None
    for col in ["module_label", "published_module"]:
        if col in assign.columns:
            color_col = col
            break
    out = assign.groupby("published_module").agg(module_gene_count=("gene_symbol", "nunique")).reset_index()
    out = out.rename(columns={"published_module": "AHBA_module_ID"})
    if color_col:
        labels = assign.groupby("published_module")[color_col].agg(lambda x: sorted(set(map(str, x)))[0]).reset_index()
        labels = labels.rename(columns={"published_module": "AHBA_module_ID", color_col: "module_label_or_color"})
        out = out.merge(labels, on="AHBA_module_ID", how="left")
    return out


def concise_interpretation(rows: pd.DataFrame) -> str:
    if rows.empty:
        return "No annotation terms available from the loaded Supplementary Table 8 file."
    names = rows["Name"].dropna().astype(str).head(5).tolist()
    categories = rows["Category"].dropna().astype(str).head(5).tolist()
    phrase = "; ".join(names)
    cats = ", ".join(dict.fromkeys(categories))
    return f"Top Hawrylycz enrichment terms point to {phrase}. Categories represented: {cats}."


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    supp_path, status = download_if_needed()
    if supp_path is None:
        return
    ann = pd.read_csv(supp_path)
    required = {"Category", "Name", "q-value FDR B&H", "Module"}
    missing = sorted(required - set(ann.columns))
    if missing:
        raise ValueError(f"Supplementary annotation table missing columns: {missing}")
    ann["AHBA_module_ID"] = ann["Module"].astype(str).str.strip()
    ann["annotation_q_value"] = pd.to_numeric(ann["q-value FDR B&H"], errors="coerce")
    ann["annotation_p_value"] = pd.to_numeric(ann.get("p-value"), errors="coerce")
    ann["category_priority"] = ann["Category"].map(CATEGORY_PRIORITY).fillna(99)
    modules = relevant_modules()
    rel = ann[ann["AHBA_module_ID"].isin(modules)].copy()
    rel = rel.sort_values(["AHBA_module_ID", "annotation_q_value", "category_priority", "annotation_p_value"])
    rel_out = rel[[
        "AHBA_module_ID", "Category", "ID", "Name", "annotation_p_value", "annotation_q_value",
        "Hit Count in Query List", "Hit Count in Genome", "Hit in Query List",
    ]].rename(columns={
        "Category": "annotation_category",
        "ID": "annotation_id",
        "Name": "annotation_name",
        "Hit Count in Query List": "hit_count_in_module",
        "Hit Count in Genome": "hit_count_in_genome",
        "Hit in Query List": "annotation_hit_genes",
    })
    rel_out["source_table"] = str(supp_path)
    rel_out["source_reference"] = SOURCE_REF
    rel_out.to_csv(OUT_DIR / "published_ahba_module_annotation_terms_relevant.csv", index=False)

    color = module_color_table()
    summary_rows = []
    for module in sorted(modules, key=lambda x: int(x[1:]) if x.startswith("M") and x[1:].isdigit() else x):
        sub = rel[rel["AHBA_module_ID"].eq(module)].copy()
        top = sub.sort_values(["annotation_q_value", "category_priority"]).head(8)
        color_row = color[color["AHBA_module_ID"].eq(module)]
        module_size = int(color_row["module_gene_count"].iloc[0]) if not color_row.empty else np.nan
        module_label = color_row["module_label_or_color"].iloc[0] if not color_row.empty and "module_label_or_color" in color_row else module
        summary_rows.append({
            "AHBA_module_ID": module,
            "module_label_or_color": module_label,
            "module_gene_count": module_size,
            "top_annotation_categories": "; ".join(top["Category"].dropna().astype(str).head(5)),
            "top_annotation_terms": "; ".join(top["Name"].dropna().astype(str).head(5)),
            "best_annotation_q_value": top["annotation_q_value"].min() if not top.empty else np.nan,
            "short_biological_interpretation": concise_interpretation(top),
            "source_table": str(supp_path),
            "source_reference": SOURCE_REF,
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_DIR / "published_ahba_module_annotation_summary.csv", index=False)

    corrected = pd.read_csv(ENR_DIR / "published_ahba_enrichment_corrected_results_table.csv")
    annotated_hits = corrected.merge(summary, on="AHBA_module_ID", how="left")
    annotated_hits.to_csv(OUT_DIR / "published_ahba_corrected_suggestive_hits_with_annotations.csv", index=False)

    md = [
        "# Published AHBA Module Annotation Summary",
        "",
        f"Source: {SOURCE_REF}",
        f"Supplement file status: {status}",
        "",
        "Modules annotated are those appearing in corrected/suggestive rows and the top nominal enrichment rows.",
        "",
        summary[["AHBA_module_ID", "module_gene_count", "top_annotation_terms", "best_annotation_q_value", "short_biological_interpretation"]].to_markdown(index=False),
        "",
        "These annotations describe healthy-brain coexpression modules from Hawrylycz et al. and should be interpreted as biological context, not disease validation.",
    ]
    (OUT_DIR / "published_ahba_module_annotation_notes.md").write_text("\n".join(md), encoding="utf-8")
    report = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "supplement_table": str(supp_path),
        "supplement_status": status,
        "modules_annotated": sorted(modules),
        "summary_table": str(OUT_DIR / "published_ahba_module_annotation_summary.csv"),
        "term_table": str(OUT_DIR / "published_ahba_module_annotation_terms_relevant.csv"),
    }
    (OUT_DIR / "published_ahba_module_annotation_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
