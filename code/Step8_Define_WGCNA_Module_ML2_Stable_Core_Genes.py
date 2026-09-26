#!/usr/bin/env python
"""Define stable core genes for future module-level ML #2 features.

This script uses saved Step 12 MDD WGCNA outputs only. It does not run
classification and does not select genes by diagnosis p-value, Tier2 status, or
ML-prioritized status. Tier2/ML flags are added only after stable genes are
selected by WGCNA module membership and kME/hubness.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(r"D:/genes/GeneDepresion")
PIPELINE_ROOT = PROJECT_ROOT / "outputs" / "length_gt_100_clean_pipeline"
WGCNA_ROOT = PIPELINE_ROOT / "12_MDD_WGCNA_modules"
OUT_DIR = WGCNA_ROOT / "module_level_ml2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FOCUS_MODULES = [
    ("Ventral Striatum", "MElavenderblush3"),
    ("Ventral Striatum", "MElightcyan1"),
    ("Ventral Striatum", "MEdarkred"),
    ("Ventral Striatum", "MElightgreen"),
    ("Insula", "MEorange"),
    ("Hippocampal Formation", "MEgreen"),
    ("Frontal Cortex", "MEdarkgreen"),
]

EXPR_FILES = {
    "Frontal Cortex": PIPELINE_ROOT / "02_normalized_expression" / "normalized_GeTMM_log2_frontal_cortex.csv",
    "Hippocampal Formation": PIPELINE_ROOT / "02_normalized_expression" / "normalized_GeTMM_log2_hippocampal_formation.csv",
    "Insula": PIPELINE_ROOT / "02_normalized_expression" / "normalized_GeTMM_log2_insula.csv",
    "Ventral Striatum": PIPELINE_ROOT / "02_normalized_expression" / "normalized_GeTMM_log2_ventral_striatum.csv",
}

INPUTS = {
    "parameter_sensitivity": WGCNA_ROOT / "mdd_wgcna_parameter_sensitivity.csv",
    "assignments": WGCNA_ROOT / "best_solution_gene_to_module_assignments.csv",
    "hubs": WGCNA_ROOT / "module_hub_genes_top50.csv",
    "top_priority": WGCNA_ROOT / "mdd_wgcna_top_priority_modules_for_report.csv",
    "trait": WGCNA_ROOT / "module_trait_association_results.csv",
    "module_enrichment": WGCNA_ROOT / "module_enrichment_for_Tier2_and_ML_genes.csv",
    "tier2": PIPELINE_ROOT / "07_Tier2_candidate_gene_sets" / "Tier2_candidate_genes_all_regions.csv",
    "ml": PIPELINE_ROOT / "10_ML_prioritized_Tier2_corrected" / "ML_prioritized_Tier2_genes_all_ranked.csv",
}

CORE_KME_THRESHOLD = 0.60
MIN_STABLE_CORE_GENES = 25
TOP_HUB_CAP = 50


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def expr_gene_ids(region: str) -> set[str]:
    # Load only the first column/index; matrices are gene-by-sample CSVs with gene ID index.
    return set(pd.read_csv(EXPR_FILES[region], usecols=[0]).iloc[:, 0].astype(str))


def classify_ml_sets(ml: pd.DataFrame) -> tuple[set[str], set[str], set[str], set[str]]:
    ml = ml.copy()
    ml["gene_id"] = ml["gene_id"].astype(str)
    ml["gene_symbol"] = ml["gene_symbol"].astype(str)
    cls = ml["ML_priority_class"].fillna("").str.lower()
    high = ml[cls.str.contains("high-confidence")]
    high_mod = ml[cls.str.contains("high-confidence|moderate-confidence")]
    return set(high["gene_id"]), set(high["gene_symbol"]), set(high_mod["gene_id"]), set(high_mod["gene_symbol"])


def bool_flag(series: pd.Series) -> str:
    return "TRUE" if bool(series.any()) else "FALSE"


def main() -> None:
    assignments = read_csv(INPUTS["assignments"])
    assignments["gene_id"] = assignments["gene_id"].astype(str)
    assignments["gene_symbol"] = assignments["gene_symbol"].astype(str)

    hubs = read_csv(INPUTS["hubs"])
    hubs["gene_id"] = hubs["gene_id"].astype(str)
    hubs["gene_symbol"] = hubs["gene_symbol"].astype(str)
    hubs["abs_kME"] = pd.to_numeric(hubs["abs_kME"], errors="coerce")
    hubs["kME"] = pd.to_numeric(hubs["kME"], errors="coerce")

    parameter_sensitivity = read_csv(INPUTS["parameter_sensitivity"])
    tier2 = read_csv(INPUTS["tier2"])
    tier2_ids = set(tier2["gene_id"].dropna().astype(str))
    tier2_symbols = set(tier2["gene_symbol"].dropna().astype(str))

    ml = read_csv(INPUTS["ml"])
    high_ids, high_symbols, high_mod_ids, high_mod_symbols = classify_ml_sets(ml)

    run_level_assignment_files = list(WGCNA_ROOT.glob("**/*gene_to_module*parameter*.csv")) + list(WGCNA_ROOT.glob("**/*assignments*parameter*.csv"))
    run_level_assignment_files = [p for p in run_level_assignment_files if p.name != "best_solution_gene_to_module_assignments.csv"]
    has_run_level_assignments = len(run_level_assignment_files) > 0
    stability_method = "run_level_assignment_frequency" if has_run_level_assignments else "best_solution_kME_hubness_fallback"
    limitation_note = "Run-level gene-to-module assignment files were not saved for each parameter-grid run; assignment_stability_frequency is unavailable and stable cores use best-solution kME/hub rank."

    expression_gene_cache = {region: expr_gene_ids(region) for region, _ in FOCUS_MODULES}

    all_rows = []
    summary_rows = []
    for region, module in FOCUS_MODULES:
        module_assign = assignments[(assignments["region"] == region) & (assignments["module"] == module) & (assignments["module"] != "MEgrey")].copy()
        module_assign = module_assign.drop_duplicates("gene_id")
        module_size = len(module_assign)
        if module_size == 0:
            continue

        module_hubs = hubs[(hubs["region"] == region) & (hubs["module"] == module)].copy()
        module_hubs = module_hubs.sort_values("abs_kME", ascending=False).drop_duplicates("gene_id")
        hub_lookup = module_hubs.set_index("gene_id")[["kME", "abs_kME"]].to_dict("index")
        hub_rank = {gid: rank for rank, gid in enumerate(module_hubs["gene_id"], start=1)}

        expr_ids = expression_gene_cache[region]
        module_assign["present_in_expression_matrix"] = module_assign["gene_id"].isin(expr_ids)
        module_assign["kME"] = module_assign["gene_id"].map(lambda gid: hub_lookup.get(gid, {}).get("kME", np.nan))
        module_assign["abs_kME"] = module_assign["gene_id"].map(lambda gid: hub_lookup.get(gid, {}).get("abs_kME", np.nan))
        module_assign["hub_rank"] = module_assign["gene_id"].map(hub_rank)

        # Stable core is based only on non-grey module membership, expression presence, and kME/hubness.
        core = module_assign[(module_assign["present_in_expression_matrix"]) & (module_assign["abs_kME"] >= CORE_KME_THRESHOLD)].copy()
        core_notes = [f"selected by |kME| >= {CORE_KME_THRESHOLD} from best-solution hub table"]
        if len(core) < min(MIN_STABLE_CORE_GENES, module_size):
            fallback_n = min(MIN_STABLE_CORE_GENES, module_size)
            fallback = module_assign[module_assign["present_in_expression_matrix"]].sort_values(["abs_kME", "hub_rank"], ascending=[False, True], na_position="last").head(fallback_n)
            core = pd.concat([core, fallback], ignore_index=True).drop_duplicates("gene_id")
            core_notes.append(f"fallback added top hub-ranked genes to reach up to {fallback_n} stable-core genes")

        top_hubs = module_assign[module_assign["present_in_expression_matrix"]].sort_values(["abs_kME", "hub_rank"], ascending=[False, True], na_position="last").head(min(TOP_HUB_CAP, module_size)).copy()

        selected_core_ids = set(core["gene_id"])
        selected_hub_ids = set(top_hubs["gene_id"])

        def base_row(row, gene_set_type: str) -> dict:
            gid = str(row["gene_id"])
            sym = str(row["gene_symbol"])
            in_core = gid in selected_core_ids
            in_hub = gid in selected_hub_ids
            notes = []
            if not row.get("present_in_expression_matrix", False):
                notes.append("not present in expression matrix used for Step6")
            if gene_set_type == "stable_core_genes":
                notes.extend(core_notes)
            elif gene_set_type == "top_hub_genes":
                notes.append(f"top {min(TOP_HUB_CAP, module_size)} by abs(kME)/hub rank")
            else:
                notes.append("all non-grey best-solution module genes")
            if not has_run_level_assignments:
                notes.append("assignment stability frequency unavailable; no run-level assignment tables found")
            return {
                "region": region,
                "module": module,
                "gene_symbol": sym,
                "gene_id": gid,
                "gene_set_type": gene_set_type,
                "module_size": module_size,
                "kME": row.get("kME", np.nan),
                "abs_kME": row.get("abs_kME", np.nan),
                "hub_rank": row.get("hub_rank", np.nan),
                "present_in_expression_matrix": bool(row.get("present_in_expression_matrix", False)),
                "assignment_stability_frequency": np.nan,
                "selected_for_stable_core": in_core,
                "selected_for_top_hub": in_hub,
                "is_Tier2": gid in tier2_ids or sym in tier2_symbols,
                "is_high_confidence_ML": gid in high_ids or sym in high_symbols,
                "is_high_plus_moderate_ML": gid in high_mod_ids or sym in high_mod_symbols,
                "notes": " | ".join(notes),
            }

        for _, row in module_assign.sort_values(["abs_kME", "hub_rank", "gene_symbol"], ascending=[False, True, True], na_position="last").iterrows():
            all_rows.append(base_row(row, "all_module_genes"))
        for _, row in core.sort_values(["abs_kME", "hub_rank", "gene_symbol"], ascending=[False, True, True], na_position="last").iterrows():
            all_rows.append(base_row(row, "stable_core_genes"))
        for _, row in top_hubs.sort_values(["abs_kME", "hub_rank", "gene_symbol"], ascending=[False, True, True], na_position="last").iterrows():
            all_rows.append(base_row(row, "top_hub_genes"))

        summary_rows.append({
            "region": region,
            "module": module,
            "module_size": module_size,
            "n_all_genes": module_size,
            "n_stable_core_genes": len(selected_core_ids),
            "n_top_hub_genes": len(selected_hub_ids),
            "n_Tier2_in_stable_core": int(core["gene_id"].isin(tier2_ids).sum() + (~core["gene_id"].isin(tier2_ids) & core["gene_symbol"].isin(tier2_symbols)).sum()),
            "n_ML_in_stable_core": int(core["gene_id"].isin(high_mod_ids).sum() + (~core["gene_id"].isin(high_mod_ids) & core["gene_symbol"].isin(high_mod_symbols)).sum()),
            "median_kME_stable_core": float(pd.to_numeric(core["abs_kME"], errors="coerce").median()) if len(core) else np.nan,
            "stability_method": stability_method,
            "notes": "; ".join(core_notes + ([limitation_note] if not has_run_level_assignments else [])),
        })

    stable_table = pd.DataFrame(all_rows)
    summary = pd.DataFrame(summary_rows)

    stable_table.to_csv(OUT_DIR / "stable_core_genes_by_module.csv", index=False)
    summary.to_csv(OUT_DIR / "stable_core_gene_summary_by_module.csv", index=False)

    enough = summary[summary["n_stable_core_genes"] >= 25].copy()
    report_lines = [
        "# Stable Core Genes for Module-Level ML #2",
        "",
        "This folder defines stable gene sets for future module-level ML #2 features/eigengenes. No classification was run in this step.",
        "",
        "## Selection Principle",
        "",
        "Stable genes were selected from WGCNA/module information only: non-grey best-solution module membership, presence in the normalized expression matrix, and module membership/hubness measured by kME or hub rank. Diagnosis p-values were not used for selection. Tier2 and ML-prioritized flags were added only after selection as annotation columns.",
        "",
        "## Stability Method",
        "",
        f"Run-level gene-to-module assignment files found: {len(run_level_assignment_files)}.",
        limitation_note if not has_run_level_assignments else "Assignment stability was estimated from run-level assignment frequencies.",
        f"Because run-level assignments were unavailable, stable_core_genes use |kME| >= {CORE_KME_THRESHOLD}; if a module had too few such genes, the top hub-ranked genes were added until up to {MIN_STABLE_CORE_GENES} genes were available.",
        "",
        "## Gene Set Types",
        "",
        "- all_module_genes: every non-grey gene assigned to the prioritized module in the best solution.",
        "- stable_core_genes: genes selected for future module score/eigengene features using kME/hubness stability rules.",
        f"- top_hub_genes: top {TOP_HUB_CAP} genes by abs(kME)/hub rank, or all genes if the module is smaller.",
        "",
        "## Module Summary",
        "",
        "| Region | Module | Module size | Stable core genes | Top hub genes | Median stable-core abs(kME) | ML #2 readiness |",
        "| --- | --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for _, row in summary.iterrows():
        readiness = "enough stable core genes" if row["n_stable_core_genes"] >= 25 else "small module; use cautiously"
        report_lines.append(
            f"| {row['region']} | {row['module']} | {int(row['module_size'])} | {int(row['n_stable_core_genes'])} | {int(row['n_top_hub_genes'])} | {row['median_kME_stable_core']:.3f} | {readiness} |"
        )
    report_lines.extend([
        "",
        "## How These Genes Can Be Used Later",
        "",
        "For ML #2, each module can be represented by module scores computed from stable_core_genes, such as an eigengene/first principal component, mean scaled expression, or robust median scaled expression. The all_module_genes set remains available for sensitivity analyses, while top_hub_genes can be used for compact interpretable feature variants.",
        "",
        "## Limitations",
        "",
        "- Run-level assignment stability could not be estimated because parameter-grid gene assignment tables were not saved.",
        "- kME values are available only for the saved top hub table, so genes outside the top hub table are retained in all_module_genes but are not selected into stable cores unless they appear in the hub/kME table.",
        "- Stable core gene selection does not imply diagnosis association and should not be interpreted as biomarker validation.",
        "- Future classification should be run in a separate step with strict cross-validation and no leakage from diagnosis labels into gene selection.",
        "",
        "## Outputs",
        "",
        "- stable_core_genes_by_module.csv",
        "- stable_core_gene_summary_by_module.csv",
        "- stable_core_gene_manifest.json",
    ])
    (OUT_DIR / "Stable_Core_Genes_for_Module_Level_ML2.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    manifest = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "project_root": str(PROJECT_ROOT),
        "inputs": {k: str(v) for k, v in INPUTS.items()},
        "outputs": [
            str(OUT_DIR / "stable_core_genes_by_module.csv"),
            str(OUT_DIR / "stable_core_gene_summary_by_module.csv"),
            str(OUT_DIR / "Stable_Core_Genes_for_Module_Level_ML2.md"),
            str(OUT_DIR / "stable_core_gene_manifest.json"),
        ],
        "focus_modules": [{"region": r, "module": m} for r, m in FOCUS_MODULES],
        "selection_rules": {
            "diagnosis_p_values_used_for_selection": False,
            "tier2_or_ml_status_used_for_selection": False,
            "stable_core_abs_kME_threshold": CORE_KME_THRESHOLD,
            "minimum_stable_core_fallback_target": MIN_STABLE_CORE_GENES,
            "top_hub_cap": TOP_HUB_CAP,
        },
        "stability_method": stability_method,
        "run_level_assignment_files_found": [str(p) for p in run_level_assignment_files],
        "limitation_note": limitation_note if not has_run_level_assignments else "Run-level assignment frequencies were available.",
    }
    (OUT_DIR / "stable_core_gene_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"Wrote {OUT_DIR / 'stable_core_genes_by_module.csv'}")
    print(f"Wrote {OUT_DIR / 'stable_core_gene_summary_by_module.csv'}")
    print(f"Wrote {OUT_DIR / 'Stable_Core_Genes_for_Module_Level_ML2.md'}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
