"""Shared helpers for the clean Length > 100 Tier2-only pipeline.

Each step below is intentionally lightweight: completed old outputs are reused,
copied, filtered, or summarized. The library avoids rerunning LMM, enrichment,
or broad ML analyses.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CONFIG_PATH = Path(__file__).with_name("00_pipeline_config.py")
spec = importlib.util.spec_from_file_location("clean_pipeline_config", CONFIG_PATH)
cfg = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(cfg)

MAIN_REGIONS = [
    "Ventral Striatum",
    "Insula",
    "Hippocampal Formation",
    "Cingulate gyrus 25",
    "dlPFC BA8/9",
    "OFC BA11",
]

STEP_NAMES = {
    1: "Prepare Length > 100 inputs",
    2: "Reuse GeTMM log2 normalized expression",
    3: "Reuse broad-region LMM outputs",
    4: "Reuse standardized effect-size outputs",
    5: "Create clean evidence tiers",
    6: "Reuse frontal subregion sensitivity outputs",
    7: "Package final Tier 2 candidate gene sets",
    8: "Reuse Tier 2 enrichment outputs",
    9: "ML v1 extended prioritization",
    10: "Extended ML-prioritized Tier 2 genes",
    11: "Collect report images",
    12: "Write integrated clean report and manifest",
}


def ensure(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def readme(path: Path, title: str, lines: list[str]) -> None:
    text = [f"# {title}", "", *lines, ""]
    path.write_text("\n".join(text), encoding="utf-8")


def copy_file(src: Path, dst: Path, skip_existing: bool = False) -> bool:
    ensure(dst.parent)
    if skip_existing and dst.exists():
        return False
    shutil.copy2(src, dst)
    return True


def copy_tree_files(src_dir: Path, dst_dir: Path, pattern: str = "*", skip_existing: bool = False) -> list[Path]:
    ensure(dst_dir)
    copied = []
    for src in src_dir.glob(pattern):
        if src.is_file():
            dst = dst_dir / src.name
            if copy_file(src, dst, skip_existing=skip_existing):
                copied.append(dst)
    return copied


def safe_ratio(value: object) -> float:
    text = str(value)
    if "/" in text:
        num, den = text.split("/", 1)
        try:
            return float(num) / float(den)
        except ValueError:
            return np.nan
    try:
        return float(text)
    except ValueError:
        return np.nan


def clean_region_from_file(name: str) -> str:
    region = name.replace("__Tier2_all.csv", "")
    region = region.replace("broad_", "")
    mapping = {
        "Hippocampal_Formation": "Hippocampal Formation",
        "Ventral_Striatum": "Ventral Striatum",
        "Cingulate_gyrus_25": "Cingulate gyrus 25",
        "dlPFC_BA8_9": "dlPFC BA8/9",
        "OFC_BA11": "OFC BA11",
        "Frontal_Cortex": "broad Frontal Cortex",
        "Insula": "Insula",
    }
    return mapping.get(region, region.replace("_", " "))


def step_output(step: int) -> Path:
    names = {
        1: "01_inputs",
        2: "02_normalized_expression",
        3: "03_broad_region_LMM",
        4: "04_standardized_effect_size",
        5: "05_evidence_tiers",
        6: "06_frontal_subregion_sensitivity",
        7: "07_Tier2_candidate_gene_sets",
        8: "08_Tier2_enrichment",
        9: "09_ML_v1_extended_prioritization",
        10: "10_ML_prioritized_Tier2_extended",
    }
    return cfg.CLEAN_OUTPUT_DIR / names[step]


def run_step_1(skip_existing: bool = False) -> dict:
    out = ensure(step_output(1))
    copy_file(cfg.GENE_FILTER_SUMMARY, out / "gene_filter_length_gt_100_summary.csv", skip_existing)
    metadata_candidates = list((cfg.OLD_ANALYSIS_DIR / "step1_from_raw" / "Length_gt_100" / "01_metadata").rglob("metadata_clean*.csv"))
    metadata_source = metadata_candidates[0] if metadata_candidates else cfg.FRONTAL_SUBREGION_DIR / "metadata_original_region_candidates.csv"
    if metadata_source.exists():
        copy_file(metadata_source, out / "sample_metadata_used.csv", skip_existing)
    lmm = pd.read_csv(cfg.LMM_ALL_GENES, usecols=["brain_region", "gene_id"])
    universe = lmm.groupby("brain_region")["gene_id"].nunique().reset_index(name="genes_tested")
    universe.to_csv(out / "tested_gene_universe_by_region.csv", index=False)
    readme(
        out / "README.md",
        "Step 1: Length > 100 Inputs",
        [
            "Reuses completed Length > 100 gene-filtering and metadata outputs.",
            f"Gene filter summary source: `{cfg.GENE_FILTER_SUMMARY}`",
            f"Sample metadata source: `{metadata_source}`",
            "The tested gene universe was summarized from the completed broad-region LMM table.",
        ],
    )
    return {"outputs": [str(out)], "reused": [str(cfg.GENE_FILTER_SUMMARY), str(metadata_source)]}


def run_step_2(skip_existing: bool = False) -> dict:
    out = ensure(step_output(2))
    copied = copy_tree_files(cfg.NORMALIZED_EXPRESSION_DIR, out, "*.csv", skip_existing)
    readme(
        out / "README.md",
        "Step 2: GeTMM + log2 Normalized Expression",
        [
            "Reuses completed normalized expression matrices.",
            "Normalization method: filtered raw counts -> GeTMM -> log2(GeTMM + 1).",
            f"Source folder: `{cfg.NORMALIZED_EXPRESSION_DIR}`",
            "These matrices are local pipeline outputs and are not required for GitHub review if file size is a concern.",
        ],
    )
    return {"outputs": [str(p) for p in copied], "reused": [str(cfg.NORMALIZED_EXPRESSION_DIR)]}


def run_step_3(skip_existing: bool = False) -> dict:
    out = ensure(step_output(3))
    copy_file(cfg.LMM_ALL_GENES, out / "LMM_all_genes_length_gt_100.csv", skip_existing)
    copy_file(cfg.LMM_SUMMARY_BY_REGION, out / "summary_by_region_length_gt_100.csv", skip_existing)
    readme(
        out / "README.md",
        "Step 3: Broad-Region LMM",
        [
            "Reuses the completed broad-region LMM results for Ventral Striatum, Insula, Hippocampal Formation, and broad Frontal Cortex.",
            "LMM means linear mixed model.",
            f"Source LMM table: `{cfg.LMM_ALL_GENES}`",
        ],
    )
    return {"outputs": [str(out)], "reused": [str(cfg.LMM_ALL_GENES), str(cfg.LMM_SUMMARY_BY_REGION)]}


def run_step_4(skip_existing: bool = False) -> dict:
    out = ensure(step_output(4))
    copy_file(cfg.STANDARDIZED_EFFECT_TABLE, out / "LMM_all_genes_with_standardized_effect_size_4groups.csv", skip_existing)
    if cfg.STANDARDIZED_EFFECT_SUMMARY.exists():
        copy_file(cfg.STANDARDIZED_EFFECT_SUMMARY, out / "standardized_effect_size_4groups_summary.csv", skip_existing)
    sign = pd.read_csv(cfg.STANDARDIZED_EFFECT_TABLE, usecols=["standardized_effect_sign_matches_beta"])
    summary = pd.DataFrame(
        {
            "rows_checked": [len(sign)],
            "sign_check_passed_rows": [int(sign["standardized_effect_sign_matches_beta"].fillna(False).sum())],
            "sign_check_failed_rows": [int((~sign["standardized_effect_sign_matches_beta"].fillna(False)).sum())],
            "definition": ["standardized_effect_size_diagnosis = beta_diagnosis / SD(expression of gene)"],
        }
    )
    summary.to_csv(out / "sign_check_summary.csv", index=False)
    readme(
        out / "README.md",
        "Step 4: Standardized Effect Sizes",
        [
            "Reuses completed standardized-effect-size outputs.",
            "Definition: standardized_effect_size_diagnosis = beta_diagnosis / SD(expression of gene).",
            "The sign-check summary confirms whether standardized effect signs match beta signs.",
        ],
    )
    return {"outputs": [str(out)], "reused": [str(cfg.STANDARDIZED_EFFECT_TABLE)]}


def run_step_5(skip_existing: bool = False) -> dict:
    out = ensure(step_output(5))
    all_genes = pd.read_csv(cfg.STANDARDIZED_EFFECT_TABLE)
    all_genes["Tier1"] = (all_genes["qvalue_Storey"] < 0.10) & (all_genes["abs_standardized_effect_size_diagnosis"] > 0.8)
    all_genes["Tier2"] = (all_genes["p_value_diagnosis"] < 0.01) & (all_genes["abs_standardized_effect_size_diagnosis"] > 0.8)
    all_genes["Tier3"] = (all_genes["p_value_diagnosis"] < 0.01) & (all_genes["abs_beta_diagnosis"] > 0.2)
    all_genes.to_csv(out / "evidence_tiers_all_genes.csv", index=False)
    all_genes[all_genes["Tier1"]].to_csv(out / "Tier1_genes.csv", index=False)
    all_genes[all_genes["Tier2"]].to_csv(out / "Tier2_genes_broad_regions.csv", index=False)
    all_genes[all_genes["Tier3"]].to_csv(out / "Tier3_genes.csv", index=False)
    summary = (
        all_genes.groupby("grouped_region")
        .agg(genes_tested=("gene_id", "nunique"), Tier1_genes=("Tier1", "sum"), Tier2_genes=("Tier2", "sum"), Tier3_genes=("Tier3", "sum"))
        .reset_index()
        .rename(columns={"grouped_region": "brain_region"})
    )
    summary.to_csv(out / "evidence_tier_summary_by_region.csv", index=False)
    readme(
        out / "README.md",
        "Step 5: Evidence Tiers",
        [
            "Tier 1: Storey q < 0.10 and abs(standardized_effect_size_diagnosis) > 0.8.",
            "Tier 2: p_value_diagnosis < 0.01 and abs(standardized_effect_size_diagnosis) > 0.8.",
            "Tier 3: p_value_diagnosis < 0.01 and abs(beta_diagnosis) > 0.2.",
            "Tier 4 is intentionally excluded from this clean Tier2-only pipeline.",
        ],
    )
    return {"outputs": [str(out)], "reused": [str(cfg.STANDARDIZED_EFFECT_TABLE)]}


def run_step_6(skip_existing: bool = False) -> dict:
    out = ensure(step_output(6))
    src_table = cfg.FRONTAL_SUBREGION_DIR / "frontal_subregion_LMM_all_genes_with_standardized_effect_size.csv"
    df = pd.read_csv(src_table)
    tier2 = df[df["Tier2_nominal_p001_absStdES08"]].copy()
    tier2.to_csv(out / "frontal_subregion_Tier2_genes.csv", index=False)
    summary_src = cfg.FRONTAL_SUBREGION_DIR / "frontal_subregion_summary.csv"
    if summary_src.exists():
        copy_file(summary_src, out / "frontal_subregion_summary.csv", skip_existing)
    feasibility = pd.DataFrame(
        {
            "model": ["expression_value ~ Diagnosis + Age + Sex + RIN"],
            "random_effect": ["BrainID random effect not estimable within subregion because each BrainID appears once"],
            "source": [str(src_table)],
        }
    )
    feasibility.to_csv(out / "frontal_subregion_model_feasibility.csv", index=False)
    readme(
        out / "README.md",
        "Step 6: Frontal Subregion Sensitivity",
        [
            "Frontal subregion candidates come from a fixed-effect sensitivity model, not the same grouped-region LMM.",
            "Documented model: Diagnosis + Age + Sex + RIN.",
            f"Source table: `{src_table}`",
        ],
    )
    return {"outputs": [str(out)], "reused": [str(src_table)]}


def run_step_7(skip_existing: bool = False) -> dict:
    out = ensure(step_output(7))
    symbols_dir = ensure(out / "Tier2_gene_symbols_by_region")
    rows = []
    broad_rows = []
    frontal_rows = []
    for path in cfg.TIER2_GENE_SET_DIR.glob("*__Tier2_all.csv"):
        region = clean_region_from_file(path.name)
        df = pd.read_csv(path)
        if region == "broad Frontal Cortex":
            source = "grouped_region_LMM_documentation_only"
        elif region in {"Ventral Striatum", "Insula", "Hippocampal Formation"}:
            source = "grouped_region_LMM"
        else:
            source = "frontal_subregion_fixed_effect_sensitivity"
        df["region_or_subregion"] = region
        df["model_source"] = source
        rows.append(df)
        if source == "grouped_region_LMM":
            broad_rows.append(df)
        if source == "frontal_subregion_fixed_effect_sensitivity":
            frontal_rows.append(df)
        if region in MAIN_REGIONS:
            df[["gene_symbol"]].dropna().drop_duplicates().to_csv(
                symbols_dir / f"{region.replace('/', '_').replace(' ', '_')}_gene_symbols.txt",
                index=False,
                header=False,
            )
    all_df = pd.concat(rows, ignore_index=True)
    main_df = all_df[all_df["region_or_subregion"].isin(MAIN_REGIONS)].copy()
    main_df.to_csv(out / "Tier2_candidate_genes_all_regions.csv", index=False)
    pd.concat(broad_rows, ignore_index=True).to_csv(out / "Tier2_candidate_genes_broad_regions.csv", index=False)
    pd.concat(frontal_rows, ignore_index=True).to_csv(out / "Tier2_candidate_genes_frontal_subregions.csv", index=False)
    summary = (
        main_df.groupby(["region_or_subregion", "model_source"])
        .agg(Tier2_genes=("gene_id", "nunique"), up_in_MDD=("direction", lambda s: (s == "up_in_MDD").sum()), down_in_MDD=("direction", lambda s: (s == "down_in_MDD").sum()))
        .reset_index()
    )
    summary.to_csv(out / "Tier2_summary_by_region.csv", index=False)
    readme(
        out / "README.md",
        "Step 7: Final Tier 2 Candidate Gene Sets",
        [
            "Combines broad-region LMM-derived Tier 2 genes with frontal-subregion fixed-effect sensitivity Tier 2 genes.",
            "Main candidate sets: Ventral Striatum, Insula, Hippocampal Formation, Cingulate gyrus 25, dlPFC BA8/9, and OFC BA11.",
            "Broad Frontal Cortex is documentation only because it has one Tier 2 gene.",
            "Tier 2 genes are exploratory candidate genes, not FDR-significant discoveries.",
        ],
    )
    return {"outputs": [str(out)], "reused": [str(cfg.TIER2_GENE_SET_DIR)]}


def run_step_8(skip_existing: bool = False) -> dict:
    out = ensure(step_output(8))
    tables_out = ensure(out / "tables")
    copied_tables = copy_tree_files(cfg.ENRICHMENT_DIR / "tables", tables_out, "*.csv", skip_existing)
    for name in ["enrichment_summary_by_gene_set.csv", "top_enrichment_terms_all_sets.csv", "README.md", "Tier2_Enrichment_Report.md"]:
        src = cfg.ENRICHMENT_DIR / name
        if src.exists():
            copy_file(src, out / name, skip_existing)
    plot_src = cfg.OLD_IMAGES_DIR / "enrichment_report_ready"
    if plot_src.exists():
        copy_tree_files(plot_src, ensure(out / "selected_report_ready_enrichment_plots"), "*.png", skip_existing)
        copy_file(plot_src / "README.md", out / "selected_report_ready_enrichment_plots" / "README.md", skip_existing)
    readme(
        out / "README_clean_pipeline.md",
        "Step 8: Tier 2 Enrichment",
        [
            "Reuses completed Tier 2 enrichment outputs.",
            "Enrichment used Tier 2 gene lists with all tested genes in the same region/subregion as background.",
            "Enrichment is exploratory biological interpretation.",
            "Broad Frontal Cortex was not interpreted because it had only one Tier 2 gene.",
        ],
    )
    return {"outputs": [str(out)], "reused": [str(cfg.ENRICHMENT_DIR)]}


def run_step_9(skip_existing: bool = False) -> dict:
    out = ensure(step_output(9))
    allowed = {"Tier2_fixed", "Tier4_top25", "Tier4_top50", "Tier4_top100", "Tier4_top250", "Tier4_top500", "strict_univariate_top100"}
    perf = pd.read_csv(cfg.ML_V1_DIR / "ML_performance_summary.csv")
    ext_perf = perf[perf["feature_set"].isin(allowed)].copy()
    ext_perf.to_csv(out / "ML_v1_extended_performance_summary.csv", index=False)
    best = ext_perf.sort_values(["region_or_subregion", "roc_auc_mean", "balanced_accuracy_mean"], ascending=[True, False, False]).groupby("region_or_subregion").head(1)
    best.to_csv(out / "ML_v1_extended_best_model_by_region.csv", index=False)
    feature_set_summary = (
        ext_perf.groupby(["feature_set", "analysis_level"])
        .agg(n_rows=("model", "size"), n_regions=("region_or_subregion", "nunique"), median_roc_auc=("roc_auc_mean", "median"), median_balanced_accuracy=("balanced_accuracy_mean", "median"))
        .reset_index()
    )
    feature_set_summary.to_csv(out / "ML_v1_extended_feature_set_summary.csv", index=False)
    fi = pd.read_csv(cfg.ML_V1_DIR / "feature_importance" / "feature_importance_summary.csv")
    ext_fi = fi[fi["feature_set"].isin(allowed)].copy()
    ext_fi.to_csv(out / "ML_v1_extended_feature_importance.csv", index=False)
    plots = ensure(out / "plots")
    old_plot_dir = cfg.ML_V1_DIR / "plots"
    for src in old_plot_dir.glob("*.png"):
        if src.is_file():
            copy_file(src, plots / src.name, skip_existing)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    best_plot = best.sort_values("roc_auc_mean")
    labels = best_plot["region_or_subregion"]
    axes[0].barh(labels, best_plot["balanced_accuracy_mean"], color="#4878a8")
    axes[0].set_title("Balanced accuracy")
    axes[1].barh(labels, best_plot["roc_auc_mean"], color="#d95f02")
    axes[1].set_title("ROC-AUC")
    for ax in axes:
        ax.axvline(0.5, color="#555", linestyle="--", linewidth=1)
        ax.set_xlim(0.45, 1.02)
        ax.grid(axis="x", color="#e5e8eb")
    fig.suptitle("Extended ML v1: exploratory performance")
    fig.tight_layout()
    fig.savefig(plots / "MLv1_extended_performance_summary_exploratory.png", dpi=200)
    plt.close(fig)
    readme(
        out / "README.md",
        "Step 9: Extended ML v1 Prioritization",
        [
            "This clean branch reuses existing ML v1 outputs and keeps Tier2_fixed plus Tier4_top-k feature sets.",
            "Tier 2 genes define the main exploratory candidate pool. Tier 4 ranked signatures are used only as an additional ML-support layer to test whether candidate-gene prioritization is stable across feature-set definitions.",
            "strict_univariate_top100 is retained only as a reference/sanity branch, not as the main candidate source.",
            "ML v1 is interpreted as gene prioritization, not validated diagnostic prediction. Classification performance is descriptive because feature sets were derived from the same dataset.",
        ],
    )
    archive = cfg.CLEAN_OUTPUT_DIR / "10_ML_prioritized_Tier2FixedOnly" / "README.md"
    if archive.exists():
        with archive.open("a", encoding="utf-8") as handle:
            handle.write("\nArchive note: Tier2-only sensitivity analysis, not the main clean workflow result.\n")
    return {"outputs": [str(out)], "reused": [str(cfg.ML_V1_DIR)], "ml_rerun": False}


def compute_lmm_evidence_score(p: pd.Series, effect: pd.Series) -> pd.Series:
    p_clip = p.clip(lower=np.nextafter(0, 1), upper=1)
    return (-np.log10(p_clip)) * effect.abs()


def plot_extended_ml_prioritized(ranked: pd.DataFrame, summary: pd.DataFrame, plots: Path) -> None:
    """Create clean report-ready plots for the restored extended ML strategy."""
    fig, ax = plt.subplots(figsize=(8, 4.8))
    s = summary.sort_values("n_high_confidence")
    ax.barh(s["region_or_subregion"], s["n_moderate_confidence"], color="#9bc2d1", label="Moderate")
    ax.barh(s["region_or_subregion"], s["n_high_confidence"], color="#2b6f93", label="High")
    ax.set_xlabel("ML-prioritized Tier 2 genes")
    ax.set_title("Extended ML-prioritized Tier 2 genes")
    ax.legend()
    ax.grid(axis="x", color="#e5e8eb")
    fig.tight_layout()
    fig.savefig(plots / "MLprioritized_counts_by_region.png", dpi=220)
    plt.close(fig)

    top = ranked.groupby("region_or_subregion").head(8).copy()
    top = top.sort_values("ML_priority_score")
    fig, ax = plt.subplots(figsize=(9.5, max(5.5, len(top) * 0.22)))
    labels = top["region_or_subregion"] + " | " + top["gene_symbol"].fillna(top["gene_id"].astype(str))
    ax.barh(labels, top["ML_priority_score"], color="#2b6f93")
    ax.set_xlabel("ML priority score")
    ax.set_title("Top extended ML-prioritized Tier 2 genes")
    ax.grid(axis="x", color="#e5e8eb")
    fig.tight_layout()
    fig.savefig(plots / "MLprioritized_top_genes_by_region.png", dpi=220)
    plt.close(fig)

    class_colors = {
        "High-confidence ML-prioritized Tier 2": "#1b6e8f",
        "Moderate-confidence ML-prioritized Tier 2": "#d95f02",
    }
    for xcol, filename, xlabel in [
        ("standardized_effect_size_diagnosis", "MLpriority_vs_standardized_effect_size.png", "Standardized effect size"),
        ("p_value_diagnosis", "MLpriority_vs_neglog10p.png", "-log10(p-value)"),
    ]:
        fig, ax = plt.subplots(figsize=(7.7, 5.2))
        x = -np.log10(ranked[xcol].clip(lower=np.nextafter(0, 1))) if xcol == "p_value_diagnosis" else ranked[xcol]
        for cls, sub in ranked.assign(_x=x).groupby("ML_priority_class"):
            ax.scatter(sub["_x"], sub["ML_priority_score"], s=30, alpha=0.75, color=class_colors.get(cls, "#7d8790"), label=cls.replace(" ML-prioritized Tier 2", ""))
        ax.set_xlabel(xlabel)
        ax.set_ylabel("ML priority score")
        ax.set_title("Extended ML support for Tier 2 genes")
        ax.grid(color="#e5e8eb")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(plots / filename, dpi=220)
        plt.close(fig)

    support_cols = ["sparse_model_support", "linear_model_support", "rf_support", "v2_support_flag", "lmm_ml_sign_agreement"]
    heat = ranked.head(30).copy()
    fig, ax = plt.subplots(figsize=(8.0, 8.8))
    ax.imshow(heat[support_cols].fillna(False).astype(float), aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_yticks(range(len(heat)), heat["gene_symbol"].fillna(heat["gene_id"].astype(str)))
    ax.set_xticks(range(len(support_cols)), ["Sparse\nmodels", "Linear\nmodels", "Random\nforest", "Strict v2\nsupport", "LMM/ML\nsign agree"])
    ax.set_title("Model support for top extended ML-prioritized Tier 2 genes")
    fig.tight_layout()
    fig.savefig(plots / "MLpriority_model_support_heatmap.png", dpi=220)
    plt.close(fig)

    agreement = ranked.groupby(["region_or_subregion", "lmm_ml_sign_agreement"]).size().unstack(fill_value=0)
    agreement = agreement.rename(columns={0.0: "Disagree", 1.0: "Agree", False: "Disagree", True: "Agree"}).reset_index()
    fig, ax = plt.subplots(figsize=(8, 4.8))
    agree = agreement.get("Agree", pd.Series([0] * len(agreement)))
    disagree = agreement.get("Disagree", pd.Series([0] * len(agreement)))
    ax.barh(agreement["region_or_subregion"], disagree, color="#c55a4d", label="Disagree")
    ax.barh(agreement["region_or_subregion"], agree, left=disagree, color="#5c8f55", label="Agree")
    ax.set_xlabel("ML-prioritized Tier 2 genes")
    ax.set_title("LMM and ML direction agreement")
    ax.legend()
    fig.tight_layout()
    fig.savefig(plots / "MLpriority_sign_agreement.png", dpi=220)
    plt.close(fig)


def run_step_10(skip_existing: bool = False) -> dict:
    out = ensure(step_output(10))
    stable_script = cfg.PROJECT_ROOT / "code" / "extract_stable_ML_prioritized_Tier2_genes.py"
    spec = importlib.util.spec_from_file_location("stable_ml_prioritized_tier2", stable_script)
    stable = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(stable)
    stable.OUT_DIR = out
    stable.PLOT_DIR = out / "plots"
    stable.GENE_LIST_DIR = out / "gene_lists"
    stable.main()

    copy_file(out / "ML_prioritized_Tier2_genes_all_ranked.csv", out / "ML_prioritized_Tier2_extended_all_ranked.csv", skip_existing=False)
    copy_file(out / "ML_prioritized_Tier2_genes_high_confidence.csv", out / "ML_prioritized_Tier2_extended_high_confidence.csv", skip_existing=False)
    copy_file(out / "ML_prioritized_Tier2_genes_moderate_confidence.csv", out / "ML_prioritized_Tier2_extended_moderate_confidence.csv", skip_existing=False)
    copy_file(out / "ML_prioritized_Tier2_gene_summary_by_region.csv", out / "ML_prioritized_Tier2_extended_summary_by_region.csv", skip_existing=False)
    heatmap = out / "plots" / "heatmap_top_stable_ML_prioritized_genes.png"

    ranked = pd.read_csv(out / "ML_prioritized_Tier2_extended_all_ranked.csv")
    high_df = pd.read_csv(out / "ML_prioritized_Tier2_extended_high_confidence.csv")
    moderate_df = pd.read_csv(out / "ML_prioritized_Tier2_extended_moderate_confidence.csv")
    summary = pd.read_csv(out / "ML_prioritized_Tier2_extended_summary_by_region.csv")
    for region, sub in pd.concat([high_df, moderate_df], ignore_index=True).groupby("region_or_subregion"):
        sub[["gene_symbol"]].dropna().drop_duplicates().to_csv(
            out / "gene_lists" / f"{region.replace('/', '_').replace(' ', '_')}_prioritized_gene_symbols.txt",
            index=False,
            header=False,
        )
    plot_extended_ml_prioritized(ranked, summary, out / "plots")
    if heatmap.exists():
        copy_file(heatmap, out / "plots" / "MLpriority_model_support_heatmap.png", skip_existing=False)
    cg25_plot = out / "plots" / "Cingulate_gyrus_25_top20_ML_priority_genes.png"
    if cg25_plot.exists():
        copy_file(cg25_plot, out / "plots" / "MLprioritized_top_genes_Cg25.png", skip_existing=False)
    manifest_path = out / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else {}
    manifest.update(
        {
            "clean_pipeline_step": "10_extract_extended_ML_prioritized_Tier2_genes",
            "source_procedure": str(stable_script),
            "source_strategy": "stable ML-prioritized Tier 2 procedure using Tier2_fixed, Tier4_top-k support, and strict v2 support",
            "ml_rerun": False,
            "ranked_count": int(len(ranked)),
            "high_confidence_count": int(len(high_df)),
            "moderate_confidence_count": int(len(moderate_df)),
            "key_heatmap": str(heatmap),
            "tier4_interpretation": "Tier 4 ranked signatures are support/sensitivity only, not final gene lists.",
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    readme(
        out / "README.md",
        "Step 10: Extended ML-Prioritized Tier 2 Genes",
        [
            "This step reuses the previous stable ML-prioritized Tier 2 output package.",
            "Tier 2 genes define the main exploratory candidate pool. Tier 4 ranked signatures are used only as an additional ML-support layer to test whether candidate-gene prioritization is stable across feature-set definitions.",
            "Tier 4 top-k feature sets are selected by absolute signed standardized log-p score, while the exported tables keep signed_score, abs_score, and direction so up-in-MDD and down-in-MDD genes remain distinguishable.",
            "The final candidates are ML-prioritized Tier 2 genes, not Tier 4 genes.",
            "The final prioritized genes should be interpreted as internally supported candidate genes, not validated biomarkers.",
        ],
    )
    return {"outputs": [str(out)], "reused": [str(stable_script), str(cfg.ML_V1_DIR), str(cfg.OLD_ANALYSIS_DIR / "ML_strict_v2")], "high_count": int(len(high_df)), "moderate_count": int(len(moderate_df)), "ranked_count": int(len(ranked))}


def plot_ml_prioritized(priority: pd.DataFrame, summary: pd.DataFrame, plots: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.8))
    s = summary.sort_values("n_high_confidence")
    ax.barh(s["region_or_subregion"], s["n_moderate_confidence"], color="#9bc2d1", label="Moderate")
    ax.barh(s["region_or_subregion"], s["n_high_confidence"], color="#2b6f93", label="High")
    ax.set_xlabel("Tier2FixedOnly ML-prioritized genes")
    ax.set_title("Tier2FixedOnly ML-prioritized genes by area")
    ax.legend()
    fig.tight_layout()
    fig.savefig(plots / "MLprioritized_Tier2FixedOnly_counts_by_region.png", dpi=200)
    plt.close(fig)
    top = priority[priority["ML_priority_class_Tier2FixedOnly"].ne("Not prioritized")].groupby("region_or_subregion").head(8)
    if len(top):
        top = top.sort_values("ML_priority_score_Tier2FixedOnly")
        fig, ax = plt.subplots(figsize=(9.5, max(5, len(top) * 0.22)))
        labels = top["region_or_subregion"] + " | " + top["gene_symbol"].fillna(top["gene_id"])
        ax.barh(labels, top["ML_priority_score_Tier2FixedOnly"], color="#2b6f93")
        ax.set_xlabel("ML priority score")
        ax.set_title("Top Tier2FixedOnly ML-prioritized genes")
        fig.tight_layout()
        fig.savefig(plots / "MLprioritized_Tier2FixedOnly_top_genes_by_region.png", dpi=200)
        plt.close(fig)
    for x, name, xlabel in [
        ("standardized_effect_size_diagnosis", "MLpriority_Tier2FixedOnly_vs_standardized_effect_size.png", "Standardized effect size"),
        ("p_value_diagnosis", "MLpriority_Tier2FixedOnly_vs_neglog10p.png", "-log10(p-value)"),
    ]:
        fig, ax = plt.subplots(figsize=(7.5, 5))
        xx = -np.log10(priority[x].clip(lower=np.nextafter(0, 1))) if x == "p_value_diagnosis" else priority[x]
        ax.scatter(xx, priority["ML_priority_score_Tier2FixedOnly"], s=25, alpha=0.75, color="#2b6f93")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Tier2FixedOnly ML priority score")
        ax.set_title("Tier2FixedOnly ML priority")
        ax.grid(color="#e5e8eb")
        fig.tight_layout()
        fig.savefig(plots / name, dpi=200)
        plt.close(fig)
    support_cols = ["sparse_model_support", "linear_model_support", "random_forest_support", "lmm_ml_sign_agreement"]
    heat = priority.head(30)
    fig, ax = plt.subplots(figsize=(7.5, 8.5))
    ax.imshow(heat[support_cols].fillna(False).astype(float), aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_yticks(range(len(heat)), heat["gene_symbol"].fillna(heat["gene_id"]))
    ax.set_xticks(range(len(support_cols)), ["Sparse", "Linear", "RF", "Sign\nagree"])
    ax.set_title("Top Tier2FixedOnly model support")
    fig.tight_layout()
    fig.savefig(plots / "MLpriority_Tier2FixedOnly_model_support_heatmap.png", dpi=200)
    plt.close(fig)


def run_step_11(skip_existing: bool = False) -> dict:
    out = ensure(cfg.CLEAN_IMAGES_DIR)
    inventory = []
    archived_ml_image_dir = out / "ml_prioritized"
    if archived_ml_image_dir.exists():
        for old_plot in archived_ml_image_dir.glob("*Tier2FixedOnly*"):
            old_plot.unlink()
    mappings = [
        ("tier2", cfg.OLD_IMAGES_DIR / "tier2", ["Tier2_counts_by_region.png", "Tier2_up_down_by_region.png", "Tier2_volcano_Cg25.png", "Tier2_volcano_VentralStriatum.png"]),
        ("enrichment", cfg.OLD_IMAGES_DIR / "enrichment_report_ready", ["*.png"]),
        ("ml_v1", step_output(9) / "plots", ["MLv1_extended_performance_summary_exploratory.png"]),
        ("ml_prioritized", step_output(10) / "plots", ["*.png"]),
    ]
    for category, src_dir, patterns in mappings:
        dst_dir = ensure(out / category)
        for pattern in patterns:
            for src in src_dir.glob(pattern):
                dst = dst_dir / src.name
                copy_file(src, dst, skip_existing)
                inventory.append({"category": category, "filename": dst.name, "relative_path": str(dst.relative_to(cfg.CLEAN_OUTPUT_DIR)).replace("\\", "/"), "source_path": str(src)})
    inv_dir = ensure(out / "inventory")
    pd.DataFrame(inventory).to_csv(inv_dir / "image_inventory.csv", index=False)
    readme(
        inv_dir / "README.md",
        "Clean Pipeline Image Inventory",
        [
            "Collects Tier 2, report-ready enrichment, extended ML v1, and extended ML-prioritized Tier 2 figures.",
            "All ML figures are descriptive and exploratory.",
        ],
    )
    return {"outputs": [str(out)], "reused": [str(cfg.OLD_IMAGES_DIR), str(step_output(9) / "plots"), str(step_output(10) / "plots")]}


def run_step_12(skip_existing: bool = False) -> dict:
    ensure(cfg.CLEAN_OUTPUT_DIR)
    ensure(cfg.FINAL_REPORT.parent)
    tier2_summary = pd.read_csv(step_output(7) / "Tier2_summary_by_region.csv")
    ml_summary = pd.read_csv(step_output(10) / "ML_prioritized_Tier2_extended_summary_by_region.csv")
    ranked_count = len(pd.read_csv(step_output(10) / "ML_prioritized_Tier2_extended_all_ranked.csv"))
    high_count = len(pd.read_csv(step_output(10) / "ML_prioritized_Tier2_extended_high_confidence.csv"))
    moderate_count = len(pd.read_csv(step_output(10) / "ML_prioritized_Tier2_extended_moderate_confidence.csv"))
    report = build_report(tier2_summary, ml_summary)
    cfg.FINAL_REPORT.write_text(report, encoding="utf-8")
    cfg.FINAL_REPORT_COPY.write_text(report, encoding="utf-8")
    manifest = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "scripts_used": [str(p) for p in sorted(Path(__file__).parent.glob("*.py"))],
        "input_files_used": [
            str(cfg.LMM_ALL_GENES),
            str(cfg.STANDARDIZED_EFFECT_TABLE),
            str(cfg.FRONTAL_SUBREGION_DIR / "frontal_subregion_LMM_all_genes_with_standardized_effect_size.csv"),
            str(cfg.ML_V1_DIR / "feature_importance" / "feature_importance_summary.csv"),
            str(cfg.ENRICHMENT_DIR),
        ],
        "reused_old_outputs": [str(cfg.OLD_ANALYSIS_DIR), str(cfg.ML_V1_DIR), str(cfg.OLD_ML_PRIORITIZED_DIR)],
        "newly_generated_outputs": [str(cfg.CLEAN_OUTPUT_DIR), str(cfg.FINAL_REPORT)],
        "skipped_steps": [],
        "ml_rerun": False,
        "enrichment_rerun": False,
        "main_ml_strategy": "extended ML v1 using Tier2_fixed plus Tier4_top-k feature sets as support",
        "tier4_role": "Tier 4 is used only as a ranked ML-support layer, not as a separate biological conclusion or final gene list.",
        "tier2_fixed_only_outputs_retained_as_sensitivity": str(cfg.CLEAN_OUTPUT_DIR / "10_ML_prioritized_Tier2FixedOnly"),
        "final_ranked_count": int(ranked_count),
        "final_high_confidence_count": int(high_count),
        "final_moderate_confidence_count": int(moderate_count),
        "final_high_moderate_gene_counts": ml_summary.to_dict(orient="records"),
        "final_report": str(cfg.FINAL_REPORT),
    }
    cfg.PIPELINE_MANIFEST.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return {"outputs": [str(cfg.FINAL_REPORT), str(cfg.FINAL_REPORT_COPY), str(cfg.PIPELINE_MANIFEST)], "reused": []}


def build_report(tier2_summary: pd.DataFrame, ml_summary: pd.DataFrame) -> str:
    tier2_md = tier2_summary.to_markdown(index=False)
    ml_md = ml_summary.to_markdown(index=False)
    return f"""# Length > 100 MDD transcriptomics: Tier 2 candidate genes, enrichment, and Tier2-only ML prioritization

## 1. Executive summary
This clean pipeline organizes the completed Length > 100 analysis into a reproducible path from filtered/normalized gene expression through LMM-derived Tier 2 candidate genes, enrichment, extended ML v1 prioritization, and stable ML-prioritized Tier 2 genes.

Tier 2 genes are exploratory candidate genes, not FDR-significant discoveries. Tier 2 genes define the main exploratory candidate pool. Tier 4 ranked signatures are used only as an additional ML-support layer to test whether candidate-gene prioritization is stable across feature-set definitions.

ML v1 is interpreted as gene prioritization, not validated diagnostic prediction. Classification performance is descriptive because feature sets were derived from the same dataset. The final prioritized genes should be interpreted as internally supported candidate genes, not validated biomarkers.

## 2. Clean pipeline overview
The clean path is: raw counts/expression inputs -> Length > 100 gene filtering -> GeTMM -> log2(GeTMM + 1) -> broad-region LMM -> standardized effect sizes -> evidence tiers -> frontal subregion fixed-effect sensitivity -> final Tier 2 candidate pools -> Tier 2 enrichment -> extended ML v1 using Tier2_fixed plus Tier4_top-k feature sets -> stable ML-prioritized Tier 2 genes -> report-ready images.

## 3. Input data and normalization
The clean pipeline reuses completed Length > 100 outputs from `{cfg.OLD_ANALYSIS_DIR}`. Normalized expression uses filtered raw counts, GeTMM normalization, and log2(GeTMM + 1).

## 4. Broad-region LMM
Broad-region Tier 2 candidates come from the grouped-region LMM for Ventral Striatum, Insula, Hippocampal Formation, and broad Frontal Cortex.

## 5. Standardized effect sizes
Standardized effect size was defined as `beta_diagnosis / SD(expression of gene)`. The standardized-effect branch was used to define Tier 2 genes.

## 6. Evidence tiers
Tier 1: Storey q < 0.10 and abs(standardized effect size) > 0.8.
Tier 2: p_value_diagnosis < 0.01 and abs(standardized effect size) > 0.8.
Tier 3: p_value_diagnosis < 0.01 and abs(beta_diagnosis) > 0.2.
Tier 4 is a ranking of all tested genes, not a separate selected gene list or biological conclusion. It helps avoid relying on one hard cutoff by providing ranked-feature support in the ML prioritization stage.

## 7. Frontal subregion fixed-effect sensitivity analysis
Frontal subregion Tier 2 candidates come from the frontal fixed-effect sensitivity model because the random effect was not estimable after splitting. The documented model is Diagnosis + Age + Sex + RIN.

## 8. Final Tier 2 candidate gene pools
{tier2_md}

## 9. Tier 2 enrichment
Tier 2 enrichment was reused from completed outputs. It used Tier 2 gene lists and region/subregion-specific tested-gene backgrounds. Enrichment is exploratory biological interpretation.

## 10. ML v1 extended prioritization
In ML, samples are rows and genes are features. ML v1 used Tier2_fixed as the main candidate set and Tier4_top-k ranked signatures as additional support/sensitivity. Tier 4 is a ranking of all tested genes, not a separate selected gene list. Training/testing was by cross-validation within the same dataset. Every sample can be held out as test in some fold, but the test data are not external.

Tier 2 defines the main exploratory candidate-gene pool using nominal diagnosis evidence and standardized effect size. Tier 4 is a ranked signature of all tested genes and is used only as an additional ML-support layer. In the extended ML v1 prioritization, models are trained using Tier2_fixed and Tier4_top-k feature sets. The final prioritization remains Tier 2-centered: genes are reported as ML-prioritized Tier 2 candidates when they belong to the Tier 2 pool and receive repeated support from the ML feature-importance outputs.

Tier2FixedOnly was tried as a simplified sensitivity analysis but produced less informative prioritization plots, so the extended ML v1 strategy is used as the main workflow.

## 11. Stable ML-prioritized Tier 2 genes
The clean prioritization reuses the previous extended ML v1 strategy. The final candidates are ML-prioritized Tier 2 genes. Tier 4 only provides additional ranked-feature support.

{ml_md}

## 12. Recommended figures
Recommended Tier 2 figures are in `outputs/length_gt_100_clean_pipeline/images/tier2`.
Recommended enrichment figures are in `outputs/length_gt_100_clean_pipeline/images/enrichment`.
Recommended ML figures are in `outputs/length_gt_100_clean_pipeline/images/ml_v1` and `outputs/length_gt_100_clean_pipeline/images/ml_prioritized`, and come from the extended ML v1 prioritization branch.

## 13. Conclusions
The clean path preserves the exploratory Tier 2 biological signal while using Tier 4 ranked signatures only as ML support for stable prioritization.

## 14. Limitations
Tier 2 selection used the full dataset, so downstream ML performance is optimistic. There is no external validation cohort in this clean pipeline. Candidate genes are not validated biomarkers.

## 15. Next steps
Use extended ML-prioritized Tier 2 genes for enrichment follow-up, AHBA/module mapping, and lower-dimensional pathway or module-level ML. Do not start deep learning until stable lower-dimensional features or external validation data are available.

## 16. Output files
Primary clean outputs are under `{cfg.CLEAN_OUTPUT_DIR}`. The reproducibility manifest is `{cfg.PIPELINE_MANIFEST}`.

## 17. Glossary
LMM: linear mixed model.
Tier 2: nominal exploratory candidate genes with p_value_diagnosis < 0.01 and abs standardized effect size > 0.8.
ML v1: exploratory classical ML using Tier 2 expression features.
"""


STEP_FUNCS: dict[int, Callable[[bool], dict]] = {
    1: run_step_1,
    2: run_step_2,
    3: run_step_3,
    4: run_step_4,
    5: run_step_5,
    6: run_step_6,
    7: run_step_7,
    8: run_step_8,
    9: run_step_9,
    10: run_step_10,
    11: run_step_11,
    12: run_step_12,
}


def run_numbered_step(step: int, skip_existing: bool = False) -> dict:
    print(f"STEP {step}: {STEP_NAMES[step]}")
    result = STEP_FUNCS[step](skip_existing)
    print(f"Completed step {step}: {STEP_NAMES[step]}")
    return result


def runner_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Run the clean Length > 100 Tier2-only pipeline.")
    parser.add_argument("--from-step", type=int, default=1)
    parser.add_argument("--to-step", type=int, default=12)
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    selected = [s for s in sorted(STEP_FUNCS) if args.from_step <= s <= args.to_step]
    print("Clean Length > 100 MDD pipeline steps:")
    for s in selected:
        print(f"  {s:02d}. {STEP_NAMES[s]}")
    if args.dry_run:
        print("Dry run only; no steps executed.")
        return
    results = {}
    for s in selected:
        results[s] = run_numbered_step(s, skip_existing=args.skip_existing)
    print("Pipeline run complete.")
    print(f"Clean output folder: {cfg.CLEAN_OUTPUT_DIR}")
