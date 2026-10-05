#!/usr/bin/env python3
"""Create Tier2 ML–RRA diagnostics, manifest finalization, report, and checksums."""

from __future__ import annotations

import argparse
import json
import sys
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "code" / "ml_rra"))
from common import MODEL_ORDER, atomic_write_csv, atomic_write_json, sha256_file  # noqa: E402


DUMMY_MODEL = "dummy_most_frequent"


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "outputs/length_gt_100_clean_pipeline/18_Tier2_ML_RRA_prioritization",
    )
    return parser.parse_args()


def markdown_table(frame: pd.DataFrame) -> str:
    values = frame.fillna("").astype(str)
    header = "| " + " | ".join(values.columns) + " |"
    rule = "| " + " | ".join(["---"] * len(values.columns)) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in values.to_numpy()]
    return "\n".join([header, rule, *rows])


def save_figure(fig: plt.Figure, base: Path) -> None:
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = arguments()
    output = args.output_dir.resolve()
    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    performance = pd.read_csv(output / "Tier2_outer_performance.csv")
    summary = pd.read_csv(output / "Tier2_model_performance_summary.csv")
    rankings = pd.read_csv(output / "Tier2_fold_gene_rankings.csv", dtype={"gene_id": str}, low_memory=False)
    level1 = pd.read_csv(output / "Tier2_model_level_rra_rankings.csv", dtype={"gene_id": str})
    consensus = pd.read_csv(output / "Tier2_consensus_rra_rankings.csv", dtype={"gene_id": str})
    accounting = pd.read_csv(output / "Tier2_candidate_accounting.csv")
    comparison_summary = pd.read_csv(output / "Tier2_vs_allgene_RRA_summary.csv")
    manifest_path = output / "Tier2_run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    consensus["rra_score_interpretation"] = "rank aggregation score; not an MDD-association p-value"
    atomic_write_csv(consensus, output / "Tier2_consensus_rra_rankings.csv")
    sns.set_theme(style="whitegrid", context="talk")

    model_palette = dict(zip(MODEL_ORDER + [DUMMY_MODEL], sns.color_palette("deep", 6)))
    fig, ax = plt.subplots(figsize=(13, 7))
    sns.boxplot(
        data=performance,
        x="region",
        y="balanced_accuracy",
        hue="model",
        hue_order=MODEL_ORDER + [DUMMY_MODEL],
        palette=model_palette,
        ax=ax,
    )
    ax.axhline(0.5, color="black", linestyle="--", linewidth=1)
    ax.set(
        title="Descriptive internal performance for Tier2-restricted models",
        xlabel="",
        ylabel="Balanced accuracy",
    )
    ax.tick_params(axis="x", rotation=28)
    ax.legend(title="Model", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    save_figure(fig, figures / "Tier2_descriptive_classifier_performance")

    top = consensus.sort_values(["region", "consensus_rank"]).groupby("region", sort=False).head(25).copy()
    rank_columns = [f"{model}_rank" for model in MODEL_ORDER]
    matrix = top.set_index(["region", "gene_symbol", "gene_id"])[rank_columns].astype(float)
    matrix = -np.log10(matrix.clip(lower=1))
    fig, ax = plt.subplots(figsize=(10, max(10, len(matrix) * 0.13)))
    sns.heatmap(matrix, cmap="viridis", ax=ax, cbar_kws={"label": "-log10(model-level rank)"})
    ax.set(
        title="Top 25 Tier2-RRA consensus genes per region across five models",
        xlabel="",
        ylabel="Region / Tier2 gene",
    )
    save_figure(fig, figures / "Tier2_top_genes_model_rank_heatmap")

    recurrence = level1.copy()
    recurrence["negative_log10_rra"] = -np.log10(recurrence["rra_score"].clip(lower=np.finfo(float).tiny))
    fig, axes = plt.subplots(2, 3, figsize=(17, 10), sharex=True, sharey=True)
    for ax, (region, part) in zip(axes.flat, recurrence.groupby("region", sort=False)):
        sns.scatterplot(
            data=part,
            x="outer_list_appearance_proportion",
            y="negative_log10_rra",
            hue="model",
            hue_order=MODEL_ORDER,
            s=28,
            alpha=0.65,
            ax=ax,
            legend=False,
        )
        ax.set_title(region)
        ax.set_xlabel("")
        ax.set_ylabel("")
    fig.suptitle("Tier2 gene recurrence across outer fits versus model-level RRA", y=1.01)
    fig.supxlabel("Outer-fit ranking appearance proportion")
    fig.supylabel("-log10(RRA score; ranking score, not MDD p-value)")
    save_figure(fig, figures / "Tier2_outer_fit_recurrence_vs_rra")

    agreement_rows: list[dict] = []
    for region, part in consensus.groupby("region", sort=False):
        for left, right in combinations(MODEL_ORDER, 2):
            pair = part[[f"{left}_rank", f"{right}_rank"]].dropna()
            rho = float(spearmanr(pair.iloc[:, 0], pair.iloc[:, 1]).statistic) if len(pair) >= 3 else np.nan
            agreement_rows.append(
                {"region": region, "model_pair": f"{left} | {right}", "n_common": len(pair), "spearman_rho": rho}
            )
    agreement = pd.DataFrame(agreement_rows)
    atomic_write_csv(agreement, output / "Tier2_model_rank_agreement.csv")
    agreement_matrix = agreement.pivot(index="region", columns="model_pair", values="spearman_rho")
    fig, ax = plt.subplots(figsize=(14, 5))
    sns.heatmap(agreement_matrix, cmap="vlag", center=0, vmin=-1, vmax=1, annot=True, fmt=".2f", ax=ax)
    ax.set(title="Spearman agreement among Tier2 model-level ranks", xlabel="Model pair", ylabel="")
    save_figure(fig, figures / "Tier2_model_rank_agreement")

    manifest.update(
        {
            "status": "complete",
            "rra": {
                "package": "RobustRankAggreg",
                "version": "1.2.1",
                "function": "RobustRankAggreg::aggregateRanks",
                "method": "RRA",
                "level1_lists_per_region_model": 25,
                "level2_lists_per_region": 5,
                "dummy_contributes_gene_evidence": False,
                "predictive_performance_weighting": False,
            },
            "result_counts": {
                "tier2_candidates": int(accounting["total_tier2_candidates"].sum()),
                "fold_ranking_rows": len(rankings),
                "model_level_rra_rows": len(level1),
                "consensus_rra_rows": len(consensus),
                "candidates_absent_from_every_fold_ranking": int(
                    accounting["candidates_absent_from_every_fold_ranking"].sum()
                ),
            },
            "frozen_allgene_comparison": {
                "role": "within-study convergence/sensitivity comparison; not independent validation",
                "allgene_outputs_modified": False,
            },
        }
    )
    atomic_write_json(manifest, manifest_path)

    sample_table = pd.DataFrame(manifest["analysis_units"])[
        ["region", "n_samples", "n_mdd", "n_control", "n_tier2_candidates"]
    ]
    sample_table.columns = ["Region", "Samples", "MDD", "Control", "Tier2 genes"]
    performance_report = summary.copy()
    for column in [c for c in performance_report if c.endswith("_mean") or c.endswith("_sd")]:
        performance_report[column] = performance_report[column].map(lambda value: f"{value:.3f}")
    accounting_report = accounting.copy()
    convergence_report = comparison_summary.copy()
    for column in [
        "proportion_present_in_allgene_rra",
        "spearman_rho_common_genes",
        "median_allgene_rank_percentile",
    ]:
        convergence_report[column] = convergence_report[column].map(lambda value: "" if pd.isna(value) else f"{value:.3f}")
    top_sections = []
    for region, part in consensus.groupby("region", sort=False):
        region_top = part.nsmallest(20, "consensus_rank")[
            ["consensus_rank", "gene_id", "gene_symbol", "consensus_rra_score", "model_appearance_count", *rank_columns]
        ].copy()
        region_top["consensus_rra_score"] = region_top["consensus_rra_score"].map(lambda value: f"{value:.4g}")
        for column in rank_columns:
            region_top[column] = region_top[column].map(lambda value: "" if pd.isna(value) else str(int(value)))
        top_sections.append(f"### {region}\n\n{markdown_table(region_top)}")

    report = f"""# Primary Tier2 → ML → RRA gene-prioritisation report

## A. Scientific objective

The gene-wise diagnosis models first define a prespecified exploratory Tier2 candidate universe. The primary ML question is then: among those candidates, which genes receive the most stable multivariate support across five model families and repeated resampling? The result is described as **ML-prioritized exploratory Tier2 candidates**, not validated biomarkers or diagnostic genes.

The canonical upstream source contains 671 region–gene pairs and was generated before historical ML scoring. Only region, gene ID, and gene symbol were imported. Historical ML scores, confidence labels, selection frequencies, Tier4, WGCNA, AHBA, and all-gene RRA ranks did not enter fitting or RRA.

## B. Difference from the frozen all-gene sensitivity analysis

The frozen all-gene analysis asks whether ML can discover useful genes from the full transcriptome using strict train-only feature selection. This primary analysis asks which genes are most consistently prioritized **within the already screened Tier2 universe**. Accordingly, it has no ANOVA/SMD selector choice and no nested top-k grid. The all-gene analysis remains frozen as a separate sensitivity analysis.

## C. Important limitation

Tier2 membership was defined from the complete cohort using Diagnosis. The candidate universe is therefore not independent of Diagnosis, and outer-fold classifier metrics may be optimistic. They are descriptive/internal sanity checks, not unbiased diagnostic-validation estimates. The primary endpoint is within-study gene prioritisation. There is no external validation, and no gene is a validated biomarker.

## D. Exact ML procedure

{markdown_table(sample_table)}

- Outer resampling: five stratified folds × five repeats; seeds 2026–2030.
- Inner tuning: three stratified folds using balanced accuracy.
- Preprocessing: training-only median imputation; training-only standardization for linear models.
- Candidate universe: all region-specific Tier2 genes, fixed throughout each regional analysis.
- L1/L2 logistic and linear SVM: `C ∈ {{0.01, 0.1, 1}}`.
- Elastic Net logistic: `C ∈ {{0.01, 0.1, 1}}`, `l1_ratio = 0.5`.
- Random Forest: 200 trees, `max_depth ∈ {{3, None}}`, `min_samples_leaf = 2`.
- L1/Elastic Net rank nonzero absolute coefficients; L2/SVM rank all Tier2 genes; Random Forest ranks nonzero feature importance and is nondirectional.
- Empty sparse-model lists are retained as empty evidence.
- A most-frequent dummy classifier is reported only as a descriptive reference and contributes no gene evidence.
- Level 1 RRA aggregates 25 outer-fit lists per region/model. Level 2 aggregates the five model-level lists per region.
- R 4.6.1; `RobustRankAggreg` 1.2.1; `aggregateRanks(method = "RRA")`.
- Performance and coefficient sign never weight RRA.

## E. Candidate accounting

{markdown_table(accounting_report)}

Tier2 membership and ML ranking support are separate concepts. A Tier2 gene remains a candidate even if a sparse model never assigns it nonzero evidence.

## Descriptive internal classifier performance

{markdown_table(performance_report)}

Because Tier2 was selected from the full cohort using Diagnosis, these values are optimistic/internal and must not be interpreted as unbiased diagnostic validation.

## Top 20 Tier2-RRA genes by region

{chr(10).join(top_sections)}

RRA scores are aggregation/ranking scores, not MDD-association p-values.

## F. Model agreement and recurrence

Model-level ranks and outer-fit recurrence are provided as transparent stability diagnostics. Agreement reflects repeated analyses of the same cohort and is not independent replication.

## Tier2-RRA versus frozen all-gene RRA

{markdown_table(convergence_report)}

This is a **within-study convergence/sensitivity comparison**. Neither ranking influenced the other, but both use the same cohort, so agreement is not independent validation. Direct top-10/top-25/top-50 overlaps and rank correlations are descriptive.

## G. Interpretation and next steps

The primary output is a continuous region-specific ranking of ML-prioritized exploratory Tier2 candidates. No High/Moderate/Low thresholds were created. Label permutations, Stability Selection, handmade-score benchmarking, LASSO-only baselines, limma/voom, DESeq2, and external validation remain later phases and were not run here.
"""
    (REPO_ROOT / "docs/TIER2_ML_RRA_PRIORITIZATION_REPORT.md").write_text(report, encoding="utf-8")

    hash_targets = sorted(
        path
        for path in output.rglob("*")
        if path.is_file()
        and path.name != "checksums.sha256"
        and "logs/checkpoints" not in path.as_posix()
        and path.name != "Tier2_validation_results.json"
    )
    lines = [f"{sha256_file(path)}  {path.relative_to(output).as_posix()}" for path in hash_targets]
    (output / "checksums.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Built Tier2 report outputs with {len(lines)} checksums")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
