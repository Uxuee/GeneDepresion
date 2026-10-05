#!/usr/bin/env python3
"""Annotate RRA outputs, create diagnostics/report, and freeze output hashes."""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from common import MODEL_ORDER, atomic_write_csv, atomic_write_json, sha256_file


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "outputs/length_gt_100_clean_pipeline/17_ML_RRA_consensus",
    )
    parser.add_argument("--top-n", type=int, default=50)
    return parser.parse_args()


def save_figure(fig: plt.Figure, base: Path) -> None:
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = arguments()
    out = args.output_dir.resolve()
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    performance = pd.read_csv(out / "outer_fold_performance.csv")
    summary = pd.read_csv(out / "model_performance_summary.csv")
    rankings = pd.read_csv(out / "fold_gene_rankings.csv", dtype={"gene_id": str}, low_memory=False)
    level1 = pd.read_csv(out / "model_level_rra_rankings.csv", dtype={"gene_id": str})
    consensus = pd.read_csv(out / "consensus_rra_rankings.csv", dtype={"gene_id": str})

    symbol_map = (
        rankings[["gene_id", "gene_symbol"]]
        .fillna("")
        .drop_duplicates("gene_id")
        .set_index("gene_id")["gene_symbol"]
        .to_dict()
    )
    level1.insert(2, "gene_symbol", level1["gene_id"].map(symbol_map).fillna(""))
    consensus.insert(2, "gene_symbol", consensus["gene_id"].map(symbol_map).fillna(""))
    consensus["rra_score_interpretation"] = "rank aggregation score; not an MDD-association p-value"
    atomic_write_csv(level1, out / "model_level_rra_rankings.csv")
    atomic_write_csv(consensus, out / "consensus_rra_rankings.csv")

    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(12, 7))
    sns.boxplot(data=performance, x="region", y="balanced_accuracy", hue="model", ax=ax)
    ax.axhline(0.5, color="black", linestyle="--", linewidth=1)
    ax.set(title="Held-out balanced accuracy by region and model", xlabel="", ylabel="Balanced accuracy")
    ax.tick_params(axis="x", rotation=30)
    ax.legend(title="Model", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    save_figure(fig, figures / "outer_fold_model_performance")

    top = consensus.sort_values(["region", "consensus_rank"]).groupby("region", sort=False).head(args.top_n).copy()
    rank_cols = [f"{model}_rank" for model in MODEL_ORDER]
    matrix = top.set_index(["region", "gene_symbol", "gene_id"])[rank_cols].copy()
    matrix = -np.log10(matrix.astype(float).clip(lower=1))
    fig, ax = plt.subplots(figsize=(10, max(8, len(matrix) * 0.10)))
    sns.heatmap(matrix, cmap="viridis", ax=ax, cbar_kws={"label": "-log10(model-level rank)"})
    ax.set(title=f"Model-level ranks for top {args.top_n} consensus genes per region", xlabel="", ylabel="Region / gene")
    save_figure(fig, figures / "top_consensus_gene_model_rank_heatmap")

    overlap_rows = []
    top_n = 100
    for region, part in level1.groupby("region"):
        sets = {
            model: set(part[part["model"] == model].nsmallest(top_n, "model_rank")["gene_id"])
            for model in MODEL_ORDER
        }
        for left, right in combinations(MODEL_ORDER, 2):
            overlap_rows.append(
                {"region": region, "model_pair": f"{left} | {right}", "jaccard": len(sets[left] & sets[right]) / len(sets[left] | sets[right])}
            )
    overlap = pd.DataFrame(overlap_rows)
    atomic_write_csv(overlap, out / "top100_model_overlap.csv")
    overlap_matrix = overlap.pivot(index="region", columns="model_pair", values="jaccard")
    fig, ax = plt.subplots(figsize=(14, 5))
    sns.heatmap(overlap_matrix, cmap="mako", annot=True, fmt=".2f", vmin=0, vmax=1, ax=ax)
    ax.set(title="Jaccard overlap among model-level top-100 gene lists", xlabel="Model pair", ylabel="")
    save_figure(fig, figures / "top100_model_overlap")

    plot_level1 = level1.copy()
    plot_level1["negative_log10_rra"] = -np.log10(plot_level1["rra_score"].clip(lower=np.finfo(float).tiny))
    fig, axes = plt.subplots(2, 3, figsize=(17, 10), sharex=True, sharey=True)
    for ax, (region, part) in zip(axes.flat, plot_level1.groupby("region", sort=False)):
        sns.scatterplot(
            data=part,
            x="outer_list_appearance_proportion",
            y="negative_log10_rra",
            hue="model",
            s=20,
            alpha=0.55,
            ax=ax,
            legend=False,
        )
        ax.set_title(region)
    fig.suptitle("Outer-list appearance versus model-level RRA score", y=1.01)
    fig.supxlabel("Outer-fold list appearance proportion (descriptive)")
    fig.supylabel("-log10(RRA score; not an MDD-association p-value)")
    save_figure(fig, figures / "appearance_vs_model_rra")

    manifest_path = out / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
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
                "predictive_performance_weighting": False,
            },
            "result_counts": {
                "fold_ranking_rows": len(rankings),
                "model_level_rra_rows": len(level1),
                "consensus_rra_rows": len(consensus),
            },
        }
    )
    atomic_write_json(manifest, manifest_path)

    report = f"""# Experimental whole-transcriptome ML–RRA prioritisation

## Status and boundary

This analysis is an experimental prioritisation redesign. It does not replace or alter the canonical corrected 671-row ML priority result. The six original analysis units are Ventral Striatum, Insula, Hippocampal Formation, CG25, dlPFC BA8/9, and OFC BA11; pooled broad Frontal Cortex is excluded.

## What was done

Each normalized whole-transcriptome matrix was analysed independently. Five repeats of stratified five-fold outer cross-validation produced 25 held-out evaluations per model family. Within every outer-training split, three-fold inner cross-validation jointly selected the feature-selection method (ANOVA F or absolute standardized mean difference), the number of retained genes (25, 50, 100, 250, or 500), and compact model hyperparameters. Median imputation, scaling for linear models, feature selection, and tuning were all fit only on training data.

The five model families were L1 logistic regression, Elastic Net logistic regression, L2 logistic regression, linear SVM, and Random Forest. Fold rankings used absolute linear coefficients or Random Forest feature importance. Signs were exported for linear models as annotation but did not enter aggregation.

## Robust rank aggregation

Level 1 used `RobustRankAggreg::aggregateRanks(method = "RRA")` to aggregate the 25 outer-fold rankings separately for every region and model. Level 2 applied the same canonical R implementation to the five model-level rankings per region. Held-out predictive metrics did not weight either RRA level. RRA scores quantify unusually consistent high ranking relative to the specified gene universe; they are not MDD-association p-values.

## Output summary

- Outer model fits: {manifest['n_completed_outer_model_fits']}
- Fold-ranking rows: {len(rankings):,}
- Model-level RRA rows: {len(level1):,}
- Consensus RRA rows: {len(consensus):,}
- Canonical outputs modified: no

No high/moderate confidence threshold was invented for the new consensus. The main result is the continuous within-region `consensus_rank` and `consensus_rra_score`, accompanied by transparent per-model ranks.

## Interpretation cautions

This is internal resampling-based prioritisation, not external validation and not biomarker validation. Small regional sample sizes limit performance estimates and ranking stability. Genes missing from sparse fold/model rankings are handled by the RRA partial-list procedure rather than being assigned fabricated coefficients. Model performance is reported descriptively and is never used to up-weight a model's gene evidence.
"""
    report_path = Path(__file__).resolve().parents[2] / "docs/ML_RRA_REDESIGN_REPORT.md"
    report_path.write_text(report, encoding="utf-8")

    hash_targets = sorted(
        path for path in out.rglob("*") if path.is_file() and path.name != "checksums.sha256" and "logs/checkpoints" not in path.as_posix()
    )
    lines = [f"{sha256_file(path)}  {path.relative_to(out).as_posix()}" for path in hash_targets]
    (out / "checksums.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Built {len(consensus):,} consensus rows and {len(lines)} checksums")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
