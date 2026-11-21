#!/usr/bin/env python3
"""Descriptively compare Tier2 RRA with the frozen all-gene sensitivity ranking."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "code" / "ml_rra"))
from common import atomic_write_csv, sha256_file  # noqa: E402


ALLGENE_RELATIVE = "outputs/length_gt_100_clean_pipeline/17_ML_RRA_consensus/consensus_rra_rankings.csv"
ALLGENE_SHA256 = "0aefec5fcfa100deafa1525c2d31521ae8c00ca890ef4e1b4145b1b861c2c87e"


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=REPO_ROOT)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "outputs/length_gt_100_clean_pipeline/18_Tier2_ML_RRA_prioritization",
    )
    return parser.parse_args()


def main() -> int:
    args = arguments()
    root = args.root.resolve()
    output = args.output_dir.resolve()
    allgene_path = root / ALLGENE_RELATIVE
    if sha256_file(allgene_path) != ALLGENE_SHA256:
        raise RuntimeError("Frozen all-gene consensus checksum differs; refusing comparison")
    tier2 = pd.read_csv(output / "Tier2_consensus_rra_rankings.csv", dtype={"gene_id": str})
    allgene = pd.read_csv(allgene_path, dtype={"gene_id": str})
    allgene_fields = allgene[
        ["region", "gene_id", "consensus_rank", "consensus_rra_score", "model_appearance_count"]
    ].rename(
        columns={
            "consensus_rank": "allgene_consensus_rank",
            "consensus_rra_score": "allgene_consensus_rra_score",
            "model_appearance_count": "allgene_model_appearance_count",
        }
    )
    comparison = tier2.merge(allgene_fields, on=["region", "gene_id"], how="left", validate="one_to_one")
    comparison = comparison.rename(
        columns={
            "consensus_rank": "tier2_consensus_rank",
            "consensus_rra_score": "tier2_consensus_rra_score",
            "model_appearance_count": "tier2_model_appearance_count",
        }
    )
    comparison["present_in_frozen_allgene_rra"] = comparison["allgene_consensus_rank"].notna()
    allgene_sizes = allgene.groupby("region").size().to_dict()
    comparison["allgene_rank_percentile"] = comparison.apply(
        lambda row: row["allgene_consensus_rank"] / allgene_sizes[row["region"]]
        if pd.notna(row["allgene_consensus_rank"])
        else np.nan,
        axis=1,
    )
    comparison["comparison_role"] = "within-study convergence/sensitivity; not independent validation"
    comparison = comparison.sort_values(["region", "tier2_consensus_rank"], kind="mergesort")
    atomic_write_csv(comparison, output / "Tier2_vs_allgene_RRA_comparison.csv")

    summary_rows: list[dict] = []
    for region, tier_part in comparison.groupby("region", sort=False):
        all_part = allgene[allgene["region"].eq(region)]
        common = tier_part.dropna(subset=["allgene_consensus_rank"])
        rho = np.nan
        p_value = np.nan
        if len(common) >= 3:
            result = spearmanr(common["tier2_consensus_rank"], common["allgene_consensus_rank"])
            rho, p_value = float(result.statistic), float(result.pvalue)
        row = {
            "region": region,
            "n_tier2_candidates": len(tier_part),
            "n_present_in_allgene_rra": len(common),
            "proportion_present_in_allgene_rra": len(common) / len(tier_part),
            "spearman_rho_common_genes": rho,
            "spearman_p_value_descriptive": p_value,
            "median_allgene_rank_percentile": common["allgene_rank_percentile"].median(),
        }
        for n in [10, 25, 50]:
            tier_set = set(tier_part.nsmallest(n, "tier2_consensus_rank")["gene_id"])
            all_set = set(all_part.nsmallest(n, "consensus_rank")["gene_id"])
            overlap = len(tier_set & all_set)
            row[f"top{n}_overlap_count"] = overlap
            row[f"top{n}_overlap_proportion_of_tier2_top"] = overlap / min(n, len(tier_set))
            row[f"top{n}_jaccard"] = overlap / len(tier_set | all_set) if tier_set | all_set else np.nan
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    atomic_write_csv(summary, output / "Tier2_vs_allgene_RRA_summary.csv")

    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="talk")
    plot = comparison[comparison["present_in_frozen_allgene_rra"]].copy()
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))
    for ax, (region, part) in zip(axes.flat, plot.groupby("region", sort=False)):
        sns.scatterplot(
            data=part,
            x="tier2_consensus_rank",
            y="allgene_rank_percentile",
            s=35,
            alpha=0.7,
            ax=ax,
        )
        ax.set_title(region)
        ax.set_xlabel("")
        ax.set_ylabel("")
    fig.suptitle("Tier2-RRA versus frozen all-gene RRA: within-study rank convergence", y=1.01)
    fig.supxlabel("Tier2 consensus rank")
    fig.supylabel("All-gene consensus rank percentile (lower is higher-ranked)")
    fig.savefig(figures / "Tier2_vs_allgene_rank_convergence.png", dpi=300, bbox_inches="tight")
    fig.savefig(figures / "Tier2_vs_allgene_rank_convergence.pdf", bbox_inches="tight")
    plt.close(fig)

    overlap_long = summary.melt(
        id_vars=["region"],
        value_vars=["top10_overlap_count", "top25_overlap_count", "top50_overlap_count"],
        var_name="top_n",
        value_name="overlap_count",
    )
    overlap_long["top_n"] = overlap_long["top_n"].str.extract(r"(\d+)").astype(int)
    fig, ax = plt.subplots(figsize=(12, 7))
    sns.barplot(data=overlap_long, x="region", y="overlap_count", hue="top_n", ax=ax)
    ax.set(
        title="Direct overlap of Tier2-RRA and all-gene-RRA top-N lists",
        xlabel="",
        ylabel="Shared genes",
    )
    ax.tick_params(axis="x", rotation=25)
    ax.legend(title="Top N")
    fig.savefig(figures / "Tier2_vs_allgene_topN_overlap.png", dpi=300, bbox_inches="tight")
    fig.savefig(figures / "Tier2_vs_allgene_topN_overlap.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Compared {len(comparison)} Tier2 pairs with frozen all-gene RRA")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
