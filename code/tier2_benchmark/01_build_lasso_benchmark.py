#!/usr/bin/env python3
"""Build the transparent LASSO-only benchmark from frozen primary L1 events."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PRIMARY = ROOT / "outputs/length_gt_100_clean_pipeline/18_Tier2_ML_RRA_prioritization"
OUTPUT = ROOT / "outputs/length_gt_100_clean_pipeline/19_Tier2_ML_benchmark"
TIER2 = ROOT / "outputs/length_gt_100_clean_pipeline/07_Tier2_candidate_gene_sets/Tier2_candidate_genes_all_regions.csv"
EXPECTED_TIER2_SHA = "1672ddf95a4f4c9e44ae42e7accc208c311e2e1da8533e1808dcf91bbf7b67a1"
EXPECTED_PRIMARY_SHA = "2e1ed5fc91aae752f57e46e61022a9077d7f1ab285cfc05a6973c3c6df23d358"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, lineterminator="\n")


def main() -> int:
    if sha256(TIER2) != EXPECTED_TIER2_SHA:
        raise RuntimeError("Frozen Tier2 source checksum mismatch")
    if sha256(PRIMARY / "Tier2_consensus_rra_rankings.csv") != EXPECTED_PRIMARY_SHA:
        raise RuntimeError("Frozen Tier2-RRA consensus checksum mismatch")
    tier2 = pd.read_csv(TIER2, dtype={"gene_id": str}, usecols=["region_or_subregion", "gene_id", "gene_symbol"])
    tier2["region"] = tier2["region_or_subregion"].replace({"Cingulate gyrus 25": "CG25"})
    fold = pd.read_csv(PRIMARY / "Tier2_fold_gene_rankings.csv", dtype={"gene_id": str}, low_memory=False)
    lasso = fold.loc[fold["model"].eq("logistic_l1")].copy()
    perf = pd.read_csv(PRIMARY / "Tier2_outer_performance.csv")
    lasso_perf = perf.loc[perf["model"].eq("logistic_l1")].copy()
    if lasso_perf.groupby("region")["fit_id"].nunique().ne(25).any():
        raise RuntimeError("Frozen LASSO events do not contain 25 fits per region")
    write_csv(lasso, OUTPUT / "LASSO_fold_rankings.csv")

    rows: list[dict] = []
    for region, candidates in tier2.groupby("region", sort=False):
        evidence = lasso.loc[lasso["region"].eq(region)]
        by_gene = evidence.groupby("gene_id")
        records = []
        for candidate in candidates.itertuples(index=False):
            if candidate.gene_id in by_gene.groups:
                part = by_gene.get_group(candidate.gene_id)
                count = int(part["fit_id"].nunique())
                signs = part["sign"].dropna()
                sign_consistency = float(signs.value_counts().max() / len(signs)) if len(signs) else np.nan
                direction = int(signs.sum() > 0) - int(signs.sum() < 0) if len(signs) else 0
                records.append({
                    "region": region, "gene_id": candidate.gene_id, "gene_symbol": candidate.gene_symbol,
                    "selection_count": count, "selection_frequency": count / 25,
                    "mean_rank_when_selected": float(part["within_fit_rank"].mean()),
                    "median_rank_when_selected": float(part["within_fit_rank"].median()),
                    "mean_abs_coefficient_when_selected": float(part["absolute_importance"].mean()),
                    "sign_consistency_when_selected": sign_consistency,
                    "majority_direction": direction,
                })
            else:
                records.append({
                    "region": region, "gene_id": candidate.gene_id, "gene_symbol": candidate.gene_symbol,
                    "selection_count": 0, "selection_frequency": 0.0,
                    "mean_rank_when_selected": np.nan, "median_rank_when_selected": np.nan,
                    "mean_abs_coefficient_when_selected": 0.0,
                    "sign_consistency_when_selected": np.nan, "majority_direction": 0,
                })
        table = pd.DataFrame(records)
        table = table.sort_values(
            ["selection_frequency", "mean_rank_when_selected", "mean_abs_coefficient_when_selected", "gene_id"],
            ascending=[False, True, False, True], na_position="last", kind="mergesort",
        ).reset_index(drop=True)
        table["lasso_export_rank"] = np.arange(1, len(table) + 1)
        supported = table["selection_count"].gt(0)
        worst = int(supported.sum()) + 1
        table["lasso_correlation_rank"] = table["lasso_export_rank"].where(supported, worst)
        rows.extend(table.to_dict("records"))
    rankings = pd.DataFrame(rows)
    write_csv(rankings, OUTPUT / "LASSO_gene_rankings.csv")
    metrics = ["balanced_accuracy", "f1", "roc_auc", "precision", "recall"]
    summary = lasso_perf.groupby("region", as_index=False).agg(
        n_outer_fits=("fit_id", "nunique"), runtime_seconds=("elapsed_seconds", "sum"),
        **{f"{metric}_mean": (metric, "mean") for metric in metrics},
        **{f"{metric}_sd": (metric, "std") for metric in metrics},
    )
    summary["source"] = "reused_exact_frozen_primary_logistic_l1_events"
    summary["performance_interpretation"] = "optimistic within-study diagnostics, not unbiased out-of-sample estimates"
    write_csv(summary, OUTPUT / "LASSO_performance_summary.csv")
    print(f"LASSO benchmark: {len(rankings)} candidates from {lasso_perf.fit_id.nunique()} frozen fits")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
