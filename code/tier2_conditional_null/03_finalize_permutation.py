#!/usr/bin/env python3
"""Build the frozen LASSO-alone ranking and final audit record for one permutation."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TIER2 = ROOT / "outputs/length_gt_100_clean_pipeline/07_Tier2_candidate_gene_sets/Tier2_candidate_genes_all_regions.csv"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, lineterminator="\n")


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("permutation_dir", type=Path)
    return parser.parse_args()


def main() -> int:
    output = arguments().permutation_dir.resolve()
    rankings = pd.read_csv(output / "Tier2_fold_gene_rankings.csv", dtype={"gene_id": str}, low_memory=False)
    performance = pd.read_csv(output / "Tier2_outer_performance.csv", low_memory=False)
    candidates = pd.read_csv(TIER2, dtype={"gene_id": str}, usecols=["region_or_subregion", "gene_id", "gene_symbol"])
    candidates["region"] = candidates["region_or_subregion"].replace({"Cingulate gyrus 25": "CG25"})
    lasso = rankings.loc[rankings["model"].eq("logistic_l1")].copy()
    if performance.loc[performance["model"].eq("logistic_l1"), "fit_id"].nunique() != 150:
        raise RuntimeError("Expected 150 permutation LASSO fits")
    rows = []
    for region, universe in candidates.groupby("region", sort=False):
        part = lasso.loc[lasso["region"].eq(region)]
        if part["fit_id"].nunique() > 25:
            raise RuntimeError(f"{region}: more than 25 LASSO fits")
        summary = part.groupby("gene_id", as_index=False).agg(
            selection_count=("fit_id", "nunique"),
            mean_rank_when_selected=("within_fit_rank", "mean"),
            median_rank_when_selected=("within_fit_rank", "median"),
            mean_absolute_coefficient=("absolute_importance", "mean"),
            positive_count=("sign", lambda values: int((values > 0).sum())),
            negative_count=("sign", lambda values: int((values < 0).sum())),
        )
        merged = universe[["gene_id", "gene_symbol"]].merge(summary, on="gene_id", how="left")
        merged["selection_count"] = merged["selection_count"].fillna(0).astype(int)
        merged["selection_frequency"] = merged["selection_count"] / 25.0
        merged["sign_consistency"] = np.where(
            merged["selection_count"].gt(0),
            merged[["positive_count", "negative_count"]].fillna(0).max(axis=1) / merged["selection_count"],
            np.nan,
        )
        merged["region"] = region
        merged = merged.sort_values(
            ["selection_frequency", "mean_rank_when_selected", "mean_absolute_coefficient", "gene_id"],
            ascending=[False, True, False, True], na_position="last", kind="mergesort",
        ).reset_index(drop=True)
        merged["lasso_export_rank"] = np.arange(1, len(merged) + 1)
        selected = merged["selection_count"].gt(0)
        worst = int(selected.sum()) + 1
        merged["lasso_correlation_rank"] = np.where(selected, merged["lasso_export_rank"], worst)
        rows.append(merged)
    result = pd.concat(rows, ignore_index=True)
    if len(result) != 671 or result.duplicated(["region", "gene_id"]).any():
        raise RuntimeError("Permutation LASSO ranking does not contain exact 671 pairs")
    write_csv(result, output / "LASSO_gene_rankings.csv")

    required = [
        "permuted_label_mapping.csv", "outer_fold_assignments.csv", "Tier2_fold_gene_rankings.csv",
        "Tier2_outer_predictions.csv", "Tier2_outer_performance.csv", "Tier2_model_performance_summary.csv",
        "Tier2_analysis_unit_manifest.csv", "Tier2_model_level_rra_rankings.csv",
        "Tier2_consensus_rra_rankings.csv", "Tier2_candidate_accounting.csv",
        "RRA_resample_consensus_rankings.csv", "LASSO_gene_rankings.csv",
    ]
    missing = [name for name in required if not (output / name).exists()]
    if missing:
        raise RuntimeError(f"Permutation outputs missing: {missing}")
    hashes = {name: sha256(output / name) for name in required}
    completion = json.loads((output / "ml_completion.json").read_text(encoding="utf-8"))
    payload = {
        "status": "complete", "permutation_index": completion["permutation_index"],
        "permutation_seed": completion["permutation_seed"],
        "wall_runtime_seconds_ml": completion["wall_runtime_seconds"],
        "n_feature_model_fits": 750, "n_lasso_outer_fits": 150,
        "files_sha256": hashes,
    }
    (output / "permutation_completion.json").write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(f"perm_{payload['permutation_index']:04d}: RRA and LASSO outputs finalized")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
