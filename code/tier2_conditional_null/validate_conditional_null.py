#!/usr/bin/env python3
"""Validate raw and summarized Tier2 conditional-null evidence."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "outputs/length_gt_100_clean_pipeline/20_Tier2_ML_conditional_null"
TIER2 = ROOT / "outputs/length_gt_100_clean_pipeline/07_Tier2_candidate_gene_sets/Tier2_candidate_genes_all_regions.csv"
FROZEN = {
    ROOT / "outputs/length_gt_100_clean_pipeline/18_Tier2_ML_RRA_prioritization/Tier2_consensus_rra_rankings.csv": "2e1ed5fc91aae752f57e46e61022a9077d7f1ab285cfc05a6973c3c6df23d358",
    ROOT / "outputs/length_gt_100_clean_pipeline/19_Tier2_ML_benchmark/checksums.sha256": "4e17dc8234317efc9bce92aa03fb7ae63ff72a128392d202f2f31a105c09cf10",
    ROOT / "outputs/length_gt_100_clean_pipeline/17_ML_RRA_consensus/consensus_rra_rankings.csv": "0aefec5fcfa100deafa1525c2d31521ae8c00ca890ef4e1b4145b1b861c2c87e",
}
MODELS = {"logistic_l1", "elastic_net", "logistic_l2", "linear_svm", "random_forest"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(condition: bool, name: str, failures: list[str]) -> None:
    print(("PASS " if condition else "FAIL ") + name)
    if not condition:
        failures.append(name)


def main() -> int:
    failures: list[str] = []
    tier2 = pd.read_csv(TIER2, dtype={"gene_id": str})
    tier2["region"] = tier2["region_or_subregion"].replace({"Cingulate gyrus 25": "CG25"})
    expected_pairs = set(zip(tier2["region"], tier2["gene_id"]))
    check(len(expected_pairs) == 671, "frozen_tier2_has_671_unique_pairs", failures)
    check(all(path.exists() and sha256(path) == digest for path, digest in FROZEN.items()), "all_frozen_hashes_unchanged", failures)
    completions = sorted((OUTPUT / "permutations").glob("perm_*/permutation_completion.json")) if (OUTPUT / "permutations").exists() else []
    check(len(completions) >= 1, "at_least_one_complete_permutation", failures)
    for marker in completions:
        directory = marker.parent
        label = directory.name
        completion = json.loads(marker.read_text(encoding="utf-8"))
        check(all((directory / name).exists() and sha256(directory / name) == digest for name, digest in completion["files_sha256"].items()), f"{label}_checksums", failures)
        mapping = pd.read_csv(directory / "permuted_label_mapping.csv", dtype={"sample_id": str})
        check(not mapping.duplicated(["region", "sample_id"]).any(), f"{label}_sample_ids_unique_within_region", failures)
        check((mapping["observed_label"] != mapping["permuted_label"]).any(), f"{label}_labels_changed", failures)
        counts = mapping.groupby(["region", "observed_label"]).size().sort_index()
        perm_counts = mapping.groupby(["region", "permuted_label"]).size().sort_index()
        check(counts.equals(perm_counts), f"{label}_class_counts_preserved", failures)
        folds = pd.read_csv(directory / "outer_fold_assignments.csv", dtype={"sample_id": str})
        fold_counts = folds.groupby(["region", "repeat", "fold", "role", "permuted_label"]).size()
        check((fold_counts > 0).all(), f"{label}_outer_folds_stratified", failures)
        perf = pd.read_csv(directory / "Tier2_outer_performance.csv")
        feature = perf.loc[perf["model"].isin(MODELS)]
        check(feature["fit_id"].nunique() == 750, f"{label}_750_feature_model_fits", failures)
        check(perf.loc[perf["model"].eq("logistic_l1"), "fit_id"].nunique() == 150, f"{label}_150_lasso_fits", failures)
        check(perf.loc[perf["model"].eq("dummy_most_frequent"), "fit_id"].nunique() == 150, f"{label}_150_dummy_references", failures)
        rankings = pd.read_csv(directory / "Tier2_fold_gene_rankings.csv", dtype={"gene_id": str}, low_memory=False)
        check(set(rankings["model"]).issubset(MODELS), f"{label}_dummy_never_gene_evidence", failures)
        check(set(zip(rankings["region"], rankings["gene_id"])).issubset(expected_pairs), f"{label}_only_tier2_pairs_ranked", failures)
        consensus = pd.read_csv(directory / "Tier2_consensus_rra_rankings.csv", dtype={"gene_id": str})
        check(len(consensus) == 671 and set(zip(consensus["region"], consensus["gene_id"])) == expected_pairs, f"{label}_consensus_exact_universe", failures)
        resample = pd.read_csv(directory / "RRA_resample_consensus_rankings.csv", dtype={"gene_id": str})
        check(len(resample) == 16775 and resample.groupby(["region", "repeat", "fold"]).ngroups == 150, f"{label}_resample_rra_complete", failures)
        lasso = pd.read_csv(directory / "LASSO_gene_rankings.csv", dtype={"gene_id": str})
        check(len(lasso) == 671 and set(zip(lasso["region"], lasso["gene_id"])) == expected_pairs, f"{label}_lasso_exact_universe", failures)
        check(((lasso["selection_frequency"] * 25).round().astype(int) == lasso["selection_count"]).all(), f"{label}_lasso_denominator_25", failures)
    if (OUTPUT / "summary/observed_vs_null_summary.csv").exists():
        summary = pd.read_csv(OUTPUT / "summary/observed_vs_null_summary.csv")
        check(
            "empirical_p_pilot" in summary
            and summary["empirical_p_pilot"].ge((1 / 21) - 1e-12).all()
            and summary["empirical_p_pilot"].le(1 + 1e-12).all(),
            "pilot_plus_one_empirical_tails",
            failures,
        )
        check(not any("weighted" in column.lower() or "composite" in column.lower() for column in summary.columns), "no_weighted_or_composite_score", failures)
    if failures:
        print(f"Validation failed: {len(failures)} checks")
        return 1
    print(f"Validation passed for {len(completions)} complete permutation(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
