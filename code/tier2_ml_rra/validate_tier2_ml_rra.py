#!/usr/bin/env python3
"""Validate Tier2 ML–RRA candidate, model, aggregation, and frozen-boundary invariants."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "code" / "ml_rra"))
from common import MODEL_ORDER, sha256_file  # noqa: E402


EXPECTED_COUNTS = {
    "Ventral Striatum": 146,
    "Insula": 70,
    "Hippocampal Formation": 127,
    "CG25": 186,
    "dlPFC BA8/9": 72,
    "OFC BA11": 70,
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "outputs/length_gt_100_clean_pipeline/18_Tier2_ML_RRA_prioritization",
    )
    args = parser.parse_args()
    output = args.output_dir.resolve()
    tier2_path = REPO_ROOT / "outputs/length_gt_100_clean_pipeline/07_Tier2_candidate_gene_sets/Tier2_candidate_genes_all_regions.csv"
    tier2 = pd.read_csv(tier2_path, dtype={"gene_id": str})
    tier2["region"] = tier2["region_or_subregion"].replace({"Cingulate gyrus 25": "CG25"})
    rankings = pd.read_csv(output / "Tier2_fold_gene_rankings.csv", dtype={"gene_id": str}, low_memory=False)
    performance = pd.read_csv(output / "Tier2_outer_performance.csv")
    predictions = pd.read_csv(output / "Tier2_outer_predictions.csv")
    level1 = pd.read_csv(output / "Tier2_model_level_rra_rankings.csv", dtype={"gene_id": str})
    consensus = pd.read_csv(output / "Tier2_consensus_rra_rankings.csv", dtype={"gene_id": str})
    accounting = pd.read_csv(output / "Tier2_candidate_accounting.csv")
    manifest = json.loads((output / "Tier2_run_manifest.json").read_text(encoding="utf-8"))
    valid_pairs = set(zip(tier2["region"], tier2["gene_id"]))
    ranking_pairs = set(zip(rankings["region"], rankings["gene_id"]))
    consensus_pairs = set(zip(consensus["region"], consensus["gene_id"]))
    regions = set(EXPECTED_COUNTS)
    checks: dict[str, bool] = {}
    checks["exact_six_regions"] = set(consensus["region"]) == regions
    checks["tier2_source_sha256"] = sha256_file(tier2_path) == "1672ddf95a4f4c9e44ae42e7accc208c311e2e1da8533e1808dcf91bbf7b67a1"
    checks["tier2_source_671_unique_pairs"] = len(tier2) == 671 and not tier2.duplicated(["region", "gene_id"]).any()
    checks["regional_candidate_counts"] = tier2.groupby("region").size().to_dict() == EXPECTED_COUNTS
    checks["only_tier2_pairs_in_fold_rankings"] = ranking_pairs.issubset(valid_pairs)
    checks["no_tier4_or_external_candidate_pair"] = consensus_pairs.issubset(valid_pairs)
    checks["all_671_candidates_in_consensus"] = len(consensus_pairs) == 671 and consensus_pairs == valid_pairs
    checks["fixed_tier2_feature_counts"] = all(
        part["n_tier2_features"].eq(EXPECTED_COUNTS[region]).all()
        for region, part in rankings.groupby("region")
    )
    checks["no_selector_or_topk_layer"] = manifest["feature_selector"] is None and manifest["top_k_grid"] is None
    feature_perf = performance[performance["model"].isin(MODEL_ORDER)]
    dummy_perf = performance[performance["model"].eq("dummy_most_frequent")]
    checks["exact_750_feature_model_fits"] = len(feature_perf) == 750 and feature_perf["fit_id"].nunique() == 750
    checks["exact_150_dummy_reference_fits"] = len(dummy_perf) == 150 and dummy_perf["fit_id"].nunique() == 150
    checks["25_fits_per_region_model"] = bool(
        (feature_perf.groupby(["region", "model"])["fit_id"].nunique() == 25).all()
    )
    checks["dummy_never_enters_rankings"] = "dummy_most_frequent" not in set(rankings["model"])
    checks["dummy_marked_noncontributing"] = not dummy_perf["contributes_gene_ranking"].astype(bool).any()
    checks["five_rra_models_only"] = set(level1["model"]) == set(MODEL_ORDER)
    checks["25_lists_per_level1_ranking"] = bool((level1["n_outer_lists"] == 25).all())
    checks["unique_gene_per_fit"] = not rankings.duplicated(["fit_id", "gene_id"]).any()
    checks["rf_is_nondirectional"] = rankings.loc[rankings["model"].eq("random_forest"), "sign"].isna().all()
    checks["linear_sign_stored"] = rankings.loc[~rankings["model"].eq("random_forest"), "sign"].notna().all()
    checks["prediction_coverage_once_per_repeat"] = bool(
        (predictions.groupby(["region", "repeat", "model", "sample_id"]).size() == 1).all()
    )
    checks["dummy_no_gene_evidence_manifest"] = manifest["rra"]["dummy_contributes_gene_evidence"] is False
    checks["no_performance_weighting"] = manifest["rra"]["predictive_performance_weighting"] is False
    checks["candidate_accounting_totals_671"] = int(accounting["total_tier2_candidates"].sum()) == 671
    checks["allgene_frozen_hash_unchanged"] = sha256_file(
        REPO_ROOT / "outputs/length_gt_100_clean_pipeline/17_ML_RRA_consensus/consensus_rra_rankings.csv"
    ) == "0aefec5fcfa100deafa1525c2d31521ae8c00ca890ef4e1b4145b1b861c2c87e"
    checks["old_corrected_hash_unchanged"] = sha256_file(
        REPO_ROOT / "outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_corrected/ML_prioritized_Tier2_genes_all_ranked.csv"
    ) == "6f59e2ed6f652061fcb81a7c5eccbc1211a8dcb9d0bafcc85a1f9fc3edf280f1"
    checksum_ok = True
    for line in (output / "checksums.sha256").read_text(encoding="utf-8").splitlines():
        expected, relative = line.split("  ", 1)
        checksum_ok &= sha256_file(output / relative) == expected
    checks["all_output_checksums_match"] = checksum_ok
    checks["rra_scores_labeled_as_ranking_not_mdd_p"] = consensus["rra_score_interpretation"].str.contains(
        "not an MDD-association p-value", regex=False
    ).all()
    checks["no_pooled_frontal"] = "Frontal Cortex" not in regions and "Frontal Cortex" not in set(consensus["region"])
    # pandas/numpy comparisons return numpy.bool_, which the standard JSON
    # encoder does not serialize.  Normalize every published check to a
    # built-in bool so the validation record is portable.
    checks = {name: bool(passed) for name, passed in checks.items()}
    result = {"passed": all(checks.values()), "checks": checks}
    (output / "Tier2_validation_results.json").write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    for name, passed in checks.items():
        print(f"{'PASS' if passed else 'FAIL'} {name}")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
