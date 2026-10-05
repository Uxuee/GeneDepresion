#!/usr/bin/env python3
"""Validate complete ML–RRA outputs and leakage/integrity invariants."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from common import ANALYSES, MODEL_ORDER, sha256_file


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "outputs/length_gt_100_clean_pipeline/17_ML_RRA_consensus",
    )
    args = parser.parse_args()
    out = args.output_dir.resolve()
    checks: dict[str, bool] = {}
    rankings = pd.read_csv(out / "fold_gene_rankings.csv", dtype={"gene_id": str}, low_memory=False)
    predictions = pd.read_csv(out / "outer_fold_predictions.csv")
    performance = pd.read_csv(out / "outer_fold_performance.csv")
    level1 = pd.read_csv(out / "model_level_rra_rankings.csv", dtype={"gene_id": str})
    consensus = pd.read_csv(out / "consensus_rra_rankings.csv", dtype={"gene_id": str})
    manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))

    regions = {spec.region for spec in ANALYSES}
    checks["six_exact_regions"] = set(consensus["region"]) == regions
    checks["no_pooled_frontal"] = "Frontal Cortex" not in set(consensus["region"])
    checks["750_unique_outer_model_fits"] = performance["fit_id"].nunique() == 6 * 25 * 5 == len(performance)
    checks["unique_gene_per_fit"] = not rankings.duplicated(["fit_id", "gene_id"]).any()
    checks["five_models_each_region"] = all(
        set(part["model"]) == set(MODEL_ORDER) for _, part in level1.groupby("region")
    )
    checks["25_lists_per_region_model"] = bool((level1["n_outer_lists"] == 25).all())
    checks["rf_sign_is_missing"] = bool(rankings.loc[rankings["model"] == "random_forest", "sign"].isna().all())
    linear = rankings[rankings["model"] != "random_forest"]
    checks["linear_sign_is_defined"] = bool(linear["sign"].notna().all())
    checks["prediction_coverage_once_per_repeat"] = bool(
        (predictions.groupby(["region", "repeat", "model", "sample_id"]).size() == 1).all()
    )
    checks["consensus_unique_region_gene"] = not consensus.duplicated(["region", "gene_id"]).any()
    checks["score_not_labeled_mdd_pvalue"] = bool(
        consensus["rra_score_interpretation"].str.contains("not an MDD-association p-value", regex=False).all()
    )
    checks["no_predictive_metric_weighting"] = manifest["rra"]["predictive_performance_weighting"] is False
    canonical_files = {
        Path(__file__).resolve().parents[2] / "outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_corrected/ML_prioritized_Tier2_genes_all_ranked.csv": "6f59e2ed6f652061fcb81a7c5eccbc1211a8dcb9d0bafcc85a1f9fc3edf280f1",
        Path(__file__).resolve().parents[2] / "outputs/length_gt_100_clean_pipeline/12_MDD_WGCNA_modules/best_solution_gene_to_module_assignments.csv": "b7c0a1bd65c07be12f1cca561ad6f6e07eb0619cc2e69067a2a1936390b8b92a",
    }
    checks["canonical_outputs_declared_unchanged"] = manifest["canonical_outputs_modified"] is False
    checks["canonical_file_hashes_unchanged"] = all(path.exists() and sha256_file(path) == expected for path, expected in canonical_files.items())

    checksum_ok = True
    for line in (out / "checksums.sha256").read_text(encoding="utf-8").splitlines():
        expected, relative = line.split("  ", 1)
        checksum_ok &= sha256_file(out / relative) == expected
    checks["all_checksums_match"] = checksum_ok
    result = {"passed": all(checks.values()), "checks": checks}
    (out / "validation_results.json").write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    for name, passed in checks.items():
        print(f"{'PASS' if passed else 'FAIL'} {name}")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
