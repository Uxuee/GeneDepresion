#!/usr/bin/env python3
"""Validate frozen boundaries and numerical invariants of the Tier2 benchmark."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "outputs/length_gt_100_clean_pipeline/19_Tier2_ML_benchmark"
TIER2 = ROOT / "outputs/length_gt_100_clean_pipeline/07_Tier2_candidate_gene_sets/Tier2_candidate_genes_all_regions.csv"
PRIMARY = ROOT / "outputs/length_gt_100_clean_pipeline/18_Tier2_ML_RRA_prioritization/Tier2_consensus_rra_rankings.csv"
ALLGENE = ROOT / "outputs/length_gt_100_clean_pipeline/17_ML_RRA_consensus/consensus_rra_rankings.csv"
EXPECTED_PRIMARY = "2e1ed5fc91aae752f57e46e61022a9077d7f1ab285cfc05a6973c3c6df23d358"
EXPECTED_ALLGENE = "0aefec5fcfa100deafa1525c2d31521ae8c00ca890ef4e1b4145b1b861c2c87e"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load() -> dict[str, pd.DataFrame]:
    tables = {
        "tier2": pd.read_csv(TIER2, dtype={"gene_id": str}),
        "lasso": pd.read_csv(OUTPUT / "LASSO_gene_rankings.csv", dtype={"gene_id": str}),
        "lasso_fold": pd.read_csv(OUTPUT / "LASSO_fold_rankings.csv", dtype={"gene_id": str}, low_memory=False),
        "stability": pd.read_csv(OUTPUT / "stability_selection_gene_rankings.csv", dtype={"gene_id": str}),
        "details": pd.read_csv(OUTPUT / "stability_selection_run_details.csv", dtype={"gene_id": str}, low_memory=False),
        "overlap": pd.read_csv(OUTPUT / "method_pairwise_topset_overlap.csv"),
        "correlations": pd.read_csv(OUTPUT / "method_rank_correlations.csv"),
        "stability_summary": pd.read_csv(OUTPUT / "method_stability_summary.csv"),
        "coverage": pd.read_csv(OUTPUT / "method_coverage_summary.csv"),
    }
    tables["tier2"]["region"] = tables["tier2"]["region_or_subregion"].replace({"Cingulate gyrus 25": "CG25"})
    return tables


def validation_checks() -> dict[str, bool]:
    t = load()
    valid = set(zip(t["tier2"]["region"], t["tier2"]["gene_id"]))
    lasso_pairs = set(zip(t["lasso"]["region"], t["lasso"]["gene_id"]))
    stability_pairs = set(zip(t["stability"]["region"], t["stability"]["gene_id"]))
    detail_pairs = set(zip(t["details"]["region"], t["details"]["gene_id"]))
    expected_regions = {"Ventral Striatum", "Insula", "Hippocampal Formation", "CG25", "dlPFC BA8/9", "OFC BA11"}
    checks: dict[str, bool] = {}
    checks["exact_671_unique_tier2_pairs"] = len(valid) == 671 == len(t["tier2"])
    checks["exact_six_regions"] = set(t["lasso"]["region"]) == expected_regions == set(t["stability"]["region"])
    checks["lasso_uses_exact_tier2_pairs"] = lasso_pairs == valid
    checks["stability_uses_exact_tier2_pairs"] = stability_pairs == valid == detail_pairs
    checks["no_tier4_or_external_pairs"] = lasso_pairs.issubset(valid) and stability_pairs.issubset(valid)
    checks["frozen_primary_hash_unchanged"] = sha256(PRIMARY) == EXPECTED_PRIMARY
    checks["frozen_allgene_hash_unchanged"] = sha256(ALLGENE) == EXPECTED_ALLGENE
    checks["lasso_25_outer_fits_per_region"] = bool(t["lasso_fold"].groupby("region")["fit_id"].nunique().eq(25).all())
    expected_lasso = t["lasso_fold"].groupby(["region", "gene_id"])["fit_id"].nunique()
    observed_lasso = t["lasso"].set_index(["region", "gene_id"])["selection_count"]
    checks["lasso_selection_counts_reproduce"] = observed_lasso.reindex(expected_lasso.index).eq(expected_lasso).all() and int(observed_lasso.sum()) == len(t["lasso_fold"])
    checks["lasso_denominator_exactly_25"] = np.allclose(t["lasso"]["selection_frequency"], t["lasso"]["selection_count"] / 25)
    checks["stability_100_subsamples_per_region"] = bool(t["details"].groupby("region")["subsample_id"].nunique().eq(100).all())
    checks["stability_one_event_per_gene_subsample"] = not t["details"].duplicated(["region", "subsample_id", "gene_id"]).any()
    checks["stability_preserves_sample_mapping"] = bool(t["details"].groupby(["region", "subsample_id"])["selected_sample_ids_sha256"].nunique().eq(1).all())
    expected_half = {
        "Ventral Striatum": (11, 13), "Insula": (11, 13),
        "Hippocampal Formation": (9, 12), "CG25": (7, 6),
        "dlPFC BA8/9": (11, 13), "OFC BA11": (11, 12),
    }
    checks["stability_balanced_half_samples"] = all(
        part["n_control_subsample"].eq(expected_half[region][0]).all()
        and part["n_mdd_subsample"].eq(expected_half[region][1]).all()
        for region, part in t["details"].groupby("region")
    )
    selected = t["details"]["selected"].astype(str).str.lower().eq("true")
    expected_stability = t["details"].loc[selected].groupby(["region", "gene_id"]).size()
    observed_stability = t["stability"].set_index(["region", "gene_id"])["selection_count"]
    checks["stability_selection_counts_reproduce"] = observed_stability.reindex(expected_stability.index).eq(expected_stability).all()
    checks["stability_denominator_exactly_100"] = np.allclose(t["stability"]["selection_probability"], t["stability"]["selection_count"] / 100)
    checks["stability_seed_grid_exact"] = set(t["details"]["base_seed"]) == set(range(31001, 31101))
    checks["never_selected_shared_worst_rank"] = all(
        part.loc[part["selection_count"].eq(0), rank].nunique() <= 1
        for table, rank in [(t["lasso"], "lasso_correlation_rank"), (t["stability"], "stability_correlation_rank")]
        for _, part in table.groupby("region")
    )
    expected_requested = []
    counts = t["tier2"].groupby("region").size().to_dict()
    for row in t["overlap"].itertuples(index=False):
        expected_requested.append({"top10": 10, "top25": 25, "top20percent": math.ceil(0.2 * counts[row.region])}[row.top_definition])
    checks["topset_definitions_exact"] = t["overlap"]["requested_n"].tolist() == expected_requested and len(t["overlap"]) == 54
    checks["rank_correlation_grid_complete"] = len(t["correlations"]) == 18
    checks["coverage_summary_reproduces"] = len(t["coverage"]) == 18 and set(t["coverage"]["method"]) == {"Tier2-RRA", "LASSO", "Stability Selection"}
    checks["resampling_summary_complete"] = set(t["stability_summary"]["method"]) == {"Tier2-RRA", "LASSO", "Stability Selection"}
    forbidden_columns = {"weighted_score", "priority_score", "balanced_accuracy", "roc_auc", "f1"}
    checks["no_metric_or_weighted_priority_score"] = forbidden_columns.isdisjoint(t["lasso"].columns) and forbidden_columns.isdisjoint(t["stability"].columns)
    manifest = json.loads((OUTPUT / "benchmark_manifest.json").read_text(encoding="utf-8"))
    checks["manifest_rejects_composite_and_metric_weighting"] = manifest["comparisons"]["composite_winner_score"] is False and manifest["comparisons"]["classifier_metrics_used_for_gene_priority"] is False
    checksum_ok = True
    for line in (OUTPUT / "checksums.sha256").read_text(encoding="utf-8").splitlines():
        expected, relative = line.split("  ", 1)
        checksum_ok &= sha256(OUTPUT / relative) == expected
    checks["all_output_checksums_match"] = checksum_ok
    return {name: bool(value) for name, value in checks.items()}


def main() -> int:
    checks = validation_checks()
    result = {"passed": all(checks.values()), "checks": checks}
    (OUTPUT / "benchmark_validation_results.json").write_bytes(json.dumps(result, indent=2, sort_keys=True).encode("utf-8"))
    for name, passed in checks.items():
        print(f"{'PASS' if passed else 'FAIL'} {name}")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
