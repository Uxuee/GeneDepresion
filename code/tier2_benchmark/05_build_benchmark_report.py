#!/usr/bin/env python3
"""Finalize manifest, checksums, and supervisor report for the Tier2 benchmark."""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import sklearn

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "outputs/length_gt_100_clean_pipeline/19_Tier2_ML_benchmark"
PRIMARY = ROOT / "outputs/length_gt_100_clean_pipeline/18_Tier2_ML_RRA_prioritization"
ALLGENE = ROOT / "outputs/length_gt_100_clean_pipeline/17_ML_RRA_consensus/consensus_rra_rankings.csv"
TIER2 = ROOT / "outputs/length_gt_100_clean_pipeline/07_Tier2_candidate_gene_sets/Tier2_candidate_genes_all_regions.csv"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def markdown(frame: pd.DataFrame) -> str:
    frame = frame.copy().fillna("")
    header = "| " + " | ".join(frame.columns.astype(str)) + " |"
    rule = "| " + " | ".join(["---"] * len(frame.columns)) + " |"
    rows = ["| " + " | ".join(map(str, row)) + " |" for row in frame.to_numpy()]
    return "\n".join([header, rule, *rows])


def rounded(frame: pd.DataFrame, digits: int = 3) -> pd.DataFrame:
    result = frame.copy()
    for column in result.select_dtypes(include=[np.number]).columns:
        result[column] = result[column].map(lambda value: "" if pd.isna(value) else f"{value:.{digits}f}")
    return result


def main() -> int:
    tier2 = pd.read_csv(TIER2, dtype={"gene_id": str})
    tier2["region"] = tier2["region_or_subregion"].replace({"Cingulate gyrus 25": "CG25"})
    rra = pd.read_csv(PRIMARY / "Tier2_consensus_rra_rankings.csv", dtype={"gene_id": str})
    lasso = pd.read_csv(OUTPUT / "LASSO_gene_rankings.csv", dtype={"gene_id": str})
    stability_genes = pd.read_csv(OUTPUT / "stability_selection_gene_rankings.csv", dtype={"gene_id": str})
    correlations = pd.read_csv(OUTPUT / "method_rank_correlations.csv")
    overlap = pd.read_csv(OUTPUT / "method_pairwise_topset_overlap.csv")
    stability = pd.read_csv(OUTPUT / "method_stability_summary.csv")
    coverage = pd.read_csv(OUTPUT / "method_coverage_summary.csv")
    runtime = pd.read_csv(OUTPUT / "method_runtime_summary.csv")
    allgene = pd.read_csv(OUTPUT / "benchmark_vs_allgene_RRA.csv")
    performance = pd.read_csv(OUTPUT / "LASSO_performance_summary.csv")

    top_rows = []
    for region in rra["region"].drop_duplicates():
        for method, frame, rank_col, evidence_col in [
            ("Tier2-RRA", rra, "consensus_rank", "consensus_rra_score"),
            ("LASSO", lasso, "lasso_export_rank", "selection_frequency"),
            ("Stability Selection", stability_genes, "stability_export_rank", "selection_probability"),
        ]:
            part = frame.loc[frame["region"].eq(region)].nsmallest(10, rank_col)
            for row in part.itertuples(index=False):
                top_rows.append({"region": region, "method": method, "rank": getattr(row, rank_col),
                                 "gene_id": row.gene_id, "gene_symbol": row.gene_symbol,
                                 "primary_evidence": getattr(row, evidence_col)})
    top = pd.DataFrame(top_rows)
    top.to_csv(OUTPUT / "benchmark_top_genes_by_method.csv", index=False, lineterminator="\n")

    rank_summary = correlations.groupby(["method_left", "method_right"], as_index=False).agg(
        spearman_min=("spearman_rho", "min"), spearman_median=("spearman_rho", "median"), spearman_max=("spearman_rho", "max"),
        kendall_min=("kendall_tau", "min"), kendall_median=("kendall_tau", "median"), kendall_max=("kendall_tau", "max"),
    )
    overlap_summary = overlap.groupby(["method_left", "method_right", "top_definition"], as_index=False).agg(
        jaccard_min=("jaccard", "min"), jaccard_median=("jaccard", "median"), jaccard_max=("jaccard", "max")
    )
    resampling_summary = stability.groupby(["method", "stability_metric", "top_definition"], as_index=False).agg(
        regional_mean_min=("mean", "min"), regional_mean_median=("mean", "median"), regional_mean_max=("mean", "max")
    )
    threshold_counts = stability_genes.groupby("region")["above_predeclared_0_60_threshold"].sum().rename("genes_probability_ge_0_60").reset_index()

    manifest = {
        "status": "complete",
        "source_commit": "df3596ed208bb7472e166e5fc0330b69453a9e31",
        "branch": "codex/tier2-prioritization-benchmark",
        "candidate_universe": {"pairs": 671, "source": str(TIER2.relative_to(ROOT)).replace("\\", "/"), "sha256": sha256(TIER2)},
        "frozen_inputs": {
            "tier2_rra_consensus_sha256": sha256(PRIMARY / "Tier2_consensus_rra_rankings.csv"),
            "allgene_rra_consensus_sha256": sha256(ALLGENE),
        },
        "lasso": {"source": "exact frozen primary logistic_l1 events", "outer_fits": 150, "denominator_per_region": 25,
                  "ranking": "selection frequency; mean selected rank; mean abs coefficient; stable gene ID"},
        "stability_selection": {"implementation": "transparent classical stratified half-sampling L1 logistic regression",
                                "subsamples_per_region": 100, "total_subsamples": 600, "inner_folds": 3,
                                "c_grid": [0.01, 0.1, 1.0], "seeds": [31001, 31100],
                                "descriptive_threshold": 0.60, "continuous_probability_primary": True},
        "comparisons": {"top_sets": ["top10", "top25", "ceil(20% of regional Tier2 universe)"],
                        "never_selected_rank_rule": "shared worst-rank block",
                        "classifier_metrics_used_for_gene_priority": False,
                        "composite_winner_score": False},
        "software": {"python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__,
                     "scipy": scipy.__version__, "scikit_learn": sklearn.__version__,
                     "R": "4.6.1", "RobustRankAggreg": "1.2.1"},
        "interpretation": "Scenario E with Scenario B nuance: RRA improves resampling consensus stability, while LASSO recovers much of the final shortlist; Stability Selection supports substantial top-gene convergence but is itself unstable under half-sampling.",
        "next_phase_not_run": "Tier2 conditional label-permutation/null testing",
    }
    (OUTPUT / "benchmark_manifest.json").write_bytes(json.dumps(manifest, indent=2, sort_keys=True).encode("utf-8"))

    region_counts = tier2.groupby("region").size().rename("Tier2 candidates").reset_index().rename(columns={"region": "Region"})
    top_five = top.loc[top["rank"].le(5)].copy()
    top_five["primary_evidence"] = top_five["primary_evidence"].map(lambda value: f"{value:.4g}")
    report = f"""# Tier2 ML prioritisation method benchmark

## 1. Scientific question

This benchmark asks whether the frozen five-model Tier2–RRA ranking adds useful prioritisation stability beyond two simpler sparse alternatives: LASSO alone and classical half-sampling Stability Selection. The endpoint is gene prioritisation—not classifier accuracy.

Tier2 was defined using Diagnosis in the complete cohort. Therefore all results are exploratory, within-study comparisons. They are not external validation, and the genes are not validated biomarkers.

## 2. Frozen primary Tier2–RRA method

The primary result was read without modification from frozen commit `df3596ed208bb7472e166e5fc0330b69453a9e31`. It combines L1, Elastic Net, L2, linear SVM, and Random Forest rankings across 25 repeated outer splits using two-level Robust Rank Aggregation. No benchmark output was allowed to change this ranking.

{markdown(region_counts)}

## 3. LASSO-alone benchmark

The exact 150 L1-logistic fits already present in the frozen primary run were reused: 25 fits per region, training-only median imputation and scaling, and training-only three-fold tuning over `C = 0.01, 0.1, 1`. The primary quantity is selection frequency across 25 fits. Tie-breakers are mean selected rank, mean absolute coefficient, and stable gene ID. No weighted score was created.

## 4. Stability Selection benchmark

No compatible Stability Selection package was installed in the frozen environment. The prospective design therefore used a transparent classical procedure: 100 deterministic stratified half-samples per region, L1 logistic regression, training-only preprocessing, and three-fold tuning over the same C grid within every half-sample. The continuous selection probability is primary. A predeclared 0.60 threshold is descriptive only.

{markdown(threshold_counts)}

The low threshold counts in some regions are a result, not a reason to change the threshold. In particular, half-sampling the small CG25 cohort is demanding and produces sparse, variable selections.

## 5. Fair comparison design

All methods received exactly the same 671 Tier2 region–gene pairs and the same regional samples. No Tier4, historical score, all-gene rank, WGCNA, or AHBA information entered fitting. Never-selected genes were assigned a shared worst-rank block for correlation; they were not silently dropped. Top-set comparisons used top 10, top 25, and `ceil(20% × regional Tier2 universe)`.

## 6. Rank agreement

{markdown(rounded(rank_summary))}

RRA–LASSO Spearman agreement ranged from {correlations.loc[(correlations.method_left.eq('Tier2-RRA')) & (correlations.method_right.eq('LASSO')), 'spearman_rho'].min():.3f} to {correlations.loc[(correlations.method_left.eq('Tier2-RRA')) & (correlations.method_right.eq('LASSO')), 'spearman_rho'].max():.3f}. RRA–Stability Selection agreement ranged from {correlations.loc[(correlations.method_left.eq('Tier2-RRA')) & (correlations.method_right.eq('Stability Selection')), 'spearman_rho'].min():.3f} to {correlations.loc[(correlations.method_left.eq('Tier2-RRA')) & (correlations.method_right.eq('Stability Selection')), 'spearman_rho'].max():.3f}. These are moderate-to-strong within-study convergence signals, not independent replication.

## 7. Top-set overlap

{markdown(rounded(overlap_summary))}

The frozen RRA top 10 shared 7–9 genes with LASSO and 7–9 genes with Stability Selection across regions. Thus the strongest candidates are not artifacts of one aggregation method, although ordering and lower-ranked candidates differ.

## 8. Resampling stability

{markdown(rounded(resampling_summary))}

Split-level five-model RRA had higher mean all-candidate rank stability than LASSO alone in every region and substantially higher stability than half-sampling Stability Selection. Stability Selection was especially variable for CG25 and OFC BA11. This does not invalidate its continuous final ranking; it shows that individual half-sample selections are sensitive to the small cohorts.

## 9. Coverage and sparsity

{markdown(rounded(coverage))}

RRA ranks the full Tier2 universe because dense L2/SVM models contribute evidence. LASSO is deliberately sparse. Stability Selection gives more genes at least occasional support, but generally low selection probabilities. Coverage is descriptive and is not scored as better or worse.

## 10. Descriptive predictive metrics

{markdown(rounded(performance))}

These are optimistic within-study diagnostics, not unbiased out-of-sample estimates. They do not decide which prioritisation method is preferable.

## 11. All-gene sensitivity comparison

{markdown(rounded(allgene))}

This comparison uses the frozen all-gene RRA only after benchmark rankings were fixed. It is within-study sensitivity/convergence, not independent validation.

## 12. Top genes by method

{markdown(top_five)}

## 13. Interpretation

The evidence best matches **Scenario E, with an important Scenario B nuance**. LASSO alone recovers much of the RRA shortlist and therefore provides a strong, simpler benchmark. However, split-level RRA rankings are consistently more stable than LASSO rankings, indicating that the ensemble contributes useful consensus beyond one sparse model. Stability Selection converges substantially with the final RRA top genes but is not more stable under the prespecified half-sampling design; it therefore supports rather than challenges retaining RRA.

This is not proof of biological truth. All methods reuse the same small cohort and the same Diagnosis-informed Tier2 screen.

## 14. Limitations

- Tier2 was defined in the complete cohort, so downstream performance metrics are optimistic.
- Resampling observations are not independent biological replications.
- CG25 has only 27 samples; half-samples are especially small.
- Stability Selection probabilities depend on the prespecified C grid and half-sampling design.
- Dense-model RRA and sparse methods answer slightly different questions about coverage.
- No external cohort was tested.

## 15. Recommendation for the next phase

Retain Tier2–RRA as the primary continuous consensus ranking, present LASSO as the transparent simpler benchmark, and present Stability Selection as a supportive sensitivity analysis. Emphasize that LASSO recovers most leading candidates, while RRA provides measurably greater resampling stability. The next prospective phase should be Tier2 conditional label-permutation/null testing. That phase was not run here. Limma/voom, DESeq2, and external validation were also not run.

## Runtime and complexity

{markdown(rounded(runtime))}
"""
    (ROOT / "docs/TIER2_ML_PRIORITIZATION_BENCHMARK_REPORT.md").write_bytes(report.encode("utf-8"))

    for path in OUTPUT.rglob("*"):
        if path.is_file() and path.suffix.lower() in {".csv", ".json", ".txt", ".sha256"}:
            path.write_bytes(path.read_bytes().replace(b"\r\n", b"\n"))
    targets = sorted(
        path for path in OUTPUT.rglob("*") if path.is_file()
        and "logs/checkpoints" not in path.as_posix()
        and path.name not in {"checksums.sha256", "benchmark_validation_results.json"}
    )
    lines = [f"{sha256(path)}  {path.relative_to(OUTPUT).as_posix()}" for path in targets]
    (OUTPUT / "checksums.sha256").write_bytes(("\n".join(lines) + "\n").encode("utf-8"))
    print(f"Built benchmark report and {len(lines)} checksums")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
