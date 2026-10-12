#!/usr/bin/env python3
"""Create the supervisor report, observed reference, and checksum ledger."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "outputs/length_gt_100_clean_pipeline/20_Tier2_ML_conditional_null"
SUMMARY = OUTPUT / "summary"
REPORT = ROOT / "docs/TIER2_ML_CONDITIONAL_NULL_REPORT.md"
OBSERVED = ROOT / "outputs/length_gt_100_clean_pipeline/18_Tier2_ML_RRA_prioritization"
BENCHMARK = ROOT / "outputs/length_gt_100_clean_pipeline/19_Tier2_ML_benchmark"
ALLGENE = ROOT / "outputs/length_gt_100_clean_pipeline/17_ML_RRA_consensus/consensus_rra_rankings.csv"
REGIONS = ["Ventral Striatum", "Insula", "Hippocampal Formation", "CG25", "dlPFC BA8/9", "OFC BA11"]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def metric(table: pd.DataFrame, family: str, region: str, name: str) -> pd.Series:
    return table.loc[table["family"].eq(family) & table["region"].eq(region) & table["metric"].eq(name)].iloc[0]


def f(value: float, digits: int = 3) -> str:
    return f"{float(value):.{digits}f}"


def p(value: float) -> str:
    return f"{float(value):.4f}"


def main() -> int:
    comparison = pd.read_csv(SUMMARY / "observed_vs_null_summary.csv")
    prediction = pd.read_csv(SUMMARY / "best_of_five_prediction_null.csv")
    top = pd.read_csv(SUMMARY / "RRA_observed_topgene_null_behavior.csv", dtype={"gene_id": str})
    run_manifest = json.loads((OUTPUT / "run_manifest.json").read_text(encoding="utf-8"))
    reference = {
        "role": "read-only observed references; no files copied or modified",
        "tier2_rra_consensus": {"path": str((OBSERVED / "Tier2_consensus_rra_rankings.csv").relative_to(ROOT)).replace("\\", "/"), "sha256": sha256(OBSERVED / "Tier2_consensus_rra_rankings.csv")},
        "tier2_rra_resample_rankings": {"path": str((BENCHMARK / "RRA_resample_consensus_rankings.csv").relative_to(ROOT)).replace("\\", "/"), "sha256": sha256(BENCHMARK / "RRA_resample_consensus_rankings.csv")},
        "lasso_ranking": {"path": str((BENCHMARK / "LASSO_gene_rankings.csv").relative_to(ROOT)).replace("\\", "/"), "sha256": sha256(BENCHMARK / "LASSO_gene_rankings.csv")},
        "allgene_rra": {"path": str(ALLGENE.relative_to(ROOT)).replace("\\", "/"), "sha256": sha256(ALLGENE)},
    }
    observed_dir = OUTPUT / "observed_reference"; observed_dir.mkdir(parents=True, exist_ok=True)
    (observed_dir / "observed_reference_manifest.json").write_text(json.dumps(reference, indent=2, sort_keys=True), encoding="utf-8")

    stability_rows = []
    lasso_rows = []
    top_rows = []
    prediction_rows = []
    for region in REGIONS:
        rra_s = metric(comparison, "RRA resample stability", region, "spearman_mean")
        rra_10 = metric(comparison, "RRA resample stability", region, "jaccard_mean_top10")
        rra_20 = metric(comparison, "RRA resample stability", region, "jaccard_mean_top20percent")
        stability_rows.append(f"| {region} | {f(rra_s.observed_value)} vs {f(rra_s.null_mean)} | {p(rra_s.empirical_p_pilot)} | {f(rra_10.observed_value)} vs {f(rra_10.null_mean)} | {p(rra_10.empirical_p_pilot)} | {f(rra_20.observed_value)} vs {f(rra_20.null_mean)} | {p(rra_20.empirical_p_pilot)} |")
        lmax = metric(comparison, "LASSO", region, "maximum_selection_frequency")
        lrank = metric(comparison, "LASSO", region, "spearman_mean")
        lj10 = metric(comparison, "LASSO", region, "jaccard_mean_top10")
        lasso_rows.append(f"| {region} | {f(lmax.observed_value)} vs {f(lmax.null_mean)} | {p(lmax.empirical_p_pilot)} | {f(lrank.observed_value)} vs {f(lrank.null_mean)} | {p(lrank.empirical_p_pilot)} | {f(lj10.observed_value)} vs {f(lj10.null_mean)} | {p(lj10.empirical_p_pilot)} |")
        lead = top.loc[top["region"].eq(region) & top["observed_rank"].eq(1)].iloc[0]
        top_rows.append(f"| {region} | {lead.gene_symbol or lead.gene_id} | {int(lead.null_top10_count)}/20 | {f(lead.null_rank_median, 1)} ({int(lead.null_rank_best)}–{int(lead.null_rank_worst)}) | {p(lead.exploratory_rank_empirical_p_pilot)} |")
        pred = prediction.loc[prediction["region"].eq(region)].iloc[0]
        prediction_rows.append(f"| {region} | {f(pred.observed_best_of_five_balanced_accuracy)} | {f(pred.null_mean)} | {f(pred.null_q025)}–{f(pred.null_q975)} | {p(pred.empirical_p_pilot)} |")

    runtime = run_manifest["runtime_seconds_ml"]
    report = f"""# Tier2 ML conditional-null pilot report

## 1. Scientific question

Given the already-fixed Tier2 candidate universe, does the original MDD/Control labeling produce more coherent and reproducible downstream ML gene prioritisation than labels shuffled within each region?

This 20-permutation pilot completed successfully for all six analysis units. It reran the full five-model Tier2-RRA pipeline and the exact LASSO-alone benchmark under each permuted label vector. No observed fit, preprocessing object, hyperparameter, coefficient, importance, ranking, or RRA result was reused.

## 2. Why the null is conditional

**This is not an end-to-end global-null test.** Tier2 was selected using the original Diagnosis labels in the full cohort, so the fixed 671-pair universe is already enriched for observed-label association. The experiment asks whether downstream ML remains structured with respect to those original labels after the candidate set is fixed. It cannot answer whether the full expression → gene-wise model → Tier2 → ML pipeline would produce the result under a global null.

This design is therefore expected to favor the observed labels. Any separation is conditional evidence for label-dependent prioritisation, not independent validation, a biomarker claim, or a resolved gene-level significance test.

## 3. Frozen Tier2 candidate universe

The exact 671 region–gene pairs were preserved: Ventral Striatum 146, Insula 70, Hippocampal Formation 127, CG25 186, dlPFC BA8/9 72, and OFC BA11 70. Expression values, sample IDs, gene identities, and class counts were unchanged. Only labels were shuffled within region using deterministic seeds 41001–41020.

## 4. Permutation procedure

For every seed the pipeline regenerated five-fold × five-repeat stratified outer CV, training-only median imputation/scaling, three-fold inner tuning on balanced accuracy, all five models, 750 feature-model fits, 150 DummyClassifier prediction references, fold gene rankings, two-level Robust Rank Aggregation, split-level RRA, and the LASSO selection-frequency ranking. The DummyClassifier never contributed gene evidence.

All 20 permutations completed; none failed. A completed seed was rerun with `--resume`, verified by checksums, and skipped without recomputation. Single-seed ML wall times ranged from {f(runtime['minimum'], 1)} to {f(runtime['maximum'], 1)} seconds (median {f(runtime['median'], 1)}; summed per-seed wall time {f(runtime['sum_permutation_wall_times'], 1)} seconds). Parallel execution affects wall time only, not method or seeds.

Pilot empirical probabilities use `(1 + null as-or-more-extreme) / 21`; the minimum possible value is 0.0476.

## 5. Tier2-RRA null results

RRA evidence was regionally heterogeneous. Ventral Striatum showed the clearest separation from shuffled labels across rank correlation and all top-set overlaps. Hippocampal Formation, dlPFC BA8/9, and CG25 showed separation for some broader top-set endpoints but not uniformly across rank and top-10 measures. Insula was suggestive but did not reach the pilot floor. OFC BA11 observed RRA resampling stability was similar to or worse than the shuffled-label reference.

| Region | Mean rank Spearman: observed vs null | Pilot p | Top-10 Jaccard: observed vs null | Pilot p | Top-20% Jaccard: observed vs null | Pilot p |
|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(stability_rows)}

Consensus-score concentration was also mixed. Several broader top-set median RRA-score summaries were stronger than the null, but the minimum RRA score alone was usually not unusual. RRA scores are aggregation scores, not MDD association p-values. Counts supported by at least three models were uninformative because the dense models supplied broad evidence to every candidate.

## 6. LASSO null results

Observed LASSO selection-frequency concentration exceeded the null in every region, especially for the leading 10–25 genes. Resample ranking stability was clearer in Ventral Striatum and dlPFC, moderate in Hippocampal Formation and Insula, and not better than null in CG25 or OFC BA11.

| Region | Max frequency: observed vs null | Pilot p | Mean rank Spearman: observed vs null | Pilot p | Top-10 Jaccard: observed vs null | Pilot p |
|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(lasso_rows)}

Thus the simple sparse model supports label-dependent concentration, but it does not establish uniformly reproducible rankings in all regions.

## 7. Ranking stability

The primary endpoint does not support one project-wide yes/no statement. Ventral Striatum provides clear conditional stability evidence. Hippocampal Formation and dlPFC provide partial evidence, depending on top-set breadth. CG25 and Insula are mixed. OFC BA11 fails the predeclared ranking-stability criterion despite strong predictive separation.

This regional heterogeneity is scientifically important and should not be hidden by pooling endpoints or constructing a weighted score.

## 8. Observed top-gene behavior

Leading observed genes were usually ranked better than their own null median, but several reappeared near the top under multiple shuffled labels. Gene-specific probabilities are exploratory only; with 20 permutations they have coarse resolution and no BH/FDR interpretation.

| Region | Observed rank-1 gene | Null top-10 recurrence | Null median rank (range) | Exploratory rank-tail p |
|---|---|---:|---:|---:|
{chr(10).join(top_rows)}

Full top-25 rank and RRA-score distributions are preserved in `RRA_observed_topgene_null_behavior.csv`; the corresponding LASSO distributions are in `LASSO_observed_topgene_null_behavior.csv`.

## 9. Secondary prediction diagnostics

Observed best-of-five balanced accuracy exceeded every one of the 20 corresponding null maxima in all six regions. This is a secondary diagnostic and **must not be used as evidence of an unbiased diagnostic classifier**, because Tier2 was selected from the same cohort using real Diagnosis.

| Region | Observed best-of-five BA | Null mean | Null 95% interval | Pilot p |
|---|---:|---:|---:|---:|
{chr(10).join(prediction_rows)}

## 10. Interpretation using the predeclared scenarios

No single scenario describes every region:

- Ventral Striatum is closest to **Scenario A**: RRA and LASSO stability/concentration both exceed the conditional null.
- Hippocampal Formation and dlPFC BA8/9 provide partial Scenario A support, with stronger separation for broader top sets and LASSO concentration than for every RRA endpoint.
- CG25 and Insula are mixed: leading-gene concentration exists, but ranking-stability endpoints are inconsistent.
- OFC BA11 is closest to **Scenario E**: prediction strongly exceeds null, while RRA and LASSO resampling stability do not. Prediction alone cannot justify OFC gene prioritisation.

Overall, the pilot provides conditional evidence for label-dependent prioritisation in some regions, not uniform validation of the entire Tier2-RRA ranking.

## 11. Limitations

1. Tier2 was selected with the observed labels, so the null is conditional and favoring of the observed result is expected.
2. Twenty permutations give a minimum empirical probability of 0.0476 and cannot resolve smaller tails.
3. Multiple separate endpoints and regions were inspected; the pilot probabilities are descriptive and are not multiplicity-adjusted claims.
4. Small regional sample sizes and correlated expression features can destabilize gene-level rankings.
5. RRA scores are rank-aggregation evidence, not disease-association p-values.
6. Prediction metrics are optimistic within-study diagnostics and not external validation.

## 12. Recommendation

**Stop at the completed 20-permutation pilot for now and proceed to the upstream statistical-method benchmark while retaining this analysis as conditional descriptive evidence.** Extending to 99 would improve tail resolution, but it would not repair the central candidate-selection bias and would require substantial additional computation. A 99-permutation extension should be reserved for a later explicit decision if precise conditional tails become necessary.

The primary report should retain cautious, region-specific language. It should not claim uniformly stable ML prioritisation, and OFC BA11 gene-level ML conclusions should be softened substantially. No validated MDD genes, biomarkers, or independent validation are established here.

## Resume command

The completed pilot can be verified/resumed without recomputing finished seeds:

```bat
.runtime\\python\\python.exe code\\tier2_conditional_null\\04_run_pilot.py --start 1 --end 20 --workers 3 --resume
```
"""
    REPORT.write_text(report, encoding="utf-8")

    # Create a complete non-log checksum inventory. The ledger does not include itself.
    paths = []
    for path in OUTPUT.rglob("*"):
        if not path.is_file() or path.name == "checksums.sha256" or "logs" in path.parts:
            continue
        paths.append(path)
    paths.append(REPORT)
    rows = [f"{sha256(path)}  {str(path.relative_to(ROOT)).replace(chr(92), '/')}" for path in sorted(set(paths))]
    (OUTPUT / "checksums.sha256").write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"Built supervisor report and checksum ledger for {len(rows)} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
