"""Create before/after tables for the corrected ML-priority scores."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OLD_DIR = ROOT / "outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_extended"
NEW_DIR = ROOT / "outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_corrected"
SCORER_PATH = ROOT / "code/extract_stable_ML_prioritized_Tier2_genes.py"

spec = importlib.util.spec_from_file_location("scorer", SCORER_PATH)
scorer = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(scorer)

KEYS = ["analysis_level", "region_or_subregion", "gene_id"]


def short_class(value: str) -> str:
    if str(value).startswith("High"):
        return "high"
    if str(value).startswith("Moderate"):
        return "moderate"
    return "low"


def classify(frame: pd.DataFrame, missing_sign_passes: bool) -> pd.Series:
    agreement = frame["lmm_ml_sign_agreement"].map(
        lambda value: True if pd.isna(value) else bool(value)
    )
    sign_gate = (
        frame["sign_consistency"].fillna(1.0) >= 0.70
        if missing_sign_passes
        else frame["sign_consistency"].notna() & (frame["sign_consistency"] >= 0.70)
    )
    high = (
        (frame["ML_priority_score"] >= 0.70)
        & (frame["selection_frequency_combined"] >= 0.50)
        & (frame["model_support_count"] >= 2)
        & frame["sparse_model_support"]
        & sign_gate
        & agreement
    )
    moderate = (
        (frame["ML_priority_score"] >= 0.45)
        & (frame["selection_frequency_combined"] >= 0.25)
        & (frame["model_support_count"] >= 1)
        & agreement
    )
    return pd.Series(
        np.select([high, moderate], ["high", "moderate"], default="low"),
        index=frame.index,
    )


def scenario(
    tier2: pd.DataFrame,
    v1_raw: pd.DataFrame,
    v2_raw: pd.DataFrame,
    signed_only: bool,
    remove_baseline: bool,
    missing_sign_passes: bool,
    legacy_weights: bool,
) -> pd.DataFrame:
    v1_events = scorer.annotate_importance_events(v1_raw, "v1")
    v2_events = scorer.annotate_importance_events(
        v2_raw, "v2" if remove_baseline else "legacy_v2"
    )
    if not signed_only:
        v1_events["signed_importance_value"] = v1_events["importance_value"]
        v2_events["signed_importance_value"] = v2_events["importance_value"]
    v1_summary = scorer.summarize_events(v1_events, tier2, "v1")
    v2_summary = scorer.summarize_events(v2_events, tier2, "v2")
    result = scorer.build_priority_table(tier2, v1_summary, v2_summary)
    if legacy_weights:
        result["ML_priority_score"] = (
            0.30 * result["selection_frequency_combined"]
            + 0.20 * result["model_support_scaled"]
            + 0.15 * result["sparse_model_support_binary"]
            + 0.15 * result["importance_rank_score"]
            + 0.10 * result["sign_consistency_score"]
            + 0.05 * result["lmm_evidence_scaled"]
            + 0.05 * result["v2_support_score"]
        ).clip(0, 1)
    result["class_simple"] = classify(result, missing_sign_passes)
    return result


def count_row(name: str, frame: pd.DataFrame) -> dict:
    counts = frame["class_simple"].value_counts()
    return {
        "scenario": name,
        "high": int(counts.get("high", 0)),
        "moderate": int(counts.get("moderate", 0)),
        "low": int(counts.get("low", 0)),
    }


def main() -> None:
    used: list[str] = []
    missing: list[str] = []
    tier2 = scorer.load_tier2_tables(used, missing)
    v1_raw = pd.read_csv(
        scorer.V1_IMPORTANCE_DIR / "fold_level_feature_importance.csv",
        dtype={"gene_id": str},
    )
    v1_raw = v1_raw[v1_raw["feature_set"].eq("Tier2_fixed")].copy()
    v2_raw = pd.read_csv(
        scorer.V2_STABILITY_DIR / "selected_genes_per_fold_all_models.csv",
        dtype={"gene_id": str},
    )

    scenarios = {
        "legacy_reconstructed": scenario(tier2, v1_raw, v2_raw, False, False, True, True),
        "verified_weights_only": scenario(tier2, v1_raw, v2_raw, False, False, True, False),
        "RF_sign_fix_only": scenario(tier2, v1_raw, v2_raw, True, False, True, True),
        "baseline_removal_only": scenario(tier2, v1_raw, v2_raw, False, True, True, True),
        "missing_sign_gate_only": scenario(tier2, v1_raw, v2_raw, False, False, False, True),
        "fully_corrected": scenario(tier2, v1_raw, v2_raw, True, True, False, False),
    }
    legacy = scenarios["legacy_reconstructed"]
    corrected = scenarios["fully_corrected"]

    old = pd.read_csv(
        OLD_DIR / "ML_prioritized_Tier2_genes_all_ranked.csv", dtype={"gene_id": str}
    )
    new = pd.read_csv(
        NEW_DIR / "ML_prioritized_Tier2_genes_all_ranked.csv", dtype={"gene_id": str}
    )
    old["class_simple"] = old["ML_priority_class"].map(short_class)
    new["class_simple"] = new["ML_priority_class"].map(short_class)

    check = new[KEYS + ["ML_priority_score", "class_simple"]].merge(
        corrected[KEYS + ["ML_priority_score", "class_simple"]],
        on=KEYS,
        suffixes=("_csv", "_scenario"),
        validate="one_to_one",
    )
    new_max_diff = float(
        (check["ML_priority_score_csv"] - check["ML_priority_score_scenario"]).abs().max()
    )
    if new_max_diff > 1e-12 or not (
        check["class_simple_csv"] == check["class_simple_scenario"]
    ).all():
        raise AssertionError("Corrected CSV failed reconstruction.")

    comparison = old[
        KEYS
        + [
            "gene_symbol",
            "ML_priority_score",
            "ML_priority_class",
            "sign_consistency",
            "selection_frequency_v2",
        ]
    ].merge(
        new[
            KEYS
            + [
                "ML_priority_score",
                "ML_priority_class",
                "sign_consistency",
                "selection_frequency_v2",
                "n_events_v1",
                "n_events_v2",
            ]
        ],
        on=KEYS,
        suffixes=("_old", "_corrected"),
        validate="one_to_one",
    )
    comparison["old_class"] = comparison["ML_priority_class_old"].map(short_class)
    comparison["corrected_class"] = comparison["ML_priority_class_corrected"].map(short_class)
    comparison["transition"] = comparison["old_class"] + "_to_" + comparison["corrected_class"]
    comparison["score_change"] = (
        comparison["ML_priority_score_corrected"] - comparison["ML_priority_score_old"]
    )
    comparison["absolute_score_change"] = comparison["score_change"].abs()
    comparison["sign_consistency_change"] = (
        comparison["sign_consistency_corrected"] - comparison["sign_consistency_old"]
    )
    comparison["v2_frequency_change"] = (
        comparison["selection_frequency_v2_corrected"]
        - comparison["selection_frequency_v2_old"]
    )
    comparison["old_region_rank"] = comparison.groupby("region_or_subregion")[
        "ML_priority_score_old"
    ].rank(method="min", ascending=False)
    comparison["corrected_region_rank"] = comparison.groupby("region_or_subregion")[
        "ML_priority_score_corrected"
    ].rank(method="min", ascending=False)
    comparison["rank_change"] = (
        comparison["corrected_region_rank"] - comparison["old_region_rank"]
    )
    comparison["absolute_rank_change"] = comparison["rank_change"].abs()
    comparison.to_csv(NEW_DIR / "before_after_all_genes.csv", index=False)

    transitions = (
        comparison.groupby("transition", as_index=False)
        .size()
        .rename(columns={"size": "n_gene_region_rows"})
        .sort_values("transition")
    )
    transitions.to_csv(NEW_DIR / "confidence_transition_summary.csv", index=False)
    comparison.nlargest(20, "absolute_score_change").to_csv(
        NEW_DIR / "top20_absolute_score_changes.csv", index=False
    )
    comparison.nlargest(20, "absolute_rank_change").to_csv(
        NEW_DIR / "top20_absolute_rank_changes.csv", index=False
    )
    comparison.assign(
        absolute_sign_change=comparison["sign_consistency_change"].abs()
    ).nlargest(20, "absolute_sign_change").to_csv(
        NEW_DIR / "top20_RF_sign_consistency_changes.csv", index=False
    )
    comparison.assign(
        absolute_v2_change=comparison["v2_frequency_change"].abs()
    ).nlargest(20, "absolute_v2_change").to_csv(
        NEW_DIR / "top20_baseline_removal_v2_frequency_changes.csv", index=False
    )

    scenario_counts = pd.DataFrame(
        [count_row(name, frame) for name, frame in scenarios.items()]
    )
    scenario_counts.to_csv(NEW_DIR / "isolated_correction_confidence_counts.csv", index=False)
    legacy_indexed = legacy.set_index(KEYS)
    effect_rows = []
    for name, frame in scenarios.items():
        if name == "legacy_reconstructed":
            continue
        joined = legacy_indexed[["ML_priority_score", "class_simple"]].join(
            frame.set_index(KEYS)[["ML_priority_score", "class_simple"]],
            lsuffix="_legacy",
            rsuffix="_scenario",
        )
        difference = (
            joined["ML_priority_score_scenario"] - joined["ML_priority_score_legacy"]
        ).abs()
        effect_rows.append(
            {
                "scenario": name,
                "n_scores_changed": int((difference > 1e-12).sum()),
                "maximum_absolute_score_change": float(difference.max()),
                "n_confidence_labels_changed": int(
                    (joined["class_simple_scenario"] != joined["class_simple_legacy"]).sum()
                ),
            }
        )
    effects = pd.DataFrame(effect_rows)
    effects.to_csv(NEW_DIR / "isolated_correction_effect_summary.csv", index=False)

    frozen_check = old[KEYS + ["ML_priority_score", "class_simple"]].merge(
        legacy[KEYS + ["ML_priority_score", "class_simple"]],
        on=KEYS,
        suffixes=("_frozen", "_reconstructed"),
        validate="one_to_one",
    )
    frozen_check["score_difference"] = (
        frozen_check["ML_priority_score_reconstructed"]
        - frozen_check["ML_priority_score_frozen"]
    )
    frozen_check["class_matches"] = (
        frozen_check["class_simple_frozen"] == frozen_check["class_simple_reconstructed"]
    )
    frozen_check.reindex(
        frozen_check["score_difference"].abs().sort_values(ascending=False).index
    ).to_csv(NEW_DIR / "frozen_vs_reconstructed_legacy_check.csv", index=False)
    legacy_max_diff = float(
        (
            frozen_check["ML_priority_score_frozen"]
            - frozen_check["ML_priority_score_reconstructed"]
        )
        .abs()
        .max()
    )
    legacy_class_matches = int(
        (
            frozen_check["class_simple_frozen"]
            == frozen_check["class_simple_reconstructed"]
        ).sum()
    )

    requested = comparison[
        (comparison["region_or_subregion"] == "Cingulate gyrus 25")
        & (comparison["gene_symbol"] == "LOC107985795")
    ]
    changed = comparison[comparison["old_class"] != comparison["corrected_class"]]
    if not changed.empty:
        requested = pd.concat(
            [requested, changed.nlargest(1, "absolute_score_change")],
            ignore_index=True,
        ).drop_duplicates(KEYS)

    new_indexed = new.set_index(KEYS)
    reconstruction_rows = []
    for _, example in requested.iterrows():
        index = tuple(example[key] for key in KEYS)
        row = new_indexed.loc[index]
        arithmetic = (
            0.40 * float(row["selection_frequency_combined"])
            + 0.15 * float(row["model_support_scaled"])
            + 0.10 * float(row["sparse_model_support_binary"])
            + 0.15 * float(row["importance_rank_score"])
            + 0.10 * float(row["sign_consistency_score"])
            + 0.05 * float(row["lmm_evidence_scaled"])
            + 0.05 * float(row["v2_support_score"])
        )
        reconstruction_rows.append(
            {
                **{key: example[key] for key in KEYS},
                "gene_symbol": example["gene_symbol"],
                "old_class": example["old_class"],
                "corrected_class": example["corrected_class"],
                "v1_selected_events": int(
                    round(float(row["selection_frequency_v1"]) * float(row["n_events_v1"]))
                ),
                "v1_denominator": int(row["n_events_v1"]),
                "selection_frequency_v1": row["selection_frequency_v1"],
                "v2_selected_events": int(
                    round(float(row["selection_frequency_v2"]) * float(row["n_events_v2"]))
                ),
                "v2_denominator": int(row["n_events_v2"]),
                "selection_frequency_v2": row["selection_frequency_v2"],
                "combined_selection_frequency": row["selection_frequency_combined"],
                "model_support_scaled": row["model_support_scaled"],
                "sparse_model_support_binary": row["sparse_model_support_binary"],
                "importance_rank_score": row["importance_rank_score"],
                "corrected_sign_consistency": row["sign_consistency"],
                "sign_consistency_score": row["sign_consistency_score"],
                "lmm_evidence_scaled": row["lmm_evidence_scaled"],
                "v2_support_score": row["v2_support_score"],
                "reconstructed_score": arithmetic,
                "csv_score": row["ML_priority_score"],
                "absolute_arithmetic_difference": abs(
                    arithmetic - float(row["ML_priority_score"])
                ),
                "confidence_label": row["ML_priority_class"],
            }
        )
    reconstruction = pd.DataFrame(reconstruction_rows)
    reconstruction.to_csv(NEW_DIR / "real_gene_score_reconstruction.csv", index=False)

    old_counts = old["class_simple"].value_counts()
    new_counts = new["class_simple"].value_counts()
    report = [
        "# Corrected ML-prioritisation regeneration report",
        "",
        "## Verification",
        "",
        "- v1 SHA-256: 406ba687a9f3f4244de0d296550e1b2fcfb043f504c175c56328b447133cde1e.",
        "- v2 SHA-256: 7b33e7d0a6eb2ff2efaffdc223aced7f4722f18236a2902cc9102257b2215ef8.",
        f"- Corrected CSV reconstruction maximum difference: {new_max_diff:.3g}.",
        f"- Frozen legacy reconstruction maximum difference: {legacy_max_diff:.3g}.",
        f"- Frozen legacy class matches: {legacy_class_matches}/{len(frozen_check)}.",
        "- The frozen CSV used legacy first-three weights 0.30/0.20/0.15; the corrected CSV uses the mandated 0.40/0.15/0.10.",
        "",
        "## Counts",
        "",
        "- Old v2 denominator: 180.",
        f"- Corrected v2 denominator values: {sorted(new['n_events_v2'].dropna().unique().astype(int).tolist())}.",
        f"- Old high/moderate/low: {int(old_counts.get('high', 0))}/{int(old_counts.get('moderate', 0))}/{int(old_counts.get('low', 0))}.",
        f"- Corrected high/moderate/low: {int(new_counts.get('high', 0))}/{int(new_counts.get('moderate', 0))}/{int(new_counts.get('low', 0))}.",
        "",
        "## Confidence transitions",
        "",
        transitions.to_markdown(index=False),
        "",
        "## Isolated correction counts",
        "",
        scenario_counts.to_markdown(index=False),
        "",
        "## Isolated correction effects",
        "",
        effects.to_markdown(index=False),
        "",
        "## Real-gene arithmetic",
        "",
        reconstruction.to_markdown(index=False),
        "",
        "These remain exploratory prioritisation outputs, not validated biomarkers.",
        "",
    ]
    (NEW_DIR / "corrected_regeneration_report.md").write_text(
        "\n".join(report), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
