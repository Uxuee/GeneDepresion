"""Regression tests for gene-level ML-priority scoring correctness fixes."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "code" / "extract_stable_ML_prioritized_Tier2_genes.py"
SPEC = importlib.util.spec_from_file_location("ml_priority", SCRIPT)
ml = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(ml)


def tier2_one() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "analysis_level": ["four_group"],
            "region_or_subregion": ["Test region"],
            "gene_id": ["1"],
            "gene_symbol": ["GENE1"],
            "beta_diagnosis": [0.2],
            "p_value_diagnosis": [0.01],
            "standardized_effect_size_diagnosis": [0.5],
            "lmm_direction": ["up_in_MDD"],
            "lmm_evidence_score": [1.0],
        }
    )


def event(model: str, coefficient: float, **extra) -> dict:
    row = {
        "analysis_level": "four_group",
        "region_or_subregion": "Test region",
        "model": model,
        "gene_id": "1",
        "coefficient": coefficient,
        "selected": True,
        "outer_repeat": 0,
        "outer_fold": 0,
        "feature_selection_method": "train_only_anova_f",
        "top_k": 25,
    }
    row.update(extra)
    return row


class MLPriorityScoringTests(unittest.TestCase):
    def test_random_forest_never_enters_sign_consistency_or_direction(self) -> None:
        raw = pd.DataFrame(
            [
                event("logistic_l2", -0.4),
                event("random_forest_exploratory", 100.0),
            ]
        )
        events = ml.annotate_importance_events(raw, "v2")
        summary = ml.summarize_events(events, tier2_one(), "v2").iloc[0]
        self.assertEqual(summary["sign_consistency_v2"], 1.0)
        self.assertEqual(summary["ml_direction_v2"], "down_in_MDD")

        rf_only = ml.annotate_importance_events(pd.DataFrame([event("random_forest_exploratory", 1.0)]), "v2")
        rf_summary = ml.summarize_events(rf_only, tier2_one(), "v2").iloc[0]
        self.assertTrue(pd.isna(rf_summary["sign_consistency_v2"]))
        self.assertEqual(rf_summary["ml_direction_v2"], "unavailable")

    def test_majority_baseline_never_enters_v2_gene_evidence(self) -> None:
        raw = pd.DataFrame(
            [
                event("majority_class_baseline", np.nan),
                event("logistic_l1", 0.3),
            ]
        )
        events = ml.annotate_importance_events(raw, "v2")
        self.assertNotIn("majority_class_baseline", set(events["model"]))
        summary = ml.summarize_events(events, tier2_one(), "v2").iloc[0]
        self.assertEqual(summary["selected_event_count"], 1)
        self.assertEqual(summary["n_events_v2"], 1)
        self.assertEqual(summary["selection_frequency_v2"], 1.0)

    def test_missing_sign_consistency_cannot_pass_high_confidence_gate(self) -> None:
        tier2 = tier2_one()
        v1 = pd.DataFrame(
            {
                "analysis_level": ["four_group"],
                "region_or_subregion": ["Test region"],
                "gene_id": ["1"],
                "selection_frequency_v1": [1.0],
                "model_support_count_v1": [2],
                "support_models_v1": ["logistic_l1;logistic_l2"],
                "sparse_model_support_v1": [True],
                "linear_model_support_v1": [True],
                "rf_support_v1": [False],
                "mean_abs_importance_v1": [1.0],
                "median_abs_importance_v1": [1.0],
                "mean_rank_v1": [1.0],
                "rank_stability_v1": [1.0],
                "sign_consistency_v1": [np.nan],
                "ml_direction_v1": ["unavailable"],
                "n_events_v1": [10],
            }
        )
        result = ml.build_priority_table(tier2, v1, None).iloc[0]
        self.assertGreaterEqual(result["ML_priority_score"], 0.70)
        self.assertNotEqual(result["ML_priority_class"], "High-confidence ML-prioritized Tier 2")
        self.assertEqual(result["sign_consistency_score"], 0.5)

    def test_score_uses_intended_weights(self) -> None:
        observed = sum(ml.SCORE_WEIGHTS.values())
        self.assertAlmostEqual(observed, 1.0)
        self.assertEqual(
            ml.SCORE_WEIGHTS,
            {
                "selection_frequency_combined": 0.40,
                "model_support_scaled": 0.15,
                "sparse_model_support_binary": 0.10,
                "importance_rank_score": 0.15,
                "sign_consistency_score": 0.10,
                "lmm_evidence_scaled": 0.05,
                "v2_support_score": 0.05,
            },
        )

    def test_v2_denominator_matches_actual_nonbaseline_event_grid(self) -> None:
        rows = []
        for method in ["anova", "effect"]:
            for top_k in [25, 50]:
                for fold in [0, 1]:
                    for model in ["logistic_l1", "random_forest_exploratory", "majority_class_baseline"]:
                        rows.append(
                            event(
                                model,
                                0.2 if model != "majority_class_baseline" else np.nan,
                                feature_selection_method=method,
                                top_k=top_k,
                                outer_fold=fold,
                            )
                        )
        events = ml.annotate_importance_events(pd.DataFrame(rows), "v2")
        summary = ml.summarize_events(events, tier2_one(), "v2").iloc[0]
        expected = 2 * 2 * 2 * 2  # methods * top-k * folds * non-baseline models
        self.assertEqual(summary["n_events_v2"], expected)
        self.assertEqual(summary["selected_event_count"], expected)
        self.assertEqual(summary["selection_frequency_v2"], 1.0)


if __name__ == "__main__":
    unittest.main()
