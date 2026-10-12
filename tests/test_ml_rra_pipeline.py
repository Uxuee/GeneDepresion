from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code" / "ml_rra"))

from common import (  # noqa: E402
    MODEL_ORDER,
    extract_importance,
    fit_training_pipeline,
    make_outer_splits,
    ranked_feature_indices,
    stable_seed,
)


class TestMLRRAPipeline(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(44)
        self.x = rng.normal(size=(40, 20))
        self.y = np.array([0] * 20 + [1] * 20)
        self.x[self.y == 1, 3] += 2.5
        self.genes = np.array([f"g{i}" for i in range(20)])

    def test_outer_coverage_is_exactly_once_per_repeat(self):
        seen = {}
        for repeat, _, _, _, test in make_outer_splits(self.y, n_splits=5, seeds=[2026, 2027]):
            for index in test:
                seen[(repeat, int(index))] = seen.get((repeat, int(index)), 0) + 1
        self.assertEqual(set(seen.values()), {1})
        self.assertEqual(len(seen), 80)

    def test_selector_uses_supplied_training_data_only(self):
        first = ranked_feature_indices(self.x[:30], self.y[:30], "anova_f")
        changed_test = self.x.copy()
        changed_test[30:, 10] += 1_000_000
        second = ranked_feature_indices(changed_test[:30], self.y[:30], "anova_f")
        np.testing.assert_array_equal(first, second)

    def test_imputation_scaling_and_selection_fit(self):
        x = self.x.copy()
        x[0, 2] = np.nan
        pipe = fit_training_pipeline(x, self.y, "anova_f", 10, "logistic_l2", {"C": 0.1}, 5)
        self.assertEqual(pipe.transform(x).shape, (40, 10))

    def test_random_forest_is_nondirectional(self):
        pipe = fit_training_pipeline(self.x, self.y, "anova_f", 10, "random_forest", {"max_depth": 3, "min_samples_leaf": 2}, 5)
        table = extract_importance(pipe, "random_forest", self.genes, {})
        self.assertTrue(table["sign"].isna().all())

    def test_sparse_linear_zero_coefficients_are_omitted(self):
        pipe = fit_training_pipeline(self.x, self.y, "anova_f", 20, "logistic_l1", {"C": 0.01}, 5)
        table = extract_importance(pipe, "logistic_l1", self.genes, {})
        self.assertTrue((table["absolute_importance"] > 1e-12).all())

    def test_importance_mapping_preserves_gene_ids(self):
        pipe = fit_training_pipeline(self.x, self.y, "anova_f", 10, "logistic_l2", {"C": 0.1}, 5)
        table = extract_importance(pipe, "logistic_l2", self.genes, {})
        expected = set(self.genes[pipe.selected_indices])
        self.assertEqual(set(table["gene_id"]), expected)

    def test_seed_derivation_is_deterministic(self):
        self.assertEqual(stable_seed("x", 1, "y"), stable_seed("x", 1, "y"))
        self.assertNotEqual(stable_seed("x", 1, "y"), stable_seed("x", 2, "y"))

    def test_exact_five_model_families(self):
        self.assertEqual(
            MODEL_ORDER,
            ["logistic_l1", "elastic_net", "logistic_l2", "linear_svm", "random_forest"],
        )

    def test_training_code_has_no_legacy_priority_fields(self):
        source = (ROOT / "code" / "ml_rra" / "01_run_nested_ml_rankings.py").read_text(encoding="utf-8")
        common = (ROOT / "code" / "ml_rra" / "common.py").read_text(encoding="utf-8")
        forbidden = ["Tier2", "Tier4", "confidence_label", "priority_score"]
        executable_training_section = source.split("    manifest = {", 1)[0] + common
        self.assertFalse(any(term in executable_training_section for term in forbidden))

    def test_rra_does_not_use_sign_or_performance_weights(self):
        source = (ROOT / "code" / "ml_rra" / "02_run_rra.R").read_text(encoding="utf-8")
        self.assertNotIn("balanced_accuracy", source)
        self.assertNotIn("$sign", source)


if __name__ == "__main__":
    unittest.main()
