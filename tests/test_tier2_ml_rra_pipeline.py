from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code" / "ml_rra"))
RUNNER_PATH = ROOT / "code" / "tier2_ml_rra" / "01_run_tier2_repeated_ml.py"
SPEC = importlib.util.spec_from_file_location("tier2_runner", RUNNER_PATH)
runner = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)


class TestTier2MLRRAPipeline(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(71)
        self.x = rng.normal(size=(40, 12))
        self.y = np.array([0] * 20 + [1] * 20)
        self.x[self.y == 1, 2] += 2.0
        self.genes = np.array([f"g{i}" for i in range(12)])

    def test_canonical_source_has_671_expected_pairs(self):
        tier2 = runner.load_tier2_identities(ROOT)
        self.assertEqual(len(tier2), 671)

    def test_fixed_pipeline_keeps_every_feature(self):
        pipe = runner.fit_fixed_pipeline(self.x, self.y, "logistic_l2", {"C": 0.1}, 3)
        self.assertEqual(pipe.transform(self.x).shape, self.x.shape)

    def test_external_test_perturbation_cannot_change_training_fit(self):
        first = runner.fit_fixed_pipeline(self.x[:30], self.y[:30], "logistic_l2", {"C": 0.1}, 3)
        altered = self.x.copy()
        altered[30:] += 1_000_000
        second = runner.fit_fixed_pipeline(altered[:30], self.y[:30], "logistic_l2", {"C": 0.1}, 3)
        np.testing.assert_allclose(first.model.coef_, second.model.coef_)
        np.testing.assert_allclose(first.imputer.statistics_, second.imputer.statistics_)

    def test_random_forest_importance_has_no_sign(self):
        pipe = runner.fit_fixed_pipeline(
            self.x, self.y, "random_forest", {"max_depth": 3, "min_samples_leaf": 2}, 3
        )
        table = runner.importance_table(pipe, "random_forest", self.genes, {})
        self.assertTrue(table["sign"].isna().all())

    def test_sparse_zeros_are_not_ranked(self):
        pipe = runner.fit_fixed_pipeline(self.x, self.y, "logistic_l1", {"C": 0.01}, 3)
        table = runner.importance_table(pipe, "logistic_l1", self.genes, {})
        self.assertTrue((table["absolute_importance"] > runner.ZERO_TOLERANCE).all())

    def test_l2_ranks_all_fixed_tier2_features(self):
        pipe = runner.fit_fixed_pipeline(self.x, self.y, "logistic_l2", {"C": 0.1}, 3)
        table = runner.importance_table(pipe, "logistic_l2", self.genes, {})
        self.assertEqual(set(table["gene_id"]), set(self.genes))

    def test_no_selector_or_topk_logic_in_runner(self):
        source = RUNNER_PATH.read_text(encoding="utf-8")
        forbidden_calls = ["ranked_feature_indices(", "selector_scores(", "inner_tune(", "selected_k"]
        self.assertFalse(any(term in source for term in forbidden_calls))

    def test_dummy_is_not_in_model_order(self):
        self.assertNotIn(runner.DUMMY_MODEL, runner.MODEL_ORDER)

    def test_rra_source_excludes_dummy_and_performance_weighting(self):
        source = (ROOT / "code" / "tier2_ml_rra" / "02_run_tier2_rra.R").read_text(encoding="utf-8")
        self.assertNotIn("balanced_accuracy", source)
        self.assertIn('models_expected <- c("logistic_l1", "elastic_net", "logistic_l2", "linear_svm", "random_forest")', source)

    def test_empty_lists_are_built_from_performance_ledger(self):
        source = (ROOT / "code" / "tier2_ml_rra" / "02_run_tier2_rra.R").read_text(encoding="utf-8")
        self.assertIn("fit_ids <- unique", source)
        self.assertIn("performance$fit_id", source)


if __name__ == "__main__":
    unittest.main()
