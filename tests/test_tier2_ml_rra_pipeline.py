from __future__ import annotations

import importlib.util
import json
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code" / "ml_rra"))
RUNNER_PATH = ROOT / "code" / "tier2_ml_rra" / "01_run_tier2_repeated_ml.py"
OUTPUT = ROOT / "outputs" / "length_gt_100_clean_pipeline" / "18_Tier2_ML_RRA_prioritization"
SPEC = importlib.util.spec_from_file_location("tier2_runner", RUNNER_PATH)
runner = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)


class TestTier2MLRRAPipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tier2 = runner.load_tier2_identities(ROOT)
        cls.tier2["region"] = cls.tier2["tier2_region"].replace(
            {"Cingulate gyrus 25": "CG25"}
        )
        cls.rankings = pd.read_csv(OUTPUT / "Tier2_fold_gene_rankings.csv", dtype={"gene_id": str})
        cls.performance = pd.read_csv(OUTPUT / "Tier2_outer_performance.csv")
        cls.consensus = pd.read_csv(OUTPUT / "Tier2_consensus_rra_rankings.csv", dtype={"gene_id": str})
        cls.manifest = json.loads((OUTPUT / "Tier2_run_manifest.json").read_text(encoding="utf-8"))

    def setUp(self):
        rng = np.random.default_rng(71)
        self.x = rng.normal(size=(40, 12))
        self.y = np.array([0] * 20 + [1] * 20)
        self.x[self.y == 1, 2] += 2.0
        self.genes = np.array([f"g{i}" for i in range(12)])

    def test_canonical_source_has_671_expected_pairs(self):
        self.assertEqual(len(self.tier2), 671)

    def test_exact_intended_six_regions(self):
        expected = {"Ventral Striatum", "Insula", "Hippocampal Formation", "CG25", "dlPFC BA8/9", "OFC BA11"}
        self.assertEqual(set(self.consensus["region"]), expected)

    def test_regional_tier2_counts_match_audited_source(self):
        expected = {"Ventral Striatum": 146, "Insula": 70, "Hippocampal Formation": 127,
                    "CG25": 186, "dlPFC BA8/9": 72, "OFC BA11": 70}
        self.assertEqual(self.tier2.groupby("region").size().to_dict(), expected)

    def test_only_identity_columns_imported_from_tier2(self):
        self.assertEqual(
            self.manifest["tier2_identity_columns_used"],
            ["region_or_subregion", "gene_id", "gene_symbol"],
        )

    def test_no_historical_ml_fields_enter_manifest_inputs(self):
        forbidden = {"priority_score", "confidence", "selection_frequency", "model_support"}
        self.assertTrue(forbidden.isdisjoint(self.manifest["tier2_identity_columns_used"]))

    def test_no_tier4_or_outside_gene_enters_consensus(self):
        expected = set(zip(self.tier2["region"], self.tier2["gene_id"].astype(str)))
        observed = set(zip(self.consensus["region"], self.consensus["gene_id"].astype(str)))
        self.assertEqual(observed, expected)

    def test_every_ranked_expression_feature_belongs_to_tier2(self):
        expected = set(zip(self.tier2["region"], self.tier2["gene_id"].astype(str)))
        observed = set(zip(self.rankings["region"], self.rankings["gene_id"].astype(str)))
        self.assertTrue(observed.issubset(expected))

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

    def test_inner_tuning_configuration_is_training_only_three_fold(self):
        self.assertEqual(self.manifest["inner_cv"]["folds"], 3)
        self.assertEqual(self.manifest["inner_cv"]["metric"], "balanced_accuracy")
        source = RUNNER_PATH.read_text(encoding="utf-8")
        self.assertIn("x_train", source)
        self.assertNotIn("x_test, y_test, model", source)

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

    def test_exact_25_outer_fits_per_region_model(self):
        feature = self.performance[self.performance["model"].isin(runner.MODEL_ORDER)]
        counts = feature.groupby(["region", "model"])["fit_id"].nunique()
        self.assertTrue(counts.eq(25).all())

    def test_exact_total_feature_and_dummy_fits(self):
        feature = self.performance[self.performance["model"].isin(runner.MODEL_ORDER)]
        dummy = self.performance[self.performance["model"].eq(runner.DUMMY_MODEL)]
        self.assertEqual(feature["fit_id"].nunique(), 750)
        self.assertEqual(dummy["fit_id"].nunique(), 150)

    def test_dummy_never_enters_gene_rankings(self):
        self.assertNotIn(runner.DUMMY_MODEL, set(self.rankings["model"]))

    def test_rra_source_excludes_dummy_and_performance_weighting(self):
        source = (ROOT / "code" / "tier2_ml_rra" / "02_run_tier2_rra.R").read_text(encoding="utf-8")
        self.assertNotIn("balanced_accuracy", source)
        self.assertIn('models_expected <- c("logistic_l1", "elastic_net", "logistic_l2", "linear_svm", "random_forest")', source)
        self.assertFalse(self.manifest["rra"]["predictive_performance_weighting"])

    def test_rra_uses_only_five_specified_models(self):
        expected = ["logistic_l1", "elastic_net", "logistic_l2", "linear_svm", "random_forest"]
        self.assertEqual(self.manifest["models_contributing_gene_rankings"], expected)

    def test_empty_sparse_lists_are_not_filled_with_zero_coefficients(self):
        sparse = self.rankings[self.rankings["model"].isin(["logistic_l1", "elastic_net"])]
        self.assertTrue((sparse["absolute_importance"] > runner.ZERO_TOLERANCE).all())

    def test_frozen_output_hashes_recorded_unchanged(self):
        self.assertEqual(
            self.manifest["frozen_allgene_consensus_sha256"],
            "0aefec5fcfa100deafa1525c2d31521ae8c00ca890ef4e1b4145b1b861c2c87e",
        )
        self.assertEqual(
            self.manifest["frozen_corrected_ml_sha256"],
            "6f59e2ed6f652061fcb81a7c5eccbc1211a8dcb9d0bafcc85a1f9fc3edf280f1",
        )

    def test_empty_lists_are_built_from_performance_ledger(self):
        source = (ROOT / "code" / "tier2_ml_rra" / "02_run_tier2_rra.R").read_text(encoding="utf-8")
        self.assertIn("fit_ids <- unique", source)
        self.assertIn("performance$fit_id", source)

    def test_rra_coerces_csv_ranking_flag_to_logical(self):
        source = (ROOT / "code" / "tier2_ml_rra" / "02_run_tier2_rra.R").read_text(encoding="utf-8")
        self.assertIn('ranking_flag == "true"', source)
        self.assertIn('c("true", "false")', source)


if __name__ == "__main__":
    unittest.main()
