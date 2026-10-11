from __future__ import annotations

import hashlib
import importlib.util
import json
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs/length_gt_100_clean_pipeline/20_Tier2_ML_conditional_null"
TIER2 = ROOT / "outputs/length_gt_100_clean_pipeline/07_Tier2_candidate_gene_sets/Tier2_candidate_genes_all_regions.csv"
RUNNER = ROOT / "code/tier2_conditional_null/01_run_permutation_ml.py"
RRA = ROOT / "code/tier2_conditional_null/02_build_permutation_resample_rra.R"
ORCHESTRATOR = ROOT / "code/tier2_conditional_null/04_run_pilot.py"
SUMMARY_SCRIPT = ROOT / "code/tier2_conditional_null/05_summarize_conditional_null.py"
MODELS = {"logistic_l1", "elastic_net", "logistic_l2", "linear_svm", "random_forest"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


class Tier2ConditionalNullTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tier2 = pd.read_csv(TIER2, dtype={"gene_id": str})
        cls.tier2["region"] = cls.tier2["region_or_subregion"].replace({"Cingulate gyrus 25": "CG25"})
        cls.expected_pairs = set(zip(cls.tier2["region"], cls.tier2["gene_id"]))
        cls.perm_dirs = [OUTPUT / "permutations" / f"perm_{index:04d}" for index in range(1, 21)]
        cls.runner_text = RUNNER.read_text(encoding="utf-8")
        cls.rra_text = RRA.read_text(encoding="utf-8")
        cls.orchestrator_text = ORCHESTRATOR.read_text(encoding="utf-8")

    def test_01_exact_671_pairs(self):
        self.assertEqual(len(self.expected_pairs), 671)

    def test_02_no_tier2_recomputation(self):
        self.assertNotIn("f_classif", self.runner_text)
        self.assertNotIn("Tier2_candidate_genes_all_regions.csv\", index=False", self.runner_text)

    def test_03_expression_values_not_permuted(self):
        self.assertIn("rng.permutation(y)", self.runner_text)
        self.assertNotIn("rng.permutation(x", self.runner_text.lower())

    def test_04_sample_ids_unchanged(self):
        for directory in self.perm_dirs:
            mapping = pd.read_csv(directory / "permuted_label_mapping.csv", dtype={"sample_id": str})
            folds = pd.read_csv(directory / "outer_fold_assignments.csv", dtype={"sample_id": str})
            self.assertEqual(set(zip(mapping.region, mapping.sample_id)), set(zip(folds.region, folds.sample_id)))

    def test_05_labels_change(self):
        for directory in self.perm_dirs:
            mapping = pd.read_csv(directory / "permuted_label_mapping.csv")
            self.assertTrue((mapping.observed_label != mapping.permuted_label).any())

    def test_06_class_counts_preserved(self):
        for directory in self.perm_dirs:
            mapping = pd.read_csv(directory / "permuted_label_mapping.csv")
            observed = mapping.groupby(["region", "observed_label"]).size().sort_index()
            permuted = mapping.groupby(["region", "permuted_label"]).size().sort_index()
            pd.testing.assert_series_equal(observed, permuted, check_names=False)

    def test_07_seed_grid_exact_and_deterministic(self):
        manifest = pd.read_csv(OUTPUT / "permutation_manifest.csv")
        self.assertEqual(manifest.permutation_seed.tolist(), list(range(41001, 41021)))
        first = pd.read_csv(self.perm_dirs[0] / "permuted_label_mapping.csv")
        self.assertEqual(first.permutation_seed.nunique(), 1)
        self.assertEqual(int(first.permutation_seed.iloc[0]), 41001)

    def test_08_observed_outputs_never_overwritten(self):
        self.assertNotIn("18_Tier2_ML_RRA_prioritization", self.runner_text.split("EXPECTED_OBSERVED_SHA")[0])
        self.assertIn("checksum_guard()", self.runner_text)

    def test_09_no_observed_hyperparameters_reused(self):
        self.assertNotIn("Tier2_outer_performance.csv\")", self.runner_text.split("def checksum_guard")[0])
        self.assertIn("best = frozen.tune_model", self.runner_text)

    def test_10_preprocessing_training_only(self):
        self.assertIn("frozen.tune_model(x[train_idx], y[train_idx]", self.runner_text)
        self.assertIn("frozen.fit_fixed_pipeline(x[train_idx], y[train_idx]", self.runner_text)

    def test_11_inner_tuning_training_only(self):
        self.assertIn("make_inner_splits(y[train_idx]", self.runner_text)
        self.assertNotIn("tune_model(x, y", self.runner_text)

    def test_12_outer_folds_stratified(self):
        for directory in self.perm_dirs:
            folds = pd.read_csv(directory / "outer_fold_assignments.csv")
            counts = folds.groupby(["region", "repeat", "fold", "role", "permuted_label"]).size()
            self.assertTrue((counts > 0).all())

    def test_13_expected_750_rra_feature_fits(self):
        for directory in self.perm_dirs:
            perf = pd.read_csv(directory / "Tier2_outer_performance.csv")
            self.assertEqual(perf.loc[perf.model.isin(MODELS), "fit_id"].nunique(), 750)

    def test_14_expected_150_lasso_fits(self):
        for directory in self.perm_dirs:
            perf = pd.read_csv(directory / "Tier2_outer_performance.csv")
            self.assertEqual(perf.loc[perf.model.eq("logistic_l1"), "fit_id"].nunique(), 150)

    def test_15_rra_is_permutation_local(self):
        self.assertIn("input <- args[[1]]", self.rra_text)
        self.assertNotIn("18_Tier2_ML_RRA_prioritization", self.rra_text)

    def test_16_empty_sparse_lists_remain_empty(self):
        self.assertIn("as.character(part$gene_id)", self.rra_text)
        self.assertNotIn("fill", self.rra_text.lower())

    def test_17_dummy_never_gene_evidence(self):
        for directory in self.perm_dirs:
            rankings = pd.read_csv(directory / "Tier2_fold_gene_rankings.csv", usecols=["model"])
            self.assertNotIn("dummy_most_frequent", set(rankings.model))

    def test_18_observed_tier2_rra_hash_unchanged(self):
        path = ROOT / "outputs/length_gt_100_clean_pipeline/18_Tier2_ML_RRA_prioritization/Tier2_consensus_rra_rankings.csv"
        self.assertEqual(sha256(path), "2e1ed5fc91aae752f57e46e61022a9077d7f1ab285cfc05a6973c3c6df23d358")

    def test_19_benchmark_hash_unchanged(self):
        path = ROOT / "outputs/length_gt_100_clean_pipeline/19_Tier2_ML_benchmark/checksums.sha256"
        self.assertEqual(sha256(path), "4e17dc8234317efc9bce92aa03fb7ae63ff72a128392d202f2f31a105c09cf10")

    def test_20_allgene_hash_unchanged(self):
        path = ROOT / "outputs/length_gt_100_clean_pipeline/17_ML_RRA_consensus/consensus_rra_rankings.csv"
        self.assertEqual(sha256(path), "0aefec5fcfa100deafa1525c2d31521ae8c00ca890ef4e1b4145b1b861c2c87e")

    def test_21_completed_checkpoints_resume_without_refit(self):
        self.assertIn("complete checksums verified; skipped", self.orchestrator_text)
        for directory in self.perm_dirs:
            completion = json.loads((directory / "permutation_completion.json").read_text(encoding="utf-8"))
            self.assertEqual(completion["status"], "complete")

    def test_22_summary_grid_reconstructs_from_raw(self):
        values = pd.read_csv(OUTPUT / "summary/null_metric_values.csv")
        grouped = values.groupby(["family", "region", "metric", "source"]).size()
        for key, count in grouped.items():
            self.assertEqual(count, 1 if key[-1] == "observed" else 20)

    def test_23_plus_one_empirical_tail(self):
        summary = pd.read_csv(OUTPUT / "summary/observed_vs_null_summary.csv")
        expected = (1 + summary.n_null_as_or_more_extreme) / (1 + summary.n_permutations)
        np.testing.assert_allclose(summary.empirical_p_pilot, expected)

    def test_24_absent_lasso_evidence_is_deterministic(self):
        for directory in self.perm_dirs:
            lasso = pd.read_csv(directory / "LASSO_gene_rankings.csv")
            absent = lasso.loc[lasso.selection_count.eq(0)]
            for _, region_absent in absent.groupby("region"):
                self.assertEqual(region_absent.lasso_correlation_rank.nunique(), 1)

    def test_25_no_weighted_priority_score(self):
        text = SUMMARY_SCRIPT.read_text(encoding="utf-8").lower()
        self.assertNotIn("weighted_score", text)
        self.assertNotIn("composite_score", text)


if __name__ == "__main__":
    unittest.main()
