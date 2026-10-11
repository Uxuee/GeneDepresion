from __future__ import annotations

import importlib.util
import json
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
VALIDATOR_PATH = ROOT / "code/tier2_benchmark/validate_benchmark.py"
SPEC = importlib.util.spec_from_file_location("benchmark_validator", VALIDATOR_PATH)
validator = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = validator
SPEC.loader.exec_module(validator)
RUNNER = ROOT / "code/tier2_benchmark/02_run_stability_selection.py"
OUTPUT = ROOT / "outputs/length_gt_100_clean_pipeline/19_Tier2_ML_benchmark"


class Tier2PrioritizationBenchmarkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.checks = validator.validation_checks()
        cls.tables = validator.load()

    def check(self, name):
        self.assertTrue(self.checks[name], name)

    def test_exact_671_pairs(self): self.check("exact_671_unique_tier2_pairs")
    def test_same_six_regions(self): self.check("exact_six_regions")
    def test_no_tier4_lasso(self): self.check("lasso_uses_exact_tier2_pairs")
    def test_no_tier4_stability(self): self.check("stability_uses_exact_tier2_pairs")
    def test_no_external_pairs(self): self.check("no_tier4_or_external_pairs")
    def test_primary_rra_frozen(self): self.check("frozen_primary_hash_unchanged")
    def test_allgene_rra_frozen(self): self.check("frozen_allgene_hash_unchanged")
    def test_outer_coverage(self): self.check("lasso_25_outer_fits_per_region")
    def test_lasso_count_reproduction(self): self.check("lasso_selection_counts_reproduce")
    def test_lasso_denominator(self): self.check("lasso_denominator_exactly_25")
    def test_stability_subsample_coverage(self): self.check("stability_100_subsamples_per_region")
    def test_stability_event_uniqueness(self): self.check("stability_one_event_per_gene_subsample")
    def test_diagnosis_sample_mapping(self): self.check("stability_preserves_sample_mapping")
    def test_stratified_half_samples(self): self.check("stability_balanced_half_samples")
    def test_stability_count_reproduction(self): self.check("stability_selection_counts_reproduce")
    def test_stability_denominator(self): self.check("stability_denominator_exactly_100")
    def test_seed_grid(self): self.check("stability_seed_grid_exact")
    def test_never_selected_shared_rank(self): self.check("never_selected_shared_worst_rank")
    def test_topset_definitions(self): self.check("topset_definitions_exact")
    def test_saved_summary_grids(self):
        self.check("rank_correlation_grid_complete")
        self.check("coverage_summary_reproduces")
        self.check("resampling_summary_complete")
    def test_no_metric_weighted_score(self):
        self.check("no_metric_or_weighted_priority_score")
        self.check("manifest_rejects_composite_and_metric_weighting")
    def test_output_checksums(self): self.check("all_output_checksums_match")

    def test_training_only_preprocessing_and_tuning_source(self):
        source = RUNNER.read_text(encoding="utf-8")
        self.assertIn("imputer.fit_transform(x_train)", source)
        self.assertIn("StandardScaler().fit(imputed)", source)
        self.assertIn("x_sub[fit_idx]", source)
        self.assertNotIn("priority_score", source)

    def test_deterministic_half_sample_function(self):
        spec = importlib.util.spec_from_file_location("stability_runner", RUNNER)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        y = np.array([0] * 20 + [1] * 20)
        first = module.stratified_half_sample(y, 9182)
        second = module.stratified_half_sample(y, 9182)
        np.testing.assert_array_equal(first, second)
        self.assertEqual(int((y[first] == 0).sum()), 10)
        self.assertEqual(int((y[first] == 1).sum()), 10)


if __name__ == "__main__":
    unittest.main()
