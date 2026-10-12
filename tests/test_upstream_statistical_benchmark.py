import importlib.util
import sys
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
MODULE = REPO / "code" / "statistical_benchmark" / "validate_statistical_benchmark.py"
sys.path.insert(0, str(MODULE.parent))
spec = importlib.util.spec_from_file_location("validate_statistical_benchmark", MODULE)
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)


class UpstreamStatisticalBenchmarkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.checks = validator.collect_checks()


def _add_test(name):
    def test(self):
        self.assertTrue(self.checks[name], name)
    setattr(UpstreamStatisticalBenchmarkTests, f"test_{name}", test)


for _name in [
    "01_exact_six_regions", "02_tier2_671_exact_counts", "03_frozen_tier2_hash_unchanged",
    "04_frozen_LM_hash_unchanged", "05_DESeq2_raw_integer_counts", "06_limma_DGEList_from_raw_counts",
    "07_no_GeTMM_for_DESeq2", "08_no_GeTMM_for_voom", "09_same_matched_samples",
    "10_same_eligible_genes", "11_contrast_MDD_minus_CTRL", "12_positive_effect_is_higher_MDD",
    "13_age_sex_RIN_in_design", "14_no_Tier2_in_fitting", "15_no_ML_rank_in_fitting",
    "16_BH_denominator_documented", "17_DESeq2_native_and_same_BH_separate",
    "18_rank_definitions_reproduce", "19_gene_ids_one_to_one", "20_frozen_checksum_ledgers_unchanged",
    "21_summaries_reproduce_from_saved_results", "22_manifest_complete",
]:
    _add_test(_name)


if __name__ == "__main__":
    unittest.main()
