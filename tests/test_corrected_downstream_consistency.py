"""Invariants for corrected ML-dependent downstream outputs."""

from pathlib import Path
import json
import unittest

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CLEAN = ROOT / "outputs/length_gt_100_clean_pipeline"
P2 = CLEAN / "16_corrected_ML_downstream"


class CorrectedDownstreamTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rank = pd.read_csv(CLEAN / "10_ML_prioritized_Tier2_corrected/ML_prioritized_Tier2_genes_all_ranked.csv", dtype={"gene_id": str})
        cls.manifest = json.loads((P2 / "phase2_manifest.json").read_text(encoding="utf-8"))

    def test_canonical_counts(self):
        self.assertEqual(self.manifest["counts"], {"total": 671, "high": 30, "moderate": 119, "high_plus_moderate": 149})

    def test_wgcna_modules_are_from_preserved_assignments(self):
        sources = {
            "broad_regions": CLEAN / "12_MDD_WGCNA_modules/best_solution_gene_to_module_assignments.csv",
            "frontal_subregions": CLEAN / "mdd_wgcna_modules/frontal_subregion_wgcna/frontal_subregion_best_solution_gene_to_module_assignments.csv",
        }
        for unit, source in sources.items():
            expected = set(map(tuple, pd.read_csv(source, usecols=["region", "module"]).drop_duplicates().to_numpy()))
            actual = set(map(tuple, pd.read_csv(P2 / f"WGCNA_ML_annotations_corrected/{unit}_ML_module_enrichment_corrected.csv", usecols=["region", "module"]).drop_duplicates().to_numpy()))
            self.assertEqual(actual, expected)

    def test_all_flagged_nodes_exist_in_corrected_ranking(self):
        ids = set(self.rank.gene_id.astype(str))
        for path in (P2 / "WGCNA_ML_annotations_corrected").glob("*_network_node_annotations_corrected.csv"):
            data = pd.read_csv(path, dtype={"gene_id": str})
            flagged = data[data.is_corrected_ML_high_plus_moderate].gene_id.astype(str)
            self.assertTrue(set(flagged).issubset(ids))

    def test_backgrounds_are_nonempty_and_fixed(self):
        ahba = pd.read_csv(P2 / "AHBA_corrected/published_AHBA_ML_enrichment_corrected.csv")
        self.assertEqual(ahba.background_genes.nunique(), 1)
        self.assertGreater(int(ahba.background_genes.iloc[0]), 0)
        for path in (P2 / "WGCNA_ML_annotations_corrected").glob("*_ML_module_enrichment_corrected.csv"):
            data = pd.read_csv(path)
            self.assertTrue((data.module_gene_count > 0).all())

    def test_before_values_reproduce_retained_wgcna_outputs(self):
        comparisons = {
            "broad_regions": CLEAN / "12_MDD_WGCNA_modules/module_enrichment_for_Tier2_and_ML_genes.csv",
            "frontal_subregions": CLEAN / "mdd_wgcna_modules/frontal_subregion_wgcna/frontal_subregion_Tier2_ML_enrichment.csv",
        }
        for unit, retained_path in comparisons.items():
            before = pd.read_csv(P2 / f"WGCNA_ML_annotations_corrected/{unit}_ML_module_enrichment_before_after.csv")
            retained = pd.read_csv(retained_path)
            retained = retained[retained["gene_set"].isin([
                "high_confidence_ML_genes",
                "high_plus_moderate_ML_genes",
                "same_region_high_confidence_ML_genes",
                "same_region_high_plus_moderate_ML_genes",
            ])]
            keys = ["region", "module", "gene_set"]
            retained = retained[keys + ["overlap_count", "p_value", "fdr_BH"]].rename(columns={
                "overlap_count": "overlap_count_retained",
                "p_value": "p_value_retained",
                "fdr_BH": "fdr_BH_retained",
            })
            check = before.merge(retained, on=keys, how="inner")
            self.assertGreater(len(check), 0)
            self.assertTrue((check.overlap_count_old == check.overlap_count_retained).all())
            for col in ("p_value", "fdr_BH"):
                delta = (check[f"{col}_old"] - check[f"{col}_retained"]).abs().max()
                self.assertLess(float(delta), 1e-10)

    def test_frontal_networks_exclude_pooled_frontal_cortex(self):
        inventory = pd.read_csv(P2 / "WGCNA_ML_annotations_corrected/network_figures_frontal_subregions/frontal_subregion_network_inventory.csv")
        self.assertNotIn("Frontal Cortex", set(inventory.region))
        self.assertEqual(set(inventory.region), {"dlPFC BA8/9", "OFC BA11", "Cingulate gyrus 25"})


if __name__ == "__main__":
    unittest.main()
