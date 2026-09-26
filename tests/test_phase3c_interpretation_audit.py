import csv
import hashlib
import json
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
P3B = ROOT / "outputs/length_gt_100_clean_pipeline/17_WGCNA_reconstruction_validation"
P3C = ROOT / "outputs/length_gt_100_clean_pipeline/18_WGCNA_enrichment_interpretation_audit"
FUNC = P3C / "functional_annotation"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


class Phase3CInterpretationAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = pd.read_csv(P3B / "functional_annotation/functional_enrichment_GO_Reactome_Hallmark.csv")
        cls.evidence = pd.read_csv(P3B / "functional_annotation/prioritized_module_evidence_table.csv")
        cls.all_terms = pd.read_csv(FUNC / "all_prioritized_module_enrichment_terms.csv", low_memory=False)
        cls.fdr05 = pd.read_csv(FUNC / "fdr05_terms_with_overlap_genes.csv")
        cls.fdr10 = pd.read_csv(FUNC / "fdr10_suggestive_terms.csv")
        cls.categories = pd.read_csv(FUNC / "publication_interpretation_categories.csv")
        cls.background = pd.read_csv(FUNC / "background_universe_audit.csv")
        cls.drivers = pd.read_csv(FUNC / "gene_drivers_by_term.csv", dtype={"gene_id": str})
        cls.ml = pd.read_csv(
            ROOT / "outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_corrected/ML_prioritized_Tier2_genes_all_ranked.csv",
            dtype={"gene_id": str},
        )
        cls.manifest = json.loads((P3C / "manifests/phase3c_manifest.json").read_text(encoding="utf-8"))

    def test_exact_prespecified_modules(self):
        expected = set(map(tuple, self.evidence[["region", "old_module"]].to_records(index=False)))
        audited = set(map(tuple, self.categories[["region", "module"]].to_records(index=False)))
        self.assertEqual(len(expected), 10)
        self.assertEqual(audited, expected)

    def test_no_new_database_or_enrichment_run(self):
        self.assertEqual(set(self.all_terms["database"]), {"GO", "Reactome", "MSigDB Hallmark"})
        self.assertTrue(self.manifest["no_new_enrichment"])
        self.assertTrue(self.manifest["no_WGCNA_rerun"])

    def test_source_q_values_and_overlap_genes_match(self):
        keys = ["region", "old_module", "database", "ontology", "term_id", "term"]
        merged = self.all_terms.merge(
            self.source[keys + ["BH_q", "raw_p", "overlap_genes"]],
            on=keys, suffixes=("_audit", "_source"), validate="one_to_one"
        )
        self.assertEqual(len(merged), len(self.source))
        self.assertTrue(np.allclose(merged["BH_q_audit"], merged["BH_q_source"], rtol=0, atol=1e-15))
        self.assertTrue(np.allclose(merged["raw_p_audit"], merged["raw_p_source"], rtol=0, atol=1e-15))
        self.assertTrue((merged["overlap_genes_audit"].astype(str) == merged["overlap_genes_source"].astype(str)).all())

    def test_threshold_classification(self):
        self.assertEqual(len(self.fdr05), 200)
        self.assertEqual(len(self.fdr10), 51)
        self.assertTrue((self.fdr05["BH_q"] < 0.05).all())
        self.assertTrue(((self.fdr10["BH_q"] >= 0.05) & (self.fdr10["BH_q"] < 0.10)).all())

    def test_diagnosis_q_values_not_relabeled(self):
        merged = self.categories.merge(
            self.evidence[["region", "old_module", "reconstructed_q"]],
            left_on=["region", "module"], right_on=["region", "old_module"], validate="one_to_one"
        )
        self.assertTrue((merged["Diagnosis_BH_q"] == merged["reconstructed_q"]).all())
        self.assertTrue((merged["Diagnosis_FDR_status"] == "NOT_BH_FDR_SIGNIFICANT").all())
        self.assertTrue((merged["Diagnosis_BH_q"] >= 0.10).all())

    def test_key_ml_statuses_unchanged(self):
        vs = self.categories[(self.categories["region"] == "Ventral Striatum") & (self.categories["module"] == "MEdarkred")].iloc[0]
        self.assertAlmostEqual(vs["ML_high_plus_moderate_p"], 0.0021701362897413564)
        self.assertAlmostEqual(vs["ML_high_plus_moderate_q"], 0.10850681448706782)
        self.assertGreaterEqual(vs["ML_high_plus_moderate_q"], 0.10)
        ofc = self.categories[(self.categories["region"] == "OFC BA11") & (self.categories["module"] == "MEtan")].iloc[0]
        self.assertEqual(ofc["ML_high_overlap"], 2)
        self.assertAlmostEqual(ofc["ML_high_q"], 0.032051449164565426)
        self.assertEqual(ofc["ML_high_plus_moderate_overlap"], 3)
        self.assertAlmostEqual(ofc["ML_high_plus_moderate_q"], 0.06860675392368916)

    def test_no_legacy_ml_list(self):
        script = (ROOT / "code/phase3c_audit_wgcna_interpretation.py").read_text(encoding="utf-8")
        manifest_text = json.dumps(self.manifest)
        self.assertNotIn("10_ML_prioritized_Tier2_extended", script)
        self.assertNotIn("10_ML_prioritized_Tier2_extended", manifest_text)
        self.assertNotIn("35/130", script)
        self.assertNotIn("35/130", manifest_text)

    def test_backgrounds_match_phase3b_region_networks(self):
        self.assertTrue(self.background["region_network_background_confirmed"].all())
        self.assertTrue((self.background["universe_mapped"] <= self.background["universe_genes"]).all())

    def test_driver_ml_status_matches_global_and_same_region_corrected_ranking(self):
        def label(value):
            value = str(value)
            if "High-confidence" in value:
                return "high"
            if "Moderate-confidence" in value:
                return "moderate"
            if "Low-confidence" in value:
                return "low"
            return "absent"

        priority = {"high": 0, "moderate": 1, "low": 2, "absent": 3}
        ml = self.ml.copy()
        ml["expected"] = ml["ML_priority_class"].map(label)
        global_expected = (
            ml.sort_values("expected", key=lambda x: x.map(priority))
            .drop_duplicates("gene_id")
            .set_index("gene_id")["expected"]
        )
        same_expected = ml.set_index(["region_or_subregion", "gene_id"])["expected"]
        for row in self.drivers[["region", "gene_id", "ml_priority", "same_region_ml_priority"]].drop_duplicates().itertuples(index=False):
            self.assertEqual(row.ml_priority, global_expected.get(row.gene_id, "absent"))
            self.assertEqual(row.same_region_ml_priority, same_expected.get((row.region, row.gene_id), "absent"))

    def test_output_manifest_hashes(self):
        inventory = pd.read_csv(P3C / "manifests/phase3c_output_file_manifest.csv")
        self.assertEqual(len(inventory), 16)
        for row in inventory.itertuples(index=False):
            path = Path(row.path)
            self.assertTrue(path.is_file(), row.relative_path)
            self.assertEqual(sha256(path), row.sha256, row.relative_path)


if __name__ == "__main__":
    unittest.main()
