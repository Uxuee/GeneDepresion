"""Integrity tests for the noncanonical Phase-3B WGCNA validation candidate.

These tests validate provenance and accounting; they deliberately do not encode
desired biological significance outcomes.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path("D:/genes/GeneDepresion")
OUT = ROOT / "outputs/length_gt_100_clean_pipeline/17_WGCNA_reconstruction_validation"
REGIONS = {
    "Ventral Striatum": "Ventral_Striatum",
    "Insula": "Insula",
    "Hippocampal Formation": "Hippocampal_Formation",
    "Cingulate gyrus 25": "Cingulate_gyrus_25",
    "dlPFC BA8/9": "dlPFC_BA8_9",
    "OFC BA11": "OFC_BA11",
}


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


class Phase3BWGCNAValidationTests(unittest.TestCase):
    def test_frozen_input_hashes(self):
        manifest = rows(OUT / "inputs_manifest/wgcna_reconstruction_input_manifest.csv")
        self.assertEqual(len(manifest), 6)
        for item in manifest:
            path = Path(item["path"])
            self.assertTrue(path.is_file())
            self.assertEqual(sha256(path), item["sha256"])

    def test_construction_script_does_not_use_diagnosis_for_network(self):
        text = (ROOT / "code/phase3b_reconstruct_wgcna.R").read_text(encoding="utf-8")
        self.assertIn('construction_used_diagnosis = FALSE', text)
        residual_line = next(line for line in text.splitlines() if "covars <- usable_covariates" in line)
        self.assertNotIn("Diagnosis", residual_line)
        self.assertNotIn("tier2", text.lower())
        self.assertNotIn("ML_prioritized", text)

    def test_selected_parameters_reproduced_without_guessing(self):
        historical = {r["region"]: r for r in rows(OUT / "inputs_manifest/historical_selected_parameter_rows.csv")}
        for region, folder in REGIONS.items():
            selected = rows(OUT / f"reconstruction/{folder}/selected_solution_parameters.csv")[0]
            old = historical[region]
            for field in ("power", "deepSplit", "minModuleSize", "mergeCutHeight", "corType", "networkType"):
                self.assertEqual(selected[field], old[field], f"{region}: {field}")

    def test_assignment_integrity(self):
        for region, folder in REGIONS.items():
            assignment = rows(OUT / f"reconstruction/{folder}/reconstructed_gene_to_module_assignments.csv")
            genes = [r["gene_id"] for r in assignment]
            self.assertEqual(len(genes), len(set(genes)), region)
            self.assertTrue(all(r["region"] == region for r in assignment))
            self.assertTrue(all(r["reconstructed_module"].startswith("ME") for r in assignment))

    def test_label_invariant_metrics_and_predeclared_thresholds(self):
        threshold_path = OUT / "module_comparison/predeclared_equivalence_thresholds.json"
        summary_path = OUT / "module_comparison/wgcna_reconstruction_summary.csv"
        self.assertTrue(threshold_path.is_file())
        thresholds = json.loads(threshold_path.read_text(encoding="utf-8"))
        self.assertEqual(thresholds["effectively_equivalent_ari_min"], 0.9)
        summary = rows(summary_path)
        self.assertEqual(len(summary), 6)
        for item in summary:
            for field in ("adjusted_rand_index", "median_best_jaccard", "gene_assignment_agreement_from_best_mapping", "median_abs_eigengene_correlation"):
                value = float(item[field])
                self.assertGreaterEqual(value, -1 if field == "adjusted_rand_index" else 0)
                self.assertLessEqual(value, 1)

    def test_corrected_ml_counts(self):
        ranking = rows(ROOT / "outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_corrected/ML_prioritized_Tier2_genes_all_ranked.csv")
        labels = [r["ML_priority_class"] for r in ranking]
        high = sum("High-confidence" in x for x in labels)
        moderate = sum("Moderate-confidence" in x for x in labels)
        low = len(labels) - high - moderate
        self.assertEqual((len(ranking), high, moderate, low, high + moderate), (671, 30, 119, 522, 149))

    def test_enrichment_background_integrity(self):
        enrichment = rows(OUT / "ml_gene_enrichment_validation/corrected_ML_gene_enrichment_old_and_reconstructed.csv")
        self.assertTrue(enrichment)
        for item in enrichment:
            n = int(item["network_universe_size"])
            m = int(item["module_gene_count"])
            s = int(item["gene_set_count_in_background"])
            overlap = int(item["overlap_count"])
            self.assertGreaterEqual(n, m)
            self.assertGreaterEqual(n, s)
            self.assertGreaterEqual(overlap, 0)
            self.assertLessEqual(overlap, min(m, s))

    def test_annotation_mapping_accounting(self):
        coverage = rows(OUT / "functional_annotation/annotation_mapping_coverage_by_module.csv")
        self.assertEqual(len(coverage), 10)
        for item in coverage:
            total = int(item["module_genes_total"])
            mapped = int(item["module_genes_mapped"])
            unmapped = int(item["module_genes_unmapped"])
            self.assertEqual(total, mapped + unmapped)

    def test_manifests_exist_and_are_noncanonical(self):
        for _, folder in REGIONS.items():
            manifest = json.loads((OUT / f"reconstruction/{folder}/reconstruction_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["status"], "validation_candidate")
            self.assertFalse(manifest["construction_used_diagnosis"])
        combined = json.loads((OUT / "manifests/phase3b_comparison_annotation_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(combined["status"], "validation_candidate")


if __name__ == "__main__":
    unittest.main()
