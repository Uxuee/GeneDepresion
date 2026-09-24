"""Lightweight integrity tests for the canonical corrected ML-priority result."""

from __future__ import annotations

import json
from pathlib import Path
import unittest

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = ROOT / "outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_corrected"


class CanonicalMLPriorityResultTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.ranking = pd.read_csv(RESULT_DIR / "ML_prioritized_Tier2_genes_all_ranked.csv")
        cls.manifest = json.loads((RESULT_DIR / "canonical_manifest.json").read_text(encoding="utf-8"))

    def test_canonical_counts(self) -> None:
        self.assertEqual(len(self.ranking), 671)
        counts = self.ranking["ML_priority_class"].value_counts().to_dict()
        self.assertEqual(counts.get("High-confidence ML-prioritized Tier 2"), 30)
        self.assertEqual(counts.get("Moderate-confidence ML-prioritized Tier 2"), 119)
        self.assertEqual(counts.get("Exploratory ML-prioritized Tier 2"), 522)

    def test_event_denominators(self) -> None:
        self.assertEqual(set(self.ranking["n_events_v1"].dropna().astype(int)), {125})
        self.assertEqual(set(self.ranking["n_events_v2"].dropna().astype(int)), {150})

    def test_manifest_matches_result(self) -> None:
        population = self.manifest["candidate_population"]
        self.assertEqual(population["total"], len(self.ranking))
        self.assertEqual(population["high_confidence"], 30)
        self.assertEqual(population["moderate_confidence"], 119)
        self.assertEqual(population["exploratory_low_confidence"], 522)

    def test_high_confidence_has_observed_directional_evidence(self) -> None:
        high = self.ranking[
            self.ranking["ML_priority_class"].eq("High-confidence ML-prioritized Tier 2")
        ]
        self.assertTrue(high["sign_consistency"].notna().all())
        self.assertTrue((high["sign_consistency"] >= 0.70).all())

    def test_representative_arithmetic(self) -> None:
        reconstruction = pd.read_csv(RESULT_DIR / "real_gene_score_reconstruction.csv")
        self.assertIn("LOC107985795", set(reconstruction["gene_symbol"]))
        self.assertIn("PRELID2", set(reconstruction["gene_symbol"]))
        self.assertLessEqual(reconstruction["absolute_arithmetic_difference"].max(), 1e-12)


if __name__ == "__main__":
    unittest.main()
