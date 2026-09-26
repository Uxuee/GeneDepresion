import json
import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = json.loads((ROOT / "provenance/canonical_pipeline.json").read_text(encoding="utf-8"))


class Phase3AGuardrails(unittest.TestCase):
    def test_registry_stage_range(self):
        ids = {s["stage_id"] for s in REGISTRY["stages"]}
        self.assertTrue({f"{i:02d}" for i in range(1, 20)} <= ids)

    def test_corrected_counts_and_denominators(self):
        x = REGISTRY["canonical_ml_result"]
        self.assertEqual((x["total"], x["high"], x["moderate"], x["low"]), (671, 30, 119, 522))
        self.assertEqual((x["v1_denominator"], x["v2_denominator"]), (125, 150))

    def test_superseded_ranking_not_canonical(self):
        active = [s for s in REGISTRY["stages"] if s["status"].startswith("canonical")]
        self.assertTrue(all("10_ML_prioritized_Tier2_extended" not in str(s["canonical_outputs"]) for s in active))

    def test_missing_generator_stage_is_excluded(self):
        stage = next(s for s in REGISTRY["stages"] if s["stage_id"] == "17")
        self.assertEqual(stage["status"], "excluded_from_publication")
        self.assertIsNone(stage["canonical_script"])

    def test_ahba_post_selection_sets_historical(self):
        x = json.loads((ROOT / "provenance/ahba_module_score_ml_status.json").read_text(encoding="utf-8"))
        self.assertEqual(x["active_feature_set"], "all_32_published_AHBA_modules")
        self.assertEqual(set(x["historical_post_selection_sensitivity_sets"]), {"FDR005", "FDR010", "M26_only"})

    def test_legacy_runner_has_opt_in_gate(self):
        source = (ROOT / "code/length_gt_100_clean_pipeline/pipeline_lib.py").read_text(encoding="utf-8")
        self.assertIn("--allow-legacy", source)
        self.assertIn("legacy_selected", source)

    def test_active_corrected_scripts_do_not_reference_legacy_ranking(self):
        for rel in ("code/run_canonical_pipeline.py", "code/regenerate_corrected_ml_downstream.py", "code/summarize_corrected_ml_downstream.py", "code/canonical_downstream_paths.py"):
            self.assertNotIn("10_ML_prioritized_Tier2_extended", (ROOT / rel).read_text(encoding="utf-8"))

    def test_current_readme_has_no_stale_counts(self):
        current = (ROOT / "README.md").read_text(encoding="utf-8").split("<!-- HISTORICAL_CONTENT_START -->", 1)[0].lower()
        for stale in ("35 high", "130 moderate", "165 high+moderate", "v2 denominator = 180"):
            self.assertNotIn(stale, current)

    def test_current_navigation_has_no_legacy_corrected_claims(self):
        text = ((ROOT / "README.md").read_text(encoding="utf-8").split("<!-- HISTORICAL_CONTENT_START -->", 1)[0]
                + (ROOT / "docs/CANONICAL_PIPELINE_MAP.md").read_text(encoding="utf-8")).lower()
        self.assertNotIn("m26 q<0.10", text)
        self.assertNotIn("m26 fdr<0.10", text)
        self.assertNotIn("medarkred q<0.10", text)
        self.assertNotIn("externally validated", text)

    def test_every_registry_stage_has_provenance_fields(self):
        required = {"stage_id", "stage_name", "status", "canonical_script", "canonical_inputs", "canonical_outputs", "input_sha256", "output_sha256", "scientific_dependency", "git_milestone", "random_seed", "environment_reference", "rerun_required", "historical_replacements", "notes"}
        for stage in REGISTRY["stages"]:
            self.assertEqual(required - set(stage), set())

    def test_recovered_tier4_blobs_exist(self):
        x = json.loads((ROOT / "provenance/tier4_ranked_signature_recovery_manifest.json").read_text(encoding="utf-8"))
        for row in x["files"]:
            self.assertEqual(subprocess.run(["git", "cat-file", "-e", row["blob"]], cwd=ROOT).returncode, 0)


if __name__ == "__main__":
    unittest.main()
