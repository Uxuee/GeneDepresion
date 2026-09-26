import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class ReleaseReadinessTests(unittest.TestCase):
    def test_readme_has_canonical_entrypoints_without_legacy_result_link(self):
        text = (ROOT / "README.md").read_text(encoding="utf-8")
        for required in (
            "docs/CANONICAL_PIPELINE_MAP.md",
            "provenance/canonical_pipeline.json",
            "config/project_paths.example.yaml",
            "data/README.md",
            "environment/README.md",
            "code/run_canonical_pipeline.py",
            "10_ML_prioritized_Tier2_corrected",
        ):
            self.assertIn(required, text)
        self.assertNotIn("](outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_extended", text)

    def test_status_is_status_only(self):
        result = subprocess.run(
            [sys.executable, str(ROOT / "code/run_canonical_pipeline.py"), "--status"],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Corrected ML-priority scoring", result.stdout)
        self.assertNotIn("NO-OP", result.stdout)
        self.assertNotIn("REFUSE", result.stdout)

    def test_fresh_clone_validator_mode_is_explicit(self):
        result = subprocess.run(
            [
                sys.executable,
                str(ROOT / "audit/validate_canonical_pipeline.py"),
                "--allow-missing-local-artifacts",
            ],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        payload = json.loads(result.stdout)
        self.assertEqual(payload["errors"], [])
        self.assertEqual(payload["warnings"], [])

    def test_external_manifest_has_acquisition_and_destination_fields(self):
        with (ROOT / "provenance/external_data_manifest.csv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertGreaterEqual(len(rows), 10)
        required_fields = {
            "resource_name", "accession_or_version", "expected_filename",
            "local_destination", "check_kind", "required", "source",
            "download_method", "license_or_terms",
        }
        self.assertTrue(required_fields.issubset(rows[0]))
        for row in rows:
            self.assertTrue(row["source"])
            self.assertTrue(row["local_destination"])
            self.assertTrue(row["license_or_terms"])

    def test_missing_external_data_is_actionable(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            config = base / "paths.yaml"
            config.write_text(
                f"raw_data_root: {base / 'missing_raw'}\n"
                f"external_resource_root: {base / 'missing_external'}\n",
                encoding="utf-8",
            )
            result = subprocess.run(
                [sys.executable, str(ROOT / "audit/check_external_data.py"), "--config", str(config)],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
            )
        self.assertEqual(result.returncode, 2)
        self.assertIn("MISSING GSE102556 raw counts", result.stdout)
        self.assertIn("data/README.md", result.stdout)
        self.assertNotIn("Traceback", result.stderr + result.stdout)

    def test_phase3b_environment_is_real_and_documented(self):
        lock = json.loads((ROOT / "environment/phase3b/renv.lock").read_text(encoding="utf-8"))
        self.assertEqual(lock["R"]["Version"], "4.6.1")
        self.assertTrue((ROOT / "environment/phase3b/sessionInfo.txt").is_file())
        self.assertTrue((ROOT / "environment/phase3b/package_versions.csv").is_file())
        documentation = (ROOT / "environment/README.md").read_text(encoding="utf-8")
        self.assertIn("reconstructed-current", documentation)
        self.assertIn("environment/phase3b/renv.lock", documentation)

    def test_registry_excludes_unrecovered_publication_branches(self):
        registry = json.loads((ROOT / "provenance/canonical_pipeline.json").read_text(encoding="utf-8"))
        stages = {stage["stage_id"]: stage for stage in registry["stages"]}
        self.assertEqual(stages["08b"]["status"], "excluded_from_publication")
        self.assertEqual(stages["17"]["status"], "excluded_from_publication")
        active = [stage for stage in registry["stages"] if stage["status"] in {"canonical", "canonical_external_input"}]
        self.assertTrue(all("10_ML_prioritized_Tier2_extended" not in json.dumps(stage["canonical_outputs"]) for stage in active))


if __name__ == "__main__":
    unittest.main()
