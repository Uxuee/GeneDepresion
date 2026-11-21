from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
MODULE = REPO / "code" / "validate_final_evidence_synthesis.py"
spec = importlib.util.spec_from_file_location("final_synthesis_validator", MODULE)
validator = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = validator
spec.loader.exec_module(validator)


class FinalEvidenceSynthesisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.checks = validator.run_checks(check_determinism=False)


def make_test(name):
    def test(self):
        self.assertTrue(self.checks[name], name)
    return test


for check_name in validator.run_checks(check_determinism=False):
    setattr(FinalEvidenceSynthesisTests, f"test_{check_name}", make_test(check_name))


if __name__ == "__main__":
    unittest.main()
