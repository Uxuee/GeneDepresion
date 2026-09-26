"""Run corrected ML-priority extraction into a separate, versioned folder."""

from __future__ import annotations

import importlib.util
from pathlib import Path


REPO_DIR = Path(__file__).resolve().parents[1]
SCORER = REPO_DIR / "code" / "extract_stable_ML_prioritized_Tier2_genes.py"
OUT_DIR = (
    REPO_DIR
    / "outputs"
    / "length_gt_100_clean_pipeline"
    / "10_ML_prioritized_Tier2_corrected"
)


def main() -> None:
    spec = importlib.util.spec_from_file_location("corrected_ml_priority_scorer", SCORER)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)

    module.OUT_DIR = OUT_DIR
    module.PLOT_DIR = OUT_DIR / "plots"
    module.GENE_LIST_DIR = OUT_DIR / "gene_lists"
    module.main()


if __name__ == "__main__":
    main()
