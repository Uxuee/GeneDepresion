"""Canonical paths shared by corrected ML-dependent downstream regeneration."""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CLEAN_ROOT = PROJECT_ROOT / "outputs" / "length_gt_100_clean_pipeline"
CANONICAL_ML_DIR = CLEAN_ROOT / "10_ML_prioritized_Tier2_corrected"
CANONICAL_RANKING = CANONICAL_ML_DIR / "ML_prioritized_Tier2_genes_all_ranked.csv"
PHASE2_ROOT = CLEAN_ROOT / "16_corrected_ML_downstream"
CANONICAL_RANKING_SHA256 = "6f59e2ed6f652061fcb81a7c5eccbc1211a8dcb9d0bafcc85a1f9fc3edf280f1"
