"""Configuration for the Length > 100 clean Tier2-only MDD pipeline.

The clean pipeline reuses completed outputs from outputs/length_gt_100_elisa_style
and writes organized, report-friendly outputs to outputs/length_gt_100_clean_pipeline.
Expensive analyses are not rerun by default.
"""

from pathlib import Path


PROJECT_ROOT = Path(r"D:\genes\GeneDepresion")
OLD_ANALYSIS_DIR = PROJECT_ROOT / "outputs" / "length_gt_100_elisa_style"
CLEAN_OUTPUT_DIR = PROJECT_ROOT / "outputs" / "length_gt_100_clean_pipeline"
RAW_DATA_DIR = PROJECT_ROOT / "data"

NORMALIZED_EXPRESSION_DIR = OLD_ANALYSIS_DIR / "normalized_expression_length_gt_100_by_region"
FILTERED_RAW_COUNTS_DIR = OLD_ANALYSIS_DIR / "filtered_raw_counts_length_gt_100_by_region"
GENE_FILTER_SUMMARY = OLD_ANALYSIS_DIR / "gene_filter_length_gt_100_summary.csv"

LMM_ALL_GENES = OLD_ANALYSIS_DIR / "LMM_all_genes_length_gt_100.csv"
LMM_SUMMARY_BY_REGION = OLD_ANALYSIS_DIR / "summary_by_region_length_gt_100.csv"

STANDARDIZED_EFFECT_TABLE = (
    OLD_ANALYSIS_DIR
    / "standardized_effect_size_4groups"
    / "LMM_all_genes_with_standardized_effect_size_4groups.csv"
)
STANDARDIZED_EFFECT_SUMMARY = (
    OLD_ANALYSIS_DIR
    / "standardized_effect_size_4groups"
    / "standardized_effect_size_4groups_summary.csv"
)
STANDARDIZED_EFFECT_SIGN_CHECK = (
    OLD_ANALYSIS_DIR
    / "standardized_effect_size_4groups"
    / "standardized_effect_size_sign_check.csv"
)

EVIDENCE_TIER_DIR = OLD_ANALYSIS_DIR / "evidence_tiers_4groups"
FRONTAL_SUBREGION_DIR = OLD_ANALYSIS_DIR / "frontal_subregion_sensitivity"
TIER2_GENE_SET_DIR = OLD_ANALYSIS_DIR / "DE_candidate_gene_sets" / "Tier2_working_exploratory"
ENRICHMENT_DIR = OLD_ANALYSIS_DIR / "enrichment_Tier2"
ML_V1_DIR = OLD_ANALYSIS_DIR / "ML_pilot_Tier2"
OLD_ML_PRIORITIZED_DIR = OLD_ANALYSIS_DIR / "ML_prioritized_Tier2_genes"
OLD_IMAGES_DIR = OLD_ANALYSIS_DIR / "images"

CLEAN_IMAGES_DIR = CLEAN_OUTPUT_DIR / "images"
FINAL_REPORT = PROJECT_ROOT / "docs" / "MDD_Length_gt_100_clean_Tier2_to_ML_prioritized_Report.md"
FINAL_REPORT_COPY = CLEAN_OUTPUT_DIR / "MDD_Length_gt_100_clean_Tier2_to_ML_prioritized_Report.md"
PIPELINE_MANIFEST = CLEAN_OUTPUT_DIR / "pipeline_manifest.json"

