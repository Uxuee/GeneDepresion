from __future__ import annotations

import hashlib
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "outputs" / "length_gt_100_clean_pipeline" / "21_Statistical_method_benchmark"
INPUTS = OUT / "matched_inputs"
RESULTS = OUT / "region_results"
SUMMARY = OUT / "summary"
FIGURES = OUT / "figures"

LEGACY_REPO = Path(r"D:\genes\GeneDepresion")
GENEDATA = Path(r"D:\genes\genedata")
STEP1 = GENEDATA / "Step1_GSE102556_GeTMM_brain_area_groups"
LEGACY_FILTERED_STEP1 = LEGACY_REPO / "outputs" / "step1_length_sensitivity" / "Length_gt_100"
RAW_COUNTS = GENEDATA / "GSE102556_raw_counts_GRCh38.p13_NCBI.tsv"
ANNOTATION = GENEDATA / "Human.GRCh38.p13.annot.tsv"
BROAD_LM = LEGACY_REPO / "outputs" / "length_gt_100_clean_pipeline" / "04_standardized_effect_size" / "LMM_all_genes_with_standardized_effect_size_4groups.csv"
TIER2 = REPO / "outputs" / "length_gt_100_clean_pipeline" / "07_Tier2_candidate_gene_sets" / "Tier2_candidate_genes_all_regions.csv"

FRONTAL_BLOB = "f700cb323a47b24825b11bc1bba33cc07e6430bb"

REGIONS = [
    "Ventral Striatum",
    "Insula",
    "Hippocampal Formation",
    "CG25",
    "dlPFC BA8/9",
    "OFC BA11",
]

SLUGS = {
    "Ventral Striatum": "ventral_striatum",
    "Insula": "insula",
    "Hippocampal Formation": "hippocampal_formation",
    "CG25": "cg25",
    "dlPFC BA8/9": "dlpfc_ba8_9",
    "OFC BA11": "ofc_ba11",
}

EXPECTED_SAMPLES = {
    "Ventral Striatum": 48,
    "Insula": 48,
    "Hippocampal Formation": 43,
    "CG25": 27,
    "dlPFC BA8/9": 48,
    "OFC BA11": 47,
}

EXPECTED_GENES = {
    "Ventral Striatum": 23917,
    "Insula": 23090,
    "Hippocampal Formation": 24371,
    "CG25": 23941,
    "dlPFC BA8/9": 23941,
    "OFC BA11": 23941,
}

EXPECTED_TIER2 = {
    "Ventral Striatum": 146,
    "Insula": 70,
    "Hippocampal Formation": 127,
    "CG25": 186,
    "dlPFC BA8/9": 72,
    "OFC BA11": 70,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ensure_dirs() -> None:
    for directory in [OUT, INPUTS, RESULTS, SUMMARY, FIGURES]:
        directory.mkdir(parents=True, exist_ok=True)
