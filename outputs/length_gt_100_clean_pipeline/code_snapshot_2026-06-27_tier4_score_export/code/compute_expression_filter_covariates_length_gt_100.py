"""Compute independent filtering covariates from Length > 100 normalized matrices.

The covariates are computed from log2(GeTMM + 1) expression values only and do
not use diagnosis, p-values, beta estimates, or any MDD/control difference.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


REPO_DIR = Path(__file__).resolve().parents[1]
BASE_DIR = REPO_DIR / "outputs" / "length_gt_100_elisa_style"
AUDIT_DIR = BASE_DIR / "pvalue_correction_audit"
LMM_FILE = BASE_DIR / "LMM_all_genes_length_gt_100.csv"
OUT_FILE = AUDIT_DIR / "expression_filter_covariates_length_gt_100.csv"

REGION_STUBS = {
    "Frontal_Cortex": "frontal_cortex",
    "Insula": "insula",
    "Ventral_Striatum": "ventral_striatum",
    "Hippocampal_Formation": "hippocampal_formation",
}

SEARCH_DIRS = [
    BASE_DIR / "normalized_expression_length_gt_100_by_region",
    BASE_DIR,
    REPO_DIR / "outputs" / "step3_length_sensitivity",
    REPO_DIR / "outputs",
]


def find_matrix_file(region: str, stub: str) -> Path:
    expected_names = [
        f"normalized_GeTMM_log2_{stub}.csv",
        f"normalized_getmm_log2_{stub}.csv",
        f"{stub}.csv",
    ]
    matches: list[Path] = []

    for folder in SEARCH_DIRS:
        if not folder.exists():
            continue
        for expected in expected_names:
            matches.extend(folder.rglob(expected))

    matches = sorted(set(path.resolve() for path in matches))
    if not matches:
        raise FileNotFoundError(f"No normalized expression matrix found for {region} ({stub}).")

    preferred_root = (BASE_DIR / "normalized_expression_length_gt_100_by_region").resolve()
    preferred = [path for path in matches if preferred_root in path.parents or path.parent == preferred_root]
    chosen = preferred[0] if preferred else matches[0]
    print(f"Using normalized matrix for {region}: {chosen}")
    return chosen


def detect_orientation(matrix: pd.DataFrame, known_gene_ids: set[str], known_sample_prefixes: tuple[str, ...] = ("GSM",)) -> str:
    index_as_str = matrix.index.astype(str)
    cols_as_str = pd.Index(matrix.columns.astype(str))

    index_gene_hits = index_as_str.isin(known_gene_ids).mean()
    column_gene_hits = cols_as_str.isin(known_gene_ids).mean()
    index_sample_hits = index_as_str.str.startswith(known_sample_prefixes).mean()
    column_sample_hits = cols_as_str.str.startswith(known_sample_prefixes).mean()

    if index_gene_hits > 0.50 and column_sample_hits > 0.50:
        return "genes_rows"
    if column_gene_hits > 0.50 and index_sample_hits > 0.50:
        return "genes_columns"

    details = {
        "index_gene_hit_fraction": index_gene_hits,
        "column_gene_hit_fraction": column_gene_hits,
        "index_sample_hit_fraction": index_sample_hits,
        "column_sample_hit_fraction": column_sample_hits,
    }
    raise ValueError(f"Could not determine matrix orientation. Diagnostics: {details}")


def load_expression_matrix(path: Path, known_gene_ids: set[str]) -> pd.DataFrame:
    matrix = pd.read_csv(path, index_col=0)
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    orientation = detect_orientation(matrix, known_gene_ids)
    print(f"Detected orientation for {path.name}: {orientation}")
    if orientation == "genes_columns":
        matrix = matrix.T
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    if matrix.isna().all(axis=None):
        raise ValueError(f"Matrix contains no numeric expression values: {path}")
    return matrix


def main() -> None:
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    lmm = pd.read_csv(LMM_FILE, dtype={"gene_id": str})
    lmm["gene_id"] = lmm["gene_id"].astype(str)
    gene_symbol_map = (
        lmm[["brain_region", "gene_id", "gene_symbol"]]
        .drop_duplicates(["brain_region", "gene_id"])
        .set_index(["brain_region", "gene_id"])["gene_symbol"]
        .to_dict()
    )
    known_gene_ids = set(lmm["gene_id"])

    rows = []
    for region, stub in REGION_STUBS.items():
        path = find_matrix_file(region, stub)
        matrix = load_expression_matrix(path, known_gene_ids)

        genes = matrix.index.astype(str)
        rows.append(
            pd.DataFrame(
                {
                    "brain_region": region,
                    "gene_id": genes,
                    "gene_symbol": [gene_symbol_map.get((region, gene_id), "") for gene_id in genes],
                    "mean_expression_overall": matrix.mean(axis=1, skipna=True).to_numpy(),
                    "expression_variance_overall": matrix.var(axis=1, skipna=True, ddof=1).to_numpy(),
                    "detection_rate_overall": (matrix > 0).mean(axis=1, skipna=True).to_numpy(),
                    "n_samples_for_covariate": matrix.notna().sum(axis=1).to_numpy(),
                }
            )
        )

    covariates = pd.concat(rows, ignore_index=True)
    covariates.to_csv(OUT_FILE, index=False)
    print(f"Wrote {len(covariates):,} expression covariate rows to: {OUT_FILE}")
    print(covariates.groupby("brain_region").size().to_string())


if __name__ == "__main__":
    main()
