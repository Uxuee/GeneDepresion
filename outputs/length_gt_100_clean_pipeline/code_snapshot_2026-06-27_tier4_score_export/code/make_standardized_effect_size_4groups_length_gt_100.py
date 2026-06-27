"""Compute Elisa-style standardized effect sizes for Length > 100 results.

This script does not rerun normalization or LMM models and does not overwrite
the current exploratory DE result. It creates a separate output folder:

  outputs/length_gt_100_elisa_style/standardized_effect_size_4groups/

Standardized effect size:
  beta_diagnosis / sd_expression_gene

where sd_expression_gene is computed from the normalized log2(GeTMM + 1)
matrix across all samples in the grouped brain region.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO_DIR = Path(__file__).resolve().parents[1]
BASE_DIR = REPO_DIR / "outputs" / "length_gt_100_elisa_style"
FIXED_AUDIT_FILE = (
    BASE_DIR
    / "pvalue_correction_audit"
    / "LMM_all_genes_with_pvalue_corrections_FIXED.csv"
)
LMM_FILE = BASE_DIR / "LMM_all_genes_length_gt_100.csv"
OUT_DIR = BASE_DIR / "standardized_effect_size_4groups"

GROUPED_REGIONS = {
    "Frontal_Cortex": "frontal_cortex",
    "Insula": "insula",
    "Ventral_Striatum": "ventral_striatum",
    "Hippocampal_Formation": "hippocampal_formation",
}

REGION_DISPLAY = {
    "Frontal_Cortex": "Frontal Cortex",
    "Insula": "Insula",
    "Ventral_Striatum": "Ventral Striatum",
    "Hippocampal_Formation": "Hippocampal Formation",
}

SEARCH_DIRS = [
    BASE_DIR / "normalized_expression_length_gt_100_by_region",
    BASE_DIR / "step1_from_raw" / "Length_gt_100" / "02_normalized_expression" / "keep_all_samples",
    BASE_DIR,
    REPO_DIR / "outputs" / "step1_length_sensitivity" / "Length_gt_100" / "02_normalized_expression" / "keep_all_samples",
    REPO_DIR / "outputs",
]


def safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "_-" else "_" for ch in str(value))


def reset_output_subfolders() -> None:
    """Clear regenerated subfolders so stale filenames from older runs disappear."""
    for folder_name in ["candidate_gene_lists", "top_gene_tables", "overlaps", "plots"]:
        folder = OUT_DIR / folder_name
        if folder.exists():
            for path in folder.iterdir():
                if path.is_file():
                    path.unlink()
        folder.mkdir(parents=True, exist_ok=True)


def read_lmm_table() -> tuple[pd.DataFrame, Path]:
    path = FIXED_AUDIT_FILE if FIXED_AUDIT_FILE.exists() else LMM_FILE
    print(f"Reading LMM/effect table: {path}")
    df = pd.read_csv(path, dtype={"gene_id": str}, low_memory=False)
    required = ["brain_region", "gene_id", "gene_symbol", "beta_diagnosis", "p_value_diagnosis"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"LMM table is missing required columns: {missing}")
    df["gene_id"] = df["gene_id"].astype(str)
    df["beta_diagnosis"] = pd.to_numeric(df["beta_diagnosis"], errors="coerce")
    df["p_value_diagnosis"] = pd.to_numeric(df["p_value_diagnosis"], errors="coerce")
    return df, path


def validate_grouped_regions(df: pd.DataFrame) -> None:
    observed = set(df["brain_region"].dropna().unique())
    expected = set(GROUPED_REGIONS)
    if observed != expected:
        print("Observed regions:", sorted(observed))
        print("Expected grouped regions:", sorted(expected))
        raise ValueError("LMM table contains regions outside the expected 4 grouped regions.")


def find_matrix_file(region: str, stub: str) -> Path:
    expected_names = [
        f"normalized_GeTMM_log2_{stub}.csv",
        f"normalized_getmm_log2_{stub}.csv",
    ]
    matches: list[Path] = []
    for folder in SEARCH_DIRS:
        if not folder.exists():
            continue
        for name in expected_names:
            matches.extend(folder.rglob(name))

    matches = sorted(set(path.resolve() for path in matches))
    if not matches:
        raise FileNotFoundError(f"No normalized expression matrix found for {region}.")

    preferred = (BASE_DIR / "normalized_expression_length_gt_100_by_region").resolve()
    preferred_matches = [path for path in matches if path.parent == preferred]
    chosen = preferred_matches[0] if preferred_matches else matches[0]
    print(f"Using normalized expression matrix for {region}: {chosen}")
    return chosen


def detect_orientation(matrix: pd.DataFrame, known_gene_ids: set[str]) -> str:
    row_ids = matrix.index.astype(str)
    col_ids = pd.Index(matrix.columns.astype(str))
    row_gene_fraction = row_ids.isin(known_gene_ids).mean()
    col_gene_fraction = col_ids.isin(known_gene_ids).mean()
    row_sample_fraction = row_ids.str.startswith("GSM").mean()
    col_sample_fraction = col_ids.str.startswith("GSM").mean()

    if row_gene_fraction > 0.50 and col_sample_fraction > 0.50:
        return "genes_rows"
    if col_gene_fraction > 0.50 and row_sample_fraction > 0.50:
        return "genes_columns"

    raise ValueError(
        "Could not determine expression matrix orientation. "
        f"row_gene_fraction={row_gene_fraction:.3f}, col_gene_fraction={col_gene_fraction:.3f}, "
        f"row_sample_fraction={row_sample_fraction:.3f}, col_sample_fraction={col_sample_fraction:.3f}"
    )


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
        raise ValueError(f"Expression matrix contains no numeric values: {path}")
    return matrix


def compute_expression_summaries(lmm: pd.DataFrame) -> pd.DataFrame:
    known_gene_ids = set(lmm["gene_id"])
    rows = []
    used_files = []

    for region, stub in GROUPED_REGIONS.items():
        path = find_matrix_file(region, stub)
        matrix = load_expression_matrix(path, known_gene_ids)
        genes = matrix.index.astype(str)
        used_files.append({"grouped_region": region, "matrix_file_used": str(path)})
        rows.append(
            pd.DataFrame(
                {
                    "grouped_region": region,
                    "gene_id": genes,
                    "sd_expression_gene": matrix.std(axis=1, skipna=True, ddof=1).to_numpy(),
                    "mean_expression_gene": matrix.mean(axis=1, skipna=True).to_numpy(),
                    "n_samples_used_for_sd": matrix.notna().sum(axis=1).to_numpy(),
                    "min_expression_gene": matrix.min(axis=1, skipna=True).to_numpy(),
                    "max_expression_gene": matrix.max(axis=1, skipna=True).to_numpy(),
                }
            )
        )

    used_files_df = pd.DataFrame(used_files)
    used_files_df.to_csv(OUT_DIR / "normalized_expression_files_used.csv", index=False)
    return pd.concat(rows, ignore_index=True)


def add_direction(values: pd.Series) -> pd.Series:
    return pd.Series(
        np.select([values > 0, values < 0, values == 0], ["up_in_MDD", "down_in_MDD", "no_change"], default="no_change"),
        index=values.index,
    )


def merge_and_compute(lmm: pd.DataFrame, expr_summary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    lmm = lmm.copy()
    lmm["grouped_region"] = lmm["brain_region"]
    before = len(lmm)
    merged = lmm.merge(expr_summary, on=["grouped_region", "gene_id"], how="left", sort=False)
    merged["sd_expression_gene"] = pd.to_numeric(merged["sd_expression_gene"], errors="coerce")
    valid_sd = merged["sd_expression_gene"].notna() & (merged["sd_expression_gene"] != 0)
    merged["standardized_effect_size_diagnosis"] = np.nan
    merged.loc[valid_sd, "standardized_effect_size_diagnosis"] = (
        merged.loc[valid_sd, "beta_diagnosis"] / merged.loc[valid_sd, "sd_expression_gene"]
    )
    merged["abs_standardized_effect_size_diagnosis"] = merged["standardized_effect_size_diagnosis"].abs()
    merged["abs_beta_diagnosis"] = merged["beta_diagnosis"].abs()
    merged["direction_beta"] = add_direction(merged["beta_diagnosis"])
    merged["direction_standardized_effect_size"] = add_direction(merged["standardized_effect_size_diagnosis"])
    sign_check_mask = valid_sd & (merged["beta_diagnosis"] != 0) & merged["standardized_effect_size_diagnosis"].notna()
    merged["standardized_effect_sign_matches_beta"] = True
    merged.loc[sign_check_mask, "standardized_effect_sign_matches_beta"] = (
        np.sign(merged.loc[sign_check_mask, "beta_diagnosis"])
        == np.sign(merged.loc[sign_check_mask, "standardized_effect_size_diagnosis"])
    )
    merged["signed_standardized_logp_score"] = (
        np.sign(merged["standardized_effect_size_diagnosis"]) * -np.log10(merged["p_value_diagnosis"].clip(lower=np.finfo(float).tiny))
    )

    merge_report = pd.DataFrame(
        [
            {
                "rows_in_LMM_before_merge": before,
                "rows_after_merge": len(merged),
                "genes_with_valid_sd_expression_gene": int(valid_sd.sum()),
                "genes_with_missing_sd_expression_gene": int(merged["sd_expression_gene"].isna().sum()),
                "genes_with_sd_expression_gene_eq_0": int((merged["sd_expression_gene"] == 0).sum()),
                "genes_excluded_from_standardized_effect_size": int((~valid_sd).sum()),
                "merge_key_used": "grouped_region + gene_id",
                "gene_symbol_fallback_used": False,
                "rows_checked_for_beta_standardized_effect_sign_match": int(sign_check_mask.sum()),
                "sign_mismatch_count": int((~merged.loc[sign_check_mask, "standardized_effect_sign_matches_beta"]).sum()),
                "sign_match_confirmed": bool(merged.loc[sign_check_mask, "standardized_effect_sign_matches_beta"].all()),
            }
        ]
    )
    merge_report.to_csv(OUT_DIR / "standardized_effect_size_merge_summary.csv", index=False)
    merge_report[
        [
            "rows_checked_for_beta_standardized_effect_sign_match",
            "sign_mismatch_count",
            "sign_match_confirmed",
        ]
    ].to_csv(OUT_DIR / "standardized_effect_size_sign_check_summary.csv", index=False)
    return merged, merge_report


def mask_col(df: pd.DataFrame, col: str, cutoff: float) -> pd.Series:
    if col not in df.columns:
        return pd.Series(False, index=df.index)
    return pd.to_numeric(df[col], errors="coerce") < cutoff


def add_filter_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["current_baseline_p001_absBeta02"] = (out["p_value_diagnosis"] < 0.01) & (out["abs_beta_diagnosis"] > 0.2)
    out["nominal_p001_absStdES08"] = (out["p_value_diagnosis"] < 0.01) & (out["abs_standardized_effect_size_diagnosis"] > 0.8)
    out["BH005_absStdES08"] = mask_col(out, "fdr_BH", 0.05) & (out["abs_standardized_effect_size_diagnosis"] > 0.8)
    out["BH010_absStdES08"] = mask_col(out, "fdr_BH", 0.10) & (out["abs_standardized_effect_size_diagnosis"] > 0.8)
    if "qvalue_Storey" in out.columns:
        out["Storey_q010_absStdES08"] = mask_col(out, "qvalue_Storey", 0.10) & (out["abs_standardized_effect_size_diagnosis"] > 0.8)
    if "ihw_adj_pvalue" in out.columns:
        out["IHW_q010_absStdES08"] = mask_col(out, "ihw_adj_pvalue", 0.10) & (out["abs_standardized_effect_size_diagnosis"] > 0.8)
    return out


def save_candidate_lists(df: pd.DataFrame) -> None:
    out_dir = OUT_DIR / "candidate_gene_lists"
    out_dir.mkdir(parents=True, exist_ok=True)
    filter_labels = {
        "BH005_absStdES08": "Elisa_strict_BH005_absStdES08",
        "BH010_absStdES08": "Elisa_relaxed_BH010_absStdES08",
        "Storey_q010_absStdES08": "Storey_q010_absStdES08",
        "IHW_q010_absStdES08": "IHW_q010_absStdES08",
        "nominal_p001_absStdES08": "Nominal_p001_absStdES08",
        "current_baseline_p001_absBeta02": "Current_baseline_p001_absBeta02",
    }
    for region, group in df.groupby("grouped_region", sort=False):
        for col, label in filter_labels.items():
            if col not in df.columns:
                continue
            selected = group[group[col]].sort_values("p_value_diagnosis", na_position="last")
            selected.to_csv(out_dir / f"{safe_name(region)}_{label}.csv", index=False)


def make_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for region, group in df.groupby("grouped_region", sort=False):
        nominal = group[group["nominal_p001_absStdES08"]]
        row = {
            "grouped_region": region,
            "genes_tested": len(group),
            "genes_with_valid_sd": int(group["standardized_effect_size_diagnosis"].notna().sum()),
            "current_baseline_p001_absBeta02": int(group["current_baseline_p001_absBeta02"].sum()),
            "nominal_p001_absStdES08": int(group["nominal_p001_absStdES08"].sum()),
            "BH005_absStdES08": int(group["BH005_absStdES08"].sum()),
            "BH010_absStdES08": int(group["BH010_absStdES08"].sum()),
            "up_in_MDD_nominal_absStdES08": int((nominal["direction_standardized_effect_size"] == "up_in_MDD").sum()),
            "down_in_MDD_nominal_absStdES08": int((nominal["direction_standardized_effect_size"] == "down_in_MDD").sum()),
            "median_abs_beta": group["abs_beta_diagnosis"].median(),
            "max_abs_beta": group["abs_beta_diagnosis"].max(),
            "median_abs_standardized_effect_size": group["abs_standardized_effect_size_diagnosis"].median(),
            "max_abs_standardized_effect_size": group["abs_standardized_effect_size_diagnosis"].max(),
        }
        if "Storey_q010_absStdES08" in group.columns:
            row["Storey_q010_absStdES08"] = int(group["Storey_q010_absStdES08"].sum())
        if "IHW_q010_absStdES08" in group.columns:
            row["IHW_q010_absStdES08"] = int(group["IHW_q010_absStdES08"].sum())
        rows.append(row)
    summary = pd.DataFrame(rows)
    ordered = [
        "grouped_region",
        "genes_tested",
        "genes_with_valid_sd",
        "current_baseline_p001_absBeta02",
        "nominal_p001_absStdES08",
        "BH005_absStdES08",
        "BH010_absStdES08",
        "Storey_q010_absStdES08",
        "IHW_q010_absStdES08",
        "up_in_MDD_nominal_absStdES08",
        "down_in_MDD_nominal_absStdES08",
        "median_abs_beta",
        "max_abs_beta",
        "median_abs_standardized_effect_size",
        "max_abs_standardized_effect_size",
    ]
    return summary[[col for col in ordered if col in summary.columns]]


def save_top_gene_tables(df: pd.DataFrame) -> None:
    out_dir = OUT_DIR / "top_gene_tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    for region, group in df.groupby("grouped_region", sort=False):
        prefix = safe_name(region)
        group.sort_values("abs_standardized_effect_size_diagnosis", ascending=False, na_position="last").head(100).to_csv(
            out_dir / f"{prefix}_top100_by_abs_standardized_effect_size.csv", index=False
        )
        group[group["abs_standardized_effect_size_diagnosis"] > 0.8].sort_values("p_value_diagnosis", na_position="last").head(100).to_csv(
            out_dir / f"{prefix}_top100_by_pvalue_among_absStdES08.csv", index=False
        )
        group.assign(abs_signed_score=group["signed_standardized_logp_score"].abs()).sort_values(
            "abs_signed_score", ascending=False, na_position="last"
        ).drop(columns=["abs_signed_score"]).head(100).to_csv(
            out_dir / f"{prefix}_top100_by_signed_standardized_logp_score.csv", index=False
        )


def save_overlap_matrices(df: pd.DataFrame) -> None:
    out_dir = OUT_DIR / "overlaps"
    out_dir.mkdir(parents=True, exist_ok=True)
    base_cols = [
        "current_baseline_p001_absBeta02",
        "nominal_p001_absStdES08",
        "BH010_absStdES08",
        "Storey_q010_absStdES08",
        "IHW_q010_absStdES08",
    ]
    for region, group in df.groupby("grouped_region", sort=False):
        sets = {col: set(group.index[group[col]]) for col in base_cols if col in group.columns}
        names = list(sets)
        matrix = pd.DataFrame(index=names, columns=names, dtype=int)
        for a in names:
            for b in names:
                matrix.loc[a, b] = len(sets[a] & sets[b])
        matrix.to_csv(out_dir / f"{safe_name(region)}_overlap_matrix.csv")


def save_plots(df: pd.DataFrame) -> None:
    out_dir = OUT_DIR / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    for region, group in df.groupby("grouped_region", sort=False):
        prefix = safe_name(region)
        p = group["p_value_diagnosis"].clip(lower=np.finfo(float).tiny)
        minus_logp = -np.log10(p)

        def scatter(x, y, xlabel, ylabel, title, filename, color="#3b6ea8", mask=None):
            plt.figure(figsize=(7, 6))
            if mask is None:
                plt.scatter(x, y, s=7, alpha=0.45, color=color)
            else:
                plt.scatter(x[~mask], y[~mask], s=7, alpha=0.35, color="#747b85")
                plt.scatter(x[mask], y[mask], s=12, alpha=0.8, color="#b42318")
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title(title)
            plt.tight_layout()
            plt.savefig(out_dir / filename, dpi=170)
            plt.close()

        scatter(
            group["beta_diagnosis"],
            group["standardized_effect_size_diagnosis"],
            "beta_diagnosis",
            "standardized_effect_size_diagnosis",
            f"Beta vs standardized effect size: {region}",
            f"{prefix}_beta_vs_standardized_effect_size.png",
        )

        for col, label in [
            ("beta_diagnosis", "beta_diagnosis"),
            ("standardized_effect_size_diagnosis", "standardized_effect_size_diagnosis"),
            ("abs_standardized_effect_size_diagnosis", "abs_standardized_effect_size_diagnosis"),
        ]:
            plt.figure(figsize=(7, 5))
            plt.hist(group[col].dropna(), bins=50, color="#667085", edgecolor="white")
            plt.xlabel(label)
            plt.ylabel("Genes")
            plt.title(f"Histogram: {label} {region}")
            plt.tight_layout()
            plt.savefig(out_dir / f"{prefix}_histogram_{label}.png", dpi=160)
            plt.close()

        scatter(
            group["beta_diagnosis"],
            minus_logp,
            "beta_diagnosis",
            "-log10(p_value_diagnosis)",
            f"Beta volcano, current baseline highlighted: {region}",
            f"{prefix}_volcano_beta_current_baseline.png",
            mask=group["current_baseline_p001_absBeta02"],
        )
        scatter(
            group["standardized_effect_size_diagnosis"],
            minus_logp,
            "standardized_effect_size_diagnosis",
            "-log10(p_value_diagnosis)",
            f"Standardized-effect volcano, nominal Elisa highlighted: {region}",
            f"{prefix}_volcano_standardized_effect_nominal.png",
            mask=group["nominal_p001_absStdES08"],
        )
        scatter(
            group["sd_expression_gene"],
            group["abs_beta_diagnosis"],
            "sd_expression_gene",
            "abs(beta_diagnosis)",
            f"Expression SD vs abs beta: {region}",
            f"{prefix}_sd_expression_vs_abs_beta.png",
        )
        scatter(
            group["sd_expression_gene"],
            group["abs_standardized_effect_size_diagnosis"],
            "sd_expression_gene",
            "abs_standardized_effect_size_diagnosis",
            f"Expression SD vs abs standardized effect size: {region}",
            f"{prefix}_sd_expression_vs_abs_standardized_effect_size.png",
        )


def top20_markdown(df: pd.DataFrame, sort_col: str, ascending: bool, filter_mask: pd.Series | None = None) -> str:
    lines = []
    data = df if filter_mask is None else df[filter_mask]
    for region, group in data.groupby("grouped_region", sort=False):
        lines.append(f"### {REGION_DISPLAY.get(region, region)}")
        top = group.sort_values(sort_col, ascending=ascending, na_position="last").head(20)
        if top.empty:
            lines.extend(["", "No genes in this category.", ""])
            continue
        lines.append("")
        lines.append("| rank | gene | beta | std_effect_size | p_value | direction |")
        lines.append("|---:|---|---:|---:|---:|---|")
        for rank, (_, row) in enumerate(top.iterrows(), start=1):
            gene = row["gene_symbol"] if isinstance(row.get("gene_symbol"), str) and row["gene_symbol"] else row["gene_id"]
            lines.append(
                f"| {rank} | {gene} | {row['beta_diagnosis']:.4g} | {row['standardized_effect_size_diagnosis']:.4g} | "
                f"{row['p_value_diagnosis']:.4g} | {row['direction_standardized_effect_size']} |"
            )
        lines.append("")
    return "\n".join(lines)


def write_report(df: pd.DataFrame, summary: pd.DataFrame, merge_report: pd.DataFrame, source_file: Path) -> None:
    report = f"""# Standardized Effect Size Analysis: Length > 100, 4 Grouped Regions

## Scope

No previous result was overwritten. This analysis writes only inside `standardized_effect_size_4groups/` and uses the existing Length > 100 LMM and normalized-expression outputs.

Grouped brain regions used:

- Frontal Cortex
- Insula
- Ventral Striatum
- Hippocampal Formation

Input effect table used:

`{source_file}`

## Preserved Current Baseline

`p_value_diagnosis < 0.01 and abs(beta_diagnosis) > 0.2`

## New Standardized Effect Size

`standardized_effect_size_diagnosis = beta_diagnosis / sd_expression_gene`

`beta_diagnosis` is the raw LMM diagnosis coefficient. `sd_expression_gene` is the standard deviation of the same gene's normalized log2(GeTMM + 1) expression across all samples in the same grouped region. `standardized_effect_size_diagnosis` is beta scaled by each gene's expression variability. The sign is preserved and should match `beta_diagnosis`; this was checked below. FDR, q-value, and IHW adjusted p-values are p-value corrections, not effect sizes.

Methods note: the standardized effect size is computed as the LMM diagnosis coefficient divided by the gene-wise standard deviation of normalized expression across the same samples used for the model. This is intended to reproduce the Elisa-style interpretation of the case/control coefficient as an effect size scaled by the data distribution. A residual-standardized version, beta divided by model residual SD, may be computed only as an optional sensitivity analysis if residual SD is available, but it should not replace the primary Elisa-style effect size.

## Merge Check

{merge_report.to_markdown(index=False)}

## Summary By Grouped Region

{summary.to_markdown(index=False)}

## Current Baseline Vs Nominal Standardized-Effect List

The current baseline uses raw beta magnitude. The nominal standardized-effect list uses the same raw p-value threshold but replaces raw beta magnitude with `abs(standardized_effect_size_diagnosis) > 0.8`.

## Top 20 By Absolute Standardized Effect Size

{top20_markdown(df, "abs_standardized_effect_size_diagnosis", False)}

## Top 20 By P-Value Among Genes With Abs Standardized Effect Size > 0.8

{top20_markdown(df, "p_value_diagnosis", True, df["abs_standardized_effect_size_diagnosis"] > 0.8)}

## Recommended Use

1. Strict reporting: use corrected DE genes with BH/Storey/IHW support plus `abs_standardized_effect_size_diagnosis > 0.8`.
2. Exploratory enrichment: use nominal standardized-effect genes, clearly labeled as exploratory.
3. ML feature selection: use ranked signatures from `abs_standardized_effect_size_diagnosis` and `signed_standardized_logp_score`.
4. AHBA module mapping: use exploratory candidates and ranked standardized signatures, clearly labeled as exploratory.

Do not call nominal genes final significant DE genes.
"""
    (OUT_DIR / "Standardized_Effect_Size_4Groups_Report.md").write_text(report, encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    reset_output_subfolders()
    lmm, source_file = read_lmm_table()
    print("Column choices:")
    for col in ["brain_region", "gene_id", "gene_symbol", "beta_diagnosis", "p_value_diagnosis", "fdr_BH", "fdr_BY", "qvalue_Storey", "ihw_adj_pvalue"]:
        print(f"  {col}: {'present' if col in lmm.columns else 'not present'}")

    validate_grouped_regions(lmm)
    expr_summary = compute_expression_summaries(lmm)
    merged, merge_report = merge_and_compute(lmm, expr_summary)
    merged = add_filter_columns(merged)

    merged.to_csv(OUT_DIR / "LMM_all_genes_with_standardized_effect_size_4groups.csv", index=False)
    save_candidate_lists(merged)
    summary = make_summary(merged)
    summary.to_csv(OUT_DIR / "standardized_effect_size_4groups_summary.csv", index=False)
    save_top_gene_tables(merged)
    save_overlap_matrices(merged)
    save_plots(merged)
    write_report(merged, summary, merge_report, source_file)
    print(f"Wrote standardized effect size outputs to: {OUT_DIR}")


if __name__ == "__main__":
    main()
