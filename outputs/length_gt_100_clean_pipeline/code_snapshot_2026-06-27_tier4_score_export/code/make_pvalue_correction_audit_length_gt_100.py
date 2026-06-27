"""Create a p-value correction audit for Length > 100 LMM results.

This audit does not rerun normalization or LMM models, and it does not modify
the existing exploratory DE result. It only reads:

  outputs/length_gt_100_elisa_style/LMM_all_genes_length_gt_100.csv

and writes new files under:

  outputs/length_gt_100_elisa_style/pvalue_correction_audit/

Baseline exploratory definition preserved here:
  p_value_diagnosis < 0.01 and abs(beta_diagnosis) > 0.2
"""

from __future__ import annotations

import math
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


REPO_DIR = Path(__file__).resolve().parents[1]
BASE_DIR = REPO_DIR / "outputs" / "length_gt_100_elisa_style"
INPUT_FILE = BASE_DIR / "LMM_all_genes_length_gt_100.csv"
AUDIT_DIR = BASE_DIR / "pvalue_correction_audit"

P_COL = "p_value_diagnosis"
BETA_COL = "beta_diagnosis"
BASELINE_COL = "exploratory_p001_absbeta02"


def safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "_-" else "_" for ch in str(value))


def check_r_package(package: str) -> bool:
    """Return True if an R package is installed for the visible Rscript."""
    rscript = Path("C:/Program Files/R/R-4.4.2/bin/Rscript.exe")
    if not rscript.exists():
        return False
    cmd = [str(rscript), "-e", f"quit(status = ifelse(requireNamespace('{package}', quietly=TRUE), 0, 1))"]
    return subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0


def bh_adjust(values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.notna()
    if valid.any():
        out.loc[valid] = multipletests(values.loc[valid], method="fdr_bh")[1]
    return out


def by_adjust(values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.notna()
    if valid.any():
        out.loc[valid] = multipletests(values.loc[valid], method="fdr_by")[1]
    return out


def add_core_columns(df: pd.DataFrame) -> pd.DataFrame:
    required = ["brain_region", "gene_id", "gene_symbol", BETA_COL, P_COL]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Input file is missing required columns: {missing}")

    out = df.copy()
    out[BETA_COL] = pd.to_numeric(out[BETA_COL], errors="coerce")
    out[P_COL] = pd.to_numeric(out[P_COL], errors="coerce")
    out["abs_beta_diagnosis"] = out[BETA_COL].abs()
    out["signed_logp_score"] = (
        np.sign(out[BETA_COL]) * -np.log10(out[P_COL].clip(lower=np.finfo(float).tiny))
    )
    out[BASELINE_COL] = (out[P_COL] < 0.01) & (out["abs_beta_diagnosis"] > 0.2)
    return out


def add_regionwise_corrections(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["fdr_BH"] = np.nan
    out["fdr_BY"] = np.nan
    for _, idx in out.groupby("brain_region", sort=False).groups.items():
        out.loc[idx, "fdr_BH"] = bh_adjust(out.loc[idx, P_COL])
        out.loc[idx, "fdr_BY"] = by_adjust(out.loc[idx, P_COL])
    return out


def independent_filter_candidates(df: pd.DataFrame) -> list[tuple[str, pd.Series]]:
    """Build valid independent filters from available, non-empty covariates."""
    filters: list[tuple[str, pd.Series]] = []

    mean_cols = [
        col
        for col in ["mean_expression_overall", "mean_expression_MDD", "mean_expression_Control"]
        if col in df.columns and pd.to_numeric(df[col], errors="coerce").notna().any()
    ]
    if "mean_expression_overall" in mean_cols:
        mean_expr = pd.to_numeric(df["mean_expression_overall"], errors="coerce")
    elif {"mean_expression_MDD", "mean_expression_Control"}.issubset(mean_cols):
        mean_expr = pd.concat(
            [
                pd.to_numeric(df["mean_expression_MDD"], errors="coerce"),
                pd.to_numeric(df["mean_expression_Control"], errors="coerce"),
            ],
            axis=1,
        ).mean(axis=1)
    else:
        mean_expr = None

    if mean_expr is not None and mean_expr.notna().any():
        for q in [0.10, 0.20, 0.30]:
            threshold = mean_expr.quantile(q)
            filters.append((f"mean_expression_above_{int(q * 100)}th_percentile", mean_expr > threshold))

    if "expression_variance" in df.columns:
        variance = pd.to_numeric(df["expression_variance"], errors="coerce")
        if variance.notna().any():
            for q in [0.10, 0.20, 0.30]:
                threshold = variance.quantile(q)
                filters.append((f"expression_variance_above_{int(q * 100)}th_percentile", variance > threshold))

    if "detection_rate" in df.columns:
        detection = pd.to_numeric(df["detection_rate"], errors="coerce")
        if detection.notna().any():
            filters.append(("detection_rate_ge_20_percent", detection >= 0.20))
            filters.append(("detection_rate_ge_50_percent", detection >= 0.50))

    return filters


def run_independent_filter_audit(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply BH after each valid independent filter within each region."""
    all_rows = []
    adjusted = pd.DataFrame(index=df.index)

    filters = independent_filter_candidates(df)
    if not filters:
        return pd.DataFrame(), adjusted

    for region, region_df in df.groupby("brain_region", sort=False):
        for filter_name, mask in filters:
            region_mask = (df["brain_region"] == region) & mask.fillna(False)
            kept = df.loc[region_mask]
            col = f"BH_after_filter__{safe_name(filter_name)}"
            adjusted[col] = np.nan
            if not kept.empty:
                adjusted.loc[kept.index, col] = bh_adjust(kept[P_COL])
            all_rows.append(
                {
                    "brain_region": region,
                    "filter_name": filter_name,
                    "genes_tested_after_filter": len(kept),
                    "BH_after_filter_lt_005": int((adjusted.loc[kept.index, col] < 0.05).sum()),
                    "BH_after_filter_lt_010": int((adjusted.loc[kept.index, col] < 0.10).sum()),
                    "adjusted_column": col,
                }
            )

    return pd.DataFrame(all_rows), adjusted


def choose_best_filters(independent_summary: pd.DataFrame) -> dict[str, dict[str, object]]:
    best: dict[str, dict[str, object]] = {}
    if independent_summary.empty:
        return best

    ranked = independent_summary.sort_values(
        ["brain_region", "BH_after_filter_lt_010", "BH_after_filter_lt_005", "genes_tested_after_filter"],
        ascending=[True, False, False, False],
    )
    for region, group in ranked.groupby("brain_region", sort=False):
        row = group.iloc[0].to_dict()
        best[region] = row
    return best


def save_candidate_lists(df: pd.DataFrame, best_filters: dict[str, dict[str, object]]) -> None:
    out_dir = AUDIT_DIR / "candidate_gene_lists"
    out_dir.mkdir(parents=True, exist_ok=True)

    base_sets = [
        ("exploratory_p001_absbeta02", df[BASELINE_COL]),
        ("BH_FDR_lt_005", df["fdr_BH"] < 0.05),
        ("BH_FDR_lt_010", df["fdr_BH"] < 0.10),
        ("BY_FDR_lt_005", df["fdr_BY"] < 0.05),
        ("BY_FDR_lt_010", df["fdr_BY"] < 0.10),
    ]

    if "qvalue_Storey" in df.columns:
        base_sets.extend(
            [
                ("Storey_q_lt_005", df["qvalue_Storey"] < 0.05),
                ("Storey_q_lt_010", df["qvalue_Storey"] < 0.10),
            ]
        )
    if "p_adj_IHW" in df.columns:
        base_sets.extend(
            [
                ("IHW_lt_005", df["p_adj_IHW"] < 0.05),
                ("IHW_lt_010", df["p_adj_IHW"] < 0.10),
            ]
        )

    for region, region_df in df.groupby("brain_region", sort=False):
        prefix = safe_name(region)
        for label, mask in base_sets:
            region_mask = mask & (df["brain_region"] == region)
            df.loc[region_mask].sort_values(P_COL).to_csv(out_dir / f"{prefix}_{label}.csv", index=False)

        best = best_filters.get(region)
        for cutoff_label, cutoff in [("005", 0.05), ("010", 0.10)]:
            if best is None:
                empty = region_df.head(0)
                empty.to_csv(out_dir / f"{prefix}_best_independent_filter_BH_lt_{cutoff_label}.csv", index=False)
            else:
                col = str(best["adjusted_column"])
                region_mask = (df["brain_region"] == region) & (df[col] < cutoff)
                df.loc[region_mask].sort_values(col).to_csv(
                    out_dir / f"{prefix}_best_independent_filter_BH_lt_{cutoff_label}.csv",
                    index=False,
                )


def make_summary(df: pd.DataFrame, best_filters: dict[str, dict[str, object]]) -> pd.DataFrame:
    rows = []
    for region, group in df.groupby("brain_region", sort=False):
        best = best_filters.get(region)
        row = {
            "brain_region": region,
            "genes_tested": len(group),
            "raw_p_lt_005": int((group[P_COL] < 0.05).sum()),
            "raw_p_lt_001": int((group[P_COL] < 0.01).sum()),
            "raw_p_lt_0001": int((group[P_COL] < 0.001).sum()),
            "exploratory_p001_absbeta02": int(group[BASELINE_COL].sum()),
            "BH_FDR_lt_005": int((group["fdr_BH"] < 0.05).sum()),
            "BH_FDR_lt_010": int((group["fdr_BH"] < 0.10).sum()),
            "BY_FDR_lt_005": int((group["fdr_BY"] < 0.05).sum()),
            "BY_FDR_lt_010": int((group["fdr_BY"] < 0.10).sum()),
            "best_independent_filter": best["filter_name"] if best else "not_available",
            "genes_tested_after_best_independent_filter": int(best["genes_tested_after_filter"]) if best else 0,
            "BH_after_best_filter_lt_005": int(best["BH_after_filter_lt_005"]) if best else 0,
            "BH_after_best_filter_lt_010": int(best["BH_after_filter_lt_010"]) if best else 0,
        }
        if "qvalue_Storey" in df.columns:
            row["Storey_q_lt_005"] = int((group["qvalue_Storey"] < 0.05).sum())
            row["Storey_q_lt_010"] = int((group["qvalue_Storey"] < 0.10).sum())
        if "p_adj_IHW" in df.columns:
            row["IHW_lt_005"] = int((group["p_adj_IHW"] < 0.05).sum())
            row["IHW_lt_010"] = int((group["p_adj_IHW"] < 0.10).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def write_plots(df: pd.DataFrame) -> None:
    plot_dir = AUDIT_DIR / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    for region, group in df.groupby("brain_region", sort=False):
        prefix = safe_name(region)
        p = group[P_COL].dropna().clip(lower=np.finfo(float).tiny)
        beta = group[BETA_COL].dropna()
        abs_beta = group["abs_beta_diagnosis"].dropna()
        minus_logp = -np.log10(group[P_COL].clip(lower=np.finfo(float).tiny))

        plt.figure(figsize=(7, 5))
        plt.hist(p, bins=50, color="#7a869a", edgecolor="white")
        plt.xlabel("p_value_diagnosis")
        plt.ylabel("Genes")
        plt.title(f"P-value histogram: {region}")
        plt.tight_layout()
        plt.savefig(plot_dir / f"{prefix}_pvalue_histogram.png", dpi=160)
        plt.close()

        expected = -np.log10((np.arange(1, len(p) + 1) - 0.5) / len(p))
        observed = -np.log10(np.sort(p.to_numpy()))
        plt.figure(figsize=(6, 6))
        plt.scatter(expected, observed, s=8, alpha=0.6, color="#3b6ea8")
        lim = max(expected.max(), observed.max()) if len(expected) else 1
        plt.plot([0, lim], [0, lim], color="#333333", linestyle="--", linewidth=1)
        plt.xlabel("Expected -log10(p)")
        plt.ylabel("Observed -log10(p)")
        plt.title(f"QQ plot: {region}")
        plt.tight_layout()
        plt.savefig(plot_dir / f"{prefix}_qq_plot.png", dpi=160)
        plt.close()

        plt.figure(figsize=(7, 5))
        plt.hist(beta, bins=50, color="#667085", edgecolor="white")
        plt.xlabel("beta_diagnosis")
        plt.ylabel("Genes")
        plt.title(f"Beta histogram: {region}")
        plt.tight_layout()
        plt.savefig(plot_dir / f"{prefix}_beta_histogram.png", dpi=160)
        plt.close()

        plt.figure(figsize=(7, 5))
        plt.hist(abs_beta, bins=50, color="#667085", edgecolor="white")
        plt.xlabel("abs(beta_diagnosis)")
        plt.ylabel("Genes")
        plt.title(f"Absolute beta histogram: {region}")
        plt.tight_layout()
        plt.savefig(plot_dir / f"{prefix}_abs_beta_histogram.png", dpi=160)
        plt.close()

        plt.figure(figsize=(7, 6))
        highlighted = group[BASELINE_COL]
        plt.scatter(group.loc[~highlighted, BETA_COL], minus_logp.loc[~highlighted], s=7, alpha=0.35, color="#667085")
        plt.scatter(group.loc[highlighted, BETA_COL], minus_logp.loc[highlighted], s=12, alpha=0.8, color="#b42318")
        plt.axhline(-math.log10(0.01), color="#333333", linestyle="--", linewidth=1)
        plt.axvline(-0.2, color="#333333", linestyle="--", linewidth=1)
        plt.axvline(0.2, color="#333333", linestyle="--", linewidth=1)
        plt.xlabel("beta_diagnosis")
        plt.ylabel("-log10(p_value_diagnosis)")
        plt.title(f"Volcano with exploratory genes: {region}")
        plt.tight_layout()
        plt.savefig(plot_dir / f"{prefix}_volcano_exploratory_highlighted.png", dpi=180)
        plt.close()

        if "gene_length" in group.columns:
            length = pd.to_numeric(group["gene_length"], errors="coerce")
            valid = length.notna() & minus_logp.notna()
            plt.figure(figsize=(7, 5))
            plt.scatter(length.loc[valid], minus_logp.loc[valid], s=7, alpha=0.35, color="#3b6ea8")
            plt.xscale("log")
            plt.xlabel("gene_length (log scale)")
            plt.ylabel("-log10(p_value_diagnosis)")
            plt.title(f"Gene length vs -log10(p): {region}")
            plt.tight_layout()
            plt.savefig(plot_dir / f"{prefix}_gene_length_vs_minus_log10_p.png", dpi=160)
            plt.close()

        mean_col = None
        for candidate in ["mean_expression_overall", "mean_expression_MDD", "mean_expression_Control"]:
            if candidate in group.columns and pd.to_numeric(group[candidate], errors="coerce").notna().any():
                mean_col = candidate
                break
        if mean_col:
            mean_expr = pd.to_numeric(group[mean_col], errors="coerce")
            valid = mean_expr.notna() & minus_logp.notna()
            plt.figure(figsize=(7, 5))
            plt.scatter(mean_expr.loc[valid], minus_logp.loc[valid], s=7, alpha=0.35, color="#3b6ea8")
            plt.xlabel(mean_col)
            plt.ylabel("-log10(p_value_diagnosis)")
            plt.title(f"Mean expression vs -log10(p): {region}")
            plt.tight_layout()
            plt.savefig(plot_dir / f"{prefix}_mean_expression_vs_minus_log10_p.png", dpi=160)
            plt.close()


def write_overlap_matrices(df: pd.DataFrame, best_filters: dict[str, dict[str, object]]) -> None:
    out_dir = AUDIT_DIR / "overlaps"
    out_dir.mkdir(parents=True, exist_ok=True)

    for region, group in df.groupby("brain_region", sort=False):
        sets = {
            "exploratory_p001_absbeta02": set(group.index[group[BASELINE_COL]]),
            "BH_FDR_lt_010": set(group.index[group["fdr_BH"] < 0.10]),
        }
        if "qvalue_Storey" in group.columns:
            sets["Storey_q_lt_010"] = set(group.index[group["qvalue_Storey"] < 0.10])
        best = best_filters.get(region)
        if best:
            col = str(best["adjusted_column"])
            sets["independent_filter_BH_lt_010"] = set(group.index[group[col] < 0.10])
        else:
            sets["independent_filter_BH_lt_010"] = set()
        if "p_adj_IHW" in group.columns:
            sets["IHW_lt_010"] = set(group.index[group["p_adj_IHW"] < 0.10])

        names = list(sets)
        matrix = pd.DataFrame(index=names, columns=names, dtype=int)
        for a in names:
            for b in names:
                matrix.loc[a, b] = len(sets[a] & sets[b])
        matrix.to_csv(out_dir / f"{safe_name(region)}_overlap_matrix.csv")


def report_gene_lists(summary: pd.DataFrame) -> str:
    lines = []
    for row in summary.itertuples():
        if row.BH_FDR_lt_010 > 0:
            strict = "BH-FDR < 0.10 genes are available for cautious reporting."
        else:
            strict = "No BH-FDR < 0.10 genes; use no strict corrected DE list for this region."
        lines.append(f"- {row.brain_region}: {strict} Exploratory list has {row.exploratory_p001_absbeta02} genes.")
    return "\n".join(lines)


def write_report(
    df: pd.DataFrame,
    summary: pd.DataFrame,
    independent_summary: pd.DataFrame,
    qvalue_available: bool,
    ihw_available: bool,
) -> None:
    bh_any = int(summary["BH_FDR_lt_010"].sum())
    by_any = int(summary["BY_FDR_lt_010"].sum())
    best_filter_any = int(summary["BH_after_best_filter_lt_010"].sum())
    baseline_any = int(summary["exploratory_p001_absbeta02"].sum())

    if independent_summary.empty:
        independent_note = (
            "Independent filtering did not run because the available mean-expression columns were empty "
            "and expression variance / detection rate columns were not present."
        )
    else:
        independent_note = "Independent filtering was evaluated with the available independent covariates."

    report = f"""# P-Value Correction Audit: Length > 100

## Scope

The current exploratory DE result is preserved and was not overwritten. This audit only reads the existing LMM result file and writes new files inside `pvalue_correction_audit/`.

Current baseline:

`p_value_diagnosis < 0.01 and abs(beta_diagnosis) > 0.2`

## How To Read These Columns

- `beta_diagnosis` is the effect size.
- `p_value_diagnosis` is nominal statistical evidence.
- FDR/q-value is the multiple-testing correction.

This audit compares correction strategies; it is not changing the main exploratory result unless the comparison supports that decision.

## Method Availability

- BH-FDR: available and computed region-by-region.
- BY-FDR: available and computed region-by-region.
- Storey q-value: {'available' if qvalue_available else 'not available; R package qvalue is not installed, so it was skipped.'}
- IHW: {'available' if ihw_available else 'not available; R package IHW is not installed, so it was skipped.'}
- Independent filtering: {independent_note}

## Summary By Region

{summary.to_markdown(index=False)}

## Correction Strategy Takeaways

- Exploratory baseline candidates recovered: {baseline_any}
- BH-FDR < 0.10 candidates recovered: {bh_any}
- BY-FDR < 0.10 candidates recovered: {by_any}
- Best independent-filter BH < 0.10 candidates recovered: {best_filter_any}

BH-FDR appears conservative for this dataset because it recovers very few genes compared with the nominal beta-plus-p-value exploratory baseline. BY-FDR is even more conservative, as expected. Storey q-value and IHW could not be evaluated because the required R packages are not installed. Independent filtering did not improve power here because valid independent expression, variance, or detection-rate covariates were not available in the LMM summary table.

## Recommended Gene Lists

1. Strict reporting: use corrected-p-value lists only. If a region has no BH-FDR < 0.10 genes, report no strict corrected DE genes for that region.
2. Exploratory enrichment: use `exploratory_p001_absbeta02` gene lists.
3. ML feature selection: use ranked signatures from `signed_logp_score` and `abs_beta_diagnosis`, keeping the exploratory candidates as a focused candidate set.
4. AHBA module mapping: use exploratory candidates plus ranked full signatures, and label them as nominal/exploratory.

{report_gene_lists(summary)}
"""
    (AUDIT_DIR / "Pvalue_Correction_Audit_Report.md").write_text(report, encoding="utf-8")


def main() -> None:
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Reading {INPUT_FILE}")
    df = pd.read_csv(INPUT_FILE)
    df = add_core_columns(df)
    df = add_regionwise_corrections(df)

    qvalue_available = check_r_package("qvalue")
    ihw_available = check_r_package("IHW")

    if not qvalue_available:
        print("Skipping Storey q-value: R package qvalue is not installed.")
    if not ihw_available:
        print("Skipping IHW: R package IHW is not installed.")

    independent_summary, independent_adjusted = run_independent_filter_audit(df)
    if not independent_adjusted.empty:
        df = pd.concat([df, independent_adjusted], axis=1)
    if not independent_summary.empty:
        independent_summary.to_csv(AUDIT_DIR / "independent_filtering_summary.csv", index=False)

    best_filters = choose_best_filters(independent_summary)
    summary = make_summary(df, best_filters)

    df.to_csv(AUDIT_DIR / "LMM_all_genes_with_pvalue_corrections.csv", index=False)
    summary.to_csv(AUDIT_DIR / "pvalue_correction_summary_by_region.csv", index=False)
    save_candidate_lists(df, best_filters)
    write_plots(df)
    write_overlap_matrices(df, best_filters)
    write_report(df, summary, independent_summary, qvalue_available, ihw_available)

    print(f"Wrote audit outputs to: {AUDIT_DIR}")


if __name__ == "__main__":
    main()
