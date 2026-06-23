"""Create Elisa-style exploratory DE tables for Length > 100 results.

Input:
  outputs/length_gt_100_elisa_style/LMM_all_genes_length_gt_100.csv

Column choices used by this script:
  - diagnosis beta/effect column: beta_diagnosis
  - diagnosis p-value column: p_value_diagnosis
  - optional FDR reference column: fdr_BH

Exploratory DE definition:
  p_value_diagnosis < 0.01 and abs(beta_diagnosis) > 0.2

These outputs are nominal/exploratory only. FDR is preserved as a reference
column when present, but it is not used as the main filter.
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
OUTPUT_DIR = REPO_DIR / "outputs" / "length_gt_100_elisa_style"
INPUT_FILE = OUTPUT_DIR / "LMM_all_genes_length_gt_100.csv"

P_VALUE_CUTOFF = 0.01
ABS_BETA_CUTOFF = 0.2

BETA_CANDIDATES = ["beta_diagnosis", "beta_Diagnosis", "estimate", "logFC"]
P_VALUE_CANDIDATES = ["p_value_diagnosis", "p_value_Diagnosis", "p.value", "P.Value"]
FDR_CANDIDATES = ["fdr_BH", "FDR", "adj.P.Val", "diagnosis_FDR_BH"]


def find_column(columns: list[str], candidates: list[str], label: str) -> str:
    """Return the first matching column name from a list of acceptable names."""
    for candidate in candidates:
        if candidate in columns:
            return candidate
    raise ValueError(f"Could not find {label} column. Tried: {', '.join(candidates)}")


def clean_for_filename(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "_-" else "_" for ch in str(value))


def add_exploratory_columns(df: pd.DataFrame, beta_col: str, p_col: str) -> pd.DataFrame:
    """Add beta magnitude, signed score, direction, and pass/fail columns."""
    out = df.copy()
    out[beta_col] = pd.to_numeric(out[beta_col], errors="coerce")
    out[p_col] = pd.to_numeric(out[p_col], errors="coerce")

    out["abs_beta_diagnosis"] = out[beta_col].abs()
    clipped_p = out[p_col].clip(lower=np.finfo(float).tiny)
    out["signed_logp_score"] = np.sign(out[beta_col]) * -np.log10(clipped_p)
    out["direction"] = np.select(
        [out[beta_col] > 0, out[beta_col] < 0, out[beta_col] == 0],
        ["up_in_MDD", "down_in_MDD", "no_change"],
        default="no_change",
    )
    out["exploratory_DE_p001_absbeta02"] = (
        (out[p_col] < P_VALUE_CUTOFF) & (out["abs_beta_diagnosis"] > ABS_BETA_CUTOFF)
    )
    return out


def top100_by_region(df: pd.DataFrame, sort_col: str, ascending: bool, out_dir: Path, suffix: str) -> pd.DataFrame:
    """Save top 100 exploratory rows per region and return the combined table."""
    all_top = []
    for region, group in df.groupby("brain_region", sort=False):
        top = group.sort_values(sort_col, ascending=ascending, na_position="last").head(100)
        top.to_csv(out_dir / f"{clean_for_filename(region)}_{suffix}.csv", index=False)
        all_top.append(top)
    return pd.concat(all_top, ignore_index=True) if all_top else df.head(0)


def write_volcano_plots(df: pd.DataFrame, beta_col: str, p_col: str, out_dir: Path) -> None:
    """Make one volcano plot per brain region and highlight exploratory candidates."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for region, group in df.groupby("brain_region", sort=False):
        plot_data = group.dropna(subset=[beta_col, p_col]).copy()
        if plot_data.empty:
            continue

        plot_data["minus_log10_p"] = -np.log10(plot_data[p_col].clip(lower=np.finfo(float).tiny))
        is_candidate = plot_data["exploratory_DE_p001_absbeta02"]

        plt.figure(figsize=(7, 6))
        plt.scatter(
            plot_data.loc[~is_candidate, beta_col],
            plot_data.loc[~is_candidate, "minus_log10_p"],
            s=7,
            c="#747b85",
            alpha=0.45,
            label="Not passing exploratory filter",
        )
        plt.scatter(
            plot_data.loc[is_candidate, beta_col],
            plot_data.loc[is_candidate, "minus_log10_p"],
            s=12,
            c="#b42318",
            alpha=0.8,
            label="p < 0.01 and |beta| > 0.2",
        )
        plt.axhline(-math.log10(P_VALUE_CUTOFF), color="#333333", linestyle="--", linewidth=1)
        plt.axvline(-ABS_BETA_CUTOFF, color="#333333", linestyle="--", linewidth=1)
        plt.axvline(ABS_BETA_CUTOFF, color="#333333", linestyle="--", linewidth=1)
        plt.xlabel("beta_diagnosis")
        plt.ylabel("-log10(p_value_diagnosis)")
        plt.title(f"Exploratory DE volcano: {region}")
        plt.legend(frameon=False, fontsize=8)
        plt.tight_layout()
        plt.savefig(out_dir / f"{clean_for_filename(region)}_volcano_exploratory_p001_absbeta02.png", dpi=180)
        plt.close()


def make_summary(df: pd.DataFrame, candidates: pd.DataFrame, fdr_col: str | None) -> pd.DataFrame:
    rows = []
    for region, all_region in df.groupby("brain_region", sort=False):
        region_candidates = candidates[candidates["brain_region"] == region]
        row = {
            "brain_region": region,
            "genes_tested": len(all_region),
            "exploratory_DE_genes": len(region_candidates),
            "up_in_MDD": int((region_candidates["direction"] == "up_in_MDD").sum()),
            "down_in_MDD": int((region_candidates["direction"] == "down_in_MDD").sum()),
            "median_abs_beta": region_candidates["abs_beta_diagnosis"].median(),
            "max_abs_beta": region_candidates["abs_beta_diagnosis"].max(),
            "min_p_value": region_candidates["p_value_diagnosis"].min(),
        }
        if fdr_col is not None:
            row["number_with_FDR_lt_005"] = int((region_candidates[fdr_col] < 0.05).sum())
            row["number_with_FDR_lt_010"] = int((region_candidates[fdr_col] < 0.10).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def format_gene_label(row: pd.Series) -> str:
    symbol = row.get("gene_symbol")
    if isinstance(symbol, str) and symbol.strip():
        return symbol
    return str(row.get("gene_id", "unknown"))


def markdown_top20(candidates: pd.DataFrame, sort_col: str, ascending: bool) -> str:
    lines = []
    for region, group in candidates.groupby("brain_region", sort=False):
        lines.append(f"### {region}")
        top = group.sort_values(sort_col, ascending=ascending, na_position="last").head(20)
        if top.empty:
            lines.append("")
            lines.append("No exploratory DE candidates passed this filter.")
            lines.append("")
            continue
        lines.append("")
        lines.append("| rank | gene | beta_diagnosis | p_value_diagnosis | abs_beta_diagnosis | direction |")
        lines.append("|---:|---|---:|---:|---:|---|")
        for rank, (_, row) in enumerate(top.iterrows(), start=1):
            lines.append(
                f"| {rank} | {format_gene_label(row)} | {row['beta_diagnosis']:.4g} | "
                f"{row['p_value_diagnosis']:.4g} | {row['abs_beta_diagnosis']:.4g} | {row['direction']} |"
            )
        lines.append("")
    return "\n".join(lines)


def write_report(summary: pd.DataFrame, candidates: pd.DataFrame, beta_col: str, p_col: str, fdr_col: str | None) -> None:
    fdr_note = f"`{fdr_col}` was retained as a reference column." if fdr_col else "No FDR column was found in the input file."
    report = f"""# Exploratory DE Filter: p < 0.01 and |beta| > 0.2

## Definition

Exploratory DE candidates are genes with:

`p_value_diagnosis < 0.01` and `abs(beta_diagnosis) > 0.2`

The script used `{beta_col}` as the diagnosis beta/effect-size column and `{p_col}` as the diagnosis p-value column. {fdr_note}

## Important Interpretation

These are nominal/exploratory DE candidates, not final FDR-significant DE genes.

FDR is not the effect size. `beta_diagnosis` is the effect size, `p_value_diagnosis` is the nominal statistical evidence, and FDR is the multiple-testing-corrected p-value.

## Summary By Region

{summary.to_markdown(index=False)}

## Top 20 By P-Value

{markdown_top20(candidates, "p_value_diagnosis", True)}

## Top 20 By Absolute Beta

{markdown_top20(candidates, "abs_beta_diagnosis", False)}

## Recommended Next Step

Use these exploratory DE candidates and ranked signatures for enrichment analysis, ML feature selection, and AHBA module mapping.
"""
    (OUTPUT_DIR / "Exploratory_DE_p001_absbeta02_report.md").write_text(report, encoding="utf-8")


def main() -> None:
    print(f"Reading {INPUT_FILE}")
    df = pd.read_csv(INPUT_FILE)
    beta_col = find_column(df.columns.tolist(), BETA_CANDIDATES, "diagnosis beta")
    p_col = find_column(df.columns.tolist(), P_VALUE_CANDIDATES, "diagnosis p-value")
    fdr_col = next((col for col in FDR_CANDIDATES if col in df.columns), None)

    if beta_col != "beta_diagnosis":
        df["beta_diagnosis"] = df[beta_col]
        beta_col = "beta_diagnosis"
    if p_col != "p_value_diagnosis":
        df["p_value_diagnosis"] = df[p_col]
        p_col = "p_value_diagnosis"

    df = add_exploratory_columns(df, beta_col, p_col)
    if fdr_col is not None:
        df[fdr_col] = pd.to_numeric(df[fdr_col], errors="coerce")

    candidates = df[df["exploratory_DE_p001_absbeta02"]].copy()
    candidates = candidates.sort_values(["brain_region", "p_value_diagnosis"], na_position="last")

    by_region_dir = OUTPUT_DIR / "exploratory_DE_by_region"
    top_by_region_dir = OUTPUT_DIR / "top100_exploratory_DE_by_region"
    volcano_dir = OUTPUT_DIR / "volcano_exploratory_p001_absbeta02"
    by_region_dir.mkdir(parents=True, exist_ok=True)
    top_by_region_dir.mkdir(parents=True, exist_ok=True)

    candidates.to_csv(
        OUTPUT_DIR / "LMM_exploratory_DE_p001_absbeta02_length_gt_100.csv",
        index=False,
    )

    for region, group in candidates.groupby("brain_region", sort=False):
        group.to_csv(
            by_region_dir / f"{clean_for_filename(region)}_exploratory_DE_p001_absbeta02.csv",
            index=False,
        )

    summary = make_summary(df, candidates, fdr_col)
    summary.to_csv(
        OUTPUT_DIR / "summary_exploratory_DE_p001_absbeta02_by_region.csv",
        index=False,
    )

    combined_top_p = top100_by_region(
        candidates,
        "p_value_diagnosis",
        True,
        top_by_region_dir,
        "top100_exploratory_by_pvalue_p001_absbeta02",
    )
    combined_top_beta = top100_by_region(
        candidates,
        "abs_beta_diagnosis",
        False,
        top_by_region_dir,
        "top100_exploratory_by_abs_beta_p001_absbeta02",
    )
    combined_top_score = top100_by_region(
        candidates,
        "signed_logp_score",
        False,
        top_by_region_dir,
        "top100_exploratory_by_signed_logp_score_p001_absbeta02",
    )
    combined_top_p.to_csv(
        OUTPUT_DIR / "top100_exploratory_DE_by_pvalue_p001_absbeta02.csv",
        index=False,
    )
    combined_top_beta.to_csv(
        OUTPUT_DIR / "top100_exploratory_DE_by_abs_beta_p001_absbeta02.csv",
        index=False,
    )
    combined_top_score.to_csv(
        OUTPUT_DIR / "top100_exploratory_DE_by_signed_logp_score_p001_absbeta02.csv",
        index=False,
    )

    write_volcano_plots(df, beta_col, p_col, volcano_dir)
    write_report(summary, candidates, beta_col, p_col, fdr_col)
    print(f"Exploratory candidates written: {len(candidates)}")
    print(f"Outputs saved under: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
