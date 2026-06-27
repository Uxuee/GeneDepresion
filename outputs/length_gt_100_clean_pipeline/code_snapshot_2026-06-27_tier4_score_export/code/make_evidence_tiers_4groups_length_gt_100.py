"""Create evidence-tier summaries for the Length > 100 four-group MDD analysis.

This script reads existing Length > 100 outputs and writes only to:

  outputs/length_gt_100_elisa_style/evidence_tiers_4groups/

It does not overwrite the preserved exploratory DE result or the
standardized-effect-size branch.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


REPO_DIR = Path(__file__).resolve().parents[1]
BASE_DIR = REPO_DIR / "outputs" / "length_gt_100_elisa_style"
BASELINE_FILE = BASE_DIR / "LMM_exploratory_DE_p001_absbeta02_length_gt_100.csv"
STD_FILE = (
    BASE_DIR
    / "standardized_effect_size_4groups"
    / "LMM_all_genes_with_standardized_effect_size_4groups.csv"
)
STD_SUMMARY_FILE = (
    BASE_DIR
    / "standardized_effect_size_4groups"
    / "standardized_effect_size_4groups_summary.csv"
)
OUT_DIR = BASE_DIR / "evidence_tiers_4groups"

REGIONS = [
    "Frontal_Cortex",
    "Insula",
    "Ventral_Striatum",
    "Hippocampal_Formation",
]

DISPLAY = {
    "Frontal_Cortex": "Frontal Cortex",
    "Insula": "Insula",
    "Ventral_Striatum": "Ventral Striatum",
    "Hippocampal_Formation": "Hippocampal Formation",
}


def gene_key_frame(df: pd.DataFrame) -> pd.DataFrame:
    return df[["grouped_region", "gene_id", "gene_symbol"]].drop_duplicates()


def sort_gene_table(df: pd.DataFrame) -> pd.DataFrame:
    sort_cols = [col for col in ["grouped_region", "p_value_diagnosis", "abs_standardized_effect_size_diagnosis"] if col in df.columns]
    ascending = [True, True, False][: len(sort_cols)]
    return df.sort_values(sort_cols, ascending=ascending, na_position="last")


def selected_columns(df: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "gene_id",
        "gene_symbol",
        "grouped_region",
        "beta_diagnosis",
        "standardized_effect_size_diagnosis",
        "abs_standardized_effect_size_diagnosis",
        "p_value_diagnosis",
        "qvalue_Storey",
        "fdr_BH",
        "direction_standardized_effect_size",
        "direction_beta",
        "signed_standardized_logp_score",
        "signed_logp_score",
    ]
    return df[[col for col in cols if col in df.columns]]


def overlap_table(region: str, left: pd.DataFrame, right: pd.DataFrame, left_name: str, right_name: str) -> pd.DataFrame:
    left_ids = set(left["gene_id"].astype(str))
    right_ids = set(right["gene_id"].astype(str))
    all_ids = sorted(left_ids | right_ids)
    lookup = (
        pd.concat([gene_key_frame(left), gene_key_frame(right)], ignore_index=True)
        .drop_duplicates(["grouped_region", "gene_id"])
        .set_index("gene_id")
    )
    rows = []
    for gene_id in all_ids:
        rows.append(
            {
                "grouped_region": region,
                "gene_id": gene_id,
                "gene_symbol": lookup.loc[gene_id, "gene_symbol"] if gene_id in lookup.index else "",
                left_name: gene_id in left_ids,
                right_name: gene_id in right_ids,
                "overlap": gene_id in left_ids and gene_id in right_ids,
            }
        )
    return pd.DataFrame(rows)


def write_region_gene_list_report(title: str, df: pd.DataFrame) -> str:
    lines = [f"## {title}", ""]
    if df.empty:
        lines.append("No genes.")
        lines.append("")
        return "\n".join(lines)
    for region, group in df.groupby("grouped_region", sort=False):
        lines.append(f"### {DISPLAY.get(region, region)}")
        lines.append("")
        lines.append("| gene_id | gene_symbol | beta | std_effect_size | p_value | qvalue_Storey | fdr_BH | direction |")
        lines.append("|---|---|---:|---:|---:|---:|---:|---|")
        for row in group.itertuples():
            lines.append(
                f"| {row.gene_id} | {row.gene_symbol} | {row.beta_diagnosis:.4g} | "
                f"{row.standardized_effect_size_diagnosis:.4g} | {row.p_value_diagnosis:.4g} | "
                f"{getattr(row, 'qvalue_Storey', np.nan):.4g} | {getattr(row, 'fdr_BH', np.nan):.4g} | "
                f"{row.direction_standardized_effect_size} |"
            )
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    baseline = pd.read_csv(BASELINE_FILE, dtype={"gene_id": str})
    std = pd.read_csv(STD_FILE, dtype={"gene_id": str}, low_memory=False)
    std_summary = pd.read_csv(STD_SUMMARY_FILE)

    std["gene_id"] = std["gene_id"].astype(str)
    std["grouped_region"] = std["grouped_region"].fillna(std["brain_region"])
    observed = set(std["grouped_region"].dropna().unique())
    expected = set(REGIONS)
    if observed != expected:
        raise ValueError(f"Unexpected grouped regions. Observed={sorted(observed)}, expected={sorted(expected)}")

    std["tier1_corrected_standardized"] = (
        pd.to_numeric(std["qvalue_Storey"], errors="coerce") < 0.10
    ) & (std["abs_standardized_effect_size_diagnosis"] > 0.8)
    std["tier2_nominal_standardized"] = (
        std["p_value_diagnosis"] < 0.01
    ) & (std["abs_standardized_effect_size_diagnosis"] > 0.8)
    std["tier3_current_baseline"] = (
        std["p_value_diagnosis"] < 0.01
    ) & (std["abs_beta_diagnosis"] > 0.2)

    tier1 = sort_gene_table(std[std["tier1_corrected_standardized"]].copy())
    tier2 = sort_gene_table(std[std["tier2_nominal_standardized"]].copy())
    tier3 = sort_gene_table(std[std["tier3_current_baseline"]].copy())

    selected_columns(tier1).to_csv(OUT_DIR / "Tier1_corrected_standardized_genes.csv", index=False)
    selected_columns(tier2).to_csv(OUT_DIR / "Tier2_nominal_standardized_genes.csv", index=False)
    selected_columns(tier3).to_csv(OUT_DIR / "Tier3_current_baseline_genes.csv", index=False)

    summary_rows = []
    overlap_23 = []
    overlap_12 = []
    overlap_13 = []
    tier2_unique = []
    tier3_unique = []

    for region in REGIONS:
        r1 = tier1[tier1["grouped_region"] == region]
        r2 = tier2[tier2["grouped_region"] == region]
        r3 = tier3[tier3["grouped_region"] == region]
        ids1 = set(r1["gene_id"])
        ids2 = set(r2["gene_id"])
        ids3 = set(r3["gene_id"])

        unique2_ids = ids2 - ids3
        unique3_ids = ids3 - ids2
        tier2_unique.append(selected_columns(r2[r2["gene_id"].isin(unique2_ids)]))
        tier3_unique.append(selected_columns(r3[r3["gene_id"].isin(unique3_ids)]))

        summary_rows.append(
            {
                "grouped_region": region,
                "display_region": DISPLAY[region],
                "tier1_corrected_standardized_genes": len(ids1),
                "tier2_nominal_standardized_genes": len(ids2),
                "tier3_current_baseline_genes": len(ids3),
                "tier2_tier3_overlap": len(ids2 & ids3),
                "tier2_unique_genes": len(unique2_ids),
                "tier3_unique_genes": len(unique3_ids),
            }
        )
        overlap_23.append(overlap_table(region, r2, r3, "in_tier2_nominal_standardized", "in_tier3_current_baseline"))
        overlap_12.append(overlap_table(region, r1, r2, "in_tier1_corrected_standardized", "in_tier2_nominal_standardized"))
        overlap_13.append(overlap_table(region, r1, r3, "in_tier1_corrected_standardized", "in_tier3_current_baseline"))

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_DIR / "evidence_tier_summary_by_region.csv", index=False)
    pd.concat(overlap_23, ignore_index=True).to_csv(OUT_DIR / "Tier2_Tier3_overlap_by_region.csv", index=False)
    pd.concat(overlap_12, ignore_index=True).to_csv(OUT_DIR / "Tier1_Tier2_overlap_by_region.csv", index=False)
    pd.concat(overlap_13, ignore_index=True).to_csv(OUT_DIR / "Tier1_Tier3_overlap_by_region.csv", index=False)
    pd.concat(tier2_unique, ignore_index=True).to_csv(OUT_DIR / "Tier2_unique_genes_by_region.csv", index=False)
    pd.concat(tier3_unique, ignore_index=True).to_csv(OUT_DIR / "Tier3_unique_genes_by_region.csv", index=False)

    ranked_cols = [
        "gene_id",
        "gene_symbol",
        "grouped_region",
        "beta_diagnosis",
        "standardized_effect_size_diagnosis",
        "p_value_diagnosis",
        "signed_standardized_logp_score",
        "signed_logp_score",
    ]
    std.sort_values(
        ["grouped_region", "signed_standardized_logp_score"],
        ascending=[True, False],
        na_position="last",
    )[ranked_cols].to_csv(OUT_DIR / "Tier4_ranked_signature_by_signed_standardized_logp_score.csv", index=False)
    std.sort_values(
        ["grouped_region", "signed_logp_score"],
        ascending=[True, False],
        na_position="last",
    )[ranked_cols].to_csv(OUT_DIR / "Tier4_ranked_signature_by_signed_logp_score.csv", index=False)

    sign_check_file = BASE_DIR / "standardized_effect_size_4groups" / "standardized_effect_size_sign_check_summary.csv"
    sign_check = pd.read_csv(sign_check_file) if sign_check_file.exists() else pd.DataFrame()

    ventral_tier1 = tier1[tier1["grouped_region"] == "Ventral_Striatum"]
    report = f"""# Evidence Tiers: Length > 100 Four-Group MDD Analysis

## Scope

The current baseline was preserved. This evidence-tier summary is a derived comparison and does not overwrite previous exploratory DE, p-value audit, or standardized-effect-size outputs.

The standardized-effect-size analysis was added as an Elisa-style branch using:

`standardized_effect_size_diagnosis = beta_diagnosis / sd_expression_gene`

## Tier Definitions

- Tier 1: corrected standardized-effect-size genes, `Storey q < 0.10` and `abs_standardized_effect_size_diagnosis > 0.8`.
- Tier 2: nominal Elisa-style standardized-effect-size genes, `p_value_diagnosis < 0.01` and `abs_standardized_effect_size_diagnosis > 0.8`.
- Tier 3: current baseline raw-beta genes, `p_value_diagnosis < 0.01` and `abs(beta_diagnosis) > 0.2`.
- Tier 4: all genes ranked by `signed_standardized_logp_score` and `signed_logp_score`.

Do not call Tier 2 or Tier 3 final significant DE genes. They are exploratory, nominal, candidate signature, and ranked signature resources.

## Sign Check

{sign_check.to_markdown(index=False) if not sign_check.empty else 'Sign-check file was not found.'}

## Evidence Tier Summary

{summary.to_markdown(index=False)}

## Tier 1 Genes

{write_region_gene_list_report('All Tier 1 corrected standardized-effect-size genes', selected_columns(tier1))}

## Ventral Striatum Tier 1 Genes

{selected_columns(ventral_tier1).to_markdown(index=False) if not ventral_tier1.empty else 'No Ventral Striatum Tier 1 genes.'}

## Interpretation

Ventral Striatum is the most robust region because it has high counts in both exploratory definitions and the only Storey q < 0.10 standardized-effect-size corrected genes. Frontal Cortex is sensitive to the effect-size definition: it has many current-baseline raw-beta candidates but only one nominal standardized-effect-size candidate. Hippocampal Formation gains genes under standardized effect size, suggesting raw beta may undercount genes with consistent but smaller absolute changes.

## Recommended Use

- Tier 1: strict reporting / high-confidence anchors.
- Tier 2: exploratory enrichment.
- Tier 3: comparison with the previous baseline.
- Tier 4: ML and AHBA ranked-signature input.

## Standardized Effect Summary Reference

{std_summary.to_markdown(index=False)}
"""
    (OUT_DIR / "Evidence_Tiers_4Groups_Report.md").write_text(report, encoding="utf-8")
    print(f"Wrote evidence-tier outputs to: {OUT_DIR}")
    print(summary.to_string(index=False))
    print("\nTier 1 genes:")
    print(selected_columns(tier1).to_string(index=False))


if __name__ == "__main__":
    main()
