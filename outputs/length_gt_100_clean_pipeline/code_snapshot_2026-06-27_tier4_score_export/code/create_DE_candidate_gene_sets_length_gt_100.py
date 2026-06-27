"""Create clean DE candidate gene-set package for Length > 100 MDD analysis.

This script does not rerun any LMMs. It packages existing Length > 100
standardized-effect-size, evidence-tier, p-value audit, and frontal subregion
outputs into practical gene sets for enrichment, ML preparation, and AHBA
module mapping.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


REPO_DIR = Path(__file__).resolve().parents[1]
ROOT = REPO_DIR / "outputs" / "length_gt_100_elisa_style"
OUT = ROOT / "DE_candidate_gene_sets"
TIER2_DIR = OUT / "Tier2_working_exploratory"
TIER3_DIR = OUT / "Tier3_raw_beta_comparison"
RANKED_DIR = OUT / "ranked_signatures"
ENRICHMENT_DIR = OUT / "enrichment_ready"

FOUR_GROUP_FILE = (
    ROOT
    / "standardized_effect_size_4groups"
    / "LMM_all_genes_with_standardized_effect_size_4groups.csv"
)
TIER1_FILE = ROOT / "evidence_tiers_4groups" / "Tier1_corrected_standardized_genes.csv"
FRONTAL_SUBREGION_FILE = (
    ROOT
    / "frontal_subregion_sensitivity"
    / "frontal_subregion_LMM_all_genes_with_standardized_effect_size.csv"
)
FRONTAL_SUBREGION_SUMMARY_FILE = (
    ROOT / "frontal_subregion_sensitivity" / "frontal_subregion_summary.csv"
)
PVALUE_AUDIT_REPORT = (
    ROOT / "pvalue_correction_sensitivity" / "Pvalue_Correction_Sensitivity_Report.md"
)


FOUR_GROUP_REGIONS = {
    "Frontal_Cortex": "broad Frontal Cortex",
    "Insula": "Insula",
    "Ventral_Striatum": "Ventral Striatum",
    "Hippocampal_Formation": "Hippocampal Formation",
}

FRONTAL_SUBREGION_NAMES = {
    "Cingulate gyrus 25 (Cg25)": "Cingulate gyrus 25",
    "Dorsolateral prefrontal cortex (dlPFC; BA8/9)": "dlPFC BA8/9",
    "Orbitofrontal (OFC; BA11)": "OFC BA11",
}

RECOMMENDED_USE = {
    "Ventral Striatum": "main robust broad-region set",
    "Insula": "exploratory/supportive broad-region set",
    "Hippocampal Formation": "exploratory/supportive broad-region set",
    "broad Frontal Cortex": "record only; not recommended as primary frontal set because Tier 2 = 1",
    "Cingulate gyrus 25": "main frontal exploratory set",
    "dlPFC BA8/9": "frontal supportive set",
    "OFC BA11": "frontal supportive set",
}


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required input file not found: {path}")


def safe_name(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9]+", "_", str(value))
    return re.sub(r"^_+|_+$", "", value)


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.astype(str).str.lower().str.strip().isin(["true", "t", "1", "yes"])


def direction_from_standardized_effect(df: pd.DataFrame) -> pd.Series:
    effect = pd.to_numeric(df["standardized_effect_size_diagnosis"], errors="coerce")
    return np.select(
        [effect > 0, effect < 0],
        ["up_in_MDD", "down_in_MDD"],
        default="no_change",
    )


def add_scores(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    p = pd.to_numeric(out["p_value_diagnosis"], errors="coerce").clip(lower=np.finfo(float).tiny)
    effect = pd.to_numeric(out["standardized_effect_size_diagnosis"], errors="coerce")
    out["signed_standardized_logp_score"] = np.sign(effect) * (-np.log10(p))
    out["abs_standardized_logp_score"] = out["signed_standardized_logp_score"].abs()
    return out


def clean_gene_table(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "gene_id",
        "gene_symbol",
        "analysis_level",
        "region_or_subregion",
        "beta_diagnosis",
        "standardized_effect_size_diagnosis",
        "abs_standardized_effect_size_diagnosis",
        "p_value_diagnosis",
        "qvalue_Storey",
        "fdr_BH",
        "direction",
        "recommended_use",
    ]
    return df[[col for col in keep if col in df.columns]].copy()


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def write_gene_symbols(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    symbols = (
        df["gene_symbol"]
        .dropna()
        .astype(str)
        .loc[lambda x: x.str.len() > 0]
        .drop_duplicates()
        .sort_values()
    )
    path.write_text("\n".join(symbols.to_list()) + ("\n" if len(symbols) else ""), encoding="utf-8")


def make_tier2_sets(df: pd.DataFrame, label: str) -> dict[str, pd.DataFrame]:
    tier2 = df[df["Tier2_flag"]].copy()
    return {
        "Tier2_all": tier2,
        "Tier2_up_in_MDD": tier2[tier2["direction"] == "up_in_MDD"].copy(),
        "Tier2_down_in_MDD": tier2[tier2["direction"] == "down_in_MDD"].copy(),
    }


def make_tier3_set(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["Tier3_flag"]].copy()


def prepare_four_group(four_group: pd.DataFrame) -> pd.DataFrame:
    df = four_group.copy()
    df["group_key"] = df["grouped_region"].where(df["grouped_region"].notna(), df["brain_region"])
    df = df[df["group_key"].isin(FOUR_GROUP_REGIONS.keys())].copy()
    df["analysis_level"] = "four_group"
    df["region_or_subregion"] = df["group_key"].map(FOUR_GROUP_REGIONS)
    df["direction"] = direction_from_standardized_effect(df)
    df["Tier2_flag"] = (
        (pd.to_numeric(df["p_value_diagnosis"], errors="coerce") < 0.01)
        & (pd.to_numeric(df["abs_standardized_effect_size_diagnosis"], errors="coerce") > 0.8)
    )
    if "current_baseline_p001_absBeta02" in df.columns:
        df["Tier3_flag"] = as_bool(df["current_baseline_p001_absBeta02"])
    else:
        df["Tier3_flag"] = (
            (pd.to_numeric(df["p_value_diagnosis"], errors="coerce") < 0.01)
            & (pd.to_numeric(df["beta_diagnosis"], errors="coerce").abs() > 0.2)
        )
    df["recommended_use"] = df["region_or_subregion"].map(RECOMMENDED_USE)
    return add_scores(df)


def prepare_frontal_subregions(frontal: pd.DataFrame) -> pd.DataFrame:
    df = frontal.copy()
    df = df[df["original_frontal_subregion"].isin(FRONTAL_SUBREGION_NAMES.keys())].copy()
    df["analysis_level"] = "frontal_subregion"
    df["region_or_subregion"] = df["original_frontal_subregion"].map(FRONTAL_SUBREGION_NAMES)
    df["direction"] = direction_from_standardized_effect(df)
    if "Tier2_nominal_p001_absStdES08" in df.columns:
        df["Tier2_flag"] = as_bool(df["Tier2_nominal_p001_absStdES08"])
    else:
        df["Tier2_flag"] = (
            (pd.to_numeric(df["p_value_diagnosis"], errors="coerce") < 0.01)
            & (pd.to_numeric(df["abs_standardized_effect_size_diagnosis"], errors="coerce") > 0.8)
        )
    if "Tier3_current_baseline_p001_absBeta02" in df.columns:
        df["Tier3_flag"] = as_bool(df["Tier3_current_baseline_p001_absBeta02"])
    else:
        df["Tier3_flag"] = (
            (pd.to_numeric(df["p_value_diagnosis"], errors="coerce") < 0.01)
            & (pd.to_numeric(df["beta_diagnosis"], errors="coerce").abs() > 0.2)
        )
    df["recommended_use"] = df["region_or_subregion"].map(RECOMMENDED_USE)
    return add_scores(df)


def prepare_tier1(tier1: pd.DataFrame) -> pd.DataFrame:
    df = tier1.copy()
    df["analysis_level"] = "four_group"
    df["region_or_subregion"] = df["grouped_region"].map(FOUR_GROUP_REGIONS)
    df["direction"] = direction_from_standardized_effect(df)
    df["recommended_use"] = df["region_or_subregion"].map(RECOMMENDED_USE)
    df = add_scores(df)
    return clean_gene_table(df)


def write_ranked_signature(df: pd.DataFrame, label: str) -> None:
    ranked = df.copy()
    ranked = ranked.sort_values(
        ["abs_standardized_logp_score", "p_value_diagnosis"],
        ascending=[False, True],
        na_position="last",
    )
    keep = [
        "gene_id",
        "gene_symbol",
        "analysis_level",
        "region_or_subregion",
        "beta_diagnosis",
        "standardized_effect_size_diagnosis",
        "abs_standardized_effect_size_diagnosis",
        "p_value_diagnosis",
        "signed_standardized_logp_score",
        "abs_standardized_logp_score",
        "direction",
        "recommended_use",
    ]
    write_csv(ranked[[col for col in keep if col in ranked.columns]], RANKED_DIR / f"{safe_name(label)}__Tier4_ranked_signature.csv")


def main() -> None:
    for path in [
        FOUR_GROUP_FILE,
        TIER1_FILE,
        FRONTAL_SUBREGION_FILE,
        FRONTAL_SUBREGION_SUMMARY_FILE,
        PVALUE_AUDIT_REPORT,
    ]:
        require_file(path)

    for path in [OUT, TIER2_DIR, TIER3_DIR, RANKED_DIR, ENRICHMENT_DIR]:
        path.mkdir(parents=True, exist_ok=True)

    four_group = prepare_four_group(pd.read_csv(FOUR_GROUP_FILE))
    frontal = prepare_frontal_subregions(pd.read_csv(FRONTAL_SUBREGION_FILE))
    tier1 = prepare_tier1(pd.read_csv(TIER1_FILE))

    write_csv(tier1, OUT / "Tier1_corrected_standardized_genes.csv")

    all_prepared = pd.concat([four_group, frontal], ignore_index=True, sort=False)
    summary_rows: list[dict[str, object]] = []

    for label, df in all_prepared.groupby("region_or_subregion", sort=False):
        analysis_level = df["analysis_level"].iloc[0]
        safe = safe_name(label)

        tier2_sets = make_tier2_sets(df, label)
        for set_name, set_df in tier2_sets.items():
            clean = clean_gene_table(set_df)
            write_csv(clean, TIER2_DIR / f"{safe}__{set_name}.csv")
            write_gene_symbols(clean, ENRICHMENT_DIR / f"{safe}__{set_name}__gene_symbols.txt")

        tier3 = clean_gene_table(make_tier3_set(df))
        write_csv(tier3, TIER3_DIR / f"{safe}__Tier3_raw_beta_comparison.csv")

        write_ranked_signature(df, label)

        tier1_for_label = tier1[tier1["region_or_subregion"] == label].copy()
        if not tier1_for_label.empty:
            write_gene_symbols(tier1_for_label, ENRICHMENT_DIR / f"{safe}__Tier1_strict_all__gene_symbols.txt")
            write_gene_symbols(
                tier1_for_label[tier1_for_label["direction"] == "up_in_MDD"],
                ENRICHMENT_DIR / f"{safe}__Tier1_strict_up_in_MDD__gene_symbols.txt",
            )
            write_gene_symbols(
                tier1_for_label[tier1_for_label["direction"] == "down_in_MDD"],
                ENRICHMENT_DIR / f"{safe}__Tier1_strict_down_in_MDD__gene_symbols.txt",
            )

        summary_rows.append(
            {
                "analysis_level": analysis_level,
                "region_or_subregion": label,
                "n_Tier1": int(len(tier1_for_label)),
                "n_Tier2_all": int(len(tier2_sets["Tier2_all"])),
                "n_Tier2_up_in_MDD": int(len(tier2_sets["Tier2_up_in_MDD"])),
                "n_Tier2_down_in_MDD": int(len(tier2_sets["Tier2_down_in_MDD"])),
                "n_Tier3": int(len(tier3)),
                "n_background_genes": int(df["gene_id"].nunique()),
                "recommended_use": RECOMMENDED_USE[label],
            }
        )

    four_background = four_group[
        ["gene_id", "gene_symbol", "region_or_subregion", "analysis_level"]
    ].drop_duplicates()
    four_background["genes_tested_flag"] = True
    write_csv(four_background, OUT / "background_genes_four_group_by_region.csv")

    frontal_background = frontal[
        ["gene_id", "gene_symbol", "region_or_subregion", "analysis_level"]
    ].drop_duplicates()
    frontal_background["genes_tested_flag"] = True
    write_csv(frontal_background, OUT / "background_genes_frontal_subregion_by_subregion.csv")

    summary = pd.DataFrame(summary_rows)
    order = [
        "Ventral Striatum",
        "Insula",
        "Hippocampal Formation",
        "broad Frontal Cortex",
        "Cingulate gyrus 25",
        "dlPFC BA8/9",
        "OFC BA11",
    ]
    summary["_order"] = summary["region_or_subregion"].map({name: i for i, name in enumerate(order)})
    summary = summary.sort_values("_order").drop(columns="_order")
    write_csv(summary, OUT / "DE_candidate_gene_set_summary.csv")

    write_ml_note()
    write_report(summary)
    write_readme()

    print(f"Saved DE candidate gene-set package to: {OUT}")
    print(summary.to_string(index=False))


def write_ml_note() -> None:
    text = """# ML Preparation Note

This folder packages candidate gene sets for later ML preparation. No ML model was run in this step.

Tier 2 exploratory DE candidate genes can be used as biologically motivated candidate features. Ranked Tier 4 signatures can be used for feature ranking or for defining ordered gene lists.

For unbiased predictive ML, feature selection must be performed inside cross-validation. If Tier 2 genes were selected using the full dataset, then classifier performance using those genes on the same dataset is exploratory and may be optimistic.

Recommended first ML models:

1. Elastic Net logistic regression
2. LASSO logistic regression
3. Linear SVM
4. Random Forest / XGBoost as exploratory models

Deep learning is not recommended as the first model because sample size is small relative to gene count. Module/eigengene-level ML after AHBA/WGCNA may be more stable than raw gene-level deep learning.
"""
    (OUT / "ML_Preparation_Note.md").write_text(text, encoding="utf-8")


def write_report(summary: pd.DataFrame) -> None:
    summary_md = summary.to_markdown(index=False)
    text = f"""# DE Candidate Gene Sets Report

## Final Decision

Tier 2 is the working exploratory DE candidate set for the Length > 100 MDD analysis.

Tier 2 definition:

```text
p_value_diagnosis < 0.01 and abs_standardized_effect_size_diagnosis > 0.8
```

Tier 2 genes are not FDR-significant DE genes. They should be described as Tier 2 exploratory DE candidate genes, nominal standardized-effect-size DE candidates, or working MDD-associated gene-expression candidates.

## Evidence Structure

- Tier 1 genes are strict corrected anchors: Storey q < 0.10 and abs standardized effect size > 0.8.
- Tier 2 genes are the working exploratory DE candidate genes.
- Tier 3 genes are preserved as the raw-beta comparison baseline.
- Tier 4 ranked signatures are full-gene ranked lists for enrichment, ML feature ranking, and AHBA/module mapping.

## Summary

{summary_md}

## Recommended Use

Ventral Striatum is the strongest broad-region signal and should be treated as the main robust broad-region set.

Insula and Hippocampal Formation provide exploratory/supportive broad-region sets.

Broad Frontal Cortex is retained for documentation, but it is not recommended as the primary frontal gene set because it has only 1 Tier 2 gene.

Frontal subregion Tier 2 sets should be used for frontal biology. Cingulate gyrus 25 is the strongest frontal exploratory signal, while dlPFC BA8/9 and OFC BA11 are supportive frontal sets.

## P-Value Correction Sensitivity

The p-value correction sensitivity audit showed that under standard correction across all genes, almost nothing survives at adjusted p < 0.10 plus abs standardized effect size > 0.8, except Ventral Striatum has 4 Storey-q genes. Larger counts mainly appear under effect-size-gated and direction-specific correction scopes, so those remain sensitivity results and do not replace Tier 2.

## Output Contents

- `Tier1_corrected_standardized_genes.csv`: strict corrected anchor genes.
- `Tier2_working_exploratory/`: main working exploratory candidate sets.
- `Tier3_raw_beta_comparison/`: raw-beta comparison sets.
- `ranked_signatures/`: full-gene Tier 4 ranked signatures.
- `enrichment_ready/`: plain text gene-symbol lists, one symbol per line.
- `background_genes_four_group_by_region.csv`: enrichment background for broad regions.
- `background_genes_frontal_subregion_by_subregion.csv`: enrichment background for frontal subregions.
- `ML_Preparation_Note.md`: ML feature-selection cautions and recommended first models.

## Source Inputs

This package was created from existing outputs only. No LMM was rerun.

- `{FOUR_GROUP_FILE.relative_to(REPO_DIR)}`
- `{TIER1_FILE.relative_to(REPO_DIR)}`
- `{FRONTAL_SUBREGION_FILE.relative_to(REPO_DIR)}`
- `{FRONTAL_SUBREGION_SUMMARY_FILE.relative_to(REPO_DIR)}`
- `{PVALUE_AUDIT_REPORT.relative_to(REPO_DIR)}`
"""
    (OUT / "DE_Candidate_Gene_Sets_Report.md").write_text(text, encoding="utf-8")


def write_readme() -> None:
    text = """# Length > 100 DE Candidate Gene Sets

This folder contains a clean package of candidate gene sets for the Length > 100 MDD analysis.

Contents:

- `Tier1_corrected_standardized_genes.csv`: strict corrected genes, used as high-confidence anchors.
- `Tier2_working_exploratory/`: working exploratory DE candidate genes using p < 0.01 and abs standardized effect size > 0.8.
- `Tier3_raw_beta_comparison/`: raw-beta comparison genes, preserved for comparison only.
- `ranked_signatures/`: full-gene Tier 4 ranked signatures using signed standardized log-p scores.
- `enrichment_ready/`: plain text gene-symbol lists for enrichment tools.
- `background_genes_four_group_by_region.csv`: all tested genes for each broad-region enrichment background.
- `background_genes_frontal_subregion_by_subregion.csv`: all tested genes for each frontal-subregion enrichment background.
- `DE_candidate_gene_set_summary.csv`: counts and recommended use for each region/subregion.
- `ML_Preparation_Note.md`: guidance for using these sets in ML without leakage.
- `DE_Candidate_Gene_Sets_Report.md`: interpretation and recommended use.

Tier 2 is exploratory and should not be described as FDR-significant.
"""
    (OUT / "README.md").write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
