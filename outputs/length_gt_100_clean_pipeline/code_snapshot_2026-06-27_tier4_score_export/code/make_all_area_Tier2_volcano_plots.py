"""Create consolidated Tier 2 volcano plots for all Length > 100 areas.

Tier 2 definition:
    p_value_diagnosis < 0.01 and abs_standardized_effect_size_diagnosis > 0.8

The script reads existing standardized-effect LMM outputs and writes a new
plot folder without overwriting older volcano plot directories.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"D:\genes\GeneDepresion")
OUT_DIR = (
    ROOT
    / "outputs"
    / "length_gt_100_elisa_style"
    / "volcano_plots_Tier2_all_areas"
)

BROAD_TABLE = (
    ROOT
    / "outputs"
    / "length_gt_100_elisa_style"
    / "standardized_effect_size_4groups"
    / "LMM_all_genes_with_standardized_effect_size_4groups.csv"
)
FRONTAL_SUBREGION_TABLE = (
    ROOT
    / "outputs"
    / "length_gt_100_elisa_style"
    / "frontal_subregion_sensitivity"
    / "frontal_subregion_LMM_all_genes_with_standardized_effect_size.csv"
)


AREA_ORDER = [
    ("broad_region", "Frontal Cortex", "Frontal_Cortex"),
    ("broad_region", "Insula", "Insula"),
    ("broad_region", "Ventral Striatum", "Ventral_Striatum"),
    ("broad_region", "Hippocampal Formation", "Hippocampal_Formation"),
    ("frontal_subregion", "Cingulate gyrus 25 (Cg25)", "Cingulate_gyrus_25_Cg25"),
    ("frontal_subregion", "dlPFC BA8/9", "dlPFC_BA8_9"),
    ("frontal_subregion", "OFC BA11", "OFC_BA11"),
]


def clean_pvalues(series: pd.Series) -> pd.Series:
    """Clip p-values away from zero so -log10 remains finite."""
    values = pd.to_numeric(series, errors="coerce")
    tiny = np.nextafter(0, 1)
    return values.clip(lower=tiny, upper=1)


def prepare_broad_regions() -> pd.DataFrame:
    df = pd.read_csv(BROAD_TABLE)
    df = df.copy()
    df["area_type"] = "broad_region"
    df["area_name"] = df["grouped_region"].astype(str).str.replace("_", " ", regex=False)
    df["area_slug"] = df["grouped_region"].astype(str)
    return df


def prepare_frontal_subregions() -> pd.DataFrame:
    df = pd.read_csv(FRONTAL_SUBREGION_TABLE)
    df = df.copy()
    df["area_type"] = "frontal_subregion"
    name_map = {
        "Cingulate gyrus 25 (Cg25)": ("Cingulate gyrus 25 (Cg25)", "Cingulate_gyrus_25_Cg25"),
        "Dorsolateral prefrontal cortex (dlPFC; BA8/9)": ("dlPFC BA8/9", "dlPFC_BA8_9"),
        "Orbitofrontal (OFC; BA11)": ("OFC BA11", "OFC_BA11"),
    }
    df["area_name"] = df["original_frontal_subregion"].map(lambda x: name_map.get(x, (x, x))[0])
    df["area_slug"] = df["original_frontal_subregion"].map(lambda x: name_map.get(x, (x, x))[1])
    return df


def add_tier2_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["p_value_diagnosis"] = clean_pvalues(df["p_value_diagnosis"])
    df["minus_log10_p"] = -np.log10(df["p_value_diagnosis"])
    df["abs_standardized_effect_size_diagnosis"] = pd.to_numeric(
        df["abs_standardized_effect_size_diagnosis"], errors="coerce"
    )
    df["standardized_effect_size_diagnosis"] = pd.to_numeric(
        df["standardized_effect_size_diagnosis"], errors="coerce"
    )
    df["Tier2"] = (
        (df["p_value_diagnosis"] < 0.01)
        & (df["abs_standardized_effect_size_diagnosis"] > 0.8)
    )
    return df


def plot_area(df: pd.DataFrame, area_type: str, area_name: str, area_slug: str) -> dict[str, object]:
    area_df = df[(df["area_type"] == area_type) & (df["area_slug"] == area_slug)].copy()
    if area_df.empty:
        raise ValueError(f"No rows found for {area_name}")

    tier2 = area_df[area_df["Tier2"]].copy()
    up = tier2[tier2["standardized_effect_size_diagnosis"] > 0]
    down = tier2[tier2["standardized_effect_size_diagnosis"] < 0]

    fig, ax = plt.subplots(figsize=(8.5, 6.2), dpi=160)
    ax.scatter(
        area_df["standardized_effect_size_diagnosis"],
        area_df["minus_log10_p"],
        s=8,
        c="#b7c0c7",
        alpha=0.42,
        linewidths=0,
        label="All tested genes",
    )
    ax.scatter(
        down["standardized_effect_size_diagnosis"],
        down["minus_log10_p"],
        s=22,
        c="#2c7fb8",
        alpha=0.88,
        linewidths=0,
        label="Tier 2 down in MDD",
    )
    ax.scatter(
        up["standardized_effect_size_diagnosis"],
        up["minus_log10_p"],
        s=22,
        c="#d95f02",
        alpha=0.88,
        linewidths=0,
        label="Tier 2 up in MDD",
    )

    ax.axvline(-0.8, color="#2f3437", linestyle="--", linewidth=1)
    ax.axvline(0.8, color="#2f3437", linestyle="--", linewidth=1)
    ax.axhline(-np.log10(0.01), color="#2f3437", linestyle="--", linewidth=1)
    ax.set_xlabel("Standardized effect size for diagnosis")
    ax.set_ylabel("-log10(p-value for diagnosis)")
    ax.set_title(f"Tier 2 volcano: {area_name}")
    ax.text(
        0.01,
        0.99,
        f"Tier 2: p < 0.01 and |standardized effect| > 0.8\n"
        f"n tested = {len(area_df):,}; Tier 2 = {len(tier2):,}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "#d8dde1", "boxstyle": "round,pad=0.35"},
    )

    if not tier2.empty:
        top_labels = tier2.nsmallest(6, "p_value_diagnosis")
        for _, row in top_labels.iterrows():
            label = str(row.get("gene_symbol") or row.get("gene_id"))
            ax.annotate(
                label,
                (row["standardized_effect_size_diagnosis"], row["minus_log10_p"]),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=7.5,
                color="#1f2528",
            )

    ax.legend(loc="lower right", frameon=True, fontsize=8)
    ax.grid(color="#e8ecef", linewidth=0.7)
    fig.tight_layout()

    path = OUT_DIR / f"{area_slug}_Tier2_volcano_standardized_effect.png"
    fig.savefig(path)
    plt.close(fig)

    return {
        "area_type": area_type,
        "area_name": area_name,
        "area_slug": area_slug,
        "genes_tested": len(area_df),
        "tier2_genes": len(tier2),
        "tier2_up_in_MDD": len(up),
        "tier2_down_in_MDD": len(down),
        "plot_file": str(path),
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    broad = prepare_broad_regions()
    frontal = prepare_frontal_subregions()
    all_rows = add_tier2_columns(pd.concat([broad, frontal], ignore_index=True, sort=False))

    summaries = [
        plot_area(all_rows, area_type, area_name, area_slug)
        for area_type, area_name, area_slug in AREA_ORDER
    ]
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT_DIR / "Tier2_volcano_all_areas_summary.csv", index=False)


if __name__ == "__main__":
    main()
