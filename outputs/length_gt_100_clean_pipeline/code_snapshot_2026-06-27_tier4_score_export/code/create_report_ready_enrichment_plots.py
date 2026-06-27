"""Create report-ready enrichment plots with larger fonts.

This script uses existing completed enrichment tables only. It does not rerun
enrichment and writes improved copies to a new report-ready image folder.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"D:\genes\GeneDepresion")
ANALYSIS = ROOT / "outputs" / "length_gt_100_elisa_style"
TABLE_DIR = ANALYSIS / "enrichment_Tier2" / "tables"
OUT_DIR = ANALYSIS / "images" / "enrichment_report_ready"


PLOTS = [
    {
        "label": "Hippocampal Formation up in MDD KEGG",
        "short_title": "HF up in MDD: KEGG",
        "source": TABLE_DIR / "four_group__Hippocampal_Formation__up_in_MDD__KEGG__full.csv",
        "dotplot": "HF_up_MDD_KEGG_dotplot_largefont.png",
        "barplot": "HF_up_MDD_KEGG_barplot_largefont.png",
    },
    {
        "label": "Insula down in MDD MSigDB C2 curated",
        "short_title": "Insula down in MDD: MSigDB C2",
        "source": TABLE_DIR / "four_group__Insula__down_in_MDD__MSigDB_C2_curated__full.csv",
        "dotplot": "Insula_down_MDD_MSigDB_C2_dotplot_largefont.png",
        "barplot": "Insula_down_MDD_MSigDB_C2_barplot_largefont.png",
    },
    {
        "label": "Insula down in MDD GO molecular function",
        "short_title": "Insula down in MDD: GO MF",
        "source": TABLE_DIR / "four_group__Insula__down_in_MDD__GO_MF__full.csv",
        "dotplot": "Insula_down_MDD_GO_MF_dotplot_largefont.png",
        "barplot": "Insula_down_MDD_GO_MF_barplot_largefont.png",
    },
    {
        "label": "Insula down in MDD Hallmark",
        "short_title": "Insula down in MDD: Hallmark",
        "source": TABLE_DIR / "four_group__Insula__down_in_MDD__MSigDB_Hallmark__full.csv",
        "dotplot": "Insula_down_MDD_Hallmark_dotplot_largefont.png",
        "barplot": "Insula_down_MDD_Hallmark_barplot_largefont.png",
    },
]


def parse_gene_ratio(value: object) -> float:
    """Convert enrichment ratios like '7/26' to numeric values."""
    if pd.isna(value):
        return np.nan
    text = str(value)
    if "/" in text:
        num, den = text.split("/", 1)
        try:
            return float(num) / float(den)
        except ValueError:
            return np.nan
    try:
        return float(text)
    except ValueError:
        return np.nan


def clean_description(value: object, width: int = 38) -> str:
    """Make long enrichment names readable on horizontal report plots."""
    text = str(value).replace("_", " ")
    text = " ".join(text.split())
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False))


def load_top_terms(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if df.empty:
        return df
    df = df.copy()
    df["p.adjust"] = pd.to_numeric(df["p.adjust"], errors="coerce")
    df["Count"] = pd.to_numeric(df["Count"], errors="coerce")
    df["gene_ratio_numeric"] = df["GeneRatio"].map(parse_gene_ratio)
    df["neg_log10_padj"] = -np.log10(df["p.adjust"].clip(lower=np.nextafter(0, 1)))
    df = df.sort_values(["p.adjust", "pvalue"], na_position="last").head(10)
    df["wrapped_description"] = df["Description"].map(clean_description)
    return df.iloc[::-1].reset_index(drop=True)


def draw_dotplot(df: pd.DataFrame, title: str, dest: Path) -> bool:
    if df.empty:
        return False
    height = max(6.5, 0.72 * len(df) + 2.2)
    fig, ax = plt.subplots(figsize=(12.5, height))
    scatter = ax.scatter(
        df["gene_ratio_numeric"],
        df["wrapped_description"],
        s=(df["Count"].fillna(1) * 95).clip(lower=110),
        c=df["p.adjust"],
        cmap="viridis_r",
        edgecolor="#263238",
        linewidth=0.7,
        alpha=0.92,
    )
    ax.set_title(title, fontsize=22, fontweight="bold", pad=18)
    ax.set_xlabel("Gene ratio", fontsize=17, labelpad=10)
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    ax.grid(axis="x", color="#dfe5e8", linewidth=0.9)
    ax.set_axisbelow(True)

    cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label("Adjusted p-value", fontsize=15)
    cbar.ax.tick_params(labelsize=12)

    handles, labels = scatter.legend_elements(prop="sizes", num=4, func=lambda s: s / 95)
    legend = ax.legend(
        handles,
        labels,
        title="Gene count",
        loc="lower right",
        frameon=True,
        fontsize=12,
        title_fontsize=13,
    )
    ax.add_artist(legend)
    fig.tight_layout()
    fig.savefig(dest, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def draw_barplot(df: pd.DataFrame, title: str, dest: Path) -> bool:
    if df.empty:
        return False
    height = max(6.5, 0.72 * len(df) + 2.2)
    fig, ax = plt.subplots(figsize=(12.5, height))
    colors = plt.cm.viridis_r(
        (df["p.adjust"] - df["p.adjust"].min())
        / max(df["p.adjust"].max() - df["p.adjust"].min(), np.nextafter(0, 1))
    )
    bars = ax.barh(df["wrapped_description"], df["neg_log10_padj"], color=colors, edgecolor="#263238", linewidth=0.6)
    ax.set_title(title, fontsize=22, fontweight="bold", pad=18)
    ax.set_xlabel("-log10(adjusted p-value)", fontsize=17, labelpad=10)
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    ax.grid(axis="x", color="#dfe5e8", linewidth=0.9)
    ax.set_axisbelow(True)

    for bar, count in zip(bars, df["Count"]):
        ax.text(
            bar.get_width() + 0.03,
            bar.get_y() + bar.get_height() / 2,
            f"n={int(count)}",
            va="center",
            ha="left",
            fontsize=12,
            color="#263238",
        )

    sm = plt.cm.ScalarMappable(cmap="viridis_r", norm=plt.Normalize(df["p.adjust"].min(), df["p.adjust"].max()))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("Adjusted p-value", fontsize=15)
    cbar.ax.tick_params(labelsize=12)
    fig.tight_layout()
    fig.savefig(dest, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def write_readme(results: list[dict[str, object]]) -> None:
    created_lines = []
    source_lines = []
    for result in results:
        created_lines.append(
            f"- {result['label']}: `{result['dotplot']}` and `{result['barplot']}`"
        )
        source_lines.append(f"- {result['label']}: `{result['source']}`")

    (OUT_DIR / "README.md").write_text(
        "\n".join(
            [
                "# Report-Ready Enrichment Figures",
                "",
                "These figures were regenerated from completed enrichment result tables only. Enrichment was not rerun.",
                "",
                "## Selected Main-Report Figures",
                "- Hippocampal Formation up_in_MDD KEGG",
                "- Insula down_in_MDD MSigDB_C2_curated",
                "- Insula down_in_MDD GO_MF",
                "- Insula down_in_MDD Hallmark",
                "",
                "Each selected analysis has a larger-font dotplot and a horizontal barplot. Long term names are wrapped, the canvas is wider/taller, and plots are saved at 300 dpi.",
                "",
                "## Created Files",
                *created_lines,
                "",
                "## Source Tables",
                *source_lines,
                "",
                "## Recommendation",
                "Main report:",
                "- Hippocampal Formation up_in_MDD KEGG",
                "- Insula down_in_MDD MSigDB_C2_curated",
                "- Insula down_in_MDD GO_MF or Hallmark",
                "",
                "Supplementary:",
                "- enrichment maps / network plots",
                "- sparse single-term plots",
                "- remaining database-specific plots",
                "",
                "Network/enrichment-map figures were not selected for the main report because they are harder to read at report scale.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = []
    for plot in PLOTS:
        df = load_top_terms(plot["source"])
        dot_dest = OUT_DIR / plot["dotplot"]
        bar_dest = OUT_DIR / plot["barplot"]
        dot_ok = draw_dotplot(df, plot["short_title"], dot_dest)
        bar_ok = draw_barplot(df, plot["short_title"], bar_dest)
        results.append(
            {
                **plot,
                "dotplot_created": dot_ok,
                "barplot_created": bar_ok,
                "n_terms_plotted": len(df),
            }
        )

    pd.DataFrame(results).to_csv(OUT_DIR / "report_ready_enrichment_plot_manifest.csv", index=False)
    write_readme(results)

    print(f"output_folder={OUT_DIR}")
    print("regenerated_plots:")
    for result in results:
        print(f"- {result['dotplot']}")
        print(f"- {result['barplot']}")
    print("barplot_versions_created=" + str(all(result["barplot_created"] for result in results)))


if __name__ == "__main__":
    main()
