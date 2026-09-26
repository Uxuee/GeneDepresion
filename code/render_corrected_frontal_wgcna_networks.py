"""Render original frontal-subregion WGCNA networks with canonical ML annotations.

This is a downstream-only renderer: saved module assignments and expression matrices
are read unchanged.  No WGCNA construction or module selection is performed.
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd

PROJECT_ROOT = Path(r"D:/genes/GeneDepresion")
CODE_ROOT = PROJECT_ROOT / "code"
sys.path.insert(0, str(CODE_ROOT))

import create_mdd_wgcna_full_module_networks as network  # noqa: E402

PIPELINE_ROOT = PROJECT_ROOT / "outputs" / "length_gt_100_clean_pipeline"
SOURCE_ROOT = PIPELINE_ROOT / "mdd_wgcna_modules" / "frontal_subregion_wgcna"
OUT_DIR = (
    PIPELINE_ROOT
    / "16_corrected_ML_downstream"
    / "WGCNA_ML_annotations_corrected"
    / "network_figures_frontal_subregions"
)
ASSIGNMENTS = SOURCE_ROOT / "frontal_subregion_best_solution_gene_to_module_assignments.csv"
ENRICHMENT = (
    PIPELINE_ROOT
    / "16_corrected_ML_downstream"
    / "WGCNA_ML_annotations_corrected"
    / "frontal_subregions_ML_module_enrichment_corrected.csv"
)
EXPR_ROOT = PIPELINE_ROOT / "02_normalized_expression" / "frontal_subregion_wgcna_inputs"
EXPR_FILES = {
    "Cingulate gyrus 25": EXPR_ROOT / "normalized_GeTMM_log2_CG25.csv",
    "dlPFC BA8/9": EXPR_ROOT / "normalized_GeTMM_log2_dlPFC_BA8_9.csv",
    "OFC BA11": EXPR_ROOT / "normalized_GeTMM_log2_OFC_BA11.csv",
}
FOCUS_MODULES = [
    ("dlPFC BA8/9", "MEblue"),
    ("dlPFC BA8/9", "MEorange"),
    ("OFC BA11", "MEtan"),
    ("Cingulate gyrus 25", "MEwhite"),
]


def load_expression(region: str) -> pd.DataFrame:
    expr = pd.read_csv(EXPR_FILES[region], index_col=0)
    expr.index = expr.index.astype(str)
    expr.columns = expr.columns.astype(str)
    return expr


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    network.OUT_DIR = OUT_DIR
    assignments = pd.read_csv(ASSIGNMENTS)
    enrichment = network.to_num(pd.read_csv(ENRICHMENT), ["overlap_count", "p_value", "fdr_BH"])
    tier2, high, high_mod = network.build_gene_sets()
    graphs = {}
    inventory = []
    for region, module in FOCUS_MODULES:
        genes = assignments[(assignments["region"] == region) & (assignments["module"] == module)].copy()
        if genes.empty:
            raise ValueError(f"No saved assignments found for {region} {module}")
        graph, pos, plotted_genes, median_edge = network.compute_graph(genes, load_expression(region))
        enrich_row = network.best_enrichment(region, module, enrichment)
        graphs[(region, module)] = (graph, pos, plotted_genes, enrich_row)
        outputs = network.save_one_module(
            region, module, graph, pos, plotted_genes, enrich_row, tier2, high, high_mod
        )
        inventory.append({
            "region": region,
            "module": module,
            "n_nodes": graph.number_of_nodes(),
            "n_edges": graph.number_of_edges(),
            "median_edge_abs_correlation": median_edge,
            "unlabeled_png": outputs["unlabeled"],
            "labeled_png": outputs["labeled"],
            "module_assignments_source": str(ASSIGNMENTS),
            "ML_enrichment_source": str(ENRICHMENT),
        })
        print(f"{region} {module}: nodes={graph.number_of_nodes()} edges={graph.number_of_edges()}")

    fig, axes = plt.subplots(2, 2, figsize=(15, 13))
    for ax, key in zip(axes.ravel(), FOCUS_MODULES):
        graph, pos, genes, enrich_row = graphs[key]
        network.draw_network(
            graph, pos, genes, key[0], key[1], enrich_row,
            tier2, high, high_mod, labeled=False, summary_ax=ax,
        )
    network.add_legend(fig, axes.ravel()[0])
    fig.suptitle(
        "Original frontal-subregion WGCNA networks with corrected ML annotations",
        fontsize=18, fontweight="bold", y=0.995,
    )
    fig.text(
        0.5, 0.045,
        "Saved module membership is unchanged; node colors use the canonical corrected ML ranking.",
        ha="center", fontsize=10,
    )
    fig.tight_layout(rect=[0.02, 0.08, 0.98, 0.97])
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"frontal_subregion_network_summary_corrected_ML.{ext}", dpi=network.PLOT_DPI, bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(inventory).to_csv(OUT_DIR / "frontal_subregion_network_inventory.csv", index=False)
    (OUT_DIR / "README.md").write_text(
        "# Corrected ML-annotated frontal-subregion networks\n\n"
        "These figures reuse the saved WGCNA module assignments and normalized expression matrices. "
        "Only Tier2/canonical ML node annotations and enrichment subtitles were refreshed. "
        "Broad pooled Frontal Cortex is not included.\n",
        encoding="utf-8",
    )
    print(f"Wrote corrected frontal-subregion networks to {OUT_DIR}")


if __name__ == "__main__":
    main()
