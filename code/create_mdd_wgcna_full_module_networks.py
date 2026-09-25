#!/usr/bin/env python
"""Create full-module WGCNA network figures for selected Step 12 modules.

Uses saved MDD WGCNA outputs and normalized GeTMM-log2 matrices only.
Does not rerun WGCNA.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Tuple
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
import numpy as np
import pandas as pd
from scipy.sparse.csgraph import minimum_spanning_tree
from sklearn.preprocessing import MinMaxScaler

PROJECT_ROOT = Path(r"D:/genes/GeneDepresion")
PIPELINE_ROOT = PROJECT_ROOT / "outputs" / "length_gt_100_clean_pipeline"
WGCNA_ROOT = PIPELINE_ROOT / "12_MDD_WGCNA_modules"
OUT_DIR = Path(os.environ.get(
    "MDD_WGCNA_NETWORK_OUTPUT_DIR",
    str(WGCNA_ROOT / "plots_cluster_visualization"),
))
ENRICHMENT_FILE = Path(os.environ.get(
    "MDD_WGCNA_ENRICHMENT_FILE",
    str(WGCNA_ROOT / "module_enrichment_for_Tier2_and_ML_genes.csv"),
))

EXPR_ROOT = PROJECT_ROOT / "outputs" / "step1_length_sensitivity" / "Length_gt_100" / "02_normalized_expression" / "keep_all_samples"
EXPR_FILES = {
    "Frontal Cortex": EXPR_ROOT / "normalized_GeTMM_log2_frontal_cortex.csv",
    "Hippocampal Formation": EXPR_ROOT / "normalized_GeTMM_log2_hippocampal_formation.csv",
    "Insula": EXPR_ROOT / "normalized_GeTMM_log2_insula.csv",
    "Ventral Striatum": EXPR_ROOT / "normalized_GeTMM_log2_ventral_striatum.csv",
}

FOCUS_MODULES = [
    ("Ventral Striatum", "MElavenderblush3"),
    ("Ventral Striatum", "MElightcyan1"),
    ("Ventral Striatum", "MEdarkred"),
    ("Ventral Striatum", "MElightgreen"),
    ("Insula", "MEorange"),
    ("Hippocampal Formation", "MEgreen"),
    ("Frontal Cortex", "MEdarkgreen"),
]

NODE_COLORS = {
    "default": "#b7c7d6",
    "tier2": "#f28e2b",
    "ml_high": "#d62728",
    "ml_high_moderate": "#f2c94c",
}
EDGE_COLOR = "#6b7280"
PLOT_DPI = 240
TSNE_SKIP_REASON = "t-SNE attempted but skipped because sklearn failed in this environment."
README_AUDIT_SENTENCE = (
    "An audit confirmed the presence of 32 PNG plots, 32 PDF plots, two CSV summary files, "
    "and the README. t-SNE visualization was attempted but not generated in this environment "
    "because sklearn support failed, so no t-SNE images are included."
)


def safe_name(value: str) -> str:
    return str(value).replace("/", "_").replace("\\", "_").replace(" ", "_").replace(":", "_")


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def to_num(df: pd.DataFrame, cols: Iterable[str]) -> pd.DataFrame:
    for col in cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def load_expression(region: str) -> pd.DataFrame:
    expr = pd.read_csv(EXPR_FILES[region], index_col=0)
    expr.index = expr.index.astype(str)
    expr.columns = expr.columns.astype(str)
    return expr


def build_gene_sets() -> tuple[set[str], set[str], set[str]]:
    tier2 = read_csv(PIPELINE_ROOT / "07_Tier2_candidate_gene_sets" / "Tier2_candidate_genes_all_regions.csv")
    ml = read_csv(PIPELINE_ROOT / "10_ML_prioritized_Tier2_corrected" / "ML_prioritized_Tier2_genes_all_ranked.csv")
    tier2_ids = set(tier2["gene_id"].dropna().astype(str))
    cls = ml["ML_priority_class"].fillna("").astype(str).str.lower()
    high_ids = set(ml.loc[cls.str.contains("high-confidence"), "gene_id"].dropna().astype(str))
    high_mod_ids = set(ml.loc[cls.str.contains("high-confidence|moderate-confidence"), "gene_id"].dropna().astype(str))
    return tier2_ids, high_ids, high_mod_ids


def classify_gene(gene_id: str, tier2: set[str], high: set[str], high_mod: set[str]) -> str:
    if gene_id in high:
        return "ml_high"
    if gene_id in high_mod:
        return "ml_high_moderate"
    if gene_id in tier2:
        return "tier2"
    return "default"


def sparse_edges_from_correlation(corr_abs: np.ndarray, top_k: int = 3, min_extra_corr: float = 0.45) -> set[tuple[int, int]]:
    n = corr_abs.shape[0]
    if n <= 1:
        return set()
    safe_corr = np.nan_to_num(corr_abs, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(safe_corr, 1.0)
    # Connectivity backbone: MST on distance, equivalent to maximum spanning tree on abs correlation.
    dist = 1.0 - safe_corr
    mst = minimum_spanning_tree(dist).tocoo()
    edges: set[tuple[int, int]] = set()
    for i, j in zip(mst.row, mst.col):
        a, b = sorted((int(i), int(j)))
        if a != b:
            edges.add((a, b))
    # Add local strongest edges per node to reveal module substructure.
    for i in range(n):
        order = np.argsort(safe_corr[i])[::-1]
        added = 0
        for j in order:
            j = int(j)
            if j == i:
                continue
            if safe_corr[i, j] < min_extra_corr and added >= 1:
                break
            a, b = sorted((i, j))
            edges.add((a, b))
            added += 1
            if added >= top_k:
                break
    return edges


def compute_graph(module_genes: pd.DataFrame, expr: pd.DataFrame) -> tuple[nx.Graph, Dict[str, tuple[float, float]], pd.DataFrame, float]:
    genes = module_genes.copy()
    genes["gene_id"] = genes["gene_id"].astype(str)
    genes = genes[genes["gene_id"].isin(expr.index)].drop_duplicates("gene_id")
    if len(genes) < 3:
        raise ValueError("Not enough module genes found in expression matrix")
    mat = expr.loc[genes["gene_id"].tolist()].astype(float)
    gene_ids = mat.index.tolist()
    # Correlation across samples between genes.
    corr = np.corrcoef(mat.values)
    corr_abs = np.abs(corr)
    edges_idx = sparse_edges_from_correlation(corr_abs)
    hubness = np.nanmean(np.nan_to_num(corr_abs, nan=0.0), axis=1)
    genes = genes.set_index("gene_id").loc[gene_ids].reset_index()
    genes["hubness"] = hubness
    scaler = MinMaxScaler(feature_range=(55, 440))
    genes["node_size"] = scaler.fit_transform(genes[["hubness"]]).ravel() if len(genes) > 1 else 240
    G = nx.Graph()
    for _, row in genes.iterrows():
        G.add_node(row["gene_id"], gene_symbol=row.get("gene_symbol", row["gene_id"]), hubness=row["hubness"], node_size=row["node_size"])
    for i, j in edges_idx:
        g1, g2 = gene_ids[i], gene_ids[j]
        w = float(corr_abs[i, j]) if not np.isnan(corr_abs[i, j]) else 0.0
        G.add_edge(g1, g2, weight=w)
    # Force-directed layout. More iterations for tiny modules, fewer for larger ones.
    iterations = 180 if len(G) < 80 else 95 if len(G) < 400 else 70
    k = 1.1 / np.sqrt(max(len(G), 2))
    pos = nx.spring_layout(G, seed=17, weight="weight", iterations=iterations, k=k)
    return G, pos, genes, float(np.nanmedian([d.get("weight", 0) for _, _, d in G.edges(data=True)]))


def best_enrichment(region: str, module: str, enrich: pd.DataFrame) -> pd.Series:
    sub = enrich[(enrich["region"] == region) & (enrich["module"] == module)].copy()
    if sub.empty:
        return pd.Series({"gene_set": "NA", "fdr_BH": np.nan, "p_value": np.nan, "overlap_count": np.nan})
    return sub.sort_values(["fdr_BH", "p_value", "overlap_count"], ascending=[True, True, False]).iloc[0]


def draw_network(
    G: nx.Graph,
    pos: Dict[str, tuple[float, float]],
    genes: pd.DataFrame,
    region: str,
    module: str,
    enrich_row: pd.Series,
    tier2: set[str],
    high: set[str],
    high_mod: set[str],
    labeled: bool,
    summary_ax=None,
) -> None:
    ax = summary_ax if summary_ax is not None else plt.gca()
    genes = genes.copy()
    genes["class"] = genes["gene_id"].map(lambda g: classify_gene(g, tier2, high, high_mod))
    genes["is_overlap"] = genes["class"] != "default"
    hub_cut = genes["hubness"].quantile(0.90)
    genes["hub_overlap"] = genes["is_overlap"] & (genes["hubness"] >= hub_cut)
    node_color = genes["class"].map(NODE_COLORS).fillna(NODE_COLORS["default"]).tolist()
    node_sizes = genes.set_index("gene_id").loc[list(G.nodes), "node_size"].tolist()
    gene_class = genes.set_index("gene_id")["class"].to_dict()
    hub_overlap = genes.set_index("gene_id")["hub_overlap"].to_dict()
    edgecolors = ["#111827" if hub_overlap.get(n, False) else ("#8a4b00" if gene_class.get(n) != "default" else "#f8fafc") for n in G.nodes]
    linewidths = [1.8 if hub_overlap.get(n, False) else (1.1 if gene_class.get(n) != "default" else 0.35) for n in G.nodes]
    edge_widths = [0.15 + 0.8 * G[u][v].get("weight", 0) for u, v in G.edges]
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color=EDGE_COLOR, alpha=0.16, width=edge_widths)
    nx.draw_networkx_nodes(
        G,
        pos,
        ax=ax,
        nodelist=list(G.nodes),
        node_color=node_color,
        node_size=node_sizes,
        edgecolors=edgecolors,
        linewidths=linewidths,
        alpha=0.92,
    )
    if labeled:
        overlap = genes[genes["is_overlap"]].sort_values(["class", "hubness"], ascending=[True, False])
        # Keep the labeled version readable: all overlaps for these modules are usually small; cap if needed.
        label_rows = overlap.head(55)
        labels = {r.gene_id: str(r.gene_symbol) if pd.notna(r.gene_symbol) and str(r.gene_symbol).strip() else r.gene_id for r in label_rows.itertuples()}
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=7, font_color="#111827", ax=ax)
    module_size = len(G.nodes)
    tier2_n = int(sum(1 for n in G.nodes if n in tier2))
    ml_n = int(sum(1 for n in G.nodes if n in high_mod))
    fdr = enrich_row.get("fdr_BH", np.nan)
    fdr_txt = "NA" if pd.isna(fdr) else (f"{float(fdr):.2e}" if float(fdr) < 0.001 else f"{float(fdr):.3g}")
    title = f"{region} {module}"
    subtitle = f"n={module_size} | Tier2={tier2_n} | ML={ml_n} | best enrichment FDR={fdr_txt}"
    ax.set_title(f"{title}\n{subtitle}", fontsize=10 if summary_ax is not None else 15, fontweight="bold", pad=8)
    ax.set_axis_off()


def add_legend(fig, ax) -> None:
    handles = [
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor=NODE_COLORS["default"], markeredgecolor="#f8fafc", label="default module gene", markersize=9),
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor=NODE_COLORS["tier2"], markeredgecolor="#8a4b00", label="Tier2 overlap gene", markersize=9),
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor=NODE_COLORS["ml_high"], markeredgecolor="#8a4b00", label="high-confidence ML gene", markersize=9),
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor=NODE_COLORS["ml_high_moderate"], markeredgecolor="#8a4b00", label="high + moderate ML gene", markersize=9),
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor="#d9e2ec", markeredgecolor="#111827", markeredgewidth=1.8, label="hub-like overlap gene", markersize=11),
        Line2D([0], [0], color=EDGE_COLOR, alpha=0.35, lw=1.2, label="sparse correlation backbone"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=True, bbox_to_anchor=(0.5, 0.01), fontsize=9)


def save_one_module(region: str, module: str, G: nx.Graph, pos, genes, enrich_row, tier2, high, high_mod) -> dict:
    outputs = {}
    for labeled in [False, True]:
        stem = ("full_module_network_labeled_overlaps" if labeled else "full_module_network") + f"_{safe_name(region)}_{safe_name(module)}"
        fig, ax = plt.subplots(figsize=(10.8, 9.2))
        draw_network(G, pos, genes, region, module, enrich_row, tier2, high, high_mod, labeled=labeled, summary_ax=ax)
        add_legend(fig, ax)
        fig.text(
            0.5,
            0.075,
            "All module genes are shown as nodes. Edges are sparsified using a connectivity backbone plus strongest local correlations.",
            ha="center",
            fontsize=9,
            color="#374151",
        )
        fig.tight_layout(rect=[0.02, 0.10, 0.98, 0.98])
        png = OUT_DIR / f"{stem}.png"
        pdf = OUT_DIR / f"{stem}.pdf"
        fig.savefig(png, dpi=PLOT_DPI, bbox_inches="tight")
        fig.savefig(pdf, bbox_inches="tight")
        plt.close(fig)
        outputs["labeled" if labeled else "unlabeled"] = png.name
    return outputs


def save_summary(graphs: dict[tuple[str, str], tuple[nx.Graph, dict, pd.DataFrame, pd.Series]], tier2, high, high_mod) -> None:
    fig, axes = plt.subplots(3, 3, figsize=(16, 15))
    axes = axes.ravel()
    for ax, key in zip(axes, FOCUS_MODULES):
        G, pos, genes, enrich_row = graphs[key]
        draw_network(G, pos, genes, key[0], key[1], enrich_row, tier2, high, high_mod, labeled=False, summary_ax=ax)
    for ax in axes[len(FOCUS_MODULES):]:
        ax.set_axis_off()
    add_legend(fig, axes[0])
    fig.suptitle("Full-module WGCNA networks for prioritized modules", fontsize=18, fontweight="bold", y=0.995)
    fig.text(0.5, 0.045, "Every gene in each module is shown. Node size reflects intramodular correlation hubness; colors mark Tier2/ML evidence.", ha="center", fontsize=10, color="#374151")
    fig.tight_layout(rect=[0.02, 0.08, 0.98, 0.97])
    for ext in ["png", "pdf"]:
        fig.savefig(OUT_DIR / f"full_module_network_summary_top_modules.{ext}", dpi=PLOT_DPI, bbox_inches="tight")
    plt.close(fig)


def update_readme_and_audit(new_rows: list[dict]) -> None:
    readme_path = OUT_DIR / "README_cluster_plots.md"
    audit_path = OUT_DIR / "cluster_visualization_plot_inventory.csv"
    # Update audit CSV in the current explicit generated/skipped schema.
    audit = pd.read_csv(audit_path) if audit_path.exists() else pd.DataFrame()
    for row in new_rows:
        audit = pd.concat([audit, pd.DataFrame([row])], ignore_index=True)
    audit.to_csv(audit_path, index=False)

    png_count = len(list(OUT_DIR.glob("*.png")))
    pdf_count = len(list(OUT_DIR.glob("*.pdf")))
    csv_count = len(list(OUT_DIR.glob("*.csv")))
    readme_count = len(list(OUT_DIR.glob("README*.md")))
    text = readme_path.read_text(encoding="utf-8") if readme_path.exists() else "# Corrected ML-Annotated WGCNA Network Figures\n\n"
    text = text.replace(
        "The plot audit found 32 PNG plots, 32 matching PDF plots, two CSV summary files, and one README. PCA/module-clustering plots were generated successfully. t-SNE was attempted but skipped because the required sklearn functionality failed in this environment, so no t-SNE plots are included.",
        README_AUDIT_SENTENCE,
    )
    section = f"""\n## Full-Module Network Addendum\n\nFull-module network plots were added for the seven prioritized modules requested in the follow-up analysis. Unlike the earlier hub-only networks, these figures keep every gene in the module as a node and sparsify only the edges. Each module has one unlabeled network and one companion version labeling only Tier2/ML overlap genes. A multi-panel summary figure is also included.\n\nAfter adding the full-module network figures, the folder contains {png_count} PNG plots, {pdf_count} matching PDF plots, {csv_count} CSV summary/audit files, and {readme_count} README. No t-SNE image files are included.\n\n"""
    if "## Full-Module Network Addendum" not in text:
        insert_before = "## Most Presentation-Friendly Plots\n"
        text = text.replace(insert_before, section + insert_before)
    readme_path.write_text(text, encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    assignments = read_csv(WGCNA_ROOT / "best_solution_gene_to_module_assignments.csv")
    enrich = to_num(read_csv(ENRICHMENT_FILE), ["overlap_count", "p_value", "fdr_BH"])
    tier2, high, high_mod = build_gene_sets()
    expr_cache: dict[str, pd.DataFrame] = {}
    graphs = {}
    audit_rows = []
    for region, module in FOCUS_MODULES:
        expr_cache.setdefault(region, load_expression(region))
        module_genes = assignments[(assignments["region"] == region) & (assignments["module"] == module)].copy()
        if module_genes.empty:
            raise ValueError(f"No genes found for {region} {module}")
        G, pos, genes, median_edge = compute_graph(module_genes, expr_cache[region])
        enrich_row = best_enrichment(region, module, enrich)
        graphs[(region, module)] = (G, pos, genes, enrich_row)
        save_one_module(region, module, G, pos, genes, enrich_row, tier2, high, high_mod)
        for plot_type in ["full_module_network", "full_module_network_labeled_overlaps"]:
            audit_rows.append({
                "plot_type": plot_type,
                "region": region,
                "png_created": "TRUE",
                "pdf_created": "TRUE",
                "skipped": "FALSE",
                "skip_reason": "",
            })
        print(f"{region} {module}: nodes={G.number_of_nodes()} edges={G.number_of_edges()} median_edge_abs_cor={median_edge:.3f}")
    save_summary(graphs, tier2, high, high_mod)
    audit_rows.append({
        "plot_type": "full_module_network_summary_top_modules",
        "region": "all selected modules",
        "png_created": "TRUE",
        "pdf_created": "TRUE",
        "skipped": "FALSE",
        "skip_reason": "",
    })
    update_readme_and_audit(audit_rows)
    print(f"Wrote full-module network figures to {OUT_DIR}")


if __name__ == "__main__":
    main()
