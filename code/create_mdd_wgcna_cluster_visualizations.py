#!/usr/bin/env python
"""Create report-friendly cluster/module visualizations for Step 12 MDD WGCNA outputs.

This script works from saved Step6/Step12 CSV outputs. It does not rerun WGCNA.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

PROJECT_ROOT = Path(r"D:/genes/GeneDepresion")
PIPELINE_ROOT = PROJECT_ROOT / "outputs" / "length_gt_100_clean_pipeline"
WGCNA_ROOT = PIPELINE_ROOT / "12_MDD_WGCNA_modules"
OUT_DIR = WGCNA_ROOT / "plots_cluster_visualization"
SCRIPT_PATH_FOR_README = "code/create_mdd_wgcna_cluster_visualizations.py"

REGION_TO_EXPR = {
    "Frontal Cortex": PIPELINE_ROOT / "02_normalized_expression" / "normalized_GeTMM_log2_frontal_cortex.csv",
    "Hippocampal Formation": PIPELINE_ROOT / "02_normalized_expression" / "normalized_GeTMM_log2_hippocampal_formation.csv",
    "Insula": PIPELINE_ROOT / "02_normalized_expression" / "normalized_GeTMM_log2_insula.csv",
    "Ventral Striatum": PIPELINE_ROOT / "02_normalized_expression" / "normalized_GeTMM_log2_ventral_striatum.csv",
}

HIGHLIGHT_AHBA = ["M26", "M28", "M2", "M6", "M7"]
PLOT_DPI = 220
sns.set_theme(style="whitegrid", context="talk")


def safe_name(value: str) -> str:
    return (
        str(value)
        .replace("/", "_")
        .replace("\\", "_")
        .replace(" ", "_")
        .replace(";", "_")
        .replace(":", "_")
    )


def read_csv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, **kwargs)


def numeric(df: pd.DataFrame, cols: Iterable[str]) -> pd.DataFrame:
    for col in cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def save_fig(fig, stem: str, plot_type: str, region: str, description: str, inventory: list[dict], module: str = "") -> None:
    png = OUT_DIR / f"{stem}.png"
    pdf = OUT_DIR / f"{stem}.pdf"
    fig.savefig(png, dpi=PLOT_DPI, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    inventory.append({
        "plot_file_png": png.name,
        "plot_file_pdf": pdf.name,
        "plot_type": plot_type,
        "region": region,
        "module": module,
        "description": description,
    })


def pval_label(x) -> str:
    if pd.isna(x):
        return "NA"
    x = float(x)
    if x == 0:
        return "0"
    if abs(x) < 0.001:
        return f"{x:.1e}"
    return f"{x:.3g}"


def get_module_color(module: str) -> str:
    if not isinstance(module, str):
        return "#9e9e9e"
    color = module[2:] if module.startswith("ME") else module
    color = color.lower()
    try:
        matplotlib.colors.to_rgba(color)
        return color
    except ValueError:
        return "#607d8b"


def load_expression(region: str) -> pd.DataFrame:
    path = REGION_TO_EXPR[region]
    expr = pd.read_csv(path, index_col=0)
    expr.index = expr.index.astype(str)
    expr.columns = expr.columns.astype(str)
    return expr


def module_eigengene_corr_plot(me_region: pd.DataFrame, region: str, inventory: list[dict]) -> None:
    module_cols = [c for c in me_region.columns if c.startswith("ME") and c != "MEgrey" and me_region[c].notna().any()]
    if len(module_cols) < 2:
        return
    mat = me_region[module_cols].dropna(axis=1, how="all").astype(float)
    corr = mat.corr()
    corr = corr.dropna(axis=0, how="all").dropna(axis=1, how="all")
    if corr.shape[0] < 2:
        return
    dist = 1 - corr
    np.fill_diagonal(dist.values, 0)
    condensed = squareform(dist.clip(lower=0).values, checks=False)
    row_link = linkage(condensed, method="average")
    cg = sns.clustermap(
        corr,
        row_linkage=row_link,
        col_linkage=row_link,
        cmap="vlag",
        center=0,
        vmin=-1,
        vmax=1,
        linewidths=0.1,
        figsize=(max(8, 0.35 * corr.shape[0]), max(8, 0.35 * corr.shape[0])),
        cbar_kws={"label": "Eigengene correlation"},
    )
    cg.fig.suptitle(f"{region}: module eigengene correlation", y=1.02, fontsize=18, fontweight="bold")
    stem = f"module_eigengene_correlation_heatmap_{safe_name(region)}"
    png = OUT_DIR / f"{stem}.png"
    pdf = OUT_DIR / f"{stem}.pdf"
    cg.fig.savefig(png, dpi=PLOT_DPI, bbox_inches="tight")
    cg.fig.savefig(pdf, bbox_inches="tight")
    plt.close(cg.fig)
    inventory.append({
        "plot_file_png": png.name,
        "plot_file_pdf": pdf.name,
        "plot_type": "module_eigengene_correlation_heatmap",
        "region": region,
        "module": "",
        "description": "Clustered heatmap showing correlations among non-grey module eigengenes.",
    })


def sample_projection_plots(me_region: pd.DataFrame, meta_region: pd.DataFrame, region: str, inventory: list[dict]) -> None:
    module_cols = [c for c in me_region.columns if c.startswith("ME") and c != "MEgrey" and me_region[c].notna().any()]
    dat = me_region[["SampleID"] + module_cols].copy()
    dat = dat.dropna(axis=1, how="all")
    module_cols = [c for c in dat.columns if c.startswith("ME")]
    dat = dat.dropna(subset=module_cols)
    if len(dat) < 5 or len(module_cols) < 2:
        return
    merged = dat.merge(meta_region, on="SampleID", how="left", suffixes=("", "_meta"))
    diag_col = "Diagnosis_clean" if "Diagnosis_clean" in merged.columns else ("Diagnosis" if "Diagnosis" in merged.columns else None)
    sex_col = "Sex" if "Sex" in merged.columns else ("gender" if "gender" in merged.columns else None)
    X = StandardScaler().fit_transform(merged[module_cols].astype(float))
    pca = PCA(n_components=2, random_state=1)
    pcs = pca.fit_transform(X)
    plot_df = pd.DataFrame({"PC1": pcs[:, 0], "PC2": pcs[:, 1]})
    plot_df["Diagnosis"] = merged[diag_col].fillna("Unknown") if diag_col else "Unknown"
    plot_df["Sex"] = merged[sex_col].fillna("Unknown") if sex_col else "Unknown"
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    sns.scatterplot(data=plot_df, x="PC1", y="PC2", hue="Diagnosis", style="Sex", s=80, alpha=0.88, ax=ax)
    ax.set_title(f"{region}: sample PCA using module eigengenes", fontsize=16, fontweight="bold")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")
    ax.legend(title="Sample metadata", bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
    save_fig(fig, f"sample_PCA_module_eigengenes_{safe_name(region)}", "sample_PCA_module_eigengenes", region, "PCA of samples using module eigengene values, colored by diagnosis.", inventory)

    perplexity = max(5, min(20, (len(merged) - 1) // 3))
    if len(merged) > perplexity + 2:
        try:
            tsne = TSNE(n_components=2, init="pca", learning_rate="auto", perplexity=perplexity, random_state=1)
            emb = tsne.fit_transform(X)
        except Exception as exc:
            inventory.append({
                "plot_file_png": "",
                "plot_file_pdf": "",
                "plot_type": "sample_tSNE_module_eigengenes_skipped",
                "region": region,
                "module": "",
                "description": f"t-SNE skipped because sklearn failed in this environment: {type(exc).__name__}: {str(exc)[:140]}",
            })
            return
        tsne_df = pd.DataFrame({"tSNE1": emb[:, 0], "tSNE2": emb[:, 1], "Diagnosis": plot_df["Diagnosis"], "Sex": plot_df["Sex"]})
        fig, ax = plt.subplots(figsize=(8.5, 6.5))
        sns.scatterplot(data=tsne_df, x="tSNE1", y="tSNE2", hue="Diagnosis", style="Sex", s=80, alpha=0.88, ax=ax)
        ax.set_title(f"{region}: sample t-SNE using module eigengenes", fontsize=16, fontweight="bold")
        ax.set_xlabel("t-SNE 1")
        ax.set_ylabel("t-SNE 2")
        ax.legend(title="Sample metadata", bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
        save_fig(fig, f"sample_tSNE_module_eigengenes_{safe_name(region)}", "sample_tSNE_module_eigengenes", region, "t-SNE of samples using module eigengene values, colored by diagnosis.", inventory)


def summary_bubbleplot(region_triage: pd.DataFrame, enrich: pd.DataFrame, region: str, inventory: list[dict]) -> None:
    if region_triage.empty:
        return
    modules = region_triage.sort_values(["priority_score", "module_size"], ascending=[False, False]).head(30)["module"].tolist()
    plot_rows = []
    for _, row in region_triage[region_triage["module"].isin(modules)].iterrows():
        mod = row["module"]
        plot_rows.append({"module": mod, "metric": "Diagnosis raw p", "value": row["diagnosis_p"], "size": row["module_size"]})
        plot_rows.append({"module": mod, "metric": "Best Tier2/ML FDR", "value": row["best_Tier2_or_ML_FDR"], "size": row["module_size"]})
        best_highlight = np.nanmin([row.get(f"AHBA_{m}_FDR", np.nan) for m in HIGHLIGHT_AHBA])
        plot_rows.append({"module": mod, "metric": "Best AHBA highlight FDR", "value": best_highlight, "size": row["module_size"]})
        hm = enrich[(enrich.region == region) & (enrich.module == mod) & (enrich.gene_set.isin(["high_plus_moderate_ML_genes", "same_region_high_plus_moderate_ML_genes"]))]
        hm_fdr = hm["fdr_BH"].min() if not hm.empty else np.nan
        plot_rows.append({"module": mod, "metric": "High+moderate ML FDR", "value": hm_fdr, "size": row["module_size"]})
    df = pd.DataFrame(plot_rows)
    df["neglog10"] = -np.log10(pd.to_numeric(df["value"], errors="coerce").clip(lower=1e-300))
    df["module"] = pd.Categorical(df["module"], categories=modules[::-1], ordered=True)
    fig, ax = plt.subplots(figsize=(10.5, max(7, 0.32 * len(modules))))
    sns.scatterplot(data=df, x="metric", y="module", size="size", hue="neglog10", palette="viridis", sizes=(40, 360), edgecolor="0.25", linewidth=0.4, ax=ax)
    ax.set_title(f"{region}: module evidence summary", fontsize=16, fontweight="bold")
    ax.set_xlabel("Evidence metric")
    ax.set_ylabel("Module")
    ax.tick_params(axis="x", rotation=25)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0, title="Module size / -log10")
    save_fig(fig, f"module_summary_bubbleplot_{safe_name(region)}", "module_summary_bubbleplot", region, "Bubble plot combining module size, diagnosis p-value, Tier2/ML enrichment, and highlighted AHBA overlap.", inventory)


def expression_heatmap(expr: pd.DataFrame, meta_region: pd.DataFrame, hubs: pd.DataFrame, region: str, module: str, inventory: list[dict]) -> None:
    hub = hubs[(hubs.region == region) & (hubs.module == module)].copy()
    if hub.empty:
        return
    gene_ids = hub.sort_values("abs_kME", ascending=False)["gene_id"].astype(str).drop_duplicates().head(30).tolist()
    gene_symbols = hub.drop_duplicates("gene_id").set_index(hub["gene_id"].astype(str))["gene_symbol"].to_dict()
    gene_ids = [g for g in gene_ids if g in expr.index]
    if len(gene_ids) < 3:
        return
    sample_ids = [c for c in expr.columns if c in set(meta_region["SampleID"].astype(str))]
    if len(sample_ids) < 5:
        return
    mat = expr.loc[gene_ids, sample_ids].astype(float)
    mat = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    labels = [str(gene_symbols.get(g, g)) for g in mat.index]
    mat.index = labels
    meta = meta_region.drop_duplicates("SampleID").set_index("SampleID")
    diag = meta.reindex(sample_ids)["Diagnosis_clean"].fillna(meta.reindex(sample_ids).get("Diagnosis", "Unknown")) if "Diagnosis_clean" in meta.columns else pd.Series("Unknown", index=sample_ids)
    palette = {"MDD": "#d95f02", "Control": "#1b9e77", "CTRL": "#1b9e77", "Unknown": "#999999"}
    col_colors = diag.map(lambda x: palette.get(str(x), "#999999"))
    cg = sns.clustermap(
        mat,
        cmap="vlag",
        center=0,
        col_colors=col_colors,
        xticklabels=False,
        yticklabels=True,
        figsize=(12, max(7, 0.28 * len(labels))),
        cbar_kws={"label": "Row z-score"},
    )
    cg.fig.suptitle(f"{region} {module}: top hub-gene expression", y=1.02, fontsize=18, fontweight="bold")
    handles = [Line2D([0], [0], marker="s", linestyle="", color=color, label=label) for label, color in palette.items() if label in set(diag.astype(str))]
    if handles:
        cg.ax_heatmap.legend(handles=handles, title="Diagnosis", loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=max(1, len(handles)), frameon=True)
    cg.ax_heatmap.tick_params(axis="y", labelsize=9)
    stem = f"top_module_expression_heatmap_{safe_name(region)}_{safe_name(module)}"
    png = OUT_DIR / f"{stem}.png"
    pdf = OUT_DIR / f"{stem}.pdf"
    cg.fig.savefig(png, dpi=PLOT_DPI, bbox_inches="tight")
    cg.fig.savefig(pdf, bbox_inches="tight")
    plt.close(cg.fig)
    inventory.append({
        "plot_file_png": png.name,
        "plot_file_pdf": pdf.name,
        "plot_type": "top_module_expression_heatmap",
        "region": region,
        "module": module,
        "description": "Clustered heatmap of row-standardized expression for the top hub genes in a prioritized module.",
    })


def network_plot(expr: pd.DataFrame, hubs: pd.DataFrame, enrich: pd.DataFrame, triage: pd.DataFrame, region: str, module: str, inventory: list[dict]) -> None:
    hub = hubs[(hubs.region == region) & (hubs.module == module)].copy()
    if hub.empty:
        return
    gene_ids = hub.sort_values("abs_kME", ascending=False)["gene_id"].astype(str).drop_duplicates().head(25).tolist()
    gene_ids = [g for g in gene_ids if g in expr.index]
    if len(gene_ids) < 5:
        return
    mat = expr.loc[gene_ids].astype(float).T
    corr = mat.corr().abs()
    edges = []
    for i, g1 in enumerate(gene_ids):
        for g2 in gene_ids[i + 1:]:
            w = float(corr.loc[g1, g2])
            if not np.isnan(w):
                edges.append((g1, g2, w))
    edges = sorted(edges, key=lambda x: x[2], reverse=True)[:55]
    G = nx.Graph()
    symbol_map = hub.drop_duplicates("gene_id").set_index(hub["gene_id"].astype(str))["gene_symbol"].to_dict()
    kme_map = hub.drop_duplicates("gene_id").set_index(hub["gene_id"].astype(str))["abs_kME"].to_dict()
    for g in gene_ids:
        G.add_node(g, label=str(symbol_map.get(g, g)), kme=float(kme_map.get(g, 0)))
    for g1, g2, w in edges:
        if w >= 0.45 or len(G.edges) < 20:
            G.add_edge(g1, g2, weight=w)
    overlap_gene_ids = set()
    e = enrich[(enrich.region == region) & (enrich.module == module) & (enrich.gene_set.isin(["same_region_Tier2_genes", "high_confidence_ML_genes", "high_plus_moderate_ML_genes", "same_region_high_plus_moderate_ML_genes"]))]
    for genes in e.get("overlap_genes", pd.Series(dtype=str)).fillna(""):
        overlap_gene_ids.update([x.strip() for x in str(genes).split(";") if x.strip()])
    pos = nx.spring_layout(G, seed=7, weight="weight", k=0.65)
    fig, ax = plt.subplots(figsize=(9, 8))
    widths = [0.6 + 2.8 * G[u][v]["weight"] for u, v in G.edges]
    nx.draw_networkx_edges(G, pos, width=widths, alpha=0.28, edge_color="#455a64", ax=ax)
    node_sizes = [260 + 900 * G.nodes[n].get("kme", 0) for n in G.nodes]
    node_colors = ["#ffcc80" if n in overlap_gene_ids else "#90caf9" for n in G.nodes]
    edgecolors = ["#bf360c" if n in overlap_gene_ids else "#263238" for n in G.nodes]
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, edgecolors=edgecolors, linewidths=1.8, ax=ax)
    labels = {n: G.nodes[n]["label"] for n in G.nodes}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, ax=ax)
    ax.set_title(f"{region} {module}: top hub-gene correlation network", fontsize=15, fontweight="bold")
    ax.set_axis_off()
    legend = [
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor="#90caf9", markeredgecolor="#263238", label="Hub gene"),
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor="#ffcc80", markeredgecolor="#bf360c", label="Tier2/ML overlap gene"),
    ]
    ax.legend(handles=legend, loc="lower left", frameon=True)
    save_fig(fig, f"module_network_{safe_name(region)}_{safe_name(module)}", "module_network", region, "Small correlation network for top hub genes; highlighted nodes overlap Tier2/ML sets when available.", inventory, module=module)


def build_best_modules_summary(triage: pd.DataFrame, enrich: pd.DataFrame, ahba: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in triage.iterrows():
        region, module = row["region"], row["module"]
        e = enrich[(enrich.region == region) & (enrich.module == module)]
        def ep(set_name: str, col: str):
            sub = e[e.gene_set == set_name]
            return sub[col].min() if not sub.empty else np.nan
        ah = ahba[(ahba.region == region) & (ahba.mdd_module == module)].sort_values(["fdr_BH", "p_value", "overlap_count"], ascending=[True, True, False])
        best_ah = ah.iloc[0] if not ah.empty else pd.Series(dtype=object)
        rows.append({
            "region": region,
            "module": module,
            "module_size": row["module_size"],
            "diagnosis_beta": row["diagnosis_beta"],
            "diagnosis_p": row["diagnosis_p"],
            "diagnosis_fdr": row["diagnosis_FDR"],
            "Tier2_enrichment_p": min(ep("all_Tier2_genes", "p_value"), ep("same_region_Tier2_genes", "p_value")),
            "Tier2_enrichment_FDR": min(ep("all_Tier2_genes", "fdr_BH"), ep("same_region_Tier2_genes", "fdr_BH")),
            "high_confidence_ML_enrichment_p": ep("high_confidence_ML_genes", "p_value"),
            "high_confidence_ML_enrichment_FDR": ep("high_confidence_ML_genes", "fdr_BH"),
            "high_plus_moderate_ML_enrichment_p": min(ep("high_plus_moderate_ML_genes", "p_value"), ep("same_region_high_plus_moderate_ML_genes", "p_value")),
            "high_plus_moderate_ML_enrichment_FDR": min(ep("high_plus_moderate_ML_genes", "fdr_BH"), ep("same_region_high_plus_moderate_ML_genes", "fdr_BH")),
            "best_AHBA_overlap_module": best_ah.get("ahba_module", ""),
            "best_AHBA_overlap_p": best_ah.get("p_value", np.nan),
            "best_AHBA_overlap_FDR": best_ah.get("fdr_BH", np.nan),
            "top_hub_genes": row["top_10_hub_genes"],
            "priority_score": row["priority_score"],
        })
    return pd.DataFrame(rows).sort_values(["priority_score", "diagnosis_p"], ascending=[False, True])


def write_readme(best_summary: pd.DataFrame, inventory: pd.DataFrame, triage: pd.DataFrame) -> None:
    lines = []
    lines.append("# MDD WGCNA Cluster Visualization Plots")
    lines.append("")
    lines.append("These plots were generated from the saved Step 12 MDD WGCNA outputs. The full WGCNA pipeline was not rerun.")
    lines.append("")
    lines.append("## What Was Created")
    lines.append("")
    lines.append("- Module eigengene correlation heatmaps: show which modules have similar sample-level eigengene patterns.")
    lines.append("- Sample PCA plots: show whether samples separate in module-eigengene space by Diagnosis. t-SNE was attempted but skipped because scikit-learn failed in this Windows/Anaconda environment.")
    lines.append("- Top-module expression heatmaps: show row-standardized expression of top hub genes across samples.")
    lines.append("- Top-module hub networks: show small correlation networks among hub genes, with Tier2/ML overlap genes highlighted when present.")
    lines.append("- Module summary bubbleplots: combine module size, Diagnosis evidence, Tier2/ML enrichment, and highlighted AHBA overlap.")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append(f"- Script: `{SCRIPT_PATH_FOR_README}`")
    lines.append("- Plot inventory: `cluster_visualization_plot_inventory.csv`")
    lines.append("- Best modules table: `best_modules_summary_for_plotting.csv`")
    lines.append("")
    lines.append("## Important Limitations")
    lines.append("")
    lines.append("- MEgrey is excluded from the main biological plots because it is the unassigned gene bin.")
    lines.append("- UMAP was not available in the current Python environment. t-SNE was attempted as a fallback, but scikit-learn failed with a threadpoolctl error, so only PCA sample projections are included.")
    lines.append("- A true cleaned or zoomed WGCNA gene dendrogram cannot be recreated from the saved CSVs alone because the original gene dendrogram/TOM object was not saved. The existing Step6 dendrogram PNGs remain available in the parent `plots/` folder.")
    lines.append("- No module eigengene survived BH-FDR < 0.10 for Diagnosis. Raw Diagnosis p < 0.05 should be treated as exploratory.")
    lines.append("")
    lines.append("## Most Presentation-Friendly Plots")
    lines.append("")
    lines.append("1. `module_summary_bubbleplot_Ventral_Striatum.png` for the clearest compact evidence summary.")
    lines.append("2. `top_module_expression_heatmap_Ventral_Striatum_MElavenderblush3.png` for the strongest multi-evidence module.")
    lines.append("3. `module_network_Ventral_Striatum_MElavenderblush3.png` for a simple hub-gene module view.")
    lines.append("4. `sample_PCA_module_eigengenes_*` plots to show that module eigengenes do not produce a clean corrected Diagnosis separation.")
    lines.append("5. `module_eigengene_correlation_heatmap_*` plots to explain relationships among discovered modules.")
    lines.append("")
    lines.append("## Region Notes")
    lines.append("")
    for region, sub in best_summary.groupby("region", sort=True):
        top = sub.sort_values(["priority_score", "diagnosis_p"], ascending=[False, True]).head(5)
        raw_dx = triage[(triage.region == region) & (pd.to_numeric(triage.diagnosis_p, errors="coerce") < 0.05)]
        tier_ml = sub[pd.to_numeric(sub["Tier2_enrichment_FDR"], errors="coerce") < 0.05]
        ml = sub[pd.to_numeric(sub["high_plus_moderate_ML_enrichment_FDR"], errors="coerce") < 0.10]
        ah = sub[pd.to_numeric(sub["best_AHBA_overlap_FDR"], errors="coerce") < 0.05]
        lines.append(f"### {region}")
        lines.append("")
        lines.append("Top plotted modules: " + ", ".join([f"{r.module} (score {r.priority_score})" for r in top.itertuples()]) + ".")
        lines.append(f"Diagnosis signal: {len(raw_dx)} module(s) had raw p < 0.05, but none survived Diagnosis BH-FDR < 0.10.")
        lines.append(f"Tier2 enrichment: {len(tier_ml)} module(s) had Tier2 FDR < 0.05.")
        lines.append(f"High+moderate ML enrichment: {len(ml)} module(s) had FDR < 0.10.")
        lines.append(f"AHBA overlap: {len(ah)} module(s) had best AHBA overlap FDR < 0.05; highlighted M26/M28/M2/M6/M7 evidence is shown in the bubbleplots and triage table.")
        lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("The plots support prioritized exploratory follow-up rather than a clean MDD/control module-level separation. Ventral Striatum modules, especially MElavenderblush3, are the most biologically interesting because they combine raw Diagnosis evidence, Tier2 enrichment, highlighted AHBA overlap, and interpretable hub genes.")
    (OUT_DIR / "README_cluster_plots.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    trait = numeric(read_csv(WGCNA_ROOT / "module_trait_association_results.csv"), ["n_module_genes", "beta_diagnosis", "p_diagnosis", "fdr_diagnosis_BH"])
    enrich = numeric(read_csv(WGCNA_ROOT / "module_enrichment_for_Tier2_and_ML_genes.csv"), ["overlap_count", "p_value", "fdr_BH"])
    ahba = numeric(read_csv(WGCNA_ROOT / "mdd_modules_vs_published_AHBA_modules_enrichment.csv"), ["overlap_count", "p_value", "fdr_BH"])
    hubs = numeric(read_csv(WGCNA_ROOT / "module_hub_genes_top50.csv"), ["abs_kME"])
    eig = read_csv(WGCNA_ROOT / "best_solution_module_eigengenes.csv")
    meta = read_csv(WGCNA_ROOT / "best_solution_sample_metadata_used.csv")
    triage = numeric(read_csv(WGCNA_ROOT / "mdd_wgcna_module_triage_table.csv"), ["module_size", "diagnosis_beta", "diagnosis_p", "diagnosis_FDR", "best_Tier2_or_ML_FDR", "priority_score", "best_AHBA_FDR"])
    triage = triage[triage.module != "MEgrey"].copy()
    best_summary = build_best_modules_summary(triage, enrich, ahba)
    best_summary.to_csv(OUT_DIR / "best_modules_summary_for_plotting.csv", index=False)

    inventory: list[dict] = []
    for region in sorted(triage.region.unique()):
        me_region = eig[eig.region == region].copy()
        meta_region = meta[meta.region == region].copy()
        meta_region["SampleID"] = meta_region["SampleID"].astype(str)
        module_eigengene_corr_plot(me_region, region, inventory)
        sample_projection_plots(me_region, meta_region, region, inventory)
        summary_bubbleplot(triage[triage.region == region].copy(), enrich, region, inventory)
        expr = load_expression(region)
        top_modules = triage[triage.region == region].sort_values(["priority_score", "diagnosis_p"], ascending=[False, True]).head(3)["module"].tolist()
        network_modules = top_modules[:2]
        for module in top_modules:
            expression_heatmap(expr, meta_region, hubs, region, module, inventory)
        for module in network_modules:
            network_plot(expr, hubs, enrich, triage, region, module, inventory)

    inv = pd.DataFrame(inventory)
    inv.to_csv(OUT_DIR / "cluster_visualization_plot_inventory.csv", index=False)
    write_readme(best_summary, inv, triage)
    print(f"Created {len(inv)} plot pairs in {OUT_DIR}")
    print(f"Wrote {OUT_DIR / 'best_modules_summary_for_plotting.csv'}")
    print(f"Wrote {OUT_DIR / 'cluster_visualization_plot_inventory.csv'}")
    print(f"Wrote {OUT_DIR / 'README_cluster_plots.md'}")


if __name__ == "__main__":
    main()
