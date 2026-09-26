#!/usr/bin/env python
"""Create follow-up tables and plots for frontal subregion WGCNA results."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


PROJECT_ROOT = Path("D:/genes/GeneDepresion")
PIPELINE_ROOT = PROJECT_ROOT / "outputs" / "length_gt_100_clean_pipeline"
SUB_ROOT = PIPELINE_ROOT / "mdd_wgcna_modules" / "frontal_subregion_wgcna"
BROAD_ROOT = PIPELINE_ROOT / "12_MDD_WGCNA_modules"
PLOT_DIR = SUB_ROOT / "plots"
ANNOT_DIR = SUB_ROOT / "biological_annotation"
INPUT_DIR = PIPELINE_ROOT / "02_normalized_expression" / "frontal_subregion_wgcna_inputs"
PLOT_DIR.mkdir(parents=True, exist_ok=True)
ANNOT_DIR.mkdir(parents=True, exist_ok=True)

HIGHLIGHT_AHBA = ["M26", "M28", "M2", "M6", "M7"]

EXPR_FILES = {
    "Cingulate gyrus 25": INPUT_DIR / "normalized_GeTMM_log2_CG25.csv",
    "dlPFC BA8/9": INPUT_DIR / "normalized_GeTMM_log2_dlPFC_BA8_9.csv",
    "OFC BA11": INPUT_DIR / "normalized_GeTMM_log2_OFC_BA11.csv",
}


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def num(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    df = df.copy()
    for col in cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def safe_name(text: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in text).strip("_")


def load_inputs() -> dict[str, pd.DataFrame]:
    return {
        "param": read_csv(SUB_ROOT / "frontal_subregion_wgcna_parameter_sensitivity.csv"),
        "assign": read_csv(SUB_ROOT / "frontal_subregion_best_solution_gene_to_module_assignments.csv"),
        "sizes": num(read_csv(SUB_ROOT / "frontal_subregion_module_size_table.csv"), ["n_genes"]),
        "trait": num(read_csv(SUB_ROOT / "frontal_subregion_module_trait_association_results.csv"), ["n_module_genes", "beta_diagnosis", "p_diagnosis", "fdr_diagnosis_BH"]),
        "enrich": num(read_csv(SUB_ROOT / "frontal_subregion_Tier2_ML_enrichment.csv"), ["overlap_count", "p_value", "fdr_BH"]),
        "ahba": num(read_csv(SUB_ROOT / "frontal_subregion_AHBA_overlap.csv"), ["overlap_count", "p_value", "fdr_BH"]),
        "hubs": num(read_csv(SUB_ROOT / "frontal_subregion_hub_genes_top50.csv"), ["kME", "abs_kME"]),
        "important": num(read_csv(SUB_ROOT / "frontal_subregion_important_modules_summary.csv"), ["module_size", "diagnosis_p_value", "diagnosis_BH_FDR"]),
        "tier2": read_csv(PIPELINE_ROOT / "07_Tier2_candidate_gene_sets" / "Tier2_candidate_genes_all_regions.csv"),
        "ml": read_csv(PIPELINE_ROOT / "10_ML_prioritized_Tier2_corrected" / "ML_prioritized_Tier2_genes_all_ranked.csv"),
    }


def best_row(df: pd.DataFrame, value_col: str, module_col: str = "module") -> tuple[str, float]:
    usable = df.dropna(subset=[value_col]).copy()
    usable = usable[~usable[module_col].astype(str).eq("MEgrey")]
    if usable.empty:
        return "", np.nan
    row = usable.sort_values(value_col).iloc[0]
    return str(row[module_col]), float(row[value_col])


def make_comparison(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []

    def one_level(label: str, region: str, root: Path, subregion: bool) -> dict[str, object]:
        if subregion:
            param = data["param"][data["param"]["region"].eq(region)]
            trait = data["trait"][data["trait"]["region"].eq(region)]
            enrich = data["enrich"][data["enrich"]["region"].eq(region)]
            ahba = data["ahba"][data["ahba"]["region"].eq(region)]
            sizes = data["sizes"][data["sizes"]["region"].eq(region)]
            important = data["important"][data["important"]["region"].eq(region)]
        else:
            param = read_csv(root / "mdd_wgcna_parameter_sensitivity.csv")
            trait = num(read_csv(root / "module_trait_association_results.csv"), ["p_diagnosis", "fdr_diagnosis_BH", "n_samples"])
            enrich = num(read_csv(root / "module_enrichment_for_Tier2_and_ML_genes.csv"), ["fdr_BH"])
            ahba = num(read_csv(root / "mdd_modules_vs_published_AHBA_modules_enrichment.csv"), ["fdr_BH", "overlap_count"])
            sizes = num(read_csv(root / "best_solution_module_size_table.csv"), ["n_genes"])
            important = read_csv(root / "important_modules_summary.csv")
            param = param[param["region"].eq("Frontal Cortex")]
            trait = trait[trait["region"].eq("Frontal Cortex")]
            enrich = enrich[enrich["region"].eq("Frontal Cortex")]
            ahba = ahba[ahba["region"].eq("Frontal Cortex")]
            sizes = sizes[sizes["region"].eq("Frontal Cortex")]
            important = important[important["region"].eq("Frontal Cortex")]

        best_diag_module, best_diag_p = best_row(trait, "p_diagnosis")
        best_diag_fdr = trait.loc[trait["module"].eq(best_diag_module), "fdr_diagnosis_BH"].min() if best_diag_module else np.nan
        best_enrich_module, best_enrich_fdr = best_row(enrich, "fdr_BH")
        ahba_hi = ahba[ahba["ahba_module"].isin(HIGHLIGHT_AHBA)] if "ahba_module" in ahba.columns else ahba
        best_ahba_module, best_ahba_fdr = best_row(ahba_hi.rename(columns={"mdd_module": "module"}), "fdr_BH")
        n_samples = int(pd.to_numeric(param["n_samples"], errors="coerce").max()) if not param.empty else np.nan
        n_genes = int(pd.to_numeric(param["n_genes_after_filter"], errors="coerce").max()) if not param.empty else np.nan
        n_modules = int(sizes[~sizes["module"].eq("MEgrey")].shape[0]) if not sizes.empty else 0
        theme = ""
        if not important.empty and "interpretation" in important.columns:
            theme = str(important[important.get("important_module", False).astype(str).eq("TRUE")]["interpretation"].head(1).squeeze() if any(important.get("important_module", False).astype(str).eq("TRUE")) else "")
        return {
            "analysis_level": label,
            "n_samples": n_samples,
            "n_genes_used": n_genes,
            "n_modules": n_modules,
            "best_diagnosis_module": best_diag_module,
            "best_diagnosis_p": best_diag_p,
            "best_diagnosis_FDR": best_diag_fdr,
            "best_Tier2_ML_module": best_enrich_module,
            "best_Tier2_ML_FDR": best_enrich_fdr,
            "best_AHBA_overlap_module": best_ahba_module,
            "best_AHBA_overlap_FDR": best_ahba_fdr,
            "top_biological_theme": theme,
            "interpretation_priority": "supplementary/exploratory" if label == "broad_frontal" else "primary frontal subregion",
            "notes": "Broad Frontal Cortex retained; not overwritten." if label == "broad_frontal" else "Subregion-specific WGCNA branch.",
        }

    rows.append(one_level("broad_frontal", "Frontal Cortex", BROAD_ROOT, False))
    rows.append(one_level("OFC_BA11", "OFC BA11", SUB_ROOT, True))
    rows.append(one_level("dlPFC_BA8_9", "dlPFC BA8/9", SUB_ROOT, True))
    rows.append(one_level("CG25", "Cingulate gyrus 25", SUB_ROOT, True))
    out = pd.DataFrame(rows)
    out.to_csv(SUB_ROOT / "frontal_broad_vs_subregion_comparison.csv", index=False)
    return out


def prioritize_modules(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    imp = data["important"].copy()
    imp = imp[~imp["module"].eq("MEgrey")].copy()
    imp["important_bool"] = imp.get("important_module", "").astype(str).eq("TRUE")
    imp["score"] = 0
    imp.loc[imp["diagnosis_p_value"].lt(0.05), "score"] += 2
    imp.loc[imp["diagnosis_BH_FDR"].lt(0.10), "score"] += 3
    imp.loc[imp.get("enriched_for_Tier2", "").astype(str).eq("TRUE"), "score"] += 2
    imp.loc[imp.get("enriched_for_high_moderate_ML", "").astype(str).eq("TRUE"), "score"] += 2
    imp.loc[imp.get("AHBA_overlap_modules", "").fillna("").astype(str).ne(""), "score"] += 1
    imp.loc[imp["important_bool"], "score"] += 1
    return imp.sort_values(["score", "diagnosis_p_value", "module_size"], ascending=[False, True, False]).head(8)


def make_bubble_plot(data: dict[str, pd.DataFrame], prioritized: pd.DataFrame) -> None:
    trait = data["trait"].copy()
    trait = trait[~trait["module"].eq("MEgrey")].copy()
    enrich = data["enrich"].groupby(["region", "module"], as_index=False)["fdr_BH"].min().rename(columns={"fdr_BH": "best_Tier2_ML_FDR"})
    ahba = data["ahba"][data["ahba"]["ahba_module"].isin(HIGHLIGHT_AHBA)].groupby(["region", "mdd_module"], as_index=False)["fdr_BH"].min().rename(columns={"mdd_module": "module", "fdr_BH": "best_AHBA_FDR"})
    plot = trait.merge(enrich, on=["region", "module"], how="left").merge(ahba, on=["region", "module"], how="left")
    plot["label"] = plot["region"] + " " + plot["module"]
    plot["neglog10_diag"] = -np.log10(plot["p_diagnosis"].clip(lower=1e-300))
    plot["highlight"] = plot.set_index(["region", "module"]).index.isin(prioritized.set_index(["region", "module"]).index)
    plot = plot.sort_values(["highlight", "neglog10_diag"], ascending=[False, False]).head(35)

    fig, ax = plt.subplots(figsize=(10, max(6, len(plot) * 0.22)))
    colors = np.where(plot["highlight"], "#c0392b", "#4f6f8f")
    sizes = np.sqrt(plot["n_module_genes"].fillna(20).astype(float)) * 18
    ax.scatter(plot["neglog10_diag"], np.arange(len(plot)), s=sizes, c=colors, alpha=0.82, edgecolor="white", linewidth=0.7)
    ax.set_yticks(np.arange(len(plot)))
    ax.set_yticklabels(plot["label"], fontsize=8)
    ax.invert_yaxis()
    ax.axvline(-np.log10(0.05), color="#555555", linestyle="--", linewidth=1)
    ax.set_xlabel("-log10 Diagnosis raw p")
    ax.set_title("Frontal subregion WGCNA module summary")
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "frontal_subregion_module_summary_bubble_plot.png", dpi=220)
    fig.savefig(PLOT_DIR / "frontal_subregion_module_summary_bubble_plot.pdf")
    plt.close(fig)


def expression_for_region(region: str) -> pd.DataFrame:
    df = pd.read_csv(EXPR_FILES[region])
    df = df.rename(columns={df.columns[0]: "gene_id"})
    df["gene_id"] = df["gene_id"].astype(str)
    return df.set_index("gene_id")


def make_network_plot(region: str, module: str, data: dict[str, pd.DataFrame], labeled: bool) -> None:
    genes = data["assign"][(data["assign"]["region"].eq(region)) & (data["assign"]["module"].eq(module))].copy()
    if genes.empty:
        return
    expr = expression_for_region(region)
    gene_ids = [g for g in genes["gene_id"].astype(str) if g in expr.index]
    if len(gene_ids) < 10:
        return
    mat = expr.loc[gene_ids].T.astype(float)
    corr = mat.corr().abs().fillna(0)
    hubs = data["hubs"][(data["hubs"]["region"].eq(region)) & (data["hubs"]["module"].eq(module))].copy()
    hub_map = dict(zip(hubs["gene_id"].astype(str), hubs["abs_kME"].fillna(0)))
    symbol_map = dict(zip(genes["gene_id"].astype(str), genes["gene_symbol"].fillna(genes["gene_id"]).astype(str)))

    tier2_ids = set(data["tier2"]["gene_id"].astype(str))
    ml = data["ml"].copy()
    ml_cls = ml["ML_priority_class"].fillna("").astype(str)
    high_ids = set(ml[ml_cls.str.contains("High-confidence", case=False, na=False)]["gene_id"].astype(str))
    high_mod_ids = set(ml[ml_cls.str.contains("High-confidence|Moderate-confidence", case=False, na=False)]["gene_id"].astype(str))

    graph = nx.Graph()
    for gene in gene_ids:
        graph.add_node(gene)
    top_k = 2 if len(gene_ids) > 250 else 3
    for gene in gene_ids:
        strongest = corr[gene].drop(index=gene).sort_values(ascending=False).head(top_k)
        for other, weight in strongest.items():
            if weight >= 0.25:
                graph.add_edge(gene, other, weight=float(weight))
    if graph.number_of_edges() == 0:
        return

    pos = nx.spring_layout(graph, seed=22, k=1 / np.sqrt(max(10, len(graph))))
    node_colors = []
    node_sizes = []
    edgecolors = []
    for gene in graph.nodes:
        if gene in high_ids:
            node_colors.append("#d73027")
        elif gene in high_mod_ids:
            node_colors.append("#f6c431")
        elif gene in tier2_ids:
            node_colors.append("#f28e2b")
        else:
            node_colors.append("#9eb3c7")
        hubness = hub_map.get(gene, 0.25)
        node_sizes.append(12 + 95 * float(hubness))
        edgecolors.append("#202020" if (gene in tier2_ids or gene in high_mod_ids or gene in high_ids) and hubness >= 0.6 else "#ffffff")

    fig, ax = plt.subplots(figsize=(9, 8))
    widths = [0.2 + graph[u][v]["weight"] * 0.45 for u, v in graph.edges]
    nx.draw_networkx_edges(graph, pos, ax=ax, edge_color="#8a8f96", width=widths, alpha=0.25)
    nx.draw_networkx_nodes(graph, pos, ax=ax, node_color=node_colors, node_size=node_sizes, edgecolors=edgecolors, linewidths=0.45, alpha=0.94)
    if labeled:
        label_nodes = [g for g in graph.nodes if g in tier2_ids or g in high_mod_ids or g in high_ids]
        label_nodes = sorted(label_nodes, key=lambda g: hub_map.get(g, 0), reverse=True)[:30]
        nx.draw_networkx_labels(graph, pos, labels={g: symbol_map.get(g, g) for g in label_nodes}, font_size=6, ax=ax)
    module_size = len(genes)
    n_tier2 = len(set(gene_ids) & tier2_ids)
    n_ml = len(set(gene_ids) & high_mod_ids)
    best_fdr = data["enrich"][(data["enrich"]["region"].eq(region)) & (data["enrich"]["module"].eq(module))]["fdr_BH"].min()
    ax.set_title(f"{region} {module}\nmodule size={module_size}; Tier2={n_tier2}; ML high/mod={n_ml}; best enrichment FDR={best_fdr:.3g}", fontsize=11)
    ax.axis("off")
    fig.tight_layout()
    prefix = "full_module_network_labeled_overlaps" if labeled else "full_module_network"
    out_base = PLOT_DIR / f"{prefix}_{safe_name(region)}_{safe_name(module)}"
    fig.savefig(out_base.with_suffix(".png"), dpi=240)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)


def make_biological_annotation(data: dict[str, pd.DataFrame], prioritized: pd.DataFrame) -> None:
    rows = []
    for _, row in prioritized.iterrows():
        region, module = row["region"], row["module"]
        hubs = data["hubs"][(data["hubs"]["region"].eq(region)) & (data["hubs"]["module"].eq(module))].head(10)
        enrich = data["enrich"][(data["enrich"]["region"].eq(region)) & (data["enrich"]["module"].eq(module))].sort_values("fdr_BH").head(3)
        ahba = data["ahba"][(data["ahba"]["region"].eq(region)) & (data["ahba"]["mdd_module"].eq(module)) & (data["ahba"]["ahba_module"].isin(HIGHLIGHT_AHBA))].sort_values("fdr_BH").head(3)
        rows.append({
            "region": region,
            "module": module,
            "module_size": row.get("module_size", np.nan),
            "diagnosis_p": row.get("diagnosis_p_value", np.nan),
            "diagnosis_FDR": row.get("diagnosis_BH_FDR", np.nan),
            "top_hub_genes": ";".join(hubs["gene_symbol"].fillna(hubs["gene_id"]).astype(str)),
            "best_Tier2_ML_enrichment": "; ".join(enrich["gene_set"].astype(str) + " FDR=" + enrich["fdr_BH"].round(4).astype(str)),
            "best_AHBA_overlap": "; ".join(ahba["ahba_module"].astype(str) + " FDR=" + ahba["fdr_BH"].round(4).astype(str)),
            "interpretation": row.get("interpretation", ""),
            "caution": "Exploratory subregion module; use for follow-up prioritization, not validation.",
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(ANNOT_DIR / "frontal_subregion_prioritized_modules_biological_summary.csv", index=False)

    lines = [
        "# Frontal Subregion WGCNA Biological Annotation Report",
        "",
        "This report annotates prioritized frontal subregion WGCNA modules using saved module-trait, Tier2/ML enrichment, AHBA overlap, and hub-gene outputs. No classification was run.",
        "",
        "Broad Frontal Cortex WGCNA outputs were retained and were not overwritten. The broad Frontal Cortex result remains exploratory/supplementary.",
        "",
        "## Prioritized Modules",
        "",
    ]
    for _, r in summary.iterrows():
        lines.extend([
            f"### {r['region']} {r['module']}",
            "",
            f"- Module size: {r['module_size']}",
            f"- Diagnosis p/FDR: {r['diagnosis_p']:.4g} / {r['diagnosis_FDR']:.4g}",
            f"- Top hubs: {r['top_hub_genes']}",
            f"- Tier2/ML context: {r['best_Tier2_ML_enrichment']}",
            f"- AHBA context: {r['best_AHBA_overlap']}",
            f"- Interpretation: {r['interpretation']}",
            f"- Caution: {r['caution']}",
            "",
        ])
    lines.extend([
        "## Limitations",
        "",
        "GO/ontology plots were not generated in this follow-up because the saved WGCNA outputs do not include a local GO enrichment table for the new frontal subregion modules. The annotation here is therefore based on candidate-gene enrichment, AHBA overlap, and hub genes.",
    ])
    (ANNOT_DIR / "Frontal_Subregion_WGCNA_Biological_Annotation_Report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    data = load_inputs()
    comparison = make_comparison(data)
    prioritized = prioritize_modules(data)
    prioritized.to_csv(SUB_ROOT / "frontal_subregion_prioritized_modules_for_plots.csv", index=False)
    make_bubble_plot(data, prioritized)
    for _, row in prioritized.head(6).iterrows():
        make_network_plot(row["region"], row["module"], data, labeled=False)
        make_network_plot(row["region"], row["module"], data, labeled=True)
    make_biological_annotation(data, prioritized)
    manifest = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "comparison_rows": len(comparison),
        "prioritized_modules_for_plots": len(prioritized),
        "network_modules_plotted": min(6, len(prioritized)),
    }
    (SUB_ROOT / "frontal_subregion_followup_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
