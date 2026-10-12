from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


REPO = Path(__file__).resolve().parents[1]
BASE = REPO / "outputs" / "length_gt_100_clean_pipeline"
OUT = BASE / "22_Final_evidence_synthesis"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
REGION_SUMMARIES = OUT / "region_summaries"
MANIFEST = OUT / "manifest"

R17 = BASE / "17_ML_RRA_consensus"
R18 = BASE / "18_Tier2_ML_RRA_prioritization"
R19 = BASE / "19_Tier2_ML_benchmark"
R20 = BASE / "20_Tier2_ML_conditional_null"
R21 = BASE / "21_Statistical_method_benchmark"
WGCNA = BASE / "17_WGCNA_reconstruction_validation"

REGIONS = [
    "Ventral Striatum",
    "Insula",
    "Hippocampal Formation",
    "CG25",
    "dlPFC BA8/9",
    "OFC BA11",
]
EXPECTED_COUNTS = {
    "Ventral Striatum": 146,
    "Insula": 70,
    "Hippocampal Formation": 127,
    "CG25": 186,
    "dlPFC BA8/9": 72,
    "OFC BA11": 70,
}
SLUGS = {
    "Ventral Striatum": "Ventral_Striatum",
    "Insula": "Insula",
    "Hippocampal Formation": "Hippocampal_Formation",
    "CG25": "CG25",
    "dlPFC BA8/9": "dlPFC_BA8_9",
    "OFC BA11": "OFC_BA11",
}
WGCNA_DIRS = {
    "Ventral Striatum": "Ventral_Striatum",
    "Insula": "Insula",
    "Hippocampal Formation": "Hippocampal_Formation",
    "CG25": "Cingulate_gyrus_25",
    "dlPFC BA8/9": "dlPFC_BA8_9",
    "OFC BA11": "OFC_BA11",
}
NULL_CONTEXT = {
    "Ventral Striatum": "Clearest conditional RRA/LASSO stability support.",
    "Insula": "Mixed stability evidence; suggestive but inconsistent across endpoints.",
    "Hippocampal Formation": "Partial support, strongest for broader top-set endpoints.",
    "CG25": "Mixed stability evidence; leading-gene concentration but inconsistent ranking stability.",
    "dlPFC BA8/9": "Partial support, with stronger broader-top-set and LASSO evidence.",
    "OFC BA11": "Scenario E-like: prediction exceeded null, but gene-ranking stability did not.",
}
LEDGER_HASHES = {
    R17 / "checksums.sha256": "117c65e6cc796a9c6757b9d9e92ed493cc8bb6300d0873ea3d5cd2801ba3634c",
    R18 / "checksums.sha256": "3c45b0a92e4a5a3b8c9fc126097c312d492141734ea8e3fe72534a91152e950e",
    R19 / "checksums.sha256": "4e17dc8234317efc9bce92aa03fb7ae63ff72a128392d202f2f31a105c09cf10",
    R20 / "checksums.sha256": "d95f66712ed0ba80b7acefee7c3ee992cc14923128a47a5ddcaf9cff12e0d545",
    R21 / "checksums.sha256": "5626960b44eb446c023fd621f2808916f591667666be631a89691237a1b025ed",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype={"gene_id": str, "GeneID": str})
    if "GeneID" in frame.columns and "gene_id" not in frame.columns:
        frame = frame.rename(columns={"GeneID": "gene_id"})
    if "gene_id" in frame.columns:
        frame["gene_id"] = frame["gene_id"].astype(str).str.replace(r"\.0$", "", regex=True)
    return frame


def markdown_table(frame: pd.DataFrame) -> str:
    printable = frame.fillna("NA").astype(str)
    columns = list(printable.columns)
    rows = printable.values.tolist()
    widths = [max(len(str(columns[i])), *(len(str(row[i])) for row in rows)) for i in range(len(columns))]
    render = lambda row: "| " + " | ".join(str(value).ljust(widths[i]) for i, value in enumerate(row)) + " |"
    return "\n".join([render(columns), "| " + " | ".join("-" * width for width in widths) + " |", *(render(row) for row in rows)])


def merge_one_to_one(left: pd.DataFrame, right: pd.DataFrame, label: str) -> pd.DataFrame:
    keys = ["region", "gene_id"]
    if right.duplicated(keys).any():
        raise ValueError(f"Duplicate region/gene rows in {label}")
    before = len(left)
    merged = left.merge(right, on=keys, how="left", validate="one_to_one")
    if len(merged) != before:
        raise ValueError(f"Row count changed while merging {label}")
    return merged


def add_percentiles(frame: pd.DataFrame) -> pd.DataFrame:
    gene_counts = read_csv(R21 / "gene_universe_accounting.csv").set_index("region")["frozen_LM_eligible_genes"].to_dict()
    frame["LM_percentile_rank"] = frame.apply(lambda row: 100 * (row.LM_regional_rank - 1) / (gene_counts[row.region] - 1), axis=1)
    frame["limma_percentile_rank"] = frame["limma_percentile_rank"].astype(float)
    frame["DESeq2_percentile_rank"] = frame["DESeq2_percentile_rank"].astype(float)
    frame["Tier2_RRA_percentile_rank"] = frame.apply(
        lambda row: 100 * (row.Tier2_RRA_consensus_rank - 1) / max(row.Tier2_universe_size - 1, 1), axis=1
    )
    frame["LASSO_percentile_rank"] = frame.apply(
        lambda row: 100 * (row.LASSO_rank - 1) / max(EXPECTED_COUNTS[row.region] - 1, 1), axis=1
    )
    frame["Stability_Selection_percentile_rank"] = frame.apply(
        lambda row: 100 * (row.Stability_Selection_rank - 1) / max(EXPECTED_COUNTS[row.region] - 1, 1), axis=1
    )
    frame["allgene_RRA_percentile_rank"] = frame["allgene_RRA_percentile_rank"].astype(float) * 100
    return frame


def build_matrix() -> pd.DataFrame:
    for path, expected in LEDGER_HASHES.items():
        if sha256(path) != expected:
            raise RuntimeError(f"Frozen checksum ledger changed: {path}")

    statistical = read_csv(R21 / "tier2_statistical_method_audit.csv").rename(columns={
        "gene_symbol": "gene_symbol",
        "Frozen LM_effect": "LM_diagnosis_beta",
        "LM_standardized_effect": "LM_standardized_diagnosis_effect",
        "Frozen LM_p": "LM_raw_p",
        "Frozen LM_BH": "LM_BH",
        "Frozen LM_rank": "LM_regional_rank",
        "limma-voom_effect": "limma_diagnosis_logFC",
        "limma-voom_p": "limma_raw_p",
        "limma-voom_BH": "limma_BH",
        "limma-voom_rank": "limma_regional_rank",
        "DESeq2_effect": "DESeq2_log2FoldChange",
        "DESeq2_p": "DESeq2_raw_p",
        "DESeq2_native_padj": "DESeq2_native_padj",
        "DESeq2_BH": "DESeq2_same_universe_BH",
        "DESeq2_rank": "DESeq2_regional_rank",
        "LM_limma_direction_agreement": "LM_limma_direction_agreement",
        "LM_DESeq2_direction_agreement": "LM_DESeq2_direction_agreement",
        "three_method_direction_agreement": "three_method_direction_agreement",
    })
    keep = [
        "region", "gene_id", "gene_symbol", "LM_diagnosis_beta", "LM_standardized_diagnosis_effect",
        "LM_raw_p", "LM_BH", "LM_regional_rank", "limma_diagnosis_logFC", "limma_raw_p",
        "limma_BH", "limma_regional_rank", "DESeq2_log2FoldChange", "DESeq2_raw_p",
        "DESeq2_native_padj", "DESeq2_same_universe_BH", "DESeq2_regional_rank",
        "LM_limma_direction_agreement", "LM_DESeq2_direction_agreement", "three_method_direction_agreement",
        "limma_percentile_rank", "DESeq2_percentile_rank",
    ]
    matrix = statistical[keep].copy()
    matrix.insert(3, "Tier2_membership", True)

    rra = read_csv(R18 / "Tier2_consensus_rra_rankings.csv").rename(columns={
        "consensus_rank": "Tier2_RRA_consensus_rank",
        "consensus_rra_score": "Tier2_RRA_consensus_score",
        "model_appearance_count": "Tier2_RRA_model_support_count",
        "tier2_universe_size": "Tier2_universe_size",
        "logistic_l1_rank": "Tier2_RRA_logistic_l1_rank",
        "elastic_net_rank": "Tier2_RRA_elastic_net_rank",
        "logistic_l2_rank": "Tier2_RRA_logistic_l2_rank",
        "linear_svm_rank": "Tier2_RRA_linear_svm_rank",
        "random_forest_rank": "Tier2_RRA_random_forest_rank",
    })
    rra_cols = ["region", "gene_id"] + [column for column in rra.columns if column.startswith("Tier2_")]
    matrix = merge_one_to_one(matrix, rra[rra_cols], "Tier2-RRA")

    lasso = read_csv(R19 / "LASSO_gene_rankings.csv").rename(columns={
        "lasso_export_rank": "LASSO_rank",
        "selection_count": "LASSO_selection_count",
        "selection_frequency": "LASSO_selection_frequency",
        "mean_rank_when_selected": "LASSO_mean_rank_when_selected",
        "mean_abs_coefficient_when_selected": "LASSO_mean_abs_coefficient",
        "sign_consistency_when_selected": "LASSO_sign_consistency",
    })
    matrix = merge_one_to_one(matrix, lasso[["region", "gene_id", "LASSO_rank", "LASSO_selection_count", "LASSO_selection_frequency", "LASSO_mean_rank_when_selected", "LASSO_mean_abs_coefficient", "LASSO_sign_consistency"]], "LASSO")

    stability = read_csv(R19 / "stability_selection_gene_rankings.csv").rename(columns={
        "stability_export_rank": "Stability_Selection_rank",
        "selection_probability": "Stability_Selection_probability",
        "above_predeclared_0_60_threshold": "Stability_Selection_ge_0_60_descriptive",
    })
    matrix = merge_one_to_one(matrix, stability[["region", "gene_id", "Stability_Selection_rank", "Stability_Selection_probability", "Stability_Selection_ge_0_60_descriptive"]], "Stability Selection")

    allgene = read_csv(R18 / "Tier2_vs_allgene_RRA_comparison.csv").rename(columns={
        "allgene_consensus_rank": "allgene_RRA_rank",
        "allgene_rank_percentile": "allgene_RRA_percentile_rank",
        "allgene_consensus_rra_score": "allgene_RRA_score",
        "present_in_frozen_allgene_rra": "allgene_RRA_available",
    })
    matrix = merge_one_to_one(matrix, allgene[["region", "gene_id", "allgene_RRA_rank", "allgene_RRA_percentile_rank", "allgene_RRA_score", "allgene_RRA_available"]], "all-gene RRA")

    rra_null = read_csv(R20 / "summary" / "RRA_observed_topgene_null_behavior.csv").rename(columns={
        "null_top10_count": "RRA_null_top10_count_of_20",
        "null_top25_count": "RRA_null_top25_count_of_20",
        "null_rank_median": "RRA_null_rank_median",
        "null_rank_best": "RRA_null_rank_best",
        "exploratory_rank_empirical_p_pilot": "RRA_exploratory_rank_empirical_p_pilot",
    })
    matrix = merge_one_to_one(matrix, rra_null[["region", "gene_id", "RRA_null_top10_count_of_20", "RRA_null_top25_count_of_20", "RRA_null_rank_median", "RRA_null_rank_best", "RRA_exploratory_rank_empirical_p_pilot"]], "RRA conditional null")

    lasso_null = read_csv(R20 / "summary" / "LASSO_observed_topgene_null_behavior.csv").rename(columns={
        "null_top10_count": "LASSO_null_top10_count_of_20",
        "null_top25_count": "LASSO_null_top25_count_of_20",
        "null_rank_median": "LASSO_null_rank_median",
        "null_rank_best": "LASSO_null_rank_best",
        "exploratory_rank_empirical_p_pilot": "LASSO_exploratory_rank_empirical_p_pilot",
    })
    matrix = merge_one_to_one(matrix, lasso_null[["region", "gene_id", "LASSO_null_top10_count_of_20", "LASSO_null_top25_count_of_20", "LASSO_null_rank_median", "LASSO_null_rank_best", "LASSO_exploratory_rank_empirical_p_pilot"]], "LASSO conditional null")
    matrix["conditional_null_region_context"] = matrix.region.map(NULL_CONTEXT)

    assignments = []
    traits = []
    for region, directory in WGCNA_DIRS.items():
        folder = WGCNA / "reconstruction" / directory
        assignment = read_csv(folder / "reconstructed_gene_to_module_assignments.csv")
        assignment["region"] = region
        assignment = assignment.rename(columns={"reconstructed_module": "WGCNA_module"})
        assignments.append(assignment[["region", "gene_id", "WGCNA_module"]])
        trait = read_csv(folder / "reconstructed_module_trait_associations.csv")
        trait["region"] = region
        trait = trait.rename(columns={"reconstructed_module": "WGCNA_module", "p_diagnosis": "WGCNA_module_diagnosis_raw_p", "fdr_diagnosis_BH": "WGCNA_module_diagnosis_BH"})
        traits.append(trait[["region", "WGCNA_module", "WGCNA_module_diagnosis_raw_p", "WGCNA_module_diagnosis_BH"]])
    wgcna = pd.concat(assignments, ignore_index=True).merge(pd.concat(traits, ignore_index=True), on=["region", "WGCNA_module"], how="left", validate="many_to_one")
    matrix = merge_one_to_one(matrix, wgcna, "WGCNA contextual assignment")
    matrix["WGCNA_context_caveat"] = np.where(matrix.WGCNA_module.notna(), "context only; not gene prioritisation or independent MDD replication", pd.NA)
    matrix["AHBA_gene_level_context"] = pd.NA
    matrix["AHBA_context_caveat"] = "No unambiguous gene-region replication field; left NA and excluded from prioritisation."

    matrix = add_percentiles(matrix)
    order = {region: index for index, region in enumerate(REGIONS)}
    matrix["_region_order"] = matrix.region.map(order)
    matrix = matrix.sort_values(["_region_order", "Tier2_RRA_consensus_rank", "gene_id"]).drop(columns="_region_order").reset_index(drop=True)
    return matrix


def build_fdr_table(matrix: pd.DataFrame) -> pd.DataFrame:
    rows = []
    limma = read_csv(R21 / "limma_voom_results_all_regions.csv")
    for row in limma[limma.BH < 0.05].itertuples(index=False):
        rows.append({"method": "limma-voom", "region": row.region, "gene_id": row.gene_id, "gene_symbol": row.gene_symbol, "effect": row.logFC, "raw_p": row.raw_p, "adjusted_p": row.BH, "rank": row.primary_rank})
    deseq = read_csv(R21 / "deseq2_results_all_regions.csv")
    for row in deseq[deseq.BH_same_universe < 0.05].itertuples(index=False):
        rows.append({"method": "DESeq2 same-universe BH", "region": row.region, "gene_id": row.gene_id, "gene_symbol": row.gene_symbol, "effect": row.log2FoldChange, "raw_p": row.raw_p, "adjusted_p": row.BH_same_universe, "rank": row.primary_rank})
    output = pd.DataFrame(rows)
    lookup = matrix[["region", "gene_id", "Tier2_membership", "Tier2_RRA_consensus_rank"]]
    output = output.merge(lookup, on=["region", "gene_id"], how="left")
    output["Tier2_membership"] = output.Tier2_membership.fillna(False).astype(bool)
    output["finding_label"] = "method-specific sensitivity finding"
    return output.sort_values(["region", "method", "rank"]).reset_index(drop=True)


def build_provenance(matrix: pd.DataFrame) -> pd.DataFrame:
    specs = {}
    def add(columns, analysis, source, fields, transformation, interpretation, caveat):
        for column in columns:
            specs[column] = [analysis, source, fields.get(column, column), transformation, interpretation, caveat]

    add(["region", "gene_id", "gene_symbol", "Tier2_membership", "LM_diagnosis_beta", "LM_standardized_diagnosis_effect", "LM_raw_p", "LM_BH", "LM_regional_rank"], "Frozen LM/Tier2 statistical audit", "21_Statistical_method_benchmark/tier2_statistical_method_audit.csv", {}, "renamed only; Tier2 membership fixed TRUE", "Original frozen statistical evidence", "No frozen LM gene survived regional BH FDR")
    add(["limma_diagnosis_logFC", "limma_raw_p", "limma_BH", "limma_regional_rank", "DESeq2_log2FoldChange", "DESeq2_raw_p", "DESeq2_native_padj", "DESeq2_same_universe_BH", "DESeq2_regional_rank", "LM_limma_direction_agreement", "LM_DESeq2_direction_agreement", "three_method_direction_agreement", "limma_percentile_rank", "DESeq2_percentile_rank"], "Statistical-method benchmark", "21_Statistical_method_benchmark/tier2_statistical_method_audit.csv", {}, "renamed only", "Within-study method robustness", "Same cohort; not independent replication")
    add([column for column in matrix if column.startswith("Tier2_RRA_")] + ["Tier2_universe_size"], "Primary Tier2 five-model ML-RRA", "18_Tier2_ML_RRA_prioritization/Tier2_consensus_rra_rankings.csv", {}, "renamed; percentile calculated within frozen regional Tier2 universe", "Primary exploratory Tier2 prioritisation", "RRA score is not an MDD-association p-value")
    add([column for column in matrix if column.startswith("LASSO_") and "null" not in column], "LASSO benchmark", "19_Tier2_ML_benchmark/LASSO_gene_rankings.csv", {}, "renamed; percentile calculated within frozen regional Tier2 universe", "Same-cohort benchmark evidence", "Selection frequency is descriptive, not external validation")
    add([column for column in matrix if column.startswith("Stability_Selection_")], "Stability Selection benchmark", "19_Tier2_ML_benchmark/stability_selection_gene_rankings.csv", {}, "renamed; percentile calculated within frozen regional Tier2 universe", "Same-cohort benchmark evidence", "0.60 indicator is predeclared and descriptive, not significance")
    add([column for column in matrix if column.startswith("allgene_RRA_")], "All-gene ML-RRA sensitivity", "18_Tier2_ML_RRA_prioritization/Tier2_vs_allgene_RRA_comparison.csv", {}, "renamed; stored proportion converted to percent for percentile field", "Within-study sensitivity evidence", "Not independent validation")
    add([column for column in matrix if column.startswith("RRA_null_")] + ["RRA_exploratory_rank_empirical_p_pilot"], "20-permutation conditional-null pilot", "20_Tier2_ML_conditional_null/summary/RRA_observed_topgene_null_behavior.csv", {}, "left join; unavailable genes remain NA", "Gene-level conditional-null behavior where frozen output exists", "Only 20 permutations; exploratory resolution floor 1/21")
    add([column for column in matrix if column.startswith("LASSO_null_")] + ["LASSO_exploratory_rank_empirical_p_pilot"], "20-permutation conditional-null pilot", "20_Tier2_ML_conditional_null/summary/LASSO_observed_topgene_null_behavior.csv", {}, "left join; unavailable genes remain NA", "Gene-level conditional-null behavior where frozen output exists", "Only 20 permutations; exploratory resolution floor 1/21")
    add(["conditional_null_region_context"], "Frozen conditional-null report", "docs/TIER2_ML_CONDITIONAL_NULL_REPORT.md", {}, "verbatim conservative region-level synthesis", "Region-level conditional-null interpretation", "Not a gene-level significance classification")
    add(["WGCNA_module", "WGCNA_module_diagnosis_raw_p", "WGCNA_module_diagnosis_BH", "WGCNA_context_caveat"], "Phase-3B WGCNA reconstruction validation", "17_WGCNA_reconstruction_validation/reconstruction/<region>/reconstructed_gene_to_module_assignments.csv and reconstructed_module_trait_associations.csv", {}, "exact gene/module join where present", "Biological/module context only", "Not used to rank genes; no module survived diagnosis BH FDR")
    add(["AHBA_gene_level_context", "AHBA_context_caveat"], "AHBA context audit", "No unambiguous gene-region field used", {}, "left NA", "Documents absence of reliable gene-level linkage", "AHBA is context only and not independent MDD replication")
    add(["LM_percentile_rank"], "Derived display field", "21_Statistical_method_benchmark gene universe and LM rank", {}, "100*(rank-1)/(regional eligible genes-1)", "Comparable display percentile", "Not a score")

    missing = sorted(set(matrix.columns) - set(specs))
    if missing:
        raise RuntimeError(f"Missing provenance definitions: {missing}")
    rows = []
    for column in matrix.columns:
        analysis, source, source_field, transformation, interpretation, caveat = specs[column]
        rows.append({"column_name": column, "source_analysis": analysis, "source_file": source, "source_field": source_field, "transformation": transformation, "interpretation": interpretation, "caveat": caveat})
    return pd.DataFrame(rows)


def save_figures(matrix: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid")
    top = matrix[matrix.Tier2_RRA_consensus_rank <= 10].copy()
    top["label"] = top.region + " | " + top.gene_symbol.fillna(top.gene_id)
    percentile_columns = {
        "LM_percentile_rank": "LM",
        "limma_percentile_rank": "limma",
        "DESeq2_percentile_rank": "DESeq2",
        "Tier2_RRA_percentile_rank": "Tier2-RRA",
        "LASSO_percentile_rank": "LASSO",
        "Stability_Selection_percentile_rank": "Stability",
        "allgene_RRA_percentile_rank": "All-gene RRA",
    }
    heat = top.set_index("label")[list(percentile_columns)].rename(columns=percentile_columns)
    fig, ax = plt.subplots(figsize=(12, 18))
    sns.heatmap(heat, cmap="viridis_r", vmin=0, vmax=100, ax=ax, cbar_kws={"label": "Regional percentile rank (0 = strongest)"})
    ax.set_title("Frozen evidence ranks for top-10 Tier2-RRA candidates per region")
    ax.set_xlabel("Evidence stream (not combined)")
    ax.set_ylabel("")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right")
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        fig.savefig(FIGURES / f"top10_cross_method_rank_percentiles.{suffix}", dpi=220, metadata={"CreationDate": None, "ModDate": None} if suffix == "pdf" else None)
    plt.close(fig)

    agreement = matrix.groupby("region", sort=False)[["LM_limma_direction_agreement", "LM_DESeq2_direction_agreement", "three_method_direction_agreement"]].mean().reindex(REGIONS) * 100
    agreement.columns = ["LM–limma", "LM–DESeq2", "All three"]
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.heatmap(agreement, annot=True, fmt=".1f", cmap="Blues", vmin=90, vmax=100, ax=ax, cbar_kws={"label": "Tier2 direction agreement (%)"})
    ax.set_title("Direction agreement across frozen Tier2 pairs")
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        fig.savefig(FIGURES / f"regional_statistical_direction_agreement.{suffix}", dpi=220, metadata={"CreationDate": None, "ModDate": None} if suffix == "pdf" else None)
    plt.close(fig)

    levels = {"Ventral Striatum": 3, "Hippocampal Formation": 2, "dlPFC BA8/9": 2, "Insula": 1, "CG25": 1, "OFC BA11": 0}
    context = pd.DataFrame({"region": REGIONS, "frozen_context_display": [levels[r] for r in REGIONS]}).set_index("region")
    fig, ax = plt.subplots(figsize=(8, 4))
    sns.heatmap(context, cmap=sns.color_palette(["#d9d9d9", "#fdae6b", "#9ecae1", "#3182bd"], as_cmap=True), vmin=0, vmax=3, annot=np.array([["Clearest"], ["Mixed"], ["Partial"], ["Mixed"], ["Partial"], ["Scenario E-like"]]), fmt="", cbar=False, ax=ax)
    ax.set_title("Frozen region-level conditional-null interpretation")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        fig.savefig(FIGURES / f"conditional_null_region_interpretation.{suffix}", dpi=220, metadata={"CreationDate": None, "ModDate": None} if suffix == "pdf" else None)
    plt.close(fig)

    leaders = matrix[matrix.Tier2_RRA_consensus_rank == 1].copy().set_index("region").reindex(REGIONS)
    profiles = leaders[list(percentile_columns)].rename(columns=percentile_columns)
    fig, ax = plt.subplots(figsize=(10, 5))
    for region, row in profiles.iterrows():
        ax.plot(profiles.columns, row.values, marker="o", label=f"{region}: {leaders.loc[region, 'gene_symbol']}")
    ax.set_yscale("symlog", linthresh=0.1)
    ax.invert_yaxis()
    ax.set_ylabel("Regional percentile rank (lower is stronger; inverted axis)")
    ax.set_title("Separate evidence profiles for the six regional Tier2-RRA leaders")
    ax.tick_params(axis="x", rotation=25)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        fig.savefig(FIGURES / f"regional_leader_evidence_profiles.{suffix}", dpi=220, metadata={"CreationDate": None, "ModDate": None} if suffix == "pdf" else None)
    plt.close(fig)


def format_candidate_table(frame: pd.DataFrame) -> pd.DataFrame:
    columns = ["Tier2_RRA_consensus_rank", "gene_symbol", "LM_regional_rank", "limma_regional_rank", "DESeq2_regional_rank", "LASSO_rank", "LASSO_selection_frequency", "Stability_Selection_rank", "Stability_Selection_probability", "allgene_RRA_rank"]
    result = frame[columns].copy()
    result.columns = ["RRA rank", "Gene", "LM rank", "limma rank", "DESeq2 rank", "LASSO rank", "LASSO frequency", "Stability rank", "Stability probability", "All-gene RRA rank"]
    for column in ["LASSO frequency", "Stability probability"]:
        result[column] = result[column].map(lambda value: "NA" if pd.isna(value) else f"{value:.2f}")
    return result


def write_region_summaries(matrix: pd.DataFrame) -> None:
    for region in REGIONS:
        subset = matrix[matrix.region == region].sort_values("Tier2_RRA_consensus_rank")
        leader = subset.iloc[0]
        text = f"""# {region}: final evidence synthesis

## Statistical-method concordance

All {len(subset)} frozen Tier2 pairs retained the same direction across the LM, limma-voom, and DESeq2. The regional LM–limma and LM–DESeq2 genome-wide rank correlations are documented in the frozen statistical benchmark. No new statistical consensus score was calculated.

## Tier2-RRA, LASSO, and Stability Selection

The primary order remains the frozen Tier2-RRA order. LASSO and Stability Selection are displayed as separate same-cohort benchmark evidence and do not alter that order.

## All-gene sensitivity

All-gene RRA ranks are shown only as within-study sensitivity evidence. Missing all-gene matches remain `NA`; they are not interpreted as negative findings.

## Conditional-null result

{NULL_CONTEXT[region]} This is the frozen 20-permutation pilot interpretation, not a new gene-level significance threshold.

## Leading candidates

{markdown_table(format_candidate_table(subset.head(10)))}

The regional Tier2-RRA leader is **{leader.gene_symbol}**. Its LM/limma/DESeq2 ranks are {int(leader.LM_regional_rank)}/{int(leader.limma_regional_rank)}/{int(leader.DESeq2_regional_rank)}; its LASSO and Stability Selection ranks are {int(leader.LASSO_rank)}/{int(leader.Stability_Selection_rank)}. A number-one ML rank does not establish biological causality or diagnostic validity.

## Caveats

- All statistical and ML evidence comes from the same cohort.
- Tier2 was selected using the frozen LM.
- RRA scores are rank-aggregation quantities, not MDD-association p-values.
- Conditional-null resolution is limited to 20 permutations.
- WGCNA is contextual annotation only; AHBA is not used as gene-level replication.
- No result here is a validated biomarker or external replication.
"""
        (REGION_SUMMARIES / f"{SLUGS[region]}.md").write_text(text, encoding="utf-8", newline="\n")


def write_report(matrix: pd.DataFrame, fdr: pd.DataFrame) -> None:
    leaders = matrix[matrix.Tier2_RRA_consensus_rank == 1].copy().set_index("region").reindex(REGIONS).reset_index()
    leader_table = leaders[["region", "gene_symbol", "LM_regional_rank", "limma_regional_rank", "DESeq2_regional_rank", "LASSO_rank", "LASSO_selection_frequency", "Stability_Selection_rank", "Stability_Selection_probability", "allgene_RRA_rank", "RRA_null_top10_count_of_20"]].copy()
    leader_table.columns = ["Region", "Leader", "LM rank", "limma rank", "DESeq2 rank", "LASSO rank", "LASSO frequency", "Stability rank", "Stability probability", "All-gene RRA rank", "Null top-10 / 20"]
    report = f"""# Final cross-method evidence synthesis

## 1. Purpose

This report assembles the frozen statistical, ML, sensitivity, benchmark, and conditional-null evidence for exactly 671 frozen Tier2 region–gene pairs. It does not rerun a model, change Tier2, or create a new score or confidence category.

## 2. Frozen analysis stages

The synthesis covers the original LM/Tier2 screen, all-gene ML-RRA sensitivity, primary Tier2 five-model ML-RRA, LASSO and Stability Selection benchmarks, the 20-permutation conditional-null pilot, LM/limma-voom/DESeq2 concordance, and limited WGCNA module context.

## 3. Why no additional scoring was introduced

The evidence streams answer different questions and use different scales. Combining them numerically would impose new, post-hoc weights. The frozen Tier2-RRA rank therefore remains the sole candidate order, while every other field remains visible separately.

## 4. Statistical robustness of Tier2

All 671 Tier2 pairs retained the same effect direction across the LM, limma-voom, and DESeq2. Genome-wide raw-p rank correlations were strong in all regions. Tier2 is therefore robust to statistical-model choice as an exploratory screen, but this is within-study concordance rather than replication. No frozen LM gene survived BH FDR.

## 5. ML prioritisation benchmark results

Tier2-RRA is an exploratory multivariate prioritisation. LASSO and Stability Selection often converge with its leading candidates, but neither supplies independent validation. The all-gene RRA comparison is a sensitivity analysis using the same cohort.

## 6. Conditional-null findings

- **Ventral Striatum:** clearest conditional RRA/LASSO stability support.
- **Hippocampal Formation:** partial support.
- **dlPFC BA8/9:** partial support.
- **Insula:** mixed evidence.
- **CG25:** mixed evidence.
- **OFC BA11:** Scenario E-like; predictive performance exceeded null but ranking stability did not.

The pilot used 20 permutations, so its smallest possible plus-one empirical quantity is 1/21. It supports region-specific caution, not validated gene calls.

## 7. All-gene sensitivity results

All-gene ranks are retained as separate sensitivity evidence. Strong convergence can support a narrative example, but a weak or missing all-gene rank does not automatically invalidate a Tier2 candidate.

## 8. Region-by-region evidence

Detailed summaries are in `22_Final_evidence_synthesis/region_summaries/`. They preserve the frozen regional interpretation and list the top ten candidates without reranking.

## 9. Leading candidate evidence profiles

{markdown_table(leader_table)}

The leaders are useful audit examples, not proven biological winners. MORN5 has the most favorable regional null context. TOR2A and LINC00269 retain partial regional support. TMEM161B-DT and LOC107985795 require mixed-evidence caveats. LOC105373386 requires the strongest caution because OFC ranking stability did not exceed the conditional null and its DESeq2 rank is weaker.

## 10. Method-specific FDR sensitivity findings

The separate file `method_specific_FDR_findings.csv` contains {len(fdr)} method-specific records: one limma-voom finding and nine same-universe DESeq2 findings. These include OFC `QPRT` and the eight Ventral Striatum DESeq2 findings. They are not promoted into Tier2 or a new shortlist.

## 11. Biological-context annotations

Unambiguous Phase-3B WGCNA module assignments and module diagnosis p/BH values are appended as context only. Missing module assignments remain `NA`. No module survived diagnosis BH FDR. AHBA is not represented as a gene-level replication field because no unambiguous gene–region linkage suitable for that claim was available; it remains context only.

## 12. Main limitations

1. All evidence streams use the same cohort.
2. Tier2 was defined from the original LM, so downstream concordance is not independent.
3. Regional sample sizes are small, especially CG25.
4. Statistical and ML methods answer different questions and use different scales.
5. Conditional-null resolution is limited to 20 permutations.
6. WGCNA and AHBA are contextual, not independent MDD replication.
7. No external cohort is included.

## 13. What can be claimed now

- Tier2 is robust to statistical-model choice as an exploratory candidate screen.
- Tier2-RRA supplies an exploratory multivariate prioritisation of frozen Tier2 candidates.
- RRA broadly converges with LASSO and Stability Selection for several leading candidates.
- Ventral Striatum shows the clearest conditional evidence for label-dependent prioritisation.
- Several candidates show within-study statistical and ML convergence.

## 14. What cannot be claimed now

The project has not established a validated biomarker, diagnostic gene, causal MDD gene, clinically predictive model, proven mechanism, independent replication, or external validation. Method-specific FDR findings cannot be described as canonical LM significance.

## 15. Recommended next scientific step

The next major scientific step is preregistered external validation in an independent cohort, preserving the frozen Tier2 identities and primary Tier2-RRA ordering before looking at external outcomes.
"""
    (REPO / "docs" / "FINAL_CROSS_METHOD_EVIDENCE_SYNTHESIS.md").write_text(report, encoding="utf-8", newline="\n")


def main() -> None:
    for directory in [TABLES, FIGURES, REGION_SUMMARIES, MANIFEST]:
        directory.mkdir(parents=True, exist_ok=True)

    matrix = build_matrix()
    counts = matrix.groupby("region").size().to_dict()
    if len(matrix) != 671 or counts != EXPECTED_COUNTS:
        raise RuntimeError(f"Tier2 identity/count failure: {len(matrix)}, {counts}")
    if matrix[["Tier2_RRA_consensus_rank", "LASSO_rank", "Stability_Selection_rank"]].isna().any().any():
        raise RuntimeError("Unexpected missing primary/benchmark ranks")

    matrix.to_csv(TABLES / "Tier2_final_evidence_matrix.csv", index=False, lineterminator="\n")
    matrix[matrix.Tier2_RRA_consensus_rank <= 10].to_csv(TABLES / "Tier2_top10_by_region_evidence.csv", index=False, lineterminator="\n")
    matrix[matrix.Tier2_RRA_consensus_rank <= 25].to_csv(TABLES / "Tier2_top25_by_region_evidence.csv", index=False, lineterminator="\n")
    leaders = matrix[matrix.Tier2_RRA_consensus_rank == 1].copy()
    leaders.to_csv(TABLES / "regional_Tier2_RRA_leader_evidence_profiles.csv", index=False, lineterminator="\n")
    fdr = build_fdr_table(matrix)
    fdr.to_csv(TABLES / "method_specific_FDR_findings.csv", index=False, lineterminator="\n")
    provenance = build_provenance(matrix)
    provenance.to_csv(MANIFEST / "evidence_provenance.csv", index=False, lineterminator="\n")

    save_figures(matrix)
    write_region_summaries(matrix)
    write_report(matrix, fdr)

    source_hashes = {str(path.relative_to(REPO)).replace("\\", "/"): sha256(path) for path in LEDGER_HASHES}
    manifest = {
        "frozen_start_commit": "437b2b63d3337896cf71085160d1545df687b1ed",
        "primary_rows": len(matrix),
        "regional_counts": counts,
        "primary_order": "frozen Tier2-RRA consensus rank",
        "new_weighted_score": False,
        "new_confidence_category": False,
        "allgene_role": "within-study sensitivity",
        "conditional_null_role": "20-permutation conditional pilot",
        "wgcna_ahba_role": "context only",
        "source_ledger_hashes": source_hashes,
    }
    (MANIFEST / "synthesis_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n")

    checksum_file = OUT / "checksums.sha256"
    files = sorted(path for path in OUT.rglob("*") if path.is_file() and path != checksum_file)
    checksum_file.write_text("\n".join(f"{sha256(path)}  {path.relative_to(OUT).as_posix()}" for path in files) + "\n", encoding="utf-8", newline="\n")
    print(f"Built final synthesis: {len(matrix)} Tier2 rows, {len(files)} checksummed files")


if __name__ == "__main__":
    main()
