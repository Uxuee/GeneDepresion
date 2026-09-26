from __future__ import annotations

import hashlib
import json
import math
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
P3B = ROOT / "outputs/length_gt_100_clean_pipeline/17_WGCNA_reconstruction_validation"
OUT = ROOT / "outputs/length_gt_100_clean_pipeline/18_WGCNA_enrichment_interpretation_audit"
FUNC = OUT / "functional_annotation"
MANIFESTS = OUT / "manifests"
DOC = ROOT / "docs/WGCNA_FUNCTIONAL_ENRICHMENT_INTERPRETATION_AUDIT.md"
for directory in (OUT, FUNC, MANIFESTS):
    directory.mkdir(parents=True, exist_ok=True)

BASELINE = "dd59b1dc2b64d681ec9e7b8f0092b41b005471ef"
Q_FDR = 0.05
Q_SUGGESTIVE = 0.10
SMALL_OVERLAP_MAX = 3
SMALL_MODULE_MAX = 50
ALLOWED_DATABASES = {"GO", "Reactome", "MSigDB Hallmark"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def norm_gene(value) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    return text[:-2] if text.endswith(".0") else text


def split_genes(value) -> list[str]:
    if pd.isna(value) or not str(value).strip():
        return []
    return [norm_gene(x) for x in str(value).split("/") if norm_gene(x)]


def parse_ratio(value) -> tuple[float, float]:
    try:
        numerator, denominator = str(value).split("/", 1)
        return float(numerator), float(denominator)
    except Exception:
        return math.nan, math.nan


def clean_bool(value) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes"}


def family_label(symbol: str) -> str:
    symbol = str(symbol or "").upper()
    rules = [
        (r"^TAS[12]R", "taste_receptor"),
        (r"^OR\d|^OR\d+[A-Z]", "olfactory_receptor"),
        (r"^PCDH", "protocadherin"),
        (r"^HIST|^H[1234][A-Z]", "histone"),
        (r"^(RPL|RPS|MRPL|MRPS)", "ribosomal"),
        (r"^(IGH|IGK|IGL)", "immunoglobulin"),
        (r"^(SNORD|SNORA|RNU|SCARNA)", "small_nuclear_nucleolar_RNA"),
    ]
    for pattern, label in rules:
        if re.search(pattern, symbol):
            return label
    return "other"


THEME_RULES = [
    ("NTRK_TRKA_SIGNALING", r"\bntrk|\btrka"),
    ("RNA_PROCESSING_NUCLEOLAR", r"rna processing|rna splic|spliceos|ribonucleoprotein|nucleol|cajal|rna binding|mrna"),
    ("SYNAPTIC_PRESYNAPTIC", r"synap|presynap|synaptic vesicle|exocytic vesicle|neurotransmitter|dense core granule|node of ranvier"),
    ("CHOLESTEROL_LIPID", r"cholesterol|sterol|lipid|lipoprotein"),
    ("EXTRACELLULAR_MATRIX_WOUND", r"extracellular matrix|response to wounding|wound healing|collagen|encapsulating structure|matrix disassembly"),
    ("ACTIN_CELL_CORTEX", r"actin|cell cortex|cortical cytoskeleton"),
    ("ANGIOGENESIS_VASCULAR", r"angiogen|vascular|blood vessel|endothelial"),
    ("UBIQUITINATION", r"ubiquit"),
    ("RECEPTOR_GPCR_FAMILY", r"receptor|g protein-coupled|gpcr|taste|pheromone"),
    ("SIGNALING_IMMEDIATE_EARLY", r"mapk|kinase|interleukin|cytokine|growth factor|intracellular signal|nuclear events|transcription factor|response to steroid"),
]


def theme_for(term: str) -> str:
    text = str(term).lower()
    for label, pattern in THEME_RULES:
        if re.search(pattern, text):
            return label
    return "OTHER_ANNOTATION"


def list_join(values) -> str:
    return ";".join(sorted({str(x) for x in values if pd.notna(x) and str(x) not in {"", "nan"}}))


enrichment_path = P3B / "functional_annotation/functional_enrichment_GO_Reactome_Hallmark.csv"
evidence_path = P3B / "functional_annotation/prioritized_module_evidence_table.csv"
mapping_path = P3B / "functional_annotation/prioritized_module_gene_id_mapping.csv"
ml_enrichment_path = P3B / "ml_gene_enrichment_validation/aligned_corrected_ML_enrichment_comparison.csv"
ml_path = ROOT / "outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_corrected/ML_prioritized_Tier2_genes_all_ranked.csv"
core_broad_path = ROOT / "outputs/length_gt_100_clean_pipeline/16_corrected_ML_downstream/stable_core_annotations_corrected/broad_structural_core_with_corrected_ML_annotations.csv"
core_frontal_path = ROOT / "outputs/length_gt_100_clean_pipeline/16_corrected_ML_downstream/stable_core_annotations_corrected/frontal_structural_core_with_corrected_ML_annotations.csv"
go_relations_path = FUNC / "phase3c_go_immediate_parent_relations.csv"

phase3b_source_paths = [
    enrichment_path, evidence_path, mapping_path, ml_enrichment_path, ml_path,
    core_broad_path, core_frontal_path,
]

enrichment = pd.read_csv(enrichment_path)
evidence = pd.read_csv(evidence_path)
mapping = pd.read_csv(mapping_path, dtype={"gene_id": str})
ml_enrichment = pd.read_csv(ml_enrichment_path)
ml = pd.read_csv(ml_path, dtype={"gene_id": str})
core = pd.concat([
    pd.read_csv(core_broad_path, dtype={"gene_id": str}),
    pd.read_csv(core_frontal_path, dtype={"gene_id": str}),
], ignore_index=True, sort=False)
go_relations = pd.read_csv(go_relations_path) if go_relations_path.exists() else pd.DataFrame(
    columns=["child_term_id", "parent_term_id", "relationship", "ontology"]
)

for frame in (mapping, ml, core):
    frame["gene_id"] = frame["gene_id"].map(norm_gene)

priority_keys = evidence[["region", "old_module"]].drop_duplicates()
assert len(priority_keys) == 10, "Phase 3B must contain exactly ten prespecified prioritized modules"
assert set(enrichment["database"].unique()) == ALLOWED_DATABASES, "Unexpected enrichment database"
enrichment = enrichment.merge(priority_keys, on=["region", "old_module"], how="inner", validate="many_to_one")

ratios = enrichment["gene_ratio"].map(parse_ratio)
bg_ratios = enrichment["background_ratio"].map(parse_ratio)
enrichment["gene_ratio_numerator"] = [x[0] for x in ratios]
enrichment["gene_ratio_denominator"] = [x[1] for x in ratios]
enrichment["term_size_in_background"] = [x[0] for x in bg_ratios]
enrichment["background_ratio_denominator"] = [x[1] for x in bg_ratios]
enrichment["overlap_fraction_of_module_mapped"] = enrichment["overlap_count"] / enrichment["module_genes_mapped"]
enrichment["overlap_fraction_of_term_background"] = enrichment["overlap_count"] / enrichment["term_size_in_background"]
enrichment["small_overlap_flag"] = enrichment["overlap_count"] <= SMALL_OVERLAP_MAX
enrichment["small_module_flag"] = enrichment["module_genes_mapped"] <= SMALL_MODULE_MAX
enrichment["q_category"] = np.select(
    [enrichment["BH_q"] < Q_FDR, enrichment["BH_q"] < Q_SUGGESTIVE],
    ["FDR_LT_0_05", "SUGGESTIVE_0_05_TO_0_10"], default="NOT_LT_0_10"
)
enrichment["theme_rule"] = enrichment["term"].map(theme_for)
enrichment["interpretation_caution"] = np.select(
    [
        enrichment["small_overlap_flag"] & enrichment["small_module_flag"],
        enrichment["small_overlap_flag"],
        enrichment["small_module_flag"],
    ],
    [
        "Very small overlap in a small module; interpret only with exact driver genes stated.",
        "Enrichment is driven by no more than three genes.",
        "Small module; related terms may be statistically dependent and gene-concentrated.",
    ],
    default="No predefined small-overlap or small-module flag; redundancy and gene diversity still require review."
)

all_terms_path = FUNC / "all_prioritized_module_enrichment_terms.csv"
fdr05_path = FUNC / "fdr05_terms_with_overlap_genes.csv"
fdr10_path = FUNC / "fdr10_suggestive_terms.csv"
enrichment.to_csv(all_terms_path, index=False)
fdr05 = enrichment[enrichment["q_category"] == "FDR_LT_0_05"].copy()
fdr10 = enrichment[enrichment["q_category"] == "SUGGESTIVE_0_05_TO_0_10"].copy()
fdr05.to_csv(fdr05_path, index=False)
fdr10.to_csv(fdr10_path, index=False)

# Gene-level evidence, using only retained Phase-2/3B annotations.
mapping_small = mapping[["region", "old_module", "gene_id", "mapped", "symbol"]].drop_duplicates()
core_small = core[[
    "region", "module", "gene_id", "gene_symbol", "kME", "abs_kME", "hub_rank",
    "selected_for_stable_core", "selected_for_top_hub", "is_Tier2",
    "is_high_confidence_ML", "is_high_plus_moderate_ML",
]].copy().rename(columns={"module": "old_module", "gene_symbol": "core_gene_symbol"})
core_small = core_small.drop_duplicates(["region", "old_module", "gene_id"], keep="first")

ml_same_region = ml[[
    "region_or_subregion", "gene_id", "gene_symbol", "ML_priority_class",
    "beta_diagnosis", "p_value_diagnosis", "standardized_effect_size_diagnosis",
]].copy().rename(columns={
    "region_or_subregion": "region",
    "gene_symbol": "same_region_ml_gene_symbol",
    "ML_priority_class": "same_region_ML_priority_class",
})
ml_same_region = ml_same_region.drop_duplicates(["region", "gene_id"], keep="first")

# Phase-3B module enrichment used the global union of corrected ML genes. Retain
# that global status for direct comparability, while also exposing same-region
# status so the two concepts cannot be conflated in interpretation.
ml_global = ml[["region_or_subregion", "gene_id", "gene_symbol", "ML_priority_class"]].copy()
ml_global["global_ml_priority"] = np.select(
    [
        ml_global["ML_priority_class"].fillna("").str.contains("High-confidence"),
        ml_global["ML_priority_class"].fillna("").str.contains("Moderate-confidence"),
        ml_global["ML_priority_class"].fillna("").str.contains("Low-confidence"),
    ], ["high", "moderate", "low"], default="absent"
)
priority_order = {"high": 0, "moderate": 1, "low": 2, "absent": 3}
ml_global["priority_order"] = ml_global["global_ml_priority"].map(priority_order)
global_regions = ml_global.groupby("gene_id")["region_or_subregion"].agg(list_join).rename("ML_priority_regions")
ml_global = (
    ml_global.sort_values(["gene_id", "priority_order"])
    .drop_duplicates("gene_id", keep="first")
    .rename(columns={"gene_symbol": "global_ml_gene_symbol"})
    .merge(global_regions, on="gene_id", how="left", validate="one_to_one")
)
ml_global = ml_global[["gene_id", "global_ml_gene_symbol", "global_ml_priority", "ML_priority_regions"]]

driver_rows = []
for row in enrichment[enrichment["BH_q"] < Q_SUGGESTIVE].itertuples(index=False):
    for gene_id in split_genes(row.overlap_genes):
        driver_rows.append({
            "region": row.region, "module": row.old_module,
            "database": row.database, "ontology": row.ontology,
            "term_id": row.term_id, "term": row.term,
            "raw_p": row.raw_p, "BH_q": row.BH_q,
            "q_category": row.q_category, "theme_rule": row.theme_rule,
            "overlap_count": row.overlap_count,
            "module_genes_mapped": row.module_genes_mapped,
            "term_size_in_background": row.term_size_in_background,
            "gene_id": gene_id,
        })
drivers = pd.DataFrame(driver_rows)
drivers = drivers.merge(
    mapping_small.rename(columns={"old_module": "module"}),
    on=["region", "module", "gene_id"], how="left", validate="many_to_one"
).merge(
    core_small.rename(columns={"old_module": "module"}),
    on=["region", "module", "gene_id"], how="left", validate="many_to_one"
).merge(ml_same_region, on=["region", "gene_id"], how="left", validate="many_to_one").merge(
    ml_global, on="gene_id", how="left", validate="many_to_one"
)
drivers["gene_symbol"] = drivers["symbol"].fillna(drivers["core_gene_symbol"]).fillna(drivers["same_region_ml_gene_symbol"]).fillna(drivers["global_ml_gene_symbol"]).fillna(drivers["gene_id"])
drivers["same_region_ml_priority"] = np.select(
    [
        drivers["same_region_ML_priority_class"].fillna("").str.contains("High-confidence"),
        drivers["same_region_ML_priority_class"].fillna("").str.contains("Moderate-confidence"),
        drivers["same_region_ML_priority_class"].fillna("").str.contains("Low-confidence"),
    ], ["high", "moderate", "low"], default="absent"
)
drivers["ml_priority"] = drivers["global_ml_priority"].fillna("absent")
drivers["region_specific_Tier2"] = drivers["is_Tier2"].map(clean_bool)
drivers["structural_core_member"] = drivers["selected_for_stable_core"].map(clean_bool)
drivers["top50_hub"] = drivers["selected_for_top_hub"].map(clean_bool) | drivers["hub_rank"].fillna(np.inf).le(50)
drivers["gene_family"] = drivers["gene_symbol"].map(family_label)
drivers["directional_statistics_available"] = drivers["beta_diagnosis"].notna()
driver_path = FUNC / "gene_drivers_by_term.csv"
drivers.to_csv(driver_path, index=False)

# Add readable symbols and family-concentration fields to the term tables.
term_gene_summary = drivers.groupby(["region", "module", "database", "term_id"], dropna=False).agg(
    overlap_gene_symbols=("gene_symbol", list_join),
    overlap_gene_ids_verified=("gene_id", list_join),
    high_ML_genes=("gene_symbol", lambda x: list_join(drivers.loc[x.index][drivers.loc[x.index, "ml_priority"] == "high"]["gene_symbol"])),
    moderate_ML_genes=("gene_symbol", lambda x: list_join(drivers.loc[x.index][drivers.loc[x.index, "ml_priority"] == "moderate"]["gene_symbol"])),
    top50_hub_genes=("gene_symbol", lambda x: list_join(drivers.loc[x.index][drivers.loc[x.index, "top50_hub"]]["gene_symbol"])),
    structural_core_genes=("gene_symbol", lambda x: list_join(drivers.loc[x.index][drivers.loc[x.index, "structural_core_member"]]["gene_symbol"])),
).reset_index().rename(columns={"module": "old_module"})

def family_stats(group: pd.DataFrame) -> pd.Series:
    counts = Counter(group["gene_family"])
    non_other = {k: v for k, v in counts.items() if k != "other"}
    if not non_other:
        return pd.Series({"dominant_gene_family": "none", "dominant_family_count": 0, "dominant_family_fraction": 0.0})
    family, count = max(non_other.items(), key=lambda x: x[1])
    return pd.Series({
        "dominant_gene_family": family,
        "dominant_family_count": count,
        "dominant_family_fraction": count / max(len(group), 1),
    })

family_by_term = drivers.groupby(["region", "module", "database", "term_id"], dropna=False).apply(family_stats, include_groups=False).reset_index().rename(columns={"module": "old_module"})
enrichment = enrichment.merge(term_gene_summary, on=["region", "old_module", "database", "term_id"], how="left")
enrichment = enrichment.merge(family_by_term, on=["region", "old_module", "database", "term_id"], how="left")
enrichment.to_csv(all_terms_path, index=False)
enrichment[enrichment["q_category"] == "FDR_LT_0_05"].to_csv(fdr05_path, index=False)
enrichment[enrichment["q_category"] == "SUGGESTIVE_0_05_TO_0_10"].to_csv(fdr10_path, index=False)

recurrence = drivers.groupby(["region", "module", "gene_id", "gene_symbol"], dropna=False).agg(
    fdr05_term_count=("BH_q", lambda x: int((x < Q_FDR).sum())),
    q_lt_0_10_term_count=("term_id", "nunique"),
    databases=("database", list_join),
    themes=("theme_rule", list_join),
    ml_priority=("ml_priority", "first"),
    same_region_ml_priority=("same_region_ml_priority", "first"),
    ML_priority_regions=("ML_priority_regions", "first"),
    region_specific_Tier2=("region_specific_Tier2", "first"),
    abs_kME=("abs_kME", "first"),
    hub_rank=("hub_rank", "first"),
    top50_hub=("top50_hub", "first"),
    structural_core_member=("structural_core_member", "first"),
    gene_family=("gene_family", "first"),
).reset_index().sort_values(["fdr05_term_count", "q_lt_0_10_term_count"], ascending=False)
recurrence.to_csv(FUNC / "gene_driver_recurrence.csv", index=False)

# Redundancy clusters from overlap-gene similarity and direct GO hierarchy.
audit_terms = enrichment[enrichment["BH_q"] < Q_SUGGESTIVE].reset_index(drop=True).copy()
gene_sets = [set(split_genes(x)) for x in audit_terms["overlap_genes"]]
parent_pairs = set(zip(go_relations.get("child_term_id", []), go_relations.get("parent_term_id", [])))
parent_pairs |= {(b, a) for a, b in parent_pairs}

class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x
    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra

cluster_rows = []
for (region, module), idx_values in audit_terms.groupby(["region", "old_module"]).groups.items():
    idxs = list(idx_values)
    uf = UnionFind(len(idxs))
    edge_meta = defaultdict(list)
    for ai in range(len(idxs)):
        for bi in range(ai + 1, len(idxs)):
            i, j = idxs[ai], idxs[bi]
            union = gene_sets[i] | gene_sets[j]
            jac = len(gene_sets[i] & gene_sets[j]) / len(union) if union else 0.0
            exact = gene_sets[i] == gene_sets[j]
            hierarchy = (audit_terms.loc[i, "term_id"], audit_terms.loc[j, "term_id"]) in parent_pairs
            if exact or jac >= 0.80 or hierarchy:
                uf.union(ai, bi)
                edge_meta[(ai, bi)] = [exact, jac, hierarchy]
    groups = defaultdict(list)
    for local_idx, global_idx in enumerate(idxs):
        groups[uf.find(local_idx)].append(global_idx)
    ordered = sorted(groups.values(), key=lambda xs: float(audit_terms.loc[xs, "BH_q"].min()))
    for cluster_number, members in enumerate(ordered, 1):
        sub = audit_terms.loc[members].sort_values("BH_q")
        union_genes = set().union(*(gene_sets[i] for i in members))
        local_positions = {global_i: pos for pos, global_i in enumerate(idxs)}
        edges = []
        for a_i in members:
            for b_i in members:
                if a_i >= b_i:
                    continue
                key = tuple(sorted((local_positions[a_i], local_positions[b_i])))
                if key in edge_meta:
                    edges.append(edge_meta[key])
        cluster_rows.append({
            "region": region, "module": module,
            "cluster_id": f"{re.sub('[^A-Za-z0-9]+', '_', region)}_{module}_C{cluster_number:02d}",
            "representative_term_id": sub.iloc[0]["term_id"],
            "representative_term": sub.iloc[0]["term"],
            "theme_rule": sub.iloc[0]["theme_rule"],
            "term_count": len(sub),
            "terms": list_join(sub["term"]),
            "term_ids": list_join(sub["term_id"]),
            "databases": list_join(sub["database"]),
            "strongest_BH_q": sub["BH_q"].min(),
            "union_overlap_gene_count": len(union_genes),
            "union_overlap_genes": ";".join(sorted(union_genes)),
            "exact_overlap_edge_count": sum(int(x[0]) for x in edges),
            "jaccard_ge_0_8_edge_count": sum(int(x[1] >= .8) for x in edges),
            "direct_GO_hierarchy_edge_count": sum(int(x[2]) for x in edges),
            "redundancy_ratio_terms_per_union_gene": len(sub) / max(len(union_genes), 1),
        })
clusters = pd.DataFrame(cluster_rows)
clusters.to_csv(FUNC / "phase3c_term_clusters.csv", index=False)

# Theme-level aggregation and cross-database consistency.
theme_rows = []
for (region, module, theme), terms in audit_terms.groupby(["region", "old_module", "theme_rule"]):
    term_keys = set(zip(terms["database"], terms["term_id"]))
    dg = drivers[(drivers["region"] == region) & (drivers["module"] == module)]
    dg = dg[dg.apply(lambda r: (r["database"], r["term_id"]) in term_keys, axis=1)]
    genes = sorted(set(dg["gene_symbol"]))
    families = Counter(dg.drop_duplicates("gene_id")["gene_family"])
    non_other = {k: v for k, v in families.items() if k != "other"}
    dom_family, dom_count = max(non_other.items(), key=lambda x: x[1]) if non_other else ("none", 0)
    term_family_max = enrichment[
        (enrichment["region"] == region)
        & (enrichment["old_module"] == module)
        & (enrichment["theme_rule"] == theme)
        & (enrichment["BH_q"] < Q_SUGGESTIVE)
    ]["dominant_family_fraction"].fillna(0).max()
    theme_rows.append({
        "region": region, "module": module, "theme": theme,
        "contributing_terms": list_join(terms["term"]),
        "contributing_term_ids": list_join(terms["term_id"]),
        "databases": list_join(terms["database"]),
        "database_count": terms["database"].nunique(),
        "term_count": len(terms),
        "fdr05_term_count": int((terms["BH_q"] < Q_FDR).sum()),
        "suggestive_term_count": int(((terms["BH_q"] >= Q_FDR) & (terms["BH_q"] < Q_SUGGESTIVE)).sum()),
        "contributing_genes": ";".join(genes),
        "independent_overlap_gene_count": len(genes),
        "strongest_BH_q": terms["BH_q"].min(),
        "multiple_databases_support": terms["database"].nunique() > 1,
        "dominant_gene_family": dom_family,
        "dominant_family_fraction": dom_count / max(dg["gene_id"].nunique(), 1),
        "maximum_single_term_family_fraction": term_family_max,
        "small_overlap_term_fraction": terms["small_overlap_flag"].mean(),
        "small_module": bool(terms["small_module_flag"].any()),
    })
themes = pd.DataFrame(theme_rows)
if len(themes):
    themes["annotation_consistency"] = np.select(
        [
            (themes["database_count"] >= 2) & (themes["fdr05_term_count"] >= 1) & (themes["independent_overlap_gene_count"] >= 3),
            (themes["fdr05_term_count"] >= 1) & (themes["independent_overlap_gene_count"] >= 4),
        ], ["MULTI_DATABASE_CONVERGENCE", "SINGLE_DATABASE_SUPPORT"], default="WEAK/REDUNDANT_SUPPORT"
    )
    themes["main_caveat"] = np.select(
        [
            themes["maximum_single_term_family_fraction"] >= .70,
            themes["independent_overlap_gene_count"] <= 3,
            themes["small_module"],
            themes["database_count"] == 1,
        ], [
            "Most supporting genes belong to one recognized gene family.",
            "Theme is driven by no more than three distinct genes.",
            "Module has no more than 50 mapped genes; term dependence is substantial.",
            "Support is confined to one annotation database.",
        ], default="Interpret as module annotation only; it is not MDD-association evidence."
    )
themes.to_csv(FUNC / "module_theme_summary.csv", index=False)

cross_rows = []
dominant_by_module = {}
for row in priority_keys.itertuples(index=False):
    subset = themes[(themes["region"] == row.region) & (themes["module"] == row.old_module)].copy()
    if subset.empty:
        cross_rows.append({
            "region": row.region, "module": row.old_module, "themes_q_lt_0_10": "",
            "databases_with_q_lt_0_10": "", "cross_database_class": "NO_CLEAR_THEME",
            "dominant_theme": "NO_CLEAR_THEME", "dominant_theme_strongest_q": math.nan,
        })
        dominant_by_module[(row.region, row.old_module)] = None
        continue
    named_subset = subset[(subset["theme"] != "OTHER_ANNOTATION") & (subset["fdr05_term_count"] >= 1) & (subset["independent_overlap_gene_count"] >= 2)]
    selection_pool = named_subset if not named_subset.empty else subset
    selection_pool["selection_score"] = (
        selection_pool["multiple_databases_support"].astype(int) * 100
        + selection_pool["fdr05_term_count"].clip(upper=10) * 5
        + selection_pool["independent_overlap_gene_count"].clip(upper=20) * 2
        - (selection_pool["maximum_single_term_family_fraction"] >= .70).astype(int) * 30
        - selection_pool["small_overlap_term_fraction"] * 5
    )
    dominant = selection_pool.sort_values(["selection_score", "strongest_BH_q"], ascending=[False, True]).iloc[0]
    dominant_by_module[(row.region, row.old_module)] = dominant
    if (subset["annotation_consistency"] == "MULTI_DATABASE_CONVERGENCE").any():
        overall = "MULTI_DATABASE_CONVERGENCE"
    elif (subset["annotation_consistency"] == "SINGLE_DATABASE_SUPPORT").any():
        overall = "SINGLE_DATABASE_SUPPORT"
    else:
        overall = "WEAK/REDUNDANT_SUPPORT"
    cross_rows.append({
        "region": row.region, "module": row.old_module,
        "themes_q_lt_0_10": list_join(subset["theme"]),
        "databases_with_q_lt_0_10": list_join(subset["databases"].str.split(";").explode()),
        "cross_database_class": overall,
        "dominant_theme": dominant["theme"],
        "dominant_theme_strongest_q": dominant["strongest_BH_q"],
    })
cross = pd.DataFrame(cross_rows)
cross.to_csv(FUNC / "cross_database_theme_support.csv", index=False)

# Module-level ML/function and hub/function descriptive overlaps.
ml_func_rows = []
hub_rows = []
for row in priority_keys.itertuples(index=False):
    for threshold_label, threshold in [("BH_q_lt_0_05", Q_FDR), ("BH_q_lt_0_10", Q_SUGGESTIVE)]:
        dg = drivers[(drivers["region"] == row.region) & (drivers["module"] == row.old_module) & (drivers["BH_q"] < threshold)].drop_duplicates("gene_id")
        ml_func_rows.append({
            "region": row.region, "module": row.old_module, "functional_threshold": threshold_label,
            "functional_overlap_unique_genes": len(dg),
            "high_ML_genes": int((dg["ml_priority"] == "high").sum()),
            "moderate_ML_genes": int((dg["ml_priority"] == "moderate").sum()),
            "low_ML_genes": int((dg["ml_priority"] == "low").sum()),
            "absent_from_corrected_ML_ranking": int((dg["ml_priority"] == "absent").sum()),
            "high_gene_symbols": list_join(dg.loc[dg["ml_priority"] == "high", "gene_symbol"]),
            "moderate_gene_symbols": list_join(dg.loc[dg["ml_priority"] == "moderate", "gene_symbol"]),
            "same_region_high_ML_genes": int((dg["same_region_ml_priority"] == "high").sum()),
            "same_region_moderate_ML_genes": int((dg["same_region_ml_priority"] == "moderate").sum()),
            "same_region_high_gene_symbols": list_join(dg.loc[dg["same_region_ml_priority"] == "high", "gene_symbol"]),
            "same_region_moderate_gene_symbols": list_join(dg.loc[dg["same_region_ml_priority"] == "moderate", "gene_symbol"]),
        })
    module_core = core_small[(core_small["region"] == row.region) & (core_small["old_module"] == row.old_module)]
    module_median = module_core["abs_kME"].median()
    for (database, term_id), dg in drivers[(drivers["region"] == row.region) & (drivers["module"] == row.old_module)].groupby(["database", "term_id"]):
        unique = dg.drop_duplicates("gene_id")
        available = unique["abs_kME"].notna().sum()
        hub_rows.append({
            "region": row.region, "module": row.old_module, "database": database,
            "term_id": term_id, "term": unique["term"].iloc[0], "BH_q": unique["BH_q"].iloc[0],
            "overlap_gene_count": len(unique), "kME_available_count": int(available),
            "median_abs_kME_overlap": unique["abs_kME"].median(),
            "median_abs_kME_recorded_module_genes": module_median,
            "top50_hub_count": int(unique["top50_hub"].sum()),
            "top50_hub_fraction": unique["top50_hub"].mean(),
            "structural_core_count": int(unique["structural_core_member"].sum()),
            "structural_core_fraction": unique["structural_core_member"].mean(),
            "descriptive_only": True,
        })
pd.DataFrame(ml_func_rows).to_csv(FUNC / "ml_gene_function_overlap.csv", index=False)
pd.DataFrame(hub_rows).to_csv(FUNC / "hub_function_overlap.csv", index=False)

# Mapping and background audits.
mapping_quality = evidence[[
    "region", "old_module", "module_genes_total", "module_genes_mapped", "mapping_fraction"
]].copy().rename(columns={"old_module": "module"})
unmapped = mapping[mapping["mapped"] == False].groupby(["region", "old_module"])["gene_id"].agg(list_join).reset_index().rename(columns={"old_module": "module", "gene_id": "unmapped_gene_ids"})
mapping_quality = mapping_quality.merge(unmapped, on=["region", "module"], how="left")
mapping_quality["module_genes_unmapped"] = mapping_quality["module_genes_total"] - mapping_quality["module_genes_mapped"]
mapping_quality["mapping_percent"] = 100 * mapping_quality["mapping_fraction"]
mapping_quality["mapping_bias_caution"] = np.where(
    mapping_quality["mapping_fraction"] < .90,
    "Below 90% mapping; interpretation may underrepresent unmapped identifiers.",
    "At least 90% mapped; residual unmapped identifiers remain listed."
)
mapping_quality.to_csv(FUNC / "module_mapping_quality.csv", index=False)

background = enrichment.groupby(["region", "database", "ontology"], dropna=False).agg(
    universe_genes=("universe_genes", "first"), universe_mapped=("universe_mapped", "first"),
    universe_genes_unique_values=("universe_genes", "nunique"),
    universe_mapped_unique_values=("universe_mapped", "nunique"),
    database_eligible_background_denominators=("background_ratio_denominator", lambda x: list_join(str(int(v)) for v in sorted(set(x.dropna())))),
    term_count=("term_id", "count"),
).reset_index()
expected_universe = {
    "Ventral Striatum": 15000, "Insula": 15000, "Hippocampal Formation": 15000,
    "Cingulate gyrus 25": 5000, "dlPFC BA8/9": 5000, "OFC BA11": 5000,
}
background["expected_region_network_genes"] = background["region"].map(expected_universe)
background["region_network_background_confirmed"] = (
    background["universe_genes"] == background["expected_region_network_genes"]
) & (background["universe_genes_unique_values"] == 1) & (background["universe_mapped_unique_values"] == 1)
background.to_csv(FUNC / "background_universe_audit.csv", index=False)

# Exact prioritized module inventory, including structural/core summaries and both ML set results.
core_summary = core_small.groupby(["region", "old_module"]).agg(
    structural_core_gene_count=("selected_for_stable_core", lambda x: int(sum(map(clean_bool, x)))),
    top50_hub_gene_count=("selected_for_top_hub", lambda x: int(sum(map(clean_bool, x)))),
    Tier2_gene_count=("is_Tier2", lambda x: int(sum(map(clean_bool, x)))),
    corrected_ML_high_gene_count=("is_high_confidence_ML", lambda x: int(sum(map(clean_bool, x)))),
    corrected_ML_high_plus_moderate_gene_count=("is_high_plus_moderate_ML", lambda x: int(sum(map(clean_bool, x)))),
).reset_index()
global_ml = ml_enrichment[ml_enrichment["gene_set"].isin(["global_high", "global_high_plus_moderate"])].copy()
global_ml = global_ml.pivot_table(index=["region", "old_module"], columns="gene_set", values=["old_overlap", "old_p", "old_q"], aggfunc="first")
global_ml.columns = [f"{metric}_{gene_set}" for metric, gene_set in global_ml.columns]
global_ml = global_ml.reset_index()
prioritized = evidence.rename(columns={
    "old_module": "module", "reconstructed_beta": "Diagnosis_beta",
    "reconstructed_p": "Diagnosis_raw_p", "reconstructed_q": "Diagnosis_BH_q",
}).merge(core_summary.rename(columns={"old_module": "module"}), on=["region", "module"], how="left")
prioritized = prioritized.merge(global_ml.rename(columns={"old_module": "module"}), on=["region", "module"], how="left")
prioritized.to_csv(FUNC / "prioritized_modules_audited.csv", index=False)

# Publication categories, guardrail wording, and final decision table.
category_rows = []
for row in prioritized.itertuples(index=False):
    dominant = dominant_by_module.get((row.region, row.module))
    if dominant is None:
        category = "AVOID_INTERPRETATION"
        theme = "NO_CLEAR_THEME"
        genes = ""
        dbs = ""
        strongest_q = math.nan
        caveat = "No Phase-3B term reached BH q<0.10."
    else:
        theme = dominant["theme"]
        genes = dominant["contributing_genes"]
        dbs = dominant["databases"]
        strongest_q = dominant["strongest_BH_q"]
        caveat = dominant["main_caveat"]
        module_q10_term_count = len(audit_terms[(audit_terms["region"] == row.region) & (audit_terms["old_module"] == row.module)])
        small_module_redundancy = dominant["small_module"] and module_q10_term_count > (2 * row.module_genes_mapped)
        if dominant["maximum_single_term_family_fraction"] >= .70 or small_module_redundancy or (dominant["small_module"] and dominant["term_count"] > max(dominant["independent_overlap_gene_count"] * 2, 10)):
            category = "SUPPLEMENTARY_CONTEXT"
        elif dominant["fdr05_term_count"] >= 1 and (dominant["database_count"] >= 2 or dominant["independent_overlap_gene_count"] >= 8):
            category = "MAIN_TEXT_CANDIDATE"
        elif dominant["fdr05_term_count"] >= 1 and dominant["independent_overlap_gene_count"] >= 2:
            category = "SUPPLEMENTARY_CONTEXT"
        elif dominant["suggestive_term_count"] >= 1:
            category = "SUPPLEMENTARY_CONTEXT" if dominant["independent_overlap_gene_count"] >= 4 else "DESCRIPTIVE_ONLY"
        else:
            category = "DESCRIPTIVE_ONLY"
    safe_theme = theme.replace("_", " ").lower()
    safe = f"The {row.region} {row.module} gene set was annotated to {safe_theme} terms in the fixed Phase-3B enrichment analysis."
    avoid = f"Avoid claiming that {safe_theme} is dysregulated, causal, replicated, or a validated MDD mechanism in {row.region}."
    category_rows.append({
        "region": row.region, "module": row.module, "dominant_biological_theme": theme,
        "supporting_databases": dbs, "key_driver_genes": genes,
        "strongest_theme_BH_q": strongest_q, "mapping_coverage": row.mapping_fraction,
        "ML_high_overlap": getattr(row, "old_overlap_global_high", math.nan),
        "ML_high_p": getattr(row, "old_p_global_high", math.nan),
        "ML_high_q": getattr(row, "old_q_global_high", math.nan),
        "ML_high_plus_moderate_overlap": getattr(row, "old_overlap_global_high_plus_moderate", row.ml_high_moderate_overlap),
        "ML_high_plus_moderate_p": getattr(row, "old_p_global_high_plus_moderate", row.ml_overlap_p),
        "ML_high_plus_moderate_q": getattr(row, "old_q_global_high_plus_moderate", row.ml_overlap_q),
        "Diagnosis_BH_q": row.Diagnosis_BH_q,
        "Diagnosis_FDR_status": "NOT_BH_FDR_SIGNIFICANT",
        "interpretation_caution": caveat,
        "publication_use_category": category,
        "safe_wording": safe, "wording_to_avoid": avoid,
    })
categories = pd.DataFrame(category_rows)
categories.to_csv(FUNC / "publication_interpretation_categories.csv", index=False)

decision = categories.copy()
decision["overlap_size"] = decision["key_driver_genes"].fillna("").map(lambda x: len([g for g in x.split(";") if g]))
decision[[
    "region", "module", "dominant_biological_theme", "supporting_databases",
    "key_driver_genes", "overlap_size", "mapping_coverage",
    "ML_high_plus_moderate_overlap", "ML_high_plus_moderate_q",
    "Diagnosis_FDR_status", "interpretation_caution", "publication_use_category",
]].to_csv(FUNC / "phase3c_final_decision_table.csv", index=False)

# Human-readable adversarial audit report.
def md_table(frame: pd.DataFrame, columns: list[str]) -> str:
    view = frame[columns].copy().fillna("")
    header = "| " + " | ".join(columns) + " |"
    sep = "|" + "|".join(["---"] * len(columns)) + "|"
    rows = ["| " + " | ".join(str(x).replace("|", "/") for x in row) + " |" for row in view.itertuples(index=False, name=None)]
    return "\n".join([header, sep, *rows])

module_lines = []
for row in categories.itertuples(index=False):
    term_sub = enrichment[(enrichment["region"] == row.region) & (enrichment["old_module"] == row.module) & (enrichment["BH_q"] < Q_SUGGESTIVE)].sort_values("BH_q").head(5)
    exact_terms = "; ".join(
        f"{t.database}: {t.term} (q={t.BH_q:.3g}; n={t.overlap_count}; genes={getattr(t, 'overlap_gene_symbols', '')})"
        for t in term_sub.itertuples(index=False)
    ) or "No term reached BH q<0.10."
    ml_module = core_small[(core_small["region"] == row.region) & (core_small["old_module"] == row.module)].copy()
    high_genes = list_join(ml_module.loc[ml_module["is_high_confidence_ML"].map(clean_bool), "core_gene_symbol"])
    hm_genes = list_join(ml_module.loc[ml_module["is_high_plus_moderate_ML"].map(clean_bool), "core_gene_symbol"])
    module_lines.extend([
        f"### {row.region} {row.module}", "",
        f"- Annotation: **{row.dominant_biological_theme}**; category: **{row.publication_use_category}**.",
        f"- Diagnosis: BH q={row.Diagnosis_BH_q:.3g}; not BH-FDR significant.",
        f"- Corrected ML overlap: high={int(row.ML_high_overlap) if pd.notna(row.ML_high_overlap) else 0}; high+moderate={int(row.ML_high_plus_moderate_overlap) if pd.notna(row.ML_high_plus_moderate_overlap) else 0}; high+moderate q={row.ML_high_plus_moderate_q:.3g}.",
        f"- Exact corrected-ML genes in module: high [{high_genes or 'none'}]; high+moderate [{hm_genes or 'none'}].",
        f"- Leading fixed terms: {exact_terms}",
        f"- Main caution: {row.interpretation_caution}",
        f"- Safe wording: {row.safe_wording}",
        f"- Wording to avoid: {row.wording_to_avoid}", "",
    ])

top_recurrence = recurrence.head(20)[["region", "module", "gene_symbol", "fdr05_term_count", "q_lt_0_10_term_count", "gene_family"]]
family_terms = enrichment[(enrichment["BH_q"] < Q_SUGGESTIVE) & (enrichment["dominant_family_fraction"].fillna(0) >= .70)][["region", "old_module", "term", "BH_q", "dominant_gene_family", "dominant_family_fraction", "overlap_gene_symbols"]].head(30)
main_claims = categories[categories["publication_use_category"] == "MAIN_TEXT_CANDIDATE"]
supp_claims = categories[categories["publication_use_category"] == "SUPPLEMENTARY_CONTEXT"]
avoid_claims = categories[categories["publication_use_category"].isin(["DESCRIPTIVE_ONLY", "AVOID_INTERPRETATION"])]

report = [
    "# WGCNA Functional-Enrichment and Biological-Interpretation Audit", "",
    "## Executive summary", "",
    f"Phase 3C audited the fixed Phase-3B enrichment output without rerunning WGCNA or enrichment. The ten prespecified modules contained {len(fdr05):,} terms with BH q<0.05 and {len(fdr10):,} additional terms with 0.05≤q<0.10. Functional annotation remains separate from Diagnosis association: zero modules have Diagnosis BH q<0.10.", "",
    "## Audit scope", "",
    f"Git baseline: `{BASELINE}`. Databases were restricted to GO BP/MF/CC, Reactome, and MSigDB Hallmark. The audit added gene tracing, redundancy grouping, family-concentration flags, and publication wording guardrails; it did not calculate new enrichment p-values.", "",
    "## Prioritized modules", "",
    md_table(prioritized, ["region", "module", "module_genes_total", "module_genes_mapped", "mapping_fraction", "Diagnosis_beta", "Diagnosis_raw_p", "Diagnosis_BH_q", "structural_core_gene_count", "Tier2_gene_count", "corrected_ML_high_gene_count", "corrected_ML_high_plus_moderate_gene_count"]), "",
    "## Mapping quality", "",
    "Phase-3B reported 2,611/2,783 unique mapped identifiers (93.8%) across the prioritized collection. Per-module accounting is retained below; unmapped identifiers are listed rather than discarded.", "",
    md_table(mapping_quality, ["region", "module", "module_genes_total", "module_genes_mapped", "module_genes_unmapped", "mapping_percent", "mapping_bias_caution"]), "",
    "## FDR-supported functional terms", "",
    f"There are {len(fdr05)} BH-q<0.05 terms. A separate file contains the {len(fdr10)} suggestive terms; the categories are never merged. Every row retains exact overlap genes, module and background sizes, raw p, BH q, and ratios.", "",
    "## Term redundancy", "",
    f"The {len(audit_terms)} q<0.10 terms collapse into {len(clusters)} overlap/hierarchy clusters using exact overlap, overlap-gene Jaccard ≥0.80, or an immediate GO parent-child relation. Cluster counts are descriptive; related terms are not treated as independent findings.", "",
    "## Gene drivers", "", md_table(top_recurrence, list(top_recurrence.columns)), "",
    "## Gene-family effects", "",
    "Recognized family concentration was flagged for taste/olfactory receptors, protocadherins, histones, ribosomal genes, immunoglobulins, and small nuclear/nucleolar RNA genes. A family flag does not declare an artifact; it limits how literally broad ontology labels should be interpreted.", "",
    (md_table(family_terms, list(family_terms.columns)) if len(family_terms) else "No q<0.10 term had ≥70% of drivers in a recognized family."), "",
    "## Cross-database convergence", "", md_table(cross, ["region", "module", "dominant_theme", "databases_with_q_lt_0_10", "cross_database_class"]), "",
    "## ML-gene contribution to functional themes", "",
    "The table `ml_gene_function_overlap.csv` reports both global-union and same-region corrected-ML status for unique functional-driver genes. The global columns match the gene universe used by the fixed Phase-3B module-over-representation analysis; same-region columns prevent cross-region priority from being mistaken for local evidence. ML overlap and pathway annotation are overlapping evidence layers, not independent replications.", "",
    "## Hub/core contribution", "",
    "`hub_function_overlap.csv` gives median abs(kME), top-50-hub fraction, and structural-core fraction for each q<0.10 term. Missing kME values remain missing; no inferential hub-enrichment test was introduced.", "",
    "## Module-by-module interpretation", "", *module_lines,
    "## Claims suitable for main text", "",
    (md_table(main_claims, ["region", "module", "dominant_biological_theme", "safe_wording", "interpretation_caution"]) if len(main_claims) else "None met the prespecified annotation-coherence rule."), "",
    "## Claims suitable only for supplement", "",
    (md_table(supp_claims, ["region", "module", "dominant_biological_theme", "safe_wording", "interpretation_caution"]) if len(supp_claims) else "None."), "",
    "## Claims to avoid", "",
    (md_table(avoid_claims, ["region", "module", "dominant_biological_theme", "wording_to_avoid"]) if len(avoid_claims) else "No additional module-level exclusions."), "",
    "Across all modules, avoid causal, biomarker, replicated-mechanism, or Diagnosis-significant wording. Annotation coherence does not override the absence of module-eigengene Diagnosis BH-FDR significance.", "",
    "## Remaining biological uncertainties", "",
    "- Functional databases are overlapping and hierarchically redundant.",
    "- Small modules and small overlap sets can yield very low q-values while remaining gene-concentrated.",
    "- Several gene-level Diagnosis statistics are unavailable outside the Tier2/corrected-ML table and were not reconstructed.",
    "- Mapping improved substantially, but unmapped identifiers remain and may not be biologically random.",
    "- No literature search or external biological validation was used to rescue or strengthen a theme.", "",
    "## Final decision table", "",
    md_table(decision, ["region", "module", "dominant_biological_theme", "supporting_databases", "key_driver_genes", "overlap_size", "mapping_coverage", "ML_high_plus_moderate_q", "Diagnosis_FDR_status", "interpretation_caution", "publication_use_category"]), "",
]
DOC.write_text("\n".join(report), encoding="utf-8")

# Machine-readable provenance and recursive output inventory (excluding self-reference).
output_files = sorted(p for p in OUT.rglob("*") if p.is_file() and p.name not in {"phase3c_manifest.json", "phase3c_output_file_manifest.csv"})
output_inventory = pd.DataFrame([
    {"path": str(p.resolve()).replace("\\", "/"), "relative_path": str(p.relative_to(OUT)).replace("\\", "/"), "bytes": p.stat().st_size, "sha256": sha256(p)}
    for p in output_files
])
output_inventory.to_csv(MANIFESTS / "phase3c_output_file_manifest.csv", index=False)

manifest = {
    "stage": "WGCNA functional-enrichment and biological-interpretation audit",
    "status": "interpretation_audit",
    "git_baseline": BASELINE,
    "phase3b_source_files": [{"path": str(p.resolve()).replace("\\", "/"), "sha256": sha256(p)} for p in phase3b_source_paths],
    "scripts": [
        {"path": "code/phase3c_export_go_relations.R", "sha256": sha256(ROOT / "code/phase3c_export_go_relations.R")},
        {"path": "code/phase3c_audit_wgcna_interpretation.py", "sha256": sha256(Path(__file__))},
    ],
    "commands": [
        "Rscript code/phase3c_export_go_relations.R",
        "python code/phase3c_audit_wgcna_interpretation.py",
    ],
    "no_new_enrichment": True,
    "no_WGCNA_rerun": True,
    "thresholds": {"FDR_supported": "BH_q < 0.05", "suggestive": "0.05 <= BH_q < 0.10", "redundancy_gene_set_jaccard": 0.80},
    "databases_already_present_in_phase3b": sorted(ALLOWED_DATABASES),
    "prioritized_module_count": 10,
    "term_counts": {"BH_q_lt_0_05": int(len(fdr05)), "BH_q_0_05_to_0_10": int(len(fdr10))},
    "outputs": output_inventory.to_dict(orient="records"),
    "report": {"path": str(DOC.resolve()).replace("\\", "/"), "sha256": sha256(DOC)},
}
(MANIFESTS / "phase3c_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

print(json.dumps({
    "prioritized_modules": len(priority_keys),
    "FDR_q_lt_0_05_terms": len(fdr05),
    "suggestive_q_lt_0_10_terms": len(fdr10),
    "term_clusters": len(clusters),
    "output_files_hashed": len(output_inventory),
}, indent=2))
