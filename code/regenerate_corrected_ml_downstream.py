"""Regenerate only downstream products that depend on ML-priority membership.

This script does not rebuild AHBA/WGCNA modules, fit ML models, or alter Tier2.
It writes versioned corrected outputs under 16_corrected_ML_downstream.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, hypergeom

from canonical_downstream_paths import (
    CANONICAL_RANKING,
    CANONICAL_RANKING_SHA256,
    CLEAN_ROOT,
    LEGACY_RANKING,
    PHASE2_ROOT,
    PROJECT_ROOT,
)

HIGH = "High-confidence ML-prioritized Tier 2"
MODERATE = "Moderate-confidence ML-prioritized Tier 2"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bh(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce").to_numpy(float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    if not len(pv):
        return pd.Series(out, index=values.index)
    order = np.argsort(pv)
    ranked = pv[order]
    adjusted = np.minimum.accumulate((ranked * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
    restored = np.empty(len(pv))
    restored[order] = np.minimum(adjusted, 1.0)
    out[valid] = restored
    return pd.Series(out, index=values.index)


def norm_symbol(value: object) -> str:
    return "" if pd.isna(value) else str(value).strip().upper()


def gene_sets(ranking: pd.DataFrame) -> dict[str, set[str]]:
    ids = ranking.assign(gene_id=ranking["gene_id"].astype(str))
    return {
        "high_confidence_ML_genes": set(ids.loc[ids["ML_priority_class"].eq(HIGH), "gene_id"]),
        "high_plus_moderate_ML_genes": set(ids.loc[ids["ML_priority_class"].isin([HIGH, MODERATE]), "gene_id"]),
    }


def region_key(value: object) -> str:
    text = str(value).lower().replace("/", "_").replace(" ", "_")
    if "cingulate" in text or text == "cg25":
        return "cg25"
    if "dlpfc" in text:
        return "dlpfc_ba8_9"
    if "ofc" in text:
        return "ofc_ba11"
    return text


def dependency_inventory() -> pd.DataFrame:
    rows = [
        ("AHBA preprocessing", "legacy all-ranked symbols", "all 671 Tier2 pairs unchanged", "B", "preserve", "processed AHBA matrix", "Only the confidence labels changed; the all-ranked symbol universe did not."),
        ("De-novo AHBA module construction", "processed AHBA expression", "independent of ML confidence", "B", "preserve", "AHBA module assignments", "No ML list in module construction."),
        ("Published AHBA high/high+moderate enrichment", str(LEGACY_RANKING), "confidence membership", "A", str(CANONICAL_RANKING), "AHBA_corrected", "Candidate sets changed."),
        ("Published AHBA all-Tier2 enrichment", "Tier2 candidate table", "unchanged Tier2", "B", "preserve", "legacy all-Tier2 rows", "No corrected-score dependency."),
        ("WGCNA module construction", "normalised expression", "independent of ML confidence", "B", "preserve", "saved module assignments/eigengenes", "Modules must remain identical."),
        ("WGCNA ML enrichment", str(LEGACY_RANKING), "confidence membership", "A", str(CANONICAL_RANKING), "WGCNA_ML_annotations_corrected", "ML gene sets changed."),
        ("WGCNA structural stable core", "module membership and kME", "independent of ML confidence", "B", "preserve", "saved structural core", "ML flags are post-selection annotations."),
        ("Stable-core ML annotation", str(LEGACY_RANKING), "confidence membership", "A", str(CANONICAL_RANKING), "stable_core_annotations_corrected", "Refresh flags only."),
        ("ML-labelled network nodes", str(LEGACY_RANKING), "node annotation", "A", str(CANONICAL_RANKING), "network_node_annotations_corrected", "Base networks/edges preserved."),
        ("AHBA module-score ML", "fixed published modules", "features independent of ML gene shortlist", "B", "preserve", "module_score_ml", "No shortlist-based feature selection."),
        ("MDD module-level ML2", "fixed WGCNA module scores", "features independent of ML gene shortlist", "B", "preserve", "13_Module_Level_ML2", "Structural modules are fixed; confidence flags are annotation only."),
        ("Tier2 enrichment", "unchanged Tier2", "independent", "B", "preserve", "08_Tier2_enrichment", "No rerun."),
        ("Tier4 ranked GSEA", "unchanged Tier4 ranking", "independent", "B", "preserve", "ranked GSEA outputs", "No rerun."),
        ("Collected images/manuscript package", "mixed legacy downstream outputs", "depends on corrected figures/reports", "C", str(PHASE2_ROOT), "future corrected publication package", "Refresh only after corrected calculations are reviewed."),
    ]
    return pd.DataFrame(rows, columns=["artifact_script", "old_ML_input_path", "actual_dependency", "class", "corrected_input_or_action", "expected_output", "reason"])


def ahba_candidate_sets(ranking: pd.DataFrame) -> list[dict]:
    ranking = ranking.copy()
    ranking["gene_symbol_norm"] = ranking["gene_symbol"].map(norm_symbol)
    direction = "statistical_direction" if "statistical_direction" in ranking.columns else "direction"
    high = ranking[ranking["ML_priority_class"].eq(HIGH)]
    hm = ranking[ranking["ML_priority_class"].isin([HIGH, MODERATE])]
    sets = []
    for name, scope, frame in [("high_confidence", "all_regions", high), ("high_plus_moderate", "all_regions", hm)]:
        sets.append({"candidate_set": name, "scope": scope, "pairs": len(frame), "genes": set(frame["gene_symbol_norm"]) - {""}})
    for region, sub in hm.groupby("region_or_subregion"):
        sets.append({"candidate_set": "high_plus_moderate", "scope": f"region:{region}", "pairs": len(sub), "genes": set(sub["gene_symbol_norm"]) - {""}})
        for direct, dsub in sub.groupby(direction):
            genes = set(dsub["gene_symbol_norm"]) - {""}
            if len(genes) >= 5:
                sets.append({"candidate_set": "high_plus_moderate", "scope": f"region:{region}|direction:{direct}", "pairs": len(dsub), "genes": genes})
    return sets


def run_ahba(ranking: pd.DataFrame, label: str) -> pd.DataFrame:
    base = CLEAN_ROOT / "11_AHBA_modules_enrichment" / "published_modules"
    modules = pd.read_csv(base / "ahba_published_gene_to_module_assignments.csv")
    modules["gene_symbol_norm"] = modules["gene_symbol"].map(norm_symbol)
    modules["published_module"] = modules["published_module"].astype(str).str.strip()
    background_file = base / "enrichment" / "published_ahba_mdd_background_genes.csv"
    background = set(pd.read_csv(background_file).iloc[:, 0].map(norm_symbol)) - {""}
    modules = modules[modules["gene_symbol_norm"].isin(background)].drop_duplicates(["gene_symbol_norm", "published_module"])
    groups = {m: set(g["gene_symbol_norm"]) for m, g in modules.groupby("published_module")}
    rows = []
    for cset in ahba_candidate_sets(ranking):
        genes = cset["genes"] & background
        for module, module_genes in groups.items():
            overlap = genes & module_genes
            M, N, K, x = len(background), len(genes), len(module_genes), len(overlap)
            p = hypergeom.sf(x - 1, M, K, N) if x else 1.0
            rows.append({"version": label, "candidate_set": cset["candidate_set"], "scope": cset["scope"], "original_gene_region_pair_count": cset["pairs"], "unique_candidate_genes": len(cset["genes"]), "published_module": module, "background_genes": M, "candidate_genes_in_background": N, "module_genes_in_background": K, "overlap_genes": x, "overlap_gene_symbols": ";".join(sorted(overlap)), "p_value_hypergeometric": p, "odds_ratio_approx": ((x + .5)/(N-x+.5))/((K-x+.5)/(M-N-K+x+.5))})
    out = pd.DataFrame(rows)
    out["fdr_BH"] = out.groupby(["candidate_set", "scope"], group_keys=False)["p_value_hypergeometric"].apply(bh)
    return out


def compare_ahba(old: pd.DataFrame, new: pd.DataFrame) -> pd.DataFrame:
    keys = ["candidate_set", "scope", "published_module"]
    cols = keys + ["original_gene_region_pair_count", "unique_candidate_genes", "candidate_genes_in_background", "overlap_genes", "overlap_gene_symbols", "p_value_hypergeometric", "odds_ratio_approx", "fdr_BH"]
    merged = old[cols].merge(new[cols], on=keys, suffixes=("_old", "_corrected"), how="outer")
    def genes(x): return set(str(x).split(";")) - {"", "nan"}
    merged["genes_lost"] = [";".join(sorted(genes(a)-genes(b))) for a,b in zip(merged["overlap_gene_symbols_old"], merged["overlap_gene_symbols_corrected"])]
    merged["genes_gained"] = [";".join(sorted(genes(b)-genes(a))) for a,b in zip(merged["overlap_gene_symbols_old"], merged["overlap_gene_symbols_corrected"])]
    merged["old_significance"] = np.select([merged["fdr_BH_old"] < .10, merged["p_value_hypergeometric_old"] < .05], ["FDR<0.10", "nominal_p<0.05"], "none")
    merged["corrected_significance"] = np.select([merged["fdr_BH_corrected"] < .10, merged["p_value_hypergeometric_corrected"] < .05], ["FDR<0.10", "nominal_p<0.05"], "none")
    merged["significance_changed"] = merged["old_significance"] != merged["corrected_significance"]
    return merged


def wgcna_enrichment(assignments: pd.DataFrame, ranking: pd.DataFrame, label: str) -> pd.DataFrame:
    assignments = assignments.copy()
    assignments["gene_id"] = assignments["gene_id"].astype(str)
    ranking = ranking.copy()
    ranking["gene_id"] = ranking["gene_id"].astype(str)
    global_sets = gene_sets(ranking)
    rows = []
    for region, reg in assignments.groupby("region"):
        bg = set(reg["gene_id"])
        same = ranking[ranking["region_or_subregion"].map(region_key).eq(region_key(region))]
        same_hm = set(same.loc[same["ML_priority_class"].isin([HIGH, MODERATE]), "gene_id"]) & bg
        sets = {k: v & bg for k, v in global_sets.items()}
        sets["same_region_high_plus_moderate_ML_genes"] = same_hm
        for module, mod in reg.groupby("module"):
            mg = set(mod["gene_id"])
            for name, gs in sets.items():
                overlap = mg & gs
                table = [[len(overlap), len(mg)-len(overlap)], [len(gs)-len(overlap), len(bg)-len(mg)-len(gs)+len(overlap)]]
                odds, p = fisher_exact(table, alternative="greater")
                rows.append({"version": label, "region": region, "module": module, "module_role": "unassigned_grey" if str(module).lower()=="megrey" else "coexpression_module", "gene_set": name, "gene_set_family": name.replace("_genes", ""), "module_gene_count": len(mg), "gene_set_count_in_background": len(gs), "overlap_count": len(overlap), "overlap_genes": ";".join(sorted(overlap)), "odds_ratio": odds, "p_value": p})
    out = pd.DataFrame(rows)
    out["fdr_BH"] = np.nan
    for (_, _), idx in out.groupby(["region", "gene_set"]).groups.items():
        nongrey = [i for i in idx if out.at[i, "module_role"] != "unassigned_grey"]
        out.loc[nongrey, "fdr_BH"] = bh(out.loc[nongrey, "p_value"])
    out["significance_category"] = np.select([out["fdr_BH"] < .10, out["p_value"] < .05], ["FDR<0.10", "nominal_p<0.05"], "none")
    return out


def compare_wgcna(old: pd.DataFrame, new: pd.DataFrame) -> pd.DataFrame:
    keys = ["region", "module", "gene_set"]
    cols = keys + ["gene_set_count_in_background", "overlap_count", "overlap_genes", "odds_ratio", "p_value", "fdr_BH", "significance_category"]
    merged = old[cols].merge(new[cols], on=keys, suffixes=("_old", "_corrected"), how="outer")
    def genes(x): return set(str(x).split(";")) - {"", "nan"}
    merged["genes_lost"] = [";".join(sorted(genes(a)-genes(b))) for a,b in zip(merged["overlap_genes_old"], merged["overlap_genes_corrected"])]
    merged["genes_gained"] = [";".join(sorted(genes(b)-genes(a))) for a,b in zip(merged["overlap_genes_old"], merged["overlap_genes_corrected"])]
    merged["significance_changed"] = merged["significance_category_old"] != merged["significance_category_corrected"]
    return merged


def annotate_structural_core(source: Path, ranking: pd.DataFrame, output: Path) -> pd.DataFrame:
    table = pd.read_csv(source, dtype={"gene_id": str})
    sets = gene_sets(ranking)
    table["is_high_confidence_ML"] = table["gene_id"].astype(str).isin(sets["high_confidence_ML_genes"])
    table["is_high_plus_moderate_ML"] = table["gene_id"].astype(str).isin(sets["high_plus_moderate_ML_genes"])
    table.to_csv(output, index=False)
    return table


def main() -> None:
    if sha256(CANONICAL_RANKING) != CANONICAL_RANKING_SHA256:
        raise RuntimeError("Canonical ranking checksum mismatch; refusing downstream regeneration.")
    PHASE2_ROOT.mkdir(parents=True, exist_ok=True)
    ahba_dir = PHASE2_ROOT / "AHBA_corrected"
    wgcna_dir = PHASE2_ROOT / "WGCNA_ML_annotations_corrected"
    core_dir = PHASE2_ROOT / "stable_core_annotations_corrected"
    for folder in [ahba_dir, wgcna_dir, core_dir]: folder.mkdir(parents=True, exist_ok=True)
    old = pd.read_csv(LEGACY_RANKING, dtype={"gene_id": str})
    new = pd.read_csv(CANONICAL_RANKING, dtype={"gene_id": str})
    dependency_inventory().to_csv(PHASE2_ROOT / "downstream_dependency_inventory.csv", index=False)

    old_ahba, new_ahba = run_ahba(old, "legacy"), run_ahba(new, "corrected")
    new_ahba.to_csv(ahba_dir / "published_AHBA_ML_enrichment_corrected.csv", index=False)
    ahba_cmp = compare_ahba(old_ahba, new_ahba)
    ahba_cmp.to_csv(ahba_dir / "published_AHBA_ML_enrichment_before_after.csv", index=False)

    assignment_sources = {
        "broad_regions": CLEAN_ROOT / "12_MDD_WGCNA_modules" / "best_solution_gene_to_module_assignments.csv",
        "frontal_subregions": CLEAN_ROOT / "mdd_wgcna_modules" / "frontal_subregion_wgcna" / "frontal_subregion_best_solution_gene_to_module_assignments.csv",
    }
    wgcna_comparisons = []
    assignment_hashes = {}
    for unit, path in assignment_sources.items():
        assignments = pd.read_csv(path, dtype={"gene_id": str})
        assignment_hashes[unit] = sha256(path)
        old_e = wgcna_enrichment(assignments, old, "legacy")
        new_e = wgcna_enrichment(assignments, new, "corrected")
        new_e.to_csv(wgcna_dir / f"{unit}_ML_module_enrichment_corrected.csv", index=False)
        cmp = compare_wgcna(old_e, new_e)
        cmp.insert(0, "analysis_unit", unit)
        cmp.to_csv(wgcna_dir / f"{unit}_ML_module_enrichment_before_after.csv", index=False)
        wgcna_comparisons.append(cmp)
        ann = assignments.copy()
        sets = gene_sets(new)
        ann["is_corrected_ML_high"] = ann["gene_id"].astype(str).isin(sets["high_confidence_ML_genes"])
        ann["is_corrected_ML_high_plus_moderate"] = ann["gene_id"].astype(str).isin(sets["high_plus_moderate_ML_genes"])
        ann.to_csv(wgcna_dir / f"{unit}_network_node_annotations_corrected.csv", index=False)
    pd.concat(wgcna_comparisons, ignore_index=True).to_csv(wgcna_dir / "all_WGCNA_ML_enrichment_before_after.csv", index=False)

    core_sources = {
        "broad": CLEAN_ROOT / "12_MDD_WGCNA_modules" / "module_level_ml2" / "stable_core_genes_by_module.csv",
        "frontal": CLEAN_ROOT / "12_MDD_WGCNA_modules" / "frontal_subregion_wgcna" / "module_level_ml2" / "frontal_subregion_stable_core_genes_by_module.csv",
    }
    core_summary = []
    for name, source in core_sources.items():
        if not source.exists(): continue
        annotated = annotate_structural_core(source, new, core_dir / f"{name}_structural_core_with_corrected_ML_annotations.csv")
        stable = annotated[annotated["gene_set_type"].eq("stable_core_genes")] if "gene_set_type" in annotated else annotated
        summary = stable.groupby(["region", "module"], as_index=False).agg(n_structural_core=("gene_id", "nunique"), n_corrected_ML_high=("is_high_confidence_ML", "sum"), n_corrected_ML_high_plus_moderate=("is_high_plus_moderate_ML", "sum"))
        summary.insert(0, "analysis_unit", name)
        core_summary.append(summary)
    if core_summary:
        pd.concat(core_summary, ignore_index=True).to_csv(core_dir / "structural_core_corrected_ML_annotation_summary.csv", index=False)

    manifest = {
        "upstream_canonical_commit": "fa3b53f71a9373cd9856df127991abce89bf52df",
        "canonical_ranking": str(CANONICAL_RANKING),
        "canonical_ranking_sha256": CANONICAL_RANKING_SHA256,
        "counts": {"total": len(new), "high": int(new["ML_priority_class"].eq(HIGH).sum()), "moderate": int(new["ML_priority_class"].eq(MODERATE).sum()), "high_plus_moderate": int(new["ML_priority_class"].isin([HIGH, MODERATE]).sum())},
        "assignment_input_sha256": assignment_hashes,
        "script": str(Path(__file__).resolve()),
        "random_seed": None,
        "methods": {"AHBA": "hypergeometric over preserved published-AHBA/MDD background; BH within candidate set/scope", "WGCNA": "one-sided Fisher exact over preserved region network background; BH within region/gene set", "stable_core": "structural core preserved; corrected ML flags reannotated only"},
        "interpretation": "Corrected-list-dependent downstream annotation only; no modules, tiers, expression models, or ML models rebuilt.",
    }
    (PHASE2_ROOT / "phase2_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
