"""Read-only numerical verification for the publication end-to-end audit.

This script reads canonical or baseline-committed inputs and writes one JSON
evidence file.  It does not modify scientific pipeline outputs.
"""

from __future__ import annotations

import hashlib
import io
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, hypergeom, t as student_t

ROOT = Path(r"D:/genes/GeneDepresion")
CLEAN = ROOT / "outputs/length_gt_100_clean_pipeline"
ELISA = ROOT / "outputs/length_gt_100_elisa_style"
BASELINE = "ca25ec35af359a4d409abf0f728833d2ff54761c"
OUT = ROOT / "docs/PUBLICATION_AUDIT_VERIFICATION.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def baseline_csv(path: str, **kwargs) -> pd.DataFrame:
    raw = subprocess.check_output(["git", "show", f"{BASELINE}:{path}"], cwd=ROOT)
    return pd.read_csv(io.BytesIO(raw), **kwargs)


def bh(pvalues: pd.Series) -> np.ndarray:
    p = pd.to_numeric(pvalues, errors="coerce").to_numpy(float)
    result = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    adj = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    restored = np.empty(len(ranked))
    restored[order] = np.minimum(adj, 1.0)
    result[valid] = restored
    return result


def sample_accounting(meta: pd.DataFrame) -> list[dict]:
    groups = {
        "Ventral Striatum": meta.AnalysisGroup.eq("Ventral_Striatum"),
        "Insula": meta.AnalysisGroup.eq("Insula"),
        "Hippocampal Formation": meta.AnalysisGroup.eq("Hippocampal_Formation"),
        "pooled Frontal Cortex": meta.AnalysisGroup.eq("Frontal_Cortex"),
        "Cingulate gyrus 25": meta.OriginalRegion.eq("Cingulate gyrus 25 (Cg25)"),
        "dlPFC BA8/9": meta.OriginalRegion.eq("Dorsolateral prefrontal cortex (dlPFC; BA8/9)"),
        "OFC BA11": meta.OriginalRegion.eq("Orbitofrontal (OFC; BA11)"),
    }
    rows = []
    for label, mask in groups.items():
        d = meta.loc[mask].copy()
        counts = d.BrainID.astype(str).value_counts()
        rows.append({
            "region": label,
            "n_samples": len(d),
            "n_MDD": int(d.Diagnosis.eq("MDD").sum()),
            "n_Control": int(d.Diagnosis.isin(["CTRL", "Control"]).sum()),
            "n_donors": int(d.BrainID.astype(str).nunique()),
            "repeated_donors": int((counts > 1).sum()),
            "max_samples_per_donor": int(counts.max()) if len(counts) else 0,
            "duplicate_sample_ids": int(d.SampleID.astype(str).duplicated().sum()),
            "missing_Age": int(d.Age.isna().sum()),
            "missing_Sex": int(d.Sex.isna().sum()),
            "missing_RIN": int(d.RIN.isna().sum()),
            "missing_PMI": int(d.PMI.isna().sum()),
            "missing_pH": int(d.pH.isna().sum()),
        })
    return rows


def normalization_check() -> dict:
    base = ROOT / "outputs/step1_length_sensitivity/Length_gt_100"
    counts = pd.read_csv(base / "03_filtered_raw_counts/keep_all_samples/filtered_raw_counts_ventral_striatum.csv", index_col=0)
    annot = pd.read_csv(base / "05_gene_annotation/keep_all_samples/gene_annotation_ventral_striatum.csv", dtype={"GeneID": str})
    annot = annot.drop_duplicates("GeneID").set_index("GeneID")
    factors = pd.read_csv(base / "04_normalization_factors/keep_all_samples/GeTMM_normalization_factors_ventral_striatum.csv")
    norm = pd.read_csv(base / "02_normalized_expression/keep_all_samples/normalized_GeTMM_log2_ventral_striatum.csv", index_col=0)
    counts.index = counts.index.astype(str)
    norm.index = norm.index.astype(str)
    sample = str(factors.iloc[0].SampleID)
    gene = next(g for g in counts.index if float(counts.at[g, sample]) > 0 and g in annot.index)
    raw_count = float(counts.at[gene, sample])
    length_bp = float(annot.at[gene, "Length"])
    rpk = raw_count / (length_bp / 1000.0)
    factor_row = factors.loc[factors.SampleID.astype(str).eq(sample)].iloc[0]
    effective = float(factor_row.effective_library_size)
    getmm = rpk / effective * 1e6
    reconstructed = math.log2(getmm + 1)
    stored = float(norm.at[gene, sample])
    return {
        "region": "Ventral Striatum", "gene_id": gene, "sample_id": sample,
        "raw_count": raw_count, "length_bp": length_bp, "rpk": rpk,
        "rpk_library_size": float(factor_row.RPK_library_size),
        "TMM_factor": float(factor_row.TMM_factor),
        "effective_library_size": effective, "GeTMM": getmm,
        "reconstructed_log2_GeTMM_plus_1": reconstructed,
        "canonical_value": stored, "absolute_difference": abs(reconstructed - stored),
    }


def lm_checks(meta: pd.DataFrame) -> list[dict]:
    expr = pd.read_csv(ROOT / "outputs/step1_length_sensitivity/Length_gt_100/02_normalized_expression/keep_all_samples/normalized_GeTMM_log2_insula.csv", index_col=0)
    expr.index = expr.index.astype(str)
    result = pd.read_csv(ELISA / "LMM_all_genes_length_gt_100.csv", dtype={"gene_id": str})
    result = result[result.brain_region.eq("Insula")].copy()
    targets = [result[result.beta_diagnosis > 0].sort_values("p_value_diagnosis").iloc[0], result[result.beta_diagnosis < 0].sort_values("p_value_diagnosis").iloc[0]]
    m = meta[meta.AnalysisGroup.eq("Insula") & meta.SampleID.astype(str).isin(expr.columns)].copy()
    m = m.set_index("SampleID").loc[expr.columns].reset_index()
    m = m.rename(columns={m.columns[0]: "SampleID"})
    m = m.dropna(subset=["Diagnosis", "Age", "Sex", "RIN"])
    sex_levels = sorted(m.Sex.astype(str).unique())
    xcols = [np.ones(len(m)), m.Diagnosis.eq("MDD").astype(float).to_numpy(), m.Age.astype(float).to_numpy()]
    names = ["Intercept", "DiagnosisMDD", "Age"]
    for level in sex_levels[1:]:
        xcols.append(m.Sex.astype(str).eq(level).astype(float).to_numpy())
        names.append(f"Sex{level}")
    xcols.append(m.RIN.astype(float).to_numpy()); names.append("RIN")
    X = np.column_stack(xcols)
    out = []
    for target in targets:
        y = expr.loc[str(target.gene_id), m.SampleID.astype(str)].astype(float).to_numpy()
        coef, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ coef
        df = len(y) - X.shape[1]
        sigma2 = float(resid @ resid / df)
        cov = sigma2 * np.linalg.inv(X.T @ X)
        se = float(np.sqrt(cov[1, 1]))
        stat = float(coef[1] / se)
        p = float(2 * student_t.sf(abs(stat), df))
        out.append({
            "gene_id": str(target.gene_id), "gene_symbol": target.gene_symbol,
            "n": len(y), "design_columns": names, "df_residual": df,
            "reconstructed_beta": float(coef[1]), "canonical_beta": float(target.beta_diagnosis),
            "beta_abs_difference": abs(float(coef[1]) - float(target.beta_diagnosis)),
            "reconstructed_se": se, "canonical_se": float(target.standard_error_diagnosis),
            "se_abs_difference": abs(se - float(target.standard_error_diagnosis)),
            "reconstructed_t": stat, "reconstructed_p": p,
            "canonical_p": float(target.p_value_diagnosis),
            "p_abs_difference": abs(p - float(target.p_value_diagnosis)),
        })
    return out


def tier_and_effect_checks() -> dict:
    data = pd.read_csv(CLEAN / "04_standardized_effect_size/LMM_all_genes_with_standardized_effect_size_4groups.csv", dtype={"gene_id": str})
    row = data.sort_values("p_value_diagnosis").iloc[0]
    std = float(row.beta_diagnosis) / float(row.sd_expression_gene)
    tier1 = (data.qvalue_Storey < .10) & (data.abs_standardized_effect_size_diagnosis > .8)
    tier2 = (data.p_value_diagnosis < .01) & (data.abs_standardized_effect_size_diagnosis > .8)
    tier3 = (data.p_value_diagnosis < .01) & (data.abs_beta_diagnosis > .2)
    saved = pd.read_csv(CLEAN / "05_evidence_tiers/evidence_tiers_all_genes.csv", dtype={"gene_id": str})
    merged = data[["grouped_region", "gene_id"]].copy()
    merged["t1"] = tier1.to_numpy(); merged["t2"] = tier2.to_numpy(); merged["t3"] = tier3.to_numpy()
    check = merged.merge(saved[["grouped_region", "gene_id", "Tier1", "Tier2", "Tier3"]], on=["grouped_region", "gene_id"], how="outer")
    return {
        "standardized_effect": {
            "region": row.grouped_region, "gene_id": str(row.gene_id),
            "beta": float(row.beta_diagnosis), "sd": float(row.sd_expression_gene),
            "reconstructed": std, "canonical": float(row.standardized_effect_size_diagnosis),
            "absolute_difference": abs(std - float(row.standardized_effect_size_diagnosis)),
            "tier2_p_check": bool(row.p_value_diagnosis < .01),
            "tier2_effect_check": bool(abs(std) > .8),
        },
        "tier_rows": len(check),
        "tier1_mismatches": int((check.t1 != check.Tier1).sum()),
        "tier2_mismatches": int((check.t2 != check.Tier2).sum()),
        "tier3_mismatches": int((check.t3 != check.Tier3).sum()),
        "counts": merged.groupby("grouped_region")[["t1", "t2", "t3"]].sum().astype(int).reset_index().to_dict("records"),
        "tier4_example": {
            "gene_id": str(row.gene_id),
            "formula": "sign(standardized_effect_size_diagnosis) * -log10(p_value_diagnosis)",
            "reconstructed": float(np.sign(row.standardized_effect_size_diagnosis) * -np.log10(row.p_value_diagnosis)),
            "stored_in_effect_table": float(row.signed_standardized_logp_score),
            "canonical_rank_file_present": False,
        },
    }


def multiple_testing_check() -> dict:
    data = pd.read_csv(ELISA / "LMM_all_genes_length_gt_100.csv")
    region = "Insula"
    d = data[data.brain_region.eq(region)].copy()
    recon = bh(d.p_value_diagnosis)
    return {"region": region, "n_tests": len(d), "max_absolute_difference": float(np.nanmax(abs(recon - d.fdr_BH.to_numpy(float)))), "minimum_canonical_BH": float(d.fdr_BH.min())}


def tier2_enrichment_check() -> dict:
    path = CLEAN / "08_Tier2_enrichment/tables/four_group__Hippocampal_Formation__down_in_MDD__GO_MF__full.csv"
    d = pd.read_csv(path)
    row = d.iloc[0]
    _, n = map(int, str(row.GeneRatio).split("/"))
    K, N = map(int, str(row.BgRatio).split("/"))
    k = int(row.Count)
    raw = float(hypergeom.sf(k - 1, N, K, n))
    adj = bh(d.pvalue)
    return {
        "source": str(path), "term": row.ID, "N": N, "K": K, "n": n, "k": k,
        "reconstructed_p": raw, "canonical_p": float(row.pvalue),
        "p_absolute_difference": abs(raw - float(row.pvalue)),
        "BH_over_saved_rows": float(adj[0]), "canonical_p_adjust": float(row["p.adjust"]),
        "saved_rows_in_family": len(d),
    }


def ml_checks() -> dict:
    sys.path.insert(0, str(ROOT / "code"))
    import extract_stable_ML_prioritized_Tier2_genes as scorer
    tier2 = scorer.load_tier2_tables([], [])
    v1raw = pd.read_csv(ELISA / "ML_pilot_Tier2/feature_importance/fold_level_feature_importance.csv", dtype={"gene_id": str})
    v1raw = v1raw[v1raw.feature_set.eq("Tier2_fixed")]
    v2raw = pd.read_csv(ELISA / "ML_strict_v2/feature_stability/selected_genes_per_fold_all_models.csv", dtype={"gene_id": str})
    v1 = scorer.summarize_events(scorer.annotate_importance_events(v1raw, "v1"), tier2, "v1")
    v2 = scorer.summarize_events(scorer.annotate_importance_events(v2raw, "v2"), tier2, "v2")
    rebuilt = scorer.build_priority_table(tier2, v1, v2)
    canonical = pd.read_csv(CLEAN / "10_ML_prioritized_Tier2_corrected/ML_prioritized_Tier2_genes_all_ranked.csv", dtype={"gene_id": str})
    merged = rebuilt.merge(canonical, on=["analysis_level", "region_or_subregion", "gene_id"], suffixes=("_rebuilt", "_canonical"))
    score_diff = abs(merged.ML_priority_score_rebuilt - merged.ML_priority_score_canonical)
    examples = []
    wanted = [("Cingulate gyrus 25", "LOC107985795"), ("Insula", "PRELID2")]
    high = canonical[canonical.ML_priority_class.str.startswith("High")].sort_values("ML_priority_score", ascending=False).iloc[0]
    high_boundary = canonical[canonical.ML_priority_class.str.startswith("High")].sort_values("ML_priority_score").iloc[0]
    moderate_boundary = canonical[canonical.ML_priority_class.str.startswith("Moderate")].assign(delta=lambda x: abs(x.ML_priority_score - .45)).sort_values("delta").iloc[0]
    wanted += [(high.region_or_subregion, high.gene_symbol), (high_boundary.region_or_subregion, high_boundary.gene_symbol), (moderate_boundary.region_or_subregion, moderate_boundary.gene_symbol)]
    for region, symbol in wanted:
        row = canonical[(canonical.region_or_subregion.eq(region)) & (canonical.gene_symbol.eq(symbol))].iloc[0]
        v1_denom = 125 if pd.isna(row.n_events_v1) else int(row.n_events_v1)
        v2_denom = 150 if pd.isna(row.n_events_v2) else int(row.n_events_v2)
        arithmetic = (
            .40 * row.selection_frequency_combined + .15 * row.model_support_scaled
            + .10 * row.sparse_model_support_binary + .15 * row.importance_rank_score
            + .10 * row.sign_consistency_score + .05 * row.lmm_evidence_scaled
            + .05 * row.v2_support_score
        )
        examples.append({
            "region": region, "gene_id": str(row.gene_id), "gene_symbol": row.gene_symbol,
            "v1_selected": int(round(row.selection_frequency_v1 * v1_denom)), "v1_denominator": v1_denom,
            "v2_selected": int(round(row.selection_frequency_v2 * v2_denom)), "v2_denominator": v2_denom,
            "selection_frequency_v1": float(row.selection_frequency_v1), "selection_frequency_v2": float(row.selection_frequency_v2),
            "combined_selection_frequency": float(row.selection_frequency_combined), "model_support_scaled": float(row.model_support_scaled),
            "sparse_model_support_binary": float(row.sparse_model_support_binary), "importance_rank_score": float(row.importance_rank_score),
            "sign_consistency_score": float(row.sign_consistency_score), "observed_sign_consistency": None if pd.isna(row.sign_consistency) else float(row.sign_consistency),
            "lmm_evidence_scaled": float(row.lmm_evidence_scaled), "v2_support_score": float(row.v2_support_score),
            "reconstructed_score": float(arithmetic), "canonical_score": float(row.ML_priority_score),
            "absolute_difference": abs(float(arithmetic) - float(row.ML_priority_score)), "class": row.ML_priority_class,
        })
    return {
        "rebuilt_rows": len(rebuilt), "matched_rows": len(merged), "max_score_absolute_difference": float(score_diff.max()),
        "counts": canonical.ML_priority_class.value_counts().to_dict(), "examples": examples,
        "v1_denominators": sorted(pd.to_numeric(canonical.n_events_v1).dropna().unique().astype(int).tolist()),
        "v2_denominators": sorted(pd.to_numeric(canonical.n_events_v2).dropna().unique().astype(int).tolist()),
    }


def downstream_checks() -> dict:
    ahba = pd.read_csv(CLEAN / "16_corrected_ML_downstream/AHBA_corrected/published_AHBA_ML_enrichment_corrected.csv")
    m26 = ahba[(ahba.candidate_set.eq("high_plus_moderate")) & ahba.scope.eq("region:Ventral Striatum|direction:up_in_MDD") & ahba.published_module.eq("M26")].iloc[0]
    ahba_p = float(hypergeom.sf(int(m26.overlap_genes) - 1, int(m26.background_genes), int(m26.module_genes_in_background), int(m26.candidate_genes_in_background)))
    wg = pd.read_csv(CLEAN / "16_corrected_ML_downstream/WGCNA_ML_annotations_corrected/broad_regions_ML_module_enrichment_corrected.csv")
    wr = wg[(wg.region.eq("Ventral Striatum")) & wg.module.eq("MEdarkred") & wg.gene_set.eq("high_plus_moderate_ML_genes")].iloc[0]
    assignment = pd.read_csv(CLEAN / "12_MDD_WGCNA_modules/best_solution_gene_to_module_assignments.csv")
    region_N = assignment[assignment.region.eq(wr.region)].gene_id.astype(str).nunique()
    a, rowN, colN = int(wr.overlap_count), int(wr.module_gene_count), int(wr.gene_set_count_in_background)
    table = [[a, rowN-a], [colN-a, region_N-rowN-colN+a]]
    fisher_p = float(fisher_exact(table, alternative="greater").pvalue)
    fw = pd.read_csv(CLEAN / "16_corrected_ML_downstream/WGCNA_ML_annotations_corrected/frontal_subregions_ML_module_enrichment_corrected.csv")
    oft = fw[(fw.region.eq("OFC BA11")) & fw.module.eq("MEtan") & fw.gene_set.eq("high_confidence_ML_genes")].iloc[0]
    core_checks = []
    pairs = [
        (CLEAN / "12_MDD_WGCNA_modules/module_level_ml2/stable_core_genes_by_module.csv", CLEAN / "16_corrected_ML_downstream/stable_core_annotations_corrected/broad_structural_core_with_corrected_ML_annotations.csv"),
        (CLEAN / "12_MDD_WGCNA_modules/frontal_subregion_wgcna/module_level_ml2/frontal_subregion_stable_core_genes_by_module.csv", CLEAN / "16_corrected_ML_downstream/stable_core_annotations_corrected/frontal_structural_core_with_corrected_ML_annotations.csv"),
    ]
    for original, corrected in pairs:
        o, c = pd.read_csv(original, dtype={"gene_id": str}), pd.read_csv(corrected, dtype={"gene_id": str})
        shared = [x for x in o.columns if x in c.columns and not x.startswith("is_")]
        core_checks.append({"original": str(original), "corrected": str(corrected), "rows_original": len(o), "rows_corrected": len(c), "shared_columns_exact": bool(o[shared].equals(c[shared]))})
    network = pd.read_csv(CLEAN / "16_corrected_ML_downstream/WGCNA_ML_annotations_corrected/network_figures_frontal_subregions/frontal_subregion_network_inventory.csv")
    return {
        "AHBA_M26": {"N": int(m26.background_genes), "K": int(m26.module_genes_in_background), "n": int(m26.candidate_genes_in_background), "k": int(m26.overlap_genes), "genes": m26.overlap_gene_symbols, "reconstructed_p": ahba_p, "canonical_p": float(m26.p_value_hypergeometric), "q": float(m26.fdr_BH), "difference": abs(ahba_p-float(m26.p_value_hypergeometric))},
        "AHBA_results_q_lt_0_10": int((ahba.fdr_BH < .10).sum()),
        "WGCNA_VS_MEdarkred": {"N": region_N, "K_module": rowN, "n_set": colN, "k": a, "genes": wr.overlap_genes, "reconstructed_p": fisher_p, "canonical_p": float(wr.p_value), "q": float(wr.fdr_BH), "difference": abs(fisher_p-float(wr.p_value))},
        "OFC_MEtan_high": {"overlap": int(oft.overlap_count), "p": float(oft.p_value), "q": float(oft.fdr_BH)},
        "stable_core": core_checks,
        "frontal_network_regions": sorted(network.region.unique().tolist()),
    }


def main() -> None:
    meta = baseline_csv("outputs/length_gt_100_clean_pipeline/01_inputs/sample_metadata_used.csv", dtype={"SampleID": str, "BrainID": str})
    evidence = {
        "baseline": BASELINE,
        "canonical_ranking_sha256": sha(CLEAN / "10_ML_prioritized_Tier2_corrected/ML_prioritized_Tier2_genes_all_ranked.csv"),
        "sample_accounting": sample_accounting(meta),
        "normalization": normalization_check(),
        "gene_models": lm_checks(meta),
        "multiple_testing": multiple_testing_check(),
        "tiers_and_effects": tier_and_effect_checks(),
        "tier2_enrichment": tier2_enrichment_check(),
        "ml": ml_checks(),
        "downstream": downstream_checks(),
    }
    OUT.write_text(json.dumps(evidence, indent=2, default=str), encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
