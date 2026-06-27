# Relevant excerpts from D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py
# Extracted for ML-priority score and confidence-class audit.

# Statistical evidence score, lines 108-121
    tier2["beta_diagnosis"] = pd.to_numeric(tier2["beta_diagnosis"], errors="coerce")
    tier2["p_value_diagnosis"] = pd.to_numeric(tier2["p_value_diagnosis"], errors="coerce")
    tier2["standardized_effect_size_diagnosis"] = pd.to_numeric(
        tier2["standardized_effect_size_diagnosis"], errors="coerce"
    )
    tier2["lmm_direction"] = np.select(
        [tier2["beta_diagnosis"] > 0, tier2["beta_diagnosis"] < 0],
        ["up_in_MDD", "down_in_MDD"],
        default="no_change",
    )
    tier2["lmm_evidence_score"] = (
        -np.log10(tier2["p_value_diagnosis"].clip(lower=1e-300))
        * tier2["standardized_effect_size_diagnosis"].abs()
    )

# Selection frequency/model support/importance/sign-support inputs, lines 293-322
    result["selection_frequency_v1"] = result["selection_frequency_v1"].fillna(0.0)
    if v2_available:
        result["selection_frequency_v2"] = result["selection_frequency_v2"].fillna(0.0)
        result["selection_frequency_combined"] = 0.6 * result["selection_frequency_v2"] + 0.4 * result["selection_frequency_v1"]
        result["v2_support_score"] = result["v2_support_flag"].fillna(False).astype(bool).astype(float)
    else:
        result["selection_frequency_combined"] = result["selection_frequency_v1"]
        result["v2_support_score"] = 0.5

    result["model_support_count"] = (
        result["model_support_count_v1"].fillna(0).astype(int)
        + result["model_support_count_v2"].fillna(0).astype(int)
    )
    # Count a family once even if v1 and v2 both support it.
    for idx, row in result.iterrows():
        models = set()
        for col in ["support_models_v1", "support_models_v2"]:
            if col in result.columns and pd.notna(row.get(col)):
                models.update(str(row[col]).split(";"))
        models.discard("")
        result.at[idx, "model_support_count"] = len(models & set(MODEL_FAMILIES))

    result["sparse_model_support"] = result["sparse_model_support_v1"].fillna(False) | result["sparse_model_support_v2"].fillna(False)
    result["linear_model_support"] = result["linear_model_support_v1"].fillna(False) | result["linear_model_support_v2"].fillna(False)
    result["rf_support"] = result["rf_support_v1"].fillna(False) | result["rf_support_v2"].fillna(False)
    result["mean_abs_importance"] = result[["mean_abs_importance_v1", "mean_abs_importance_v2"]].mean(axis=1, skipna=True)
    result["median_abs_importance"] = result[["median_abs_importance_v1", "median_abs_importance_v2"]].median(axis=1, skipna=True)
    result["mean_rank"] = result[["mean_rank_v1", "mean_rank_v2"]].mean(axis=1, skipna=True)
    result["rank_stability"] = result[["rank_stability_v1", "rank_stability_v2"]].mean(axis=1, skipna=True).fillna(0.0)
    result["sign_consistency"] = result[["sign_consistency_v1", "sign_consistency_v2"]].mean(axis=1, skipna=True)

# ML direction and statistical/ML sign agreement, lines 333-338
    result["ml_direction"] = result.apply(choose_direction, axis=1)
    result["lmm_ml_sign_agreement"] = np.where(
        result["ml_direction"].isin(["up_in_MDD", "down_in_MDD"]),
        result["ml_direction"].eq(result["lmm_direction"]),
        np.nan,
    )

# Composite ML_priority_score calculation, lines 340-355
    max_models_by_region = result.groupby("region_or_subregion")["model_support_count"].transform("max").replace(0, np.nan)
    result["model_support_scaled"] = (result["model_support_count"] / max_models_by_region).fillna(0.0)
    max_rank_by_region = result.groupby("region_or_subregion")["mean_rank"].transform("max")
    result["importance_rank_score"] = (1 - ((result["mean_rank"] - 1) / max_rank_by_region.replace(0, np.nan))).clip(0, 1).fillna(0.0)
    result["sign_consistency_score"] = result["sign_consistency"].fillna(0.5)
    result["lmm_evidence_scaled"] = minmax_by_region(result["lmm_evidence_score"], result["region_or_subregion"])
    result["sparse_model_support_binary"] = result["sparse_model_support"].astype(float)
    result["ML_priority_score"] = (
        0.30 * result["selection_frequency_combined"]
        + 0.20 * result["model_support_scaled"]
        + 0.15 * result["sparse_model_support_binary"]
        + 0.15 * result["importance_rank_score"]
        + 0.10 * result["sign_consistency_score"]
        + 0.05 * result["lmm_evidence_scaled"]
        + 0.05 * result["v2_support_score"]
    ).clip(0, 1)

# High/moderate confidence class assignment, lines 357-377
    agreement_allowed = result["lmm_ml_sign_agreement"].map(lambda x: True if pd.isna(x) else bool(x))
    agreement_not_false = result["lmm_ml_sign_agreement"].map(lambda x: True if pd.isna(x) else bool(x))
    high = (
        (result["ML_priority_score"] >= 0.70)
        & (result["selection_frequency_combined"] >= 0.50)
        & (result["model_support_count"] >= 2)
        & result["sparse_model_support"]
        & (result["sign_consistency"].fillna(1.0) >= 0.70)
        & agreement_allowed
    )
    moderate = (
        (result["ML_priority_score"] >= 0.45)
        & (result["selection_frequency_combined"] >= 0.25)
        & (result["model_support_count"] >= 1)
        & agreement_not_false
    )
    result["ML_priority_class"] = np.select(
        [high, moderate],
        ["High-confidence ML-prioritized Tier 2", "Moderate-confidence ML-prioritized Tier 2"],
        default="Exploratory ML-prioritized Tier 2",
    )

# Final class-specific exports, lines 594-600
    result[output_cols].to_csv(OUT_DIR / "ML_prioritized_Tier2_genes_all_ranked.csv", index=False)
    result[result["ML_priority_class"].eq("High-confidence ML-prioritized Tier 2")][output_cols].to_csv(
        OUT_DIR / "ML_prioritized_Tier2_genes_high_confidence.csv", index=False
    )
    result[result["ML_priority_class"].eq("Moderate-confidence ML-prioritized Tier 2")][output_cols].to_csv(
        OUT_DIR / "ML_prioritized_Tier2_genes_moderate_confidence.csv", index=False
    )
