# ML Priority Score Code Audit

## Files Inspected

- Primary scoring script: `D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py`
- Clean pipeline wrapper: `D:\genes\GeneDepresion\code\length_gt_100_clean_pipeline\10_extract_extended_ML_prioritized_Tier2_genes.py`
- Current final ranked table: `D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline\10_ML_prioritized_Tier2_extended\final_ranked_ML_prioritized_Tier2_genes.csv`

This audit did not rerun ML, LMM, or enrichment.

## Exact Code: Statistical Evidence Input

```python
108:     tier2["beta_diagnosis"] = pd.to_numeric(tier2["beta_diagnosis"], errors="coerce")
109:     tier2["p_value_diagnosis"] = pd.to_numeric(tier2["p_value_diagnosis"], errors="coerce")
110:     tier2["standardized_effect_size_diagnosis"] = pd.to_numeric(
111:         tier2["standardized_effect_size_diagnosis"], errors="coerce"
112:     )
113:     tier2["lmm_direction"] = np.select(
114:         [tier2["beta_diagnosis"] > 0, tier2["beta_diagnosis"] < 0],
115:         ["up_in_MDD", "down_in_MDD"],
116:         default="no_change",
117:     )
118:     tier2["lmm_evidence_score"] = (
119:         -np.log10(tier2["p_value_diagnosis"].clip(lower=1e-300))
120:         * tier2["standardized_effect_size_diagnosis"].abs()
121:     )
```

Plain language: `lmm_evidence_score` is `-log10(p_value_diagnosis) * abs(standardized_effect_size_diagnosis)`. Despite the historical column name, this field is used for both broad-region LMM rows and frontal subregion fixed-effect rows in the downstream prioritization table.

## Exact Code: Support Columns Used Before Scoring

```python
293:     result["selection_frequency_v1"] = result["selection_frequency_v1"].fillna(0.0)
294:     if v2_available:
295:         result["selection_frequency_v2"] = result["selection_frequency_v2"].fillna(0.0)
296:         result["selection_frequency_combined"] = 0.6 * result["selection_frequency_v2"] + 0.4 * result["selection_frequency_v1"]
297:         result["v2_support_score"] = result["v2_support_flag"].fillna(False).astype(bool).astype(float)
298:     else:
299:         result["selection_frequency_combined"] = result["selection_frequency_v1"]
300:         result["v2_support_score"] = 0.5
301: 
302:     result["model_support_count"] = (
303:         result["model_support_count_v1"].fillna(0).astype(int)
304:         + result["model_support_count_v2"].fillna(0).astype(int)
305:     )
306:     # Count a family once even if v1 and v2 both support it.
307:     for idx, row in result.iterrows():
308:         models = set()
309:         for col in ["support_models_v1", "support_models_v2"]:
310:             if col in result.columns and pd.notna(row.get(col)):
311:                 models.update(str(row[col]).split(";"))
312:         models.discard("")
313:         result.at[idx, "model_support_count"] = len(models & set(MODEL_FAMILIES))
314: 
315:     result["sparse_model_support"] = result["sparse_model_support_v1"].fillna(False) | result["sparse_model_support_v2"].fillna(False)
316:     result["linear_model_support"] = result["linear_model_support_v1"].fillna(False) | result["linear_model_support_v2"].fillna(False)
317:     result["rf_support"] = result["rf_support_v1"].fillna(False) | result["rf_support_v2"].fillna(False)
318:     result["mean_abs_importance"] = result[["mean_abs_importance_v1", "mean_abs_importance_v2"]].mean(axis=1, skipna=True)
319:     result["median_abs_importance"] = result[["median_abs_importance_v1", "median_abs_importance_v2"]].median(axis=1, skipna=True)
320:     result["mean_rank"] = result[["mean_rank_v1", "mean_rank_v2"]].mean(axis=1, skipna=True)
321:     result["rank_stability"] = result[["rank_stability_v1", "rank_stability_v2"]].mean(axis=1, skipna=True).fillna(0.0)
322:     result["sign_consistency"] = result[["sign_consistency_v1", "sign_consistency_v2"]].mean(axis=1, skipna=True)
```

Columns prepared here include `selection_frequency_v1`, `selection_frequency_v2`, `selection_frequency_combined`, `v2_support_score`, `model_support_count`, `sparse_model_support`, `linear_model_support`, `rf_support`, `mean_abs_importance`, `median_abs_importance`, `mean_rank`, `rank_stability`, and `sign_consistency`.

## Exact Code: Statistical/ML Direction Agreement

```python
333:     result["ml_direction"] = result.apply(choose_direction, axis=1)
334:     result["lmm_ml_sign_agreement"] = np.where(
335:         result["ml_direction"].isin(["up_in_MDD", "down_in_MDD"]),
336:         result["ml_direction"].eq(result["lmm_direction"]),
337:         np.nan,
338:     )
```

`lmm_ml_sign_agreement` is true when the ML direction is up/down in MDD and matches the statistical direction column used in the prioritization script. Missing/unavailable agreement is later allowed; explicit false agreement blocks high/moderate classification.

## Exact Code: ML_priority_score Formula

```python
340:     max_models_by_region = result.groupby("region_or_subregion")["model_support_count"].transform("max").replace(0, np.nan)
341:     result["model_support_scaled"] = (result["model_support_count"] / max_models_by_region).fillna(0.0)
342:     max_rank_by_region = result.groupby("region_or_subregion")["mean_rank"].transform("max")
343:     result["importance_rank_score"] = (1 - ((result["mean_rank"] - 1) / max_rank_by_region.replace(0, np.nan))).clip(0, 1).fillna(0.0)
344:     result["sign_consistency_score"] = result["sign_consistency"].fillna(0.5)
345:     result["lmm_evidence_scaled"] = minmax_by_region(result["lmm_evidence_score"], result["region_or_subregion"])
346:     result["sparse_model_support_binary"] = result["sparse_model_support"].astype(float)
347:     result["ML_priority_score"] = (
348:         0.30 * result["selection_frequency_combined"]
349:         + 0.20 * result["model_support_scaled"]
350:         + 0.15 * result["sparse_model_support_binary"]
351:         + 0.15 * result["importance_rank_score"]
352:         + 0.10 * result["sign_consistency_score"]
353:         + 0.05 * result["lmm_evidence_scaled"]
354:         + 0.05 * result["v2_support_score"]
355:     ).clip(0, 1)
```

Plain-language formula:

`ML_priority_score = 0.30 * selection_frequency_combined + 0.20 * model_support_scaled + 0.15 * sparse_model_support_binary + 0.15 * importance_rank_score + 0.10 * sign_consistency_score + 0.05 * lmm_evidence_scaled + 0.05 * v2_support_score`, clipped to `[0, 1]`.

## Exact Code: High-Confidence And Moderate-Confidence Thresholds

```python
357:     agreement_allowed = result["lmm_ml_sign_agreement"].map(lambda x: True if pd.isna(x) else bool(x))
358:     agreement_not_false = result["lmm_ml_sign_agreement"].map(lambda x: True if pd.isna(x) else bool(x))
359:     high = (
360:         (result["ML_priority_score"] >= 0.70)
361:         & (result["selection_frequency_combined"] >= 0.50)
362:         & (result["model_support_count"] >= 2)
363:         & result["sparse_model_support"]
364:         & (result["sign_consistency"].fillna(1.0) >= 0.70)
365:         & agreement_allowed
366:     )
367:     moderate = (
368:         (result["ML_priority_score"] >= 0.45)
369:         & (result["selection_frequency_combined"] >= 0.25)
370:         & (result["model_support_count"] >= 1)
371:         & agreement_not_false
372:     )
373:     result["ML_priority_class"] = np.select(
374:         [high, moderate],
375:         ["High-confidence ML-prioritized Tier 2", "Moderate-confidence ML-prioritized Tier 2"],
376:         default="Exploratory ML-prioritized Tier 2",
377:     )
```

## Threshold Table

| priority_class      | criterion_name                      | operator                    | threshold_value                                  | source_column                                 | code_source_file                                                         | code_line_or_block          |
|:--------------------|:------------------------------------|:----------------------------|:-------------------------------------------------|:----------------------------------------------|:-------------------------------------------------------------------------|:----------------------------|
| score_formula       | selection_frequency_combined_weight | weight                      | 0.30                                             | selection_frequency_combined                  | D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py | lines 347-355               |
| score_formula       | model_support_scaled_weight         | weight                      | 0.20                                             | model_support_count scaled within region      | D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py | lines 340-341, 347-355      |
| score_formula       | sparse_model_support_weight         | weight                      | 0.15                                             | sparse_model_support_binary                   | D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py | lines 346-355               |
| score_formula       | importance_rank_score_weight        | weight                      | 0.15                                             | mean_rank converted to importance_rank_score  | D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py | lines 342-343, 347-355      |
| score_formula       | sign_consistency_score_weight       | weight                      | 0.10                                             | sign_consistency, fillna(0.5)                 | D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py | lines 344, 347-355          |
| score_formula       | statistical_evidence_weight         | weight                      | 0.05                                             | lmm_evidence_score scaled within region       | D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py | lines 118-121, 345, 347-355 |
| score_formula       | v2_support_score_weight             | weight                      | 0.05                                             | v2_support_score                              | D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py | lines 294-300, 347-355      |
| high_confidence     | ML_priority_score                   | >=                          | 0.70                                             | ML_priority_score                             | D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py | line 360                    |
| high_confidence     | selection_frequency_combined        | >=                          | 0.50                                             | selection_frequency_combined                  | D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py | line 361                    |
| high_confidence     | model_support_count                 | >=                          | 2                                                | model_support_count                           | D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py | line 362                    |
| high_confidence     | sparse_model_support                | is True                     | True                                             | sparse_model_support                          | D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py | line 363                    |
| high_confidence     | sign_consistency                    | >=                          | 0.70 after fillna(1.0)                           | sign_consistency                              | D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py | line 364                    |
| high_confidence     | statistical/ML sign agreement       | not False                   | True if missing else bool(lmm_ml_sign_agreement) | lmm_ml_sign_agreement via agreement_allowed   | D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py | lines 357, 365              |
| moderate_confidence | ML_priority_score                   | >=                          | 0.45                                             | ML_priority_score                             | D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py | line 368                    |
| moderate_confidence | selection_frequency_combined        | >=                          | 0.25                                             | selection_frequency_combined                  | D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py | line 369                    |
| moderate_confidence | model_support_count                 | >=                          | 1                                                | model_support_count                           | D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py | line 370                    |
| moderate_confidence | statistical/ML sign agreement       | not False                   | True if missing else bool(lmm_ml_sign_agreement) | lmm_ml_sign_agreement via agreement_not_false | D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py | lines 358, 371              |
| class_assignment    | mutual exclusivity                  | np.select first-match order | high evaluated before moderate                   | high, moderate masks                          | D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py | lines 373-377               |

## Plain-Language Threshold Explanation

High-confidence requires all of the following:

- `ML_priority_score >= 0.70`
- `selection_frequency_combined >= 0.50`
- `model_support_count >= 2`
- `sparse_model_support == True`
- `sign_consistency >= 0.70`, with missing values treated as 1.0 for this threshold
- statistical/ML sign agreement is not explicitly false

Moderate-confidence requires all of the following:

- `ML_priority_score >= 0.45`
- `selection_frequency_combined >= 0.25`
- `model_support_count >= 1`
- statistical/ML sign agreement is not explicitly false

## Are High And Moderate Mutually Exclusive?

Yes. The Boolean masks can overlap mathematically because the moderate thresholds are lower, but the final label assignment uses `np.select([high, moderate], ...)` with `high` first. Therefore, any row meeting high-confidence criteria is labeled high-confidence, not moderate-confidence. The exported moderate-confidence file is created by filtering rows whose final `ML_priority_class` equals `Moderate-confidence ML-prioritized Tier 2`, so it excludes high-confidence rows.

## Exact Code: Final Class-Specific Exports

```python
594:     result[output_cols].to_csv(OUT_DIR / "ML_prioritized_Tier2_genes_all_ranked.csv", index=False)
595:     result[result["ML_priority_class"].eq("High-confidence ML-prioritized Tier 2")][output_cols].to_csv(
596:         OUT_DIR / "ML_prioritized_Tier2_genes_high_confidence.csv", index=False
597:     )
598:     result[result["ML_priority_class"].eq("Moderate-confidence ML-prioritized Tier 2")][output_cols].to_csv(
599:         OUT_DIR / "ML_prioritized_Tier2_genes_moderate_confidence.csv", index=False
600:     )
```

## Final Counts Verified From Current Files

- Total ranked Tier 2 gene-region pairs: 671
- High-confidence count: 35
- Moderate-confidence count, excluding high-confidence: 130
- High + moderate total: 165
- Exploratory/lower-priority count: 506

These match the expected current values: 671 total, 35 high-confidence, 130 additional moderate-confidence, and 165 high/moderate total.
