### ML-priority score and confidence-class thresholds

Genes were assigned to high- or moderate-confidence classes using predefined thresholds on the composite ML-priority score and supporting stability metrics. The exact thresholds were extracted from the prioritization script `D:\genes\GeneDepresion\code\extract_stable_ML_prioritized_Tier2_genes.py`.

The composite score was:

`ML_priority_score = 0.30 * selection_frequency_combined + 0.20 * model_support_scaled + 0.15 * sparse_model_support_binary + 0.15 * importance_rank_score + 0.10 * sign_consistency_score + 0.05 * lmm_evidence_scaled + 0.05 * v2_support_score`, clipped to 0-1.

High-confidence required `ML_priority_score >= 0.70`, `selection_frequency_combined >= 0.50`, `model_support_count >= 2`, sparse-model support, `sign_consistency >= 0.70`, and no explicit statistical/ML sign disagreement. Moderate-confidence required `ML_priority_score >= 0.45`, `selection_frequency_combined >= 0.25`, `model_support_count >= 1`, and no explicit statistical/ML sign disagreement. Labels are mutually exclusive because high-confidence is assigned first with `np.select([high, moderate], ...)`; moderate-confidence therefore excludes high-confidence rows.

Current final counts are 671 ranked Tier 2 gene-region pairs, 35 high-confidence, 130 additional moderate-confidence, and 165 high/moderate candidates total.
