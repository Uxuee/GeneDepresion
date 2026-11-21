# Step 3 provenance reconciliation

The prior registry stage 04 associated the frozen diagnosis tables with the current limma-voom/dream scripts. That association was too broad. The preserved registry is `provenance/consolidation/canonical_pipeline.before.json`.

The audited six-region Tier2 screen uses fixed-effect OLS on log2(GeTMM + 1): `expression ~ Diagnosis + Age + Sex + RIN`, CTRL reference. Standardized effect is diagnosis beta divided by the gene-wise sample SD, not residual SD. Broad-region evidence resides in `outputs/length_gt_100_clean_pipeline/04_standardized_effect_size/LMM_all_genes_with_standardized_effect_size_4groups.csv`. The frontal-subregion source is Git blob `f700cb323a47b24825b11bc1bba33cc07e6430bb`, documented from historical commit `ac37ec9`. Historical filenames containing LMM do not establish that every regional fit was mixed-effects.

Evidence: `docs/UPSTREAM_STATISTICAL_BENCHMARK_DESIGN.md`, `docs/UPSTREAM_STATISTICAL_METHOD_BENCHMARK_REPORT.md`, `code/statistical_benchmark/benchmark_common.py`, and the frozen benchmark manifest. The clean `03_run_broad_region_LMM.py` and `06_run_frontal_subregion_sensitivity.py` explicitly reuse completed outputs rather than fitting models. Exact end-to-end generator recovery for the frozen LM remains incomplete; the registry now exposes this gap with no asserted direct generator.

Current Step3 scripts represent count-based limma-voom/dream branches. Pooled frontal dream uses Diagnosis + Age + Sex + RIN + OriginalRegion + (1 | BrainID). The later six-region benchmark uses matched raw counts, TMM/voom/lmFit/eBayes and DESeq2, retaining the frozen LM gene/sample universe. It is separate from the pooled dream analysis. A later same-sized but non-identical count export was rejected in that benchmark.

Consolidation changes the provenance description, not the fits, scripts, Tier2 membership, or effect estimates. No Step 1–3 script was edited. This reconciliation does not claim that rerunning a current Step3 script reproduces the frozen LM.
