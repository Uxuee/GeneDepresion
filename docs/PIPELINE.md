# Current MDD workflow

Tier2-RRA is the current primary exploratory Tier 2 ranking. This is a presentation and workflow designation, not evidence of statistical or biological superiority. The corrected weighted v1 + strict-v2 ranking remains a historical canonical benchmark. Its 30 high, 119 moderate and 522 low labels do not apply to RRA.

The consolidated root is `D:/genes/GeneDepresion`. ASD repository `Uxuee/DGEAnalysis` is separate. Paths below are relative to the MDD root. Descriptive identifiers avoid ambiguous historical output numbers; no output directory is renamed.

| Descriptive stage | Existing source / product | Execution and status |
|---|---|---|
| source_data | `D:/genes/genedata/GSE102556_raw_counts_GRCh38.p13_NCBI.tsv`; annotation and matched metadata | Local external data; inspect source manifests before rerunning |
| preprocessing_getmm | `code/Step1_length_cutoff_sensitivity_GeTMM.R`; `outputs/step1_length_sensitivity/Length_gt_100/` | Script available; requires raw inputs, metadata and R packages. Frozen results preserved |
| covariate_exploration | `code/Step2_Covariate_Exploration_FrontalCortex.R` | Reads log2 GeTMM; hard-coded local input root and first-match metadata selection need explicit checking |
| frozen_regional_lm | `outputs/length_gt_100_clean_pipeline/04_standardized_effect_size/`; frontal Git blob recorded in registry | Fixed-effect LM evidence; reuse wrappers are not model generators. See Step 3 reconciliation |
| statistical_method_benchmark | `code/statistical_benchmark/`; `outputs/length_gt_100_clean_pipeline/21_Statistical_method_benchmark/` | Frozen matched limma-voom/DESeq2 benchmarks; no regeneration in consolidation |
| pooled_dream_sensitivity | `code/Step3_frontal_dream_simple.R`, `code/Step3_all_regions_combined_simple.R` | Distinct pooled/repeated-donor branch; not generator of the six-region Tier2 LM |
| tier2_candidates | `outputs/length_gt_100_clean_pipeline/07_Tier2_candidate_gene_sets/Tier2_candidate_genes_all_regions.csv` | Frozen 671 pairs; p < 0.01 AND abs(beta/SD(expression)) > 0.8 |
| corrected_weighted_benchmark | `outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_corrected/` | Unchanged v1 + strict-v2 weighted result and method-specific labels |
| tier2_rra_primary | `code/tier2_ml_rra/`; `outputs/length_gt_100_clean_pipeline/18_Tier2_ML_RRA_prioritization/` | Current primary exploratory ranking; frozen numerical outputs |
| allgene_rra_sensitivity | `code/ml_rra/`; `outputs/length_gt_100_clean_pipeline/17_ML_RRA_consensus/` | Separate full-transcriptome sensitivity |
| tier2_method_benchmark | `code/tier2_benchmark/`; `outputs/length_gt_100_clean_pipeline/19_Tier2_ML_benchmark/` | LASSO / Stability Selection, same-cohort dependent evidence |
| tier2_conditional_null | `code/tier2_conditional_null/`; `outputs/length_gt_100_clean_pipeline/20_Tier2_ML_conditional_null/` | Frozen 20-permutation pilot |
| wgcna_context | `outputs/length_gt_100_clean_pipeline/12_MDD_WGCNA_modules/`, `17_WGCNA_reconstruction_validation/` | Frozen modules; reconstruction is a validation candidate, not automatic replacement |
| ahba_context | `outputs/length_gt_100_clean_pipeline/11_AHBA_modules_enrichment/`, `16_corrected_ML_downstream/` | Healthy-brain context; corrected weighted annotations are not RRA annotations |
| final_evidence_synthesis | `outputs/length_gt_100_clean_pipeline/22_Final_evidence_synthesis/` | Frozen 671-row evidence matrix; RRA order retained |
| literature_intersections | `external_validation/reference_studies/`, `external_validation/rra_validation/tier2_branch_primary/` | Retrieval, mapping, comparisons and direct full-Tier2/Labonte report completed locally; independent-cohort validation remains unfulfilled |

## What can be executed now

Inspection: `python code/run_canonical_pipeline.py --status` and `--dry-run` describe registered evidence. No stage has an automatic execution command: this is not a turnkey raw-data-to-paper runner. `python audit/verify_consolidation.py` checks preserved hashes and path availability without fitting models. Historical scripts may be executable after environment/input setup, but execution has not been recertified here. Do not run all validation scripts indiscriminately: some regenerate derived outputs.

See [Methods](METHODS.md), [Results](RESULTS.md), [Limitations](LIMITATIONS.md), [Reproducibility](REPRODUCIBILITY.md) and [Step 3 provenance reconciliation](STEP3_PROVENANCE_RECONCILIATION.md). The machine-readable registry and `provenance/canonical_outputs_manifest.json` distinguish designation, provenance, local availability and publication scope.
