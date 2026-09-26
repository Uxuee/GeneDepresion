# Canonical Results Index

| Stage | Current authoritative output | Generating script | Status | Notes |
|---|---|---|---|---|
| Normalised expression | `outputs/length_gt_100_clean_pipeline/02_normalized_expression/` | `02_run_getmm_log2_normalization.py` and original Step 1 | Complete | log2(GeTMM + 1); clean stage reuses matrices |
| QC/covariates | `01_inputs/` plus saved QC outputs | Step 1/2 and clean input scripts | Needs documentation | A tracked metadata file is currently deleted in the working tree |
| Gene-wise Diagnosis models | `03_broad_region_LMM/` | `03_run_broad_region_LMM.py` and original Step 3 | Complete | Statistical association layer |
| Standardised effects | `04_standardized_effect_size/` | `04_compute_standardized_effect_sizes.py` | Complete | Direction and magnitude |
| Tier definitions | `05_evidence_tiers/`, `07_Tier2_candidate_gene_sets/` | Steps 05 and 07 | Complete | Tier 2 is exploratory |
| Enrichment | `08_Tier2_enrichment/` | Step 08 and original enrichment workflow | Needs repository cleanup | Tracked tables/plots currently deleted |
| ML v1 events | `outputs/length_gt_100_elisa_style/ML_pilot_Tier2/` | archived v1 training script | Complete provenance | Final score uses `Tier2_fixed` only |
| Strict ML v2 events | `outputs/length_gt_100_elisa_style/ML_strict_v2/` | archived v2 training script | Complete provenance | Baseline excluded from gene evidence |
| Corrected ML priority | `10_ML_prioritized_Tier2_corrected/` | corrected scorer/wrapper | **Canonical/current** | 671 = 30/119/522 |
| Legacy ML priority | `10_ML_prioritized_Tier2_extended/` | historical scorer | Legacy only | 35/130/506; do not use downstream |
| AHBA mapping | `11_AHBA_modules_enrichment/` | `code/ahba_modules/*` | Needs rerun | Existing scripts use legacy ML lists; context only |
| MDD WGCNA | `12_MDD_WGCNA_modules/` | Steps 6–8 | Modules complete; annotations need rerun | Preserve modules; refresh ML-dependent outputs |
| Module-level ML | `13_Module_Level_ML2/` and Step 12 subfolders | module-level scripts | Needs dependency review | Exploratory |
| External validation | `15_External_Validation/` | validation workflow | Future validation | Not complete |
| Manuscript package | `14_Preliminary_Manuscript_Package/` | packaging scripts | Needs regeneration | Replace legacy ML-dependent artifacts |
