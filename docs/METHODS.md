# Methods and ranking roles

GSE102556 counts and NCBI GRCh38.p13 annotation are local external inputs. Metadata sample identifiers connect count columns to diagnosis, region, age, sex and RIN. Do not infer BrainBank or SeqBatch when absent. Preserve original identifiers, downloaded supplements and the existing pinned gene mappings.

Frozen preprocessing sums duplicated GeneIDs, matches gene-length annotation, retains length >100 bp, applies the regional expression eligibility filter, calculates reads per kilobase and TMM factors, then GeTMM and log2(GeTMM + 1). The original frontal group is split into CG25, dlPFC BA8/9 and OFC BA11 for regional LM interpretation. Step 2 explores covariates; it is not the Tier2 selector.

The frozen regional LM uses Diagnosis + Age + Sex + RIN. Tier2 requires nominal diagnosis p <0.01 and absolute beta/SD(expression) >0.8 simultaneously. Regional BH values are separate evidence. The later matched limma-voom and DESeq2 results are model-sensitivity benchmarks, not replacements for the frozen LM or Tier2. See `STEP3_PROVENANCE_RECONCILIATION.md`.

Corrected weighted prioritization combines v1 Tier2_fixed event support (125 events/region) and strict-v2 train-only selector-grid support (150 non-baseline events/region). Combined frequency is 0.60*v2 + 0.40*v1. The final score weights combined frequency 0.40, model support 0.15, sparse support 0.10, importance rank 0.15, sign consistency 0.10, LM evidence 0.05 and any-v2 support 0.05. Full definitions and gates remain in `ML_PRIORITY_SCORE_FINAL.md`. The 30/119/522 labels are method-specific internal prioritization labels.

Primary Tier2-RRA keeps all 671 screened pairs. Five model families (L1, L2, elastic-net logistic, linear SVM, random forest) use five outer folds times five repeats, with three inner folds. Training-only imputation/scaling is applied. RRA first aggregates 25 fit lists per region/model, then the five model lists. Performance and coefficient signs do not weight RRA. RRA scores are not gene-wise diagnosis FDR values. All-gene RRA remains a separate train-only-selection sensitivity analysis.

LASSO, Stability Selection and the conditional-null pilot provide dependent within-cohort evidence. The final synthesis preserves RRA ordering without creating a composite score. WGCNA assignments, eigengene associations, gene-set enrichment and AHBA healthy-brain context answer different questions; their significance must not be interchanged.

Literature comparisons use the already downloaded Labonte, Goes, GSE80655-specific Ramaker and Lanz supplements and the existing mappings. Published significance and tested-gene coverage are distinct. Labonte uses the same GSE102556 dataset, nominal sex-specific DEG lists, and no equivalent standardized-effect cutoff. No new analysis or literature selection was performed in consolidation.
