# AHBA Gene Module Analysis Report

## Aim

The AHBA module-analysis stage adds healthy-brain coexpression context to the
clean Length > 100 MDD transcriptomics workflow. It tests whether
ML-prioritized MDD Tier 2 candidate genes are overrepresented in gene
coexpression modules built from Allen Human Brain Atlas (AHBA) microarray data.

AHBA module analysis provides healthy-brain coexpression context for
ML-prioritized MDD candidate genes. It does not validate them as disease
biomarkers.

## Data Source

The stage expects AHBA donor-level microarray files under:

```text
data/ahba/
```

The fetch script first checks for local files. If missing, it attempts:

```python
abagen.fetch_microarray(donors="all", data_dir="data/ahba")
```

Raw AHBA downloads are ignored by git.

## Clean Pipeline Candidate Inputs

The source of truth for MDD candidate genes is:

```text
outputs/length_gt_100_clean_pipeline/
```

Selected input files:

| File | Rows | Role |
|---|---:|---|
| `10_ML_prioritized_Tier2_extended/final_ranked_ML_prioritized_Tier2_genes.csv` | 671 | all ranked Tier 2 gene-region pairs |
| `10_ML_prioritized_Tier2_extended/final_high_confidence_ML_prioritized_Tier2_genes.csv` | 35 | high-confidence ML-prioritized Tier 2 gene-region pairs |
| `10_ML_prioritized_Tier2_extended/final_moderate_confidence_ML_prioritized_Tier2_genes.csv` | 130 | additional moderate-confidence ML-prioritized Tier 2 gene-region pairs |
| `07_Tier2_candidate_gene_sets/Tier2_candidate_genes_all_regions.csv` | 671 | all Tier 2 candidate gene-region pairs |

Candidate counts by region:

| Region/subregion | All Tier 2 | High-confidence | Additional moderate-confidence |
|---|---:|---:|---:|
| Cingulate gyrus 25 | 186 | 8 | 29 |
| Hippocampal Formation | 127 | 3 | 25 |
| Insula | 70 | 6 | 18 |
| OFC BA11 | 70 | 3 | 28 |
| Ventral Striatum | 146 | 8 | 18 |
| dlPFC BA8/9 | 72 | 7 | 12 |

Direction counts in the final ranked Tier 2 table:

| Direction | Gene-region pairs |
|---|---:|
| up_in_MDD | 392 |
| down_in_MDD | 279 |

## Preprocessing

`code/ahba_modules/01_preprocess_ahba_microarray.py`:

- loads donor-level AHBA microarray expression files;
- loads probe metadata and sample annotations;
- maps probes to gene symbols;
- removes probes/genes with excessive missingness (>50%);
- collapses probes to one value per gene by selecting the probe with highest
  mean expression across AHBA samples;
- harmonizes AHBA gene symbols with the clean MDD candidate-gene tables;
- writes a gene-by-sample expression matrix and gene mapping report.

## Module Construction

`code/ahba_modules/02_build_ahba_modules.R` builds modules from the processed
AHBA gene-by-sample matrix, preferably with WGCNA.

Planned outputs:

- `modules/ahba_gene_to_module_assignments.csv`
- `modules/ahba_module_size_table.csv`
- `modules/ahba_module_eigengenes.csv`
- `modules/module_parameters.json`
- `modules/AHBA_module_construction_summary.md`



## Published AHBA Canonical Modules

De novo AHBA WGCNA was attempted first using the processed AHBA matrix. The memory-safe run and a 72-combination WGCNA parameter sensitivity grid did not produce a useful multi-module structure: every completed setting produced one giant non-grey module containing nearly all 5,000 selected genes, with only three grey/unassigned genes. Because that result is not biologically useful for module enrichment, de novo AHBA modules are documented but not used for biological enrichment.

The primary AHBA module-reference branch now uses published adult human brain canonical coexpression modules from Hawrylycz et al. 2015, *Canonical genetic signatures of the adult human brain*. The expected input is Supplementary Table 6, containing consensus module/signature assignments and eigengene/module-membership correlations for genes assigned to the published network.

Published AHBA canonical modules provide healthy-brain coexpression context for ML-prioritized MDD candidate genes. They do not validate these genes as disease biomarkers.

Scripts:

```text
code/ahba_modules/04_fetch_or_load_published_ahba_modules.py
code/ahba_modules/05_test_published_ahba_module_enrichment.py
```

Manual input folder:

```text
data/ahba/published_modules/
```

Standardized published-module output:

```text
outputs/length_gt_100_clean_pipeline/ahba_modules/published_modules/ahba_published_gene_to_module_assignments.csv
```

The published-module enrichment background is not all human genes. It is:

```text
unique genes present in both the MDD clean-pipeline tested-gene universe and the published AHBA module assignment table
```

Gene-region candidate rows are deduplicated to unique gene symbols within each tested candidate set before Fisher/hypergeometric enrichment testing.

## Candidate Gene Sets

The enrichment script tests:

- high-confidence only;
- high + moderate-confidence;
- all Tier 2;
- region-specific high + moderate sets;
- direction-specific high + moderate sets when at least five genes are present.

## Enrichment Method

`code/ahba_modules/03_test_ahba_module_enrichment.py` uses one-sided
hypergeometric tests for overrepresentation of candidate genes in AHBA modules.

The background is not all human genes. It is:

```text
genes present in both the MDD tested-gene universe and the processed AHBA module universe
```

BH-FDR correction is applied within each candidate-set/scope family.



## Published Module Enrichment Run

Hawrylycz et al. 2015 Supplementary Table 6 was loaded from the official Springer/Nature supplementary CSV. The standardized table contains 15,627 unique genes across 32 published modules. The enrichment background contains 12,042 genes, defined as the intersection of the clean MDD tested-gene universe and the published AHBA module-assignment universe.

All-region candidate overlaps with that background were:

| candidate_set      |   original_gene_region_pair_count |   unique_gene_count |   unique_genes_overlapping_published_AHBA_background |   modules_tested |
|:-------------------|----------------------------------:|--------------------:|-----------------------------------------------------:|-----------------:|
| all_Tier2          |                               671 |                 640 |                                                  276 |               32 |
| high_confidence    |                                35 |                  35 |                                                   12 |               32 |
| high_plus_moderate |                               165 |                 154 |                                                   56 |               32 |

The published-module enrichment produced 41 nominal rows with p < 0.05, 6 rows with BH-FDR < 0.10, and 3 rows with BH-FDR < 0.05. These are healthy-brain module-context findings, not disease-biomarker validation.



## Published AHBA Enrichment Interpretation By Region

The report-ready enrichment tables are stored in:

```text
outputs/length_gt_100_clean_pipeline/ahba_modules/published_modules/enrichment/
```

Main corrected result: Ventral Striatum up-in-MDD all Tier 2 genes were enriched in published AHBA module M26 (9 overlapping genes; raw p = 5.2e-06; BH-FDR = 0.000166). Additional BH-FDR < 0.05 results were Cingulate gyrus 25 down-in-MDD all Tier 2 genes in M2 and Ventral Striatum down-in-MDD all Tier 2 genes in M6.

Suggestive BH-FDR < 0.10 results included Hippocampal Formation down-in-MDD all Tier 2 genes in M7 and M28, plus high + moderate Ventral Striatum up-in-MDD genes in M26.

Nominal-only patterns were present in several additional region/module combinations, including OFC BA11, Insula, Cingulate gyrus 25, Hippocampal Formation, and Ventral Striatum. Nominal results are reported as exploratory patterns and should not be interpreted as corrected discoveries.

Region-level calls from the best region-specific hit were:

| region                | best_module_hit   | best_candidate_set   | best_direction   |   best_raw_p_value |   best_BH_FDR | evidence_call   |
|:----------------------|:------------------|:---------------------|:-----------------|-------------------:|--------------:|:----------------|
| Ventral Striatum      | M26               | all_Tier2            | up_in_MDD        |        5.19607e-06 |   0.000166274 | FDR < 0.05      |
| Insula                | M10               | all_Tier2            | up_in_MDD        |        0.00910301  |   0.291296    | nominal only    |
| Hippocampal Formation | M7                | all_Tier2            | down_in_MDD      |        0.00255056  |   0.051894    | FDR < 0.10      |
| Cingulate gyrus 25    | M2                | all_Tier2            | down_in_MDD      |        0.000405302 |   0.0129697   | FDR < 0.05      |
| dlPFC BA8/9           | M27               | all_Tier2            | all_directions   |        0.0191628   |   0.613209    | nominal only    |
| OFC BA11              | M28               | all_Tier2            | down_in_MDD      |        0.005623    |   0.179936    | nominal only    |

The corrected findings are driven mainly by all Tier 2 direction-specific gene sets. The high + moderate ML-prioritized branch contributes one suggestive corrected result, and the high-confidence-only set does not show corrected enrichment in this AHBA module run.

## Main Results

Results are produced only after AHBA raw data have been downloaded or placed
locally and modules have been constructed. The stage writes outputs to:

```text
outputs/length_gt_100_clean_pipeline/ahba_modules/
```

## Figures

Planned figures:

- candidate/module enrichment dot plot;
- region-by-module candidate overlap heatmap;
- module annotation plots if module annotation is later added.



## AHBA Module Annotations

Hawrylycz et al. 2015 Supplementary Table 8 / MOESM99 was loaded from the official Springer supplementary material source and used to annotate modules appearing in corrected, suggestive, and top nominal enrichment rows. The full annotation outputs are stored in:

```text
outputs/length_gt_100_clean_pipeline/ahba_modules/published_modules/annotations/
```

Key annotated modules included M26, M2, M6, M7, M28, and additional nominal modules. Selected module annotations:

| AHBA_module_ID   |   module_gene_count | top_annotation_terms                                                                                                                          |   best_annotation_q_value | short_biological_interpretation                                                                                                                                                                                                                |
|:-----------------|--------------------:|:----------------------------------------------------------------------------------------------------------------------------------------------|--------------------------:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| M2               |                 197 | hsa-miR-448:PITA_TOP; central nervous system neuron development; hsa-miR-548l:PITA_TOP; hsa-miR-606:PITA_TOP; hsa-miR-548n:PITA_TOP           |                  2.71e-05 | Top Hawrylycz enrichment terms point to hsa-miR-448:PITA_TOP; central nervous system neuron development; hsa-miR-548l:PITA_TOP; hsa-miR-606:PITA_TOP; hsa-miR-548n:PITA_TOP. Categories represented: MicroRNA, GO: Biological Process.         |
| M6               |                 573 | Genes in module_11; Genes in module_137; hsa-miR-219-5p:TargetScan; Genes in module_100; Genes in module_66                                   |                  5.79e-09 | Top Hawrylycz enrichment terms point to Genes in module_11; Genes in module_137; hsa-miR-219-5p:TargetScan; Genes in module_100; Genes in module_66. Categories represented: Computational, MicroRNA.                                          |
| M7               |                 164 | synaptic transmission; cell-cell signaling; Calcium signaling pathway; dendritic spine; neuron spine                                          |                  6.27e-07 | Top Hawrylycz enrichment terms point to synaptic transmission; cell-cell signaling; Calcium signaling pathway; dendritic spine; neuron spine. Categories represented: GO: Biological Process, Pathway, GO: Cellular Component.                 |
| M10              |                 166 | Heterotrimeric G-protein signaling pathway-Gi alpha and Gs alpha mediated pathway; Ro 41-9067; SCH 23390; Haloperidol; Pramipexole            |                  1.29e-07 | Top Hawrylycz enrichment terms point to Heterotrimeric G-protein signaling pathway-Gi alpha and Gs alpha mediated pathway; Ro 41-9067; SCH 23390; Haloperidol; Pramipexole. Categories represented: Pathway, Drug.                             |
| M19              |                1102 | mitochondrion; mitochondrial part; mitochondrial envelope; mitochondrial inner membrane; cellular respiration                                 |                  5.5e-82  | Top Hawrylycz enrichment terms point to mitochondrion; mitochondrial part; mitochondrial envelope; mitochondrial inner membrane; cellular respiration. Categories represented: GO: Cellular Component, GO: Biological Process.                 |
| M24              |                1280 | Gene Expression; nucleoplasm; Zinc fingers, C2H2-type; ZINC_FINGER_C2H2_2; ZINC_FINGER_C2H2_1                                                 |                  9.87e-46 | Top Hawrylycz enrichment terms point to Gene Expression; nucleoplasm; Zinc fingers, C2H2-type; ZINC_FINGER_C2H2_2; ZINC_FINGER_C2H2_1. Categories represented: Pathway, GO: Cellular Component, Gene Family, Domain.                           |
| M26              |                 795 | cilium; ciliary part; cilium morphogenesis; cilium organization; cilium assembly                                                              |                  6.09e-45 | Top Hawrylycz enrichment terms point to cilium; ciliary part; cilium morphogenesis; cilium organization; cilium assembly. Categories represented: GO: Cellular Component, GO: Biological Process.                                              |
| M27              |                 855 | mitochondrion; endoplasmic reticulum part; endoplasmic reticulum; tert-Butylhydroperoxide; Diuron                                             |                  2.98e-17 | Top Hawrylycz enrichment terms point to mitochondrion; endoplasmic reticulum part; endoplasmic reticulum; tert-Butylhydroperoxide; Diuron. Categories represented: GO: Cellular Component, Drug.                                               |
| M28              |                1132 | GPCR_Rhodpsn_supfam; Homeobox_CS; transmembrane signaling receptor activity; signaling receptor activity; G-protein coupled receptor activity |                  1.69e-07 | Top Hawrylycz enrichment terms point to GPCR_Rhodpsn_supfam; Homeobox_CS; transmembrane signaling receptor activity; signaling receptor activity; G-protein coupled receptor activity. Categories represented: Domain, GO: Molecular Function. |
| M30              |                1838 | hsa-miR-506:TargetScan; hsa-miR-124:TargetScan; hsa-miR-124:PITA_TOP; GTGCCTT,MIR-506:MSigDB; hsa-miR-590-3p:PITA_TOP                         |                  1.33e-25 | Top Hawrylycz enrichment terms point to hsa-miR-506:TargetScan; hsa-miR-124:TargetScan; hsa-miR-124:PITA_TOP; GTGCCTT,MIR-506:MSigDB; hsa-miR-590-3p:PITA_TOP. Categories represented: MicroRNA.                                               |

These annotations describe healthy adult-brain coexpression-module enrichments. They provide biological context only and are not disease validation.

## Exploratory AHBA Module-Score ML

Published AHBA modules were used as external predefined biological feature sets. Module-score ML was exploratory and was used to assess whether AHBA-informed module-level expression patterns carried MDD/Control signal, not to establish a validated diagnostic classifier.

Main leakage-aware setup:

- Feature set: all 32 published AHBA modules.
- Score definition: within each CV fold, gene expression was z-scored using training samples only; held-out samples were transformed using training means and standard deviations; each module score was the mean standardized expression of available module genes.
- CV: 20 repeated seeds; grouped stratified folds were used where BrainID repeated within a region.
- Models: majority baseline, logistic regression L2/L1/elastic net, linear SVM, and exploratory random forest.
- Post-enrichment module subsets, including M26-only and FDR-selected modules, were included only as sensitivity analyses because those modules were selected after enrichment.

Best all-32-module model per broad region:

| region                | model       | cv_type                        |   n_samples |   n_MDD |   n_Control |   balanced_accuracy_mean |   balanced_accuracy_std |   roc_auc_mean |   roc_auc_std | leakage_status                    |
|:----------------------|:------------|:-------------------------------|------------:|--------:|------------:|-------------------------:|------------------------:|---------------:|--------------:|:----------------------------------|
| Frontal_Cortex        | linear_SVM  | StratifiedGroupKFold           |         122 |      63 |          59 |                 0.58668  |                0.106298 |       0.632764 |      0.128532 | predefined_external_all32_modules |
| Hippocampal_Formation | logistic_L2 | StratifiedKFold_repeated_seeds |          43 |      24 |          19 |                 0.540083 |                0.164577 |       0.483833 |      0.175832 | predefined_external_all32_modules |
| Insula                | linear_SVM  | StratifiedKFold_repeated_seeds |          48 |      26 |          22 |                 0.630917 |                0.138324 |       0.642483 |      0.176715 | predefined_external_all32_modules |
| Ventral_Striatum      | linear_SVM  | StratifiedKFold_repeated_seeds |          48 |      26 |          22 |                 0.595333 |                0.141175 |       0.630817 |      0.173329 | predefined_external_all32_modules |

Best sensitivity models by region/feature set are stored in:

```text
outputs/length_gt_100_clean_pipeline/ahba_modules/module_score_ml/module_score_ml_best_model_by_region_feature_set.csv
```

Interpretation: all-32-module performance was modest and should be described as preliminary module-level signal, not diagnostic performance. Insula showed the strongest all-32 balanced accuracy in this pilot, followed by Ventral Striatum and Frontal Cortex; Hippocampal Formation was close to baseline. Post-enrichment module subsets can look stronger in some regions, but they are explicitly exploratory because the module subsets were chosen after enrichment.

## Limitations

- AHBA is a healthy-brain reference atlas, not an MDD case-control dataset.
- AHBA microarray probes and RNA-seq-style MDD gene symbols require careful
  symbol harmonization.
- Module enrichment is sensitive to the chosen background universe.
- Any enriched module should be interpreted as biological context, not disease
  validation.

## Interpretation

AHBA module analysis provides healthy-brain coexpression context for
ML-prioritized MDD candidate genes. It does not validate them as disease
biomarkers.
