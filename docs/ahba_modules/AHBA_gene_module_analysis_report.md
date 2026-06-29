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
