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
