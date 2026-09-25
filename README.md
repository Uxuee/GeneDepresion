# GeneDepresion

## START HERE — canonical publication workflow

The authoritative pipeline definition is
[`provenance/canonical_pipeline.json`](provenance/canonical_pipeline.json), with a
plain-language map in
[`docs/CANONICAL_PIPELINE_MAP.md`](docs/CANONICAL_PIPELINE_MAP.md). Validate the
state with:

1. [Canonical pipeline map](docs/CANONICAL_PIPELINE_MAP.md)
2. [Machine-readable registry](provenance/canonical_pipeline.json)
3. [Local path configuration example](config/project_paths.example.yaml)
4. [External data instructions](data/README.md)
5. [Environment status](environment/README.md)

```text
python code/run_canonical_pipeline.py --validate
python code/run_canonical_pipeline.py --dry-run
```

The canonical corrected ML result contains **671 gene-region pairs: 30 high,
119 moderate, and 522 low confidence** (149 high+moderate). The v1 and v2 event
denominators are 125 and 150. The historical
`10_ML_prioritized_Tier2_extended` directory is superseded and must not be used
as the current result.

Scientific baseline: `ca25ec35af359a4d409abf0f728833d2ff54761c`.
The canonical output tree is `outputs/length_gt_100_clean_pipeline`, but only
artifacts named by the registry are canonical. Other folders may be historical
and must be used only for explicit provenance investigation.

Publication guardrails: there are no validated biomarkers; no MDD WGCNA module
survived Diagnosis BH-FDR; AHBA provides context only; module-level ML is
exploratory; external validation is incomplete. Tier4 GSEA and historical MDD
module-level ML2 are excluded until their exact generators are recovered.

> Historical navigation and workflow notes follow. They are retained for
> provenance and may contain superseded paths; use the canonical map above for
> current work.

<!-- HISTORICAL_CONTENT_START -->

## READ THIS FIRST: Current Clean MDD Report

The main report for the current Length > 100 MDD transcriptomics workflow is here:

[**Open the step-by-step Length > 100 MDD transcriptomics report**](docs/Step_by_step_MDD_Length_gt_100_clean_report.md)

This is the report-style version with the analysis story, formulas, tables, figures, interpretation, limitations, and next steps. It is the best starting point for reading the project.

| Start here | Link |
| --- | --- |
| Main step-by-step report | [docs/Step_by_step_MDD_Length_gt_100_clean_report.md](docs/Step_by_step_MDD_Length_gt_100_clean_report.md) |
| Clear workflow report | [docs/Main_Workflow_Length_gt_100_Report.md](docs/Main_Workflow_Length_gt_100_Report.md) |
| Clear workflow notebook | [outputs/length_gt_100_clean_pipeline/notebooks/01_clear_length_gt_100_workflow_notebook.md](outputs/length_gt_100_clean_pipeline/notebooks/01_clear_length_gt_100_workflow_notebook.md) |
| Report-ready figures | [outputs/length_gt_100_clean_pipeline/images/](outputs/length_gt_100_clean_pipeline/images/) |
| Final ML-prioritized Tier 2 results | [outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_extended/](outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_extended/) |
| Current AHBA module enrichment | [docs/ahba_modules/AHBA_gene_module_analysis_report.md](docs/ahba_modules/AHBA_gene_module_analysis_report.md) |
| Technical workflow report | [docs/MDD_Length_gt_100_clean_Tier2_to_ML_prioritized_Report.md](docs/MDD_Length_gt_100_clean_Tier2_to_ML_prioritized_Report.md) |

Important interpretation: Tier 2 and ML-prioritized genes are exploratory candidate signatures, not final FDR-significant disease biomarkers. ML v1 is used for internal gene prioritization and feature stability, not validated diagnostic classification.

---

## Older Baseline Workflow Notes

RNA-seq preprocessing and normalization workflow for GSE102556.

The current Step 1 script prepares region-specific expression matrices using one normalization path:

```text
filtered raw counts -> GeTMM -> log2(GeTMM + 1)
```


## Current Clean Length > 100 Analysis

The current report-ready workflow is the clean `Length > 100` MDD transcriptomics path. It preserves the exploratory interpretation of Tier 2 genes and organizes the analysis from normalization through LMM-derived candidate genes, enrichment, ML-prioritized Tier 2 genes, review figures, Step 11 AHBA contextualization, and Step 12 MDD WGCNA modules.

Start here:

- [Step-by-step main report](docs/Step_by_step_MDD_Length_gt_100_clean_report.md) - narrative report from LMM/statistical evidence to Tier 2, enrichment, ML-prioritized genes, figures, interpretation, and next steps.
- [Clear workflow report](docs/Main_Workflow_Length_gt_100_Report.md) - compact map of the current clean workflow, inputs, methods, result counts, and output locations.
- [Clear workflow notebook](outputs/length_gt_100_clean_pipeline/notebooks/01_clear_length_gt_100_workflow_notebook.md) - notebook-style walkthrough of the same workflow inside the output folder.
- [Technical clean workflow report](docs/MDD_Length_gt_100_clean_Tier2_to_ML_prioritized_Report.md) - detailed pipeline/archive-style report.
- [Report-ready image folder](outputs/length_gt_100_clean_pipeline/images/) - figures collected from Tier 2, enrichment, and ML-prioritization steps.
- [Final ML-prioritized Tier 2 tables](outputs/length_gt_100_clean_pipeline/10_ML_prioritized_Tier2_extended/) - high-confidence and moderate-confidence candidate gene tables.
- [AHBA module analysis report](docs/ahba_modules/AHBA_gene_module_analysis_report.md) - current healthy-brain published-module enrichment results for MDD candidate genes.

Important interpretation: Tier 2 and ML-prioritized genes are exploratory candidate signatures, not final FDR-significant disease biomarkers. ML v1 is used for internal prioritization and feature stability, not validated diagnostic classification.

## Current AHBA Module Enrichment Status

AHBA module enrichment is now part of the current clean `Length > 100` workflow as a healthy-brain coexpression context step. It tests whether MDD Tier 2 / ML-prioritized candidate genes overlap published adult human brain canonical coexpression modules from Hawrylycz et al. 2015.

Important method note: de novo AHBA WGCNA was attempted first, but the parameter-sensitivity run collapsed into one very large module, so it was not used for biological interpretation. The current reported AHBA enrichment uses the published canonical AHBA modules instead.

Main AHBA locations:

- Report: [docs/ahba_modules/AHBA_gene_module_analysis_report.md](docs/ahba_modules/AHBA_gene_module_analysis_report.md)
- Output folder: [outputs/length_gt_100_clean_pipeline/ahba_modules/published_modules/](outputs/length_gt_100_clean_pipeline/ahba_modules/published_modules/)
- Enrichment results: [outputs/length_gt_100_clean_pipeline/ahba_modules/published_modules/enrichment/published_ahba_module_enrichment_results.csv](outputs/length_gt_100_clean_pipeline/ahba_modules/published_modules/enrichment/published_ahba_module_enrichment_results.csv)
- Region summary: [outputs/length_gt_100_clean_pipeline/ahba_modules/published_modules/enrichment/published_ahba_enrichment_region_summary.csv](outputs/length_gt_100_clean_pipeline/ahba_modules/published_modules/enrichment/published_ahba_enrichment_region_summary.csv)
- Figures: [outputs/length_gt_100_clean_pipeline/ahba_modules/published_modules/figures/](outputs/length_gt_100_clean_pipeline/ahba_modules/published_modules/figures/)

Current published-module run:

- Published AHBA table: 15,627 unique genes across 32 modules.
- Enrichment background: 12,042 genes shared between the clean MDD tested-gene universe and the published AHBA module table.
- Candidate sets tested: high-confidence ML only, high + moderate ML, all Tier 2, region-specific sets, and direction-specific sets.

Main corrected AHBA enrichment signals:

| Candidate set | Region | Direction | AHBA module | Overlap | BH-FDR | Interpretation |
| --- | --- | --- | --- | ---: | ---: | --- |
| all Tier 2 | Ventral Striatum | up in MDD | M26 | 9 | 0.000166 | strongest corrected AHBA signal |
| all Tier 2 | Cingulate gyrus 25 | down in MDD | M2 | 4 | 0.012970 | corrected AHBA signal |
| all Tier 2 | Ventral Striatum | down in MDD | M6 | 8 | 0.019053 | corrected AHBA signal |
| all Tier 2 | Hippocampal Formation | down in MDD | M7 | 3 | 0.051894 | suggestive, FDR < 0.10 |
| all Tier 2 | Hippocampal Formation | down in MDD | M28 | 5 | 0.051894 | suggestive, FDR < 0.10 |
| high + moderate ML | Ventral Striatum | up in MDD | M26 | 3 | 0.094393 | suggestive ML-prioritized support |

Conservative interpretation: AHBA enrichment gives healthy-brain module context for the MDD candidate genes. It does not validate the genes as MDD biomarkers, and it should be reported as coexpression-context support rather than disease-specific proof.

## Main Scripts

- `run_step1.R` - easiest entry point; runs the current recommended Step 1 workflow.
- `code/Step1_GSE102556_GeTMM_brain_area_groups.R` - recommended Step 1; groups exact GEO tissue labels into broader editable brain areas before GeTMM-log2 normalization.
- `code/Step2_Covariate_Exploration_FrontalCortex.R` - PCA/QC and covariate exploration for Frontal Cortex.
- `code/Step3_small_groups_limma_voom_simple.R` - limma-voom differential expression for Insula, Ventral Striatum, and Hippocampal Formation.
- `code/Step3_frontal_dream_simple.R` - dream differential expression for Frontal Cortex with repeated donor handling.
- `code/Step3_all_regions_combined_simple.R` - combined all-region dream model for a shared depression effect.
- `code/Step1_length_cutoff_sensitivity_GeTMM.R` - reruns Step 1 GeTMM preprocessing for several gene-length cutoffs without overwriting the baseline Step 1 outputs.
- `code/Step3_run_all_DGE_methods_on_length_cutoffs.R` - runs limma-voom, dream, and normalized-expression comparator models across the length-cutoff outputs.

## Brain-Area Mapping

The broad grouping is controlled by:

- `config/GSE102556_brain_area_mapping_template.csv`

Default grouping:

- OFC, dlPFC, Cg25 -> `Frontal_Cortex`
- aINS -> `Insula`
- NAc -> `Ventral_Striatum`
- Subiculum -> `Hippocampal_Formation`

Edit the mapping CSV and rerun the brain-area script if a different anatomical grouping is needed.

## Documentation

- `docs/Step1_GSE102556_GeTMM_log2_Notion.md`
- `docs/Data_Treatment_GSE102556_Notion.md`
- `docs/Full_Process_Simple_Explanation.md`
- `docs/Step3_Length_Cutoff_Sensitivity_Report.md`

## Current Analysis Summary

The current workflow is:

```text
raw counts
-> metadata cleaning
-> brain-area grouping
-> GeTMM normalization
-> log2(GeTMM + 1)
-> PCA/QC inspection
-> transcriptomics differential expression with limma-voom / dream
```

The recommended Step 3 differential expression scripts now use:

- limma-voom for `Insula`, `Ventral_Striatum`, and `Hippocampal_Formation`
- dream / variancePartition for `Frontal_Cortex`
- dream / variancePartition for the combined all-region shared-effect model

The earlier naive one-gene-at-a-time LM/LMM scripts are superseded.

Latest Step 3 run summary:

| Analysis | Genes tested | Raw p < 0.05 | FDR < 0.05 | FDR < 0.10 | Minimum FDR | Top gene |
|---|---:|---:|---:|---:|---:|---|
| Insula, limma-voom | 23,375 | 941 | 0 | 0 | 0.8351945 | 10090 |
| Ventral_Striatum, limma-voom | 24,275 | 2,472 | 0 | 2 | 0.09301074 | 107987147 |
| Hippocampal_Formation, limma-voom | 24,700 | 800 | 0 | 0 | 0.9999502 | 65083 |
| Frontal_Cortex, dream | 24,286 | 1,490 | 0 | 0 | 0.1265182 | 5147 |
| All regions combined, dream | 22,575 | 707 | 0 | 0 | 0.9995043 | 5468 |

Interpretation:

- The dream models ran successfully after package-version fixes.
- Frontal Cortex had one gene/model failure, so the final dream result table has 24,286 tested genes instead of 24,287.
- The newer models show more raw p-value signal in some groups, but almost nothing survives FDR correction.
- The only FDR < 0.10 results are two Ventral Striatum genes from the limma-voom small-group analysis.
- The result CSV files are local outputs and are not committed to GitHub.

## Length-Cutoff Sensitivity Analysis

The baseline gene-length filter is `Length > 20`.

To test whether that cutoff is too permissive, the repository now includes a
conservative sensitivity analysis for:

```text
Length > 20
Length > 50
Length > 100
Length > 200
Length > 500
```

Run the scripts in this order:

```r
source("code/Step1_length_cutoff_sensitivity_GeTMM.R")
source("code/Step3_run_all_DGE_methods_on_length_cutoffs.R")
```

The Step 1 sensitivity script writes:

```text
outputs/step1_length_sensitivity
```

The Step 3 sensitivity script writes:

```text
outputs/step3_length_sensitivity
```

Important: the Step 3 sensitivity script can take a long time because it runs
dream and per-gene LMM comparator models for every length cutoff.

Do not choose a cutoff only because it gives more FDR-significant genes.
Use the summary, top-gene stability, cross-method overlap, and effect-direction
correlation tables together.

## Simple Step-by-Step Explanation

This section explains the whole data process from normalization until the current
differential expression analysis.

### 1. Start With Raw RNA-seq Counts

The starting data are raw RNA-seq gene counts from GSE102556.

The raw count table has:

- genes in rows
- samples in columns
- count values showing how many sequencing reads mapped to each gene

Raw counts cannot be compared directly across samples because samples can have
different sequencing depths, and genes have different lengths.

### 2. Clean and Standardize Metadata

The workflow loads sample metadata from GEO.

Important metadata columns are cleaned into standard names:

- `SampleID`
- `Diagnosis`
- `Age`
- `Sex`
- `RIN`
- `OriginalRegion`
- `BrainID`

`BrainID` identifies the donor/person. This matters because the same donor can
appear more than once if multiple brain regions were sampled.

### 3. Group Brain Regions

The original tissue labels are grouped into broader analysis groups.

| Original region | Analysis group |
|---|---|
| dlPFC | Frontal_Cortex |
| OFC | Frontal_Cortex |
| Cg25 | Frontal_Cortex |
| Anterior insula | Insula |
| Nucleus accumbens | Ventral_Striatum |
| Subiculum | Hippocampal_Formation |

This mapping is controlled by:

```text
config/GSE102556_brain_area_mapping_template.csv
```

### 4. Match Expression Samples to Metadata

The workflow keeps only samples that are present in both:

- the expression count matrix
- the metadata table

This avoids mismatches between expression columns and metadata rows.

### 5. Filter and Annotate Genes

Genes are matched to the annotation file so gene length can be used.

Very short or unusable genes are removed before normalization.

### 6. Normalize With GeTMM-log2

The project uses one normalization path:

```text
raw counts
-> gene-length adjustment
-> TMM library-size adjustment
-> GeTMM
-> log2(GeTMM + 1)
```

GeTMM adjusts for:

- gene length
- sequencing-depth differences between samples

The `log2(GeTMM + 1)` transform makes expression values easier to model and
keeps zero values valid.

### 7. Save Normalized Expression and Metadata

Step 1 writes the main downstream files:

```text
normalized_GeTMM_log2_<group>.csv
metadata_<group>.csv
```

The expression files have:

- genes in rows
- samples in columns

The metadata files have:

- one row per sample

### 8. Run QC and PCA

PCA/QC was used to inspect sample structure before differential expression.

The important finding was that `OriginalRegion` strongly matters inside
`Frontal_Cortex`, because that group combines dlPFC, OFC, and Cg25.

Therefore, `OriginalRegion` is included as a covariate in the Frontal Cortex
model.

### 9. Choose the Differential Expression Models

The current recommended Step 3 uses transcriptomics-specific models.

For Insula, Ventral Striatum, and Hippocampal Formation, the script uses
limma-voom:

```r
~ Diagnosis + Age + Sex + RIN
```

limma-voom starts from filtered raw counts, estimates mean-variance weights, and
uses empirical Bayes moderation.

For Frontal Cortex, the script uses dream because Frontal Cortex can contain
repeated donors:

```r
~ Diagnosis + Age + Sex + RIN + OriginalRegion + (1 | BrainID)
```

For the combined all-region analysis, the script also uses dream:

```r
~ Diagnosis + BrainArea + Age + Sex + RIN + (1 | BrainID)
```

The combined model tests whether a shared depression effect across regions is
easier to detect than separate regional analyses.

### 10. Correct for Testing Thousands of Genes

Each gene gets a diagnosis p-value.

Because thousands of genes are tested, the workflow corrects p-values using
Benjamini-Hochberg FDR.

FDR means false discovery rate.

Simple interpretation:

- p-value: evidence for one gene by itself
- FDR: corrected evidence after testing thousands of genes
- FDR < 0.05: strict significance threshold
- FDR < 0.10: exploratory threshold

### 11. What the New Result Files Will Mean

The new Step 3 scripts save result CSV files and simple p-value histograms under:

```text
outputs/step3_limma_voom_dream
```

The CSV files contain:

- gene ID
- diagnosis coefficient
- average expression
- moderated t-statistic
- raw p-value
- BH/FDR-corrected p-value

The p-value histograms help show whether the p-values look mostly uniform or
whether there is a left-skewed signal that is not strong enough to survive FDR.

The previous naive LM/LMM result summary was:

| Group | Samples used | Genes tested | FDR < 0.05 | FDR < 0.10 |
|---|---:|---:|---:|---:|
| Frontal_Cortex | 122 | 24,287 | 0 | 1 |
| Insula | 48 | 23,375 | 0 | 0 |
| Ventral_Striatum | 48 | 24,275 | 0 | 0 |
| Hippocampal_Formation | 43 | 24,700 | 0 | 0 |

Those older results are kept only as context. They are not the recommended final
Step 3 workflow.

### 12. Important Caution

The dream scripts require the Bioconductor package `variancePartition`.

If it is not installed, install it with:

```r
BiocManager::install("variancePartition")
```

In this Windows/R 4.4.2 run, dream worked with:

```text
variancePartition 1.36.3
lmerTest 3.2.1
lme4 1.1.37
```

The local newer `lme4 2.0.1` caused dream model extraction errors, so `lme4`
was changed back to `1.1.37` for this analysis.

## Inputs

The script expects local input files under:

```text
D:/genes/genedata
```

Required:

- `GSE102556_raw_counts_GRCh38.p13_NCBI.tsv`
- `Human.GRCh38.p13.annot.tsv`

Optional:

- `GSE102556_metadata_clean.csv`

If the optional metadata file is missing, the script fetches GEO metadata using:

```r
GEOquery::getGEO("GSE102556", GSEMatrix = TRUE)
```

## Outputs

The brain-area grouped Step 1 writes outputs to:

```text
D:/genes/genedata/Step1_GSE102556_GeTMM_brain_area_groups
```

Output folders are organized by file type:

```text
00_config/                 mapping used for the run
01_metadata/               matched metadata and per-group metadata
02_normalized_expression/  final log2(GeTMM + 1) matrices
03_filtered_raw_counts/    raw counts after gene filtering
04_normalization_factors/  TMM and effective-library-size tables
05_gene_annotation/        annotation for retained genes
06_repetition_reports/     repeated donor / gene-entry reports
07_qc/                     compact QC and model-decision summaries
08_rds/                    compact R object for later steps
```

Large input data and generated CSV outputs are intentionally ignored by Git.

