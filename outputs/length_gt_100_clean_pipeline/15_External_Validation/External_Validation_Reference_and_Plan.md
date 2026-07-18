# External Validation Reference and Plan

Generated for the preliminary exploratory MDD transcriptomics manuscript package.

No external data were downloaded or analyzed for this step. This document prepares references and a concrete validation plan only.

## Summary

External validation will project fixed prioritized gene/module sets from the clean GSE102556 discovery workflow into independent postmortem MDD datasets. The primary endpoint will be module-score MDD vs Control effect direction, standardized effect size, and enrichment/overlap consistency. The validation stage should not retrain a large classifier in external datasets.

## Candidate Datasets

| dataset | primary paper | GEO accession | platform | regions | why useful | proposed validation test | caveats |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Discovery dataset | Labonte et al. 2017, Nature Medicine, "Sex-specific transcriptional signatures in human depression" | GSE102556 | RNA-seq | multiple MDD-relevant brain regions including frontal, hippocampal, insula, and ventral striatal/nucleus accumbens-related regions | Current discovery dataset for the clean pipeline | Not a validation dataset; used to define prioritized genes/modules | Same-dataset derivation limits external claims |
| Primary validation target | Ramaker et al. 2017, Genome Medicine, "Post-mortem molecular profiling of three psychiatric disorders" | GSE80655 | RNA-seq | anterior cingulate cortex, dorsolateral prefrontal cortex, nucleus accumbens | Best primary validation candidate because it includes MDD and controls, relevant frontal/striatal regions, and transdiagnostic groups | Project fixed WGCNA/Tier2/ML module gene sets; compute module scores; test MDD vs Control direction/effect size; compare schizophrenia and bipolar groups separately | Region matching is partial; psychiatric comorbidity and cohort differences require cautious interpretation |
| dlPFC-specific validation target | Pantazatos et al. 2017, Molecular Psychiatry, "Whole-transcriptome brain expression and exon-usage profiling in major depression and suicide: evidence for altered glial, endothelial and ATPase activity" | GSE101521 | RNA-seq / small RNA-seq | dorsolateral prefrontal cortex, Brodmann area 9 | Useful for dlPFC BA8/9-focused validation of prioritized frontal modules | Project dlPFC module scores and test MDD vs Control direction/effect size; evaluate overlap with glial/endothelial/ATPase themes | Includes suicide-related MDD strata; sample design differs from discovery dataset |
| Cross-platform exploratory validation | Lanz et al. 2019, Translational Psychiatry, "Postmortem transcriptional profiling reveals widespread increase in inflammation in schizophrenia..." | GSE53987 | Affymetrix Human Genome U133 Plus 2.0 microarray | prefrontal cortex/BA46, hippocampus, associative striatum | Useful for cross-platform and transdiagnostic checks across regions overlapping the discovery framework | Map genes to array probes; compute fixed module scores where enough genes map; test MDD vs Control and compare with schizophrenia/bipolar patterns | Microarray/RNA-seq platform differences; region definitions not identical; MDD effects may be weaker than schizophrenia effects |
| AHBA context resource | Hawrylycz et al. 2015, Nature Neuroscience, "Transcriptional architecture of the human brain" | Allen Human Brain Atlas | Microarray atlas | adult human brain | Healthy-brain coexpression reference context | Use only as healthy-brain coexpression context for overlap/enrichment | Not a disease validation dataset |
| AHBA processing reference | Markello et al. 2021, eLife, "Standardizing workflows in imaging transcriptomics with the abagen toolbox" | abagen/AHBA tooling | software/workflow reference | AHBA processing | Relevant because AHBA fetching/processing used abagen in this project | Cite for AHBA workflow/tooling context | Not itself a disease dataset |

## Proposed Validation Steps

1. Freeze prioritized gene/module sets from the discovery workflow before loading validation labels.
2. For each external dataset, harmonize identifiers to gene symbols/IDs and keep only genes present in both the validation dataset and fixed module set.
3. Compute module scores per sample using fixed gene membership. Prefer mean/median z-score scores and PC1 scores fit within each validation cohort without using labels for feature selection.
4. Test MDD vs Control effect direction and standardized effect size for each module score.
5. Report nominal p-values and BH-FDR within each external dataset, while emphasizing effect direction and consistency over significance.
6. For GSE80655 and GSE53987, run transdiagnostic specificity checks comparing MDD, schizophrenia, bipolar disorder, and controls.
7. Avoid retraining large classifiers in validation datasets. If classification is performed later, restrict it to small predefined module-score panels and report it as exploratory.

## Transdiagnostic Specificity Analysis

Because GSE80655 and GSE53987 include schizophrenia and bipolar disorder in addition to MDD, future analyses can test whether prioritized module scores are MDD-enriched or reflect broader transdiagnostic psychiatric variation. This is especially important for inflammation, glial, synaptic, and stress-response modules that may not be specific to MDD.

## Current Status

External validation was not performed in this preliminary version and is required before generalization or biomarker claims.
