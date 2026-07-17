# MDD WGCNA Top Module Interpretation

This report prioritizes non-grey WGCNA modules from the full Step6 MDD/control run using diagnosis association, Tier 2/ML enrichment, AHBA M26/M28/M2/M6/M7 overlap, and hub-gene interpretability.

Important caution: no module eigengene survived BH-FDR < 0.10 for Diagnosis in the saved full-run module-trait association results. Diagnosis raw p < 0.05 is treated only as exploratory evidence.

## Outputs

- Full triage table: `mdd_wgcna_module_triage_table.csv`
- Report-priority subset: `mdd_wgcna_top_priority_modules_for_report.csv`

## Scoring Rule

Priority score = +3 Diagnosis BH-FDR < 0.10; +2 Diagnosis raw p < 0.05; +3 best Tier2/ML enrichment BH-FDR < 0.05; +2 high+moderate ML enrichment BH-FDR < 0.10; +2 AHBA M26/M28/M2/M6/M7 BH-FDR < 0.05; +1 interpretable hub genes.

## Top Prioritized Modules

| Rank | Region | Module | Size | Score | Diagnosis p | Diagnosis FDR | Best Tier2/ML set | Best Tier2/ML FDR | Best AHBA | Best AHBA FDR | Highlight AHBA evidence | Top hub genes | Interpretation |
| ---: | --- | --- | ---: | ---: | ---: | ---: | --- | ---: | --- | ---: | --- | --- | --- |
| 1 | Ventral Striatum | MElavenderblush3 | 35 | 8 | 0.007441 | 0.3278 | same_region_Tier2_genes | 0.00159 | M6 | 9.08e-04 | M6 FDR 9.08e-04 | SCN2A;CLTC;ITFG1;PREPL;CLCN4 | exploratory Diagnosis association; strong Tier2/ML enrichment; highlighted AHBA module overlap; interpretable hubs |
| 2 | Ventral Striatum | MElightcyan1 | 586 | 6 | 0.03009 | 0.4384 | same_region_Tier2_genes | 0.00159 | M18 | 0.001771 | none below 0.05 | RAB18;MOB4;NUDT21;FAM20A;FAM169A | exploratory Diagnosis association; strong Tier2/ML enrichment; interpretable hubs |
| 3 | Ventral Striatum | MEdarkred | 227 | 6 | 0.6318 | 0.9422 | high_plus_moderate_ML_genes | 0.00638 | M31 | 1.12e-12 | none below 0.05 | SPPL2A;PTTG1IP;ADD3;MDFIC;TMX1 | strong Tier2/ML enrichment; high+moderate ML support; interpretable hubs |
| 4 | Ventral Striatum | MElightgreen | 299 | 5 | 0.01311 | 0.3278 | same_region_Tier2_genes | 0.4782 | M9 | 4.43e-12 | M6 FDR 0.0148 | GRIA2;PNMA2;SLC2A13;SLC4A10;GDAP1 | exploratory Diagnosis association; highlighted AHBA module overlap; interpretable hubs |
| 5 | Insula | MEorange | 53 | 5 | 0.03238 | 0.4695 | same_region_Tier2_genes | 0.2188 | M2 | 1.08e-04 | M2 FDR 1.08e-04 | OLFML2B;HS3ST4;GRIK3;GALNT14;PTPRU | exploratory Diagnosis association; highlighted AHBA module overlap; interpretable hubs |
| 6 | Frontal Cortex | MEdarkgreen | 479 | 5 | 0.03799 | 0.4029 | all_Tier2_genes | 1 | M26 | 4.21e-04 | M26 FDR 4.21e-04 | LOC107987081;LOC644285;LINC01004;TNFSF4;GPR18 | exploratory Diagnosis association; highlighted AHBA module overlap; interpretable hubs |
| 7 | Frontal Cortex | MEblack | 1382 | 5 | 0.04889 | 0.4029 | all_Tier2_genes | 1 | M17 | 9.24e-06 | M28 FDR 0.01034; M6 FDR 0.005843 | LOC112268125;CELF2-AS2;LOC105370588;LOC107984615;LOC112268323 | exploratory Diagnosis association; highlighted AHBA module overlap; interpretable hubs |
| 8 | Hippocampal Formation | MEgreen | 857 | 4 | 0.07365 | 0.9328 | same_region_Tier2_genes | 1.83e-07 | M30 | 1.79e-101 | none below 0.05 | USP54;PLEKHH1;C10orf90;KCNH8;ENOSF1 | strong Tier2/ML enrichment; interpretable hubs |
| 9 | Frontal Cortex | MEsalmon | 297 | 4 | 0.8538 | 0.9249 | all_Tier2_genes | 1.88e-06 | M14 | 0.2509 | none below 0.05 | RNU4-1;RNY1;SIRT4;SCARNA4;RNU4-2 | strong Tier2/ML enrichment; interpretable hubs |
| 10 | Insula | MEdarkorange | 52 | 3 | 0.02917 | 0.4695 | all_Tier2_genes | 0.9611 | M24 | 0.1206 | none below 0.05 | KDM5D;TTTY14;TTTY15;UTY;TTTY10 | exploratory Diagnosis association; interpretable hubs |

## Interpretation Notes

- Non-grey modules evaluated: 174.
- Diagnosis BH-FDR < 0.10 modules: 0.
- Modules with raw Diagnosis p < 0.05 should be treated as exploratory because the module-trait test did not survive BH-FDR.
- Tier2/ML enrichment prioritizes modules that concentrate candidate genes, but this does not make the module disease-causal or externally validated.
- AHBA highlighted overlaps use M26/M28/M2/M6/M7 as biological-context evidence from healthy-brain modules.
- MEgrey is excluded from this triage table because it is the unassigned gene bin.

## Recommended Use

Use the top-priority rows as a short list for report discussion and follow-up visualization. Keep the caution that the strongest modules are prioritized by convergent exploratory evidence rather than corrected Diagnosis association.
