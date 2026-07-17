# Frontal Subregion WGCNA Module Report FULL

Generated: 2026-07-18 01:12:30 KST

Run scope: Frontal subregion run. All three available frontal subregion matrices were processed with the top 5,000 MAD genes, the full dynamic-tree-cut parameter grid, and unadjusted plus residualized-if-available variants.

## Limited Validation vs Full Run

Existing broad Frontal Cortex WGCNA outputs were not deleted or overwritten. The broad Frontal Cortex WGCNA should be retained as exploratory/supplementary, while frontal interpretation should prioritize subregion-specific OFC BA11, dlPFC BA8/9, and Cingulate gyrus 25 modules when those modules are usable.

Frontal subregion WGCNA was added because frontal biology may be diluted when OFC, dlPFC, and cingulate regions are pooled.

The earlier broad Frontal Cortex-only limited validation result was a smoke test run with STEP6_MAX_REGIONS=1, STEP6_GENE_FILTER_TARGETS=5000, STEP6_PARAM_GRID_LIMIT=1, and STEP6_VARIANTS=unadjusted_only. That limited run checked that the workflow could execute and should not be interpreted as the final biological result.

## What was done

This analysis built WGCNA coexpression modules directly from frontal subregion normalized MDD/control GeTMM-log2 expression matrices.

The modules were not built from DEG-only genes, Tier 2-only genes, or ML-prioritized-only genes. Genes were filtered by expression quality and variability/MAD, then existing Tier 2 and ML-prioritized gene sets were tested for enrichment after module construction.

For this frontal subregion branch, module construction used the top 5,000 genes by MAD for each subregion/variant. A default 10k/15k attempt was started but was not feasible for the subregion add-on because the larger TOM steps were too slow; the final saved run therefore prioritizes all subregions with the full dynamic-tree-cut parameter grid at 5k genes.

Signed WGCNA networks were used as the main analysis. The script tried robust bicor correlation first and used Pearson only if bicor failed.

MEgrey is the WGCNA unassigned module. MEgrey rows are retained in raw output files for auditability, but MEgrey is excluded from the main biological interpretation tables.

## Sample-size and usability overview

- Regions processed: Cingulate gyrus 25, dlPFC BA8/9, OFC BA11
- Parameter-grid runs recorded: 108
- Usable parameter-grid runs: 99
- Regions with at least one usable module solution: Cingulate gyrus 25, dlPFC BA8/9, OFC BA11

A solution was marked usable when it had at least 3 non-grey modules, the largest module was below 70% of genes, grey genes were below 50%, and the region had at least 15 samples.

## Diagnosis-associated modules

No module eigengene reached BH-FDR < 0.10 for Diagnosis in the saved best solutions.

Exploratory raw-p Diagnosis modules:

- Cingulate gyrus 25 MEwhite: beta=0.149, p=0.0467, FDR=0.897
- dlPFC BA8/9 MEorange: beta=-0.107, p=0.0154, FDR=0.272
- dlPFC BA8/9 MEblue: beta=-0.0946, p=0.029, FDR=0.272
- dlPFC BA8/9 MEgreen: beta=-0.0921, p=0.0389, FDR=0.272

## Tier 2 and ML-prioritized enrichment

- Cingulate gyrus 25 MEwhite / all_Tier2_genes: overlap=5, p=0.00568, FDR=0.0904
- Cingulate gyrus 25 MEwhite / same_region_Tier2_genes: overlap=5, p=4.58e-06, FDR=0.000124
- Cingulate gyrus 25 MEdarkgrey / same_region_Tier2_genes: overlap=3, p=0.00288, FDR=0.0259
- Cingulate gyrus 25 MEdarkolivegreen / all_Tier2_genes: overlap=4, p=0.0067, FDR=0.0904
- Cingulate gyrus 25 MEdarkolivegreen / same_region_Tier2_genes: overlap=3, p=0.000614, FDR=0.00829
- dlPFC BA8/9 MEblue / all_Tier2_genes: overlap=37, p=4.81e-05, FDR=0.00101
- dlPFC BA8/9 MEblue / same_region_Tier2_genes: overlap=21, p=1.08e-11, FDR=2.27e-10
- dlPFC BA8/9 MEorange / all_Tier2_genes: overlap=5, p=0.000336, FDR=0.00352
- OFC BA11 MEtan / high_confidence_ML_genes: overlap=2, p=0.00228, FDR=0.041
- OFC BA11 MEtan / high_plus_moderate_ML_genes: overlap=3, p=0.0041, FDR=0.0738

## Published AHBA module overlap

The AHBA overlap step tests each MDD-derived module against published AHBA modules and highlights M26, M28, M2, M6, and M7.

- Cingulate gyrus 25 MEgreen vs M6: overlap=125, p=2.43e-45, FDR=5.25e-43
- Cingulate gyrus 25 MEbrown vs M2: overlap=18, p=0.0026, FDR=0.0548
- Cingulate gyrus 25 MEbrown vs M6: overlap=76, p=1.07e-08, FDR=5.45e-07
- Cingulate gyrus 25 MEbrown vs M7: overlap=18, p=0.00199, FDR=0.0452
- Cingulate gyrus 25 MEblack vs M26: overlap=71, p=1.49e-66, FDR=4.29e-64
- Cingulate gyrus 25 MEfloralwhite vs M2: overlap=7, p=0.00116, FDR=0.0278
- Cingulate gyrus 25 MEfloralwhite vs M6: overlap=16, p=0.00273, FDR=0.0562
- Cingulate gyrus 25 MEdarkred vs M28: overlap=3, p=0.00133, FDR=0.0311
- dlPFC BA8/9 MEturquoise vs M26: overlap=39, p=4.34e-06, FDR=0.000208
- dlPFC BA8/9 MEbrown vs M6: overlap=44, p=2.7e-15, FDR=2.59e-13
- dlPFC BA8/9 MEtan vs M2: overlap=5, p=0.00149, FDR=0.0323
- dlPFC BA8/9 MElightyellow vs M6: overlap=7, p=5.96e-06, FDR=0.000267
- dlPFC BA8/9 MEorange vs M2: overlap=4, p=8.06e-05, FDR=0.00285
- OFC BA11 MEyellow vs M2: overlap=13, p=0.000157, FDR=0.00648
- OFC BA11 MEyellow vs M6: overlap=39, p=3.41e-09, FDR=1.96e-07
- OFC BA11 MEyellow vs M7: overlap=13, p=0.00188, FDR=0.0433
- OFC BA11 MElightcyan vs M6: overlap=5, p=5.46e-05, FDR=0.00242
- OFC BA11 MEgreenyellow vs M6: overlap=7, p=0.00438, FDR=0.0832
- OFC BA11 MEtan vs M6: overlap=5, p=0.00513, FDR=0.0895

## Important Module Summary

- Cingulate gyrus 25 MEwhite: Exploratory diagnosis signal with candidate-gene or AHBA support. Top hubs: RNU4-1;RNU4-2;RNU5B-1;RNY1;RNY3;SNORA10;SNORA11;SNORA21;SNORA38B;SNORA68
- Cingulate gyrus 25 MEdarkolivegreen: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: LOC105379524;LOC105379574;RNU1-1;RNU1-2;RNU1-3;RNU1-4;RNVU1-18;RNVU1-2A;RNVU1-7;SEC22B4P
- Cingulate gyrus 25 MEdarkgrey: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: SCARNA18;SNORA14B;SNORA2C;SNORA48;SNORA64;SNORA8;SNORD3B-2;SNORD8;SNORD9;WDR74
- Cingulate gyrus 25 MEfloralwhite: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: BSN;CAMTA2;DLGAP2;IQSEC1;IQSEC2;LUARIS;MRTFA;SRRM4;SYT7;UNC5A
- Cingulate gyrus 25 MEdarkred: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: CELF2-AS2;DLGAP1-AS5;LINC00710;LOC101927432;LOC102724949;LOC105369667;LOC105376405;LOC105377790;LOC105378747;LOC107985137
- Cingulate gyrus 25 MEblack: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ARMC3;C7orf57;CCDC170;CFAP54;LRRIQ1;ODAD2;SPAG17;SPATA17;TEKT1;VWA3B
- Cingulate gyrus 25 MEgreen: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: CAMK2A;DLGAP1;GABBR2;HECW1;KALRN;MICAL2;NRGN;OLFM1;RIMS3;UNC13A
- Cingulate gyrus 25 MEbrown: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ATP6V1C1;CAP2;CYRIA;FAM3C;FAXC;NEGR1;PPP3CB;SCN2A;SERPINI1;TSPYL1
- Cingulate gyrus 25 MEyellowgreen: Nominal candidate-gene or AHBA context only. Top hubs: ARL17B;ATP6;ATP8;CYTB;LRRC37A;LRRC37A2;ND1;ND2;ND4L;ND5
- Cingulate gyrus 25 MEdarkgreen: Nominal candidate-gene or AHBA context only. Top hubs: DANT1;LOC105371486;LOC105371800;LOC105375613;LOC105376328;LOC105378369;LOC105378444;LOC107985035;LOC112268101;LOC286083
- Cingulate gyrus 25 MEorange: Nominal candidate-gene or AHBA context only. Top hubs: DGKK;ECEL1;GPR39;INSYN2A;LYPD1;NRP2;SHISA6;SLC5A7;STAC;ZIC3
- Cingulate gyrus 25 MEskyblue: Nominal candidate-gene or AHBA context only. Top hubs: ANO3;DCLK3;DRD1;GNAL;GPR52;GPR88;LAMP5;MUC15;STRIP2;SYNPR
- OFC BA11 MElightcyan: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: FILIP1;ISLR;KRT17;LOC102724068;LOC105374966;LOC112268195;PDZRN4;SEMA3E;THEMIS;TMEM233
- OFC BA11 MEtan: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: CLSPN;DGAT2;FBXO40;GPR6;KCNS2;MCHR2;RASL11B;RSPO2;SOWAHB;TPBG
- OFC BA11 MEyellow: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: CAMTA2;CLIP3;DLGAP4;LARP1;MAPT;PACSIN1;PRKACA;STX1B;SYN1;SYNGAP1
- OFC BA11 MEgreenyellow: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ATP6V1G2;BEX1;CMAS;GAP43;GNG3;NEDD8;NEDD8-MDP1;SNCA;TCEAL2;YWHAH
- OFC BA11 MEpurple: Nominal candidate-gene or AHBA context only. Top hubs: CELF2-AS1;HECW1-IT1;LOC105370074;LOC105371967;LOC105376405;LOC105378369;LOC107984131;LOC107985137;LOC112268101;LOC112268151
- OFC BA11 MEblue: Nominal candidate-gene or AHBA context only. Top hubs: KCNIP4-IT1;LOC105370588;LOC105372521;LOC105374516;LOC105375634;LOC105375689;LOC105377790;LOC107984007;LOC112268125;LOC112268323
- OFC BA11 MEgreen: Nominal candidate-gene or AHBA context only. Top hubs: BMPR1B;FERMT2;GPAM;METTL7A;PLPP3;PON2;PREX2;RANBP3L;SLC15A2;SLC1A2
- dlPFC BA8/9 MEorange: Exploratory diagnosis signal with candidate-gene or AHBA support. Top hubs: ARC;DNAJB5;DUSP5;DUSP6;EGR1;EGR2;GPR3;LOC102724428;NR4A3;SPRY4
- dlPFC BA8/9 MEblue: Exploratory diagnosis signal with candidate-gene or AHBA support. Top hubs: KCNIP4-IT1;LOC105370588;LOC105374516;LOC105377694;LOC105377789;LOC107984294;LOC107984629;LOC107986348;LOC112268277;PWAR1
- dlPFC BA8/9 MEgreen: Exploratory diagnosis signal with candidate-gene or AHBA support. Top hubs: CELF2-AS1;LINC00513;LOC105371077;LOC105371475;LOC105371967;LOC105375724;LOC107984131;LOC107985137;LOC112268101;LOC399975
- dlPFC BA8/9 MEbrown: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: CAP2;G3BP2;GABRA1;GPRASP2;MDH1;NAP1L2;PREPL;RANBP6;RGS4;SUCLA2
- dlPFC BA8/9 MEtan: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ATP1A1;ATP1A1-AS1;ATP1A3;DIRAS1;HCN1;KCNQ5;MAP1B;SCN8A;TESPA1;UNC13A
- dlPFC BA8/9 MEturquoise: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: CRB1;HEPH;NOTCH2;PDLIM5;PON2;PPARA;PREX2;SOX9;TP53BP2;YAP1
- dlPFC BA8/9 MElightyellow: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: CDH6;CLSTN2;DACH1;ISLR;NPFFR2;RXFP1;SEMA3E;SMYD1;THEMIS;VAT1L
- dlPFC BA8/9 MEdarkred: Nominal candidate-gene or AHBA context only. Top hubs: CES4A;LENG8;LINC00342;LINC01000;LOC105371414;NPIPA2;NPIPB13;NPIPB3;NPIPB4;ZFHX2

## Limitations

- This is exploratory and should be interpreted as module discovery, not validation.
- Some regions have modest sample sizes, so module stability may be limited.
- Residualized analyses are sensitivity analyses and may use Age/Sex if those covariates are available.
- No external validation cohort has been used yet.

## Main output files

- `frontal_subregion_wgcna_parameter_sensitivity.csv`
- `frontal_subregion_best_solution_gene_to_module_assignments.csv`
- `frontal_subregion_module_size_table.csv`
- `frontal_subregion_module_eigengenes.csv`
- `frontal_subregion_module_trait_association_results.csv`
- `frontal_subregion_Tier2_ML_enrichment.csv`
- `frontal_subregion_AHBA_overlap.csv`
- `frontal_subregion_hub_genes_top50.csv`
- `frontal_subregion_important_modules_summary.csv`
- `frontal_subregion_wgcna_manifest.json`
## Follow-up outputs added after WGCNA

- `frontal_broad_vs_subregion_comparison.csv`: compares broad Frontal Cortex with OFC BA11, dlPFC BA8/9, and Cingulate gyrus 25.
- `frontal_subregion_prioritized_modules_for_plots.csv`: prioritized non-grey modules used for follow-up visualization.
- `plots/frontal_subregion_module_summary_bubble_plot.png` and `.pdf`: module summary bubble plot.
- `plots/full_module_network_*`: full-module network plots for prioritized frontal subregion modules, plus labeled-overlap versions.
- `biological_annotation/frontal_subregion_prioritized_modules_biological_summary.csv`: biological annotation summary.
- `biological_annotation/Frontal_Subregion_WGCNA_Biological_Annotation_Report.md`: short interpretation report.
