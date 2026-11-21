> **Historical scope notice (2026-11-21):** This document preserves its original method/phase claims. Current workflow designation and completed literature work are described in [PIPELINE.md](PIPELINE.md). Tier2-RRA is the current primary exploratory ranking; corrected weighted v1 + strict-v2 results and confidence labels remain method-specific historical benchmarks. Earlier references to planned validation or canonical weighted ranking are not current project-wide status.

# MDD WGCNA Module Report FULL

Generated: 2026-07-17 16:38:32 KST

Run scope: Full default run. No STEP6 limiting environment variables were set; all available normalized region matrices, default gene filters, full parameter grid, and unadjusted plus residualized-if-available variants were attempted.

## Limited Validation vs Full Run

The earlier Frontal Cortex-only result was a smoke test run with STEP6_MAX_REGIONS=1, STEP6_GENE_FILTER_TARGETS=5000, STEP6_PARAM_GRID_LIMIT=1, and STEP6_VARIANTS=unadjusted_only. That limited run checked that the workflow could execute and should not be interpreted as the final biological result. This full report is the intended biological summary after running without those limiting controls.

## What was done

This analysis built WGCNA coexpression modules directly from the normalized MDD/control GeTMM-log2 expression matrices.

The modules were not built from DEG-only genes, Tier 2-only genes, or ML-prioritized-only genes. Genes were filtered by expression quality and variability/MAD, then existing Tier 2 and ML-prioritized gene sets were tested for enrichment after module construction.

Signed WGCNA networks were used as the main analysis. The script tried robust bicor correlation first and used Pearson only if bicor failed.

MEgrey is the WGCNA unassigned module. MEgrey rows are retained in raw output files for auditability, but MEgrey is excluded from the main biological interpretation tables.

## Sample-size and usability overview

- Regions processed: Frontal Cortex, Hippocampal Formation, Insula, Ventral Striatum
- Parameter-grid runs recorded: 432
- Usable parameter-grid runs: 396
- Regions with at least one usable module solution: Frontal Cortex, Hippocampal Formation, Insula, Ventral Striatum

A solution was marked usable when it had at least 3 non-grey modules, the largest module was below 70% of genes, grey genes were below 50%, and the region had at least 15 samples.

## Diagnosis-associated modules

No module eigengene reached BH-FDR < 0.10 for Diagnosis in the saved best solutions.

Exploratory raw-p Diagnosis modules:

- Frontal Cortex MEdarkmagenta: beta=0.0327, p=0.0433, FDR=0.403
- Frontal Cortex MEblack: beta=-0.0335, p=0.0489, FDR=0.403
- Frontal Cortex MEdarkgreen: beta=-0.0348, p=0.038, FDR=0.403
- Insula MEorange: beta=-0.0939, p=0.0324, FDR=0.469
- Insula MEdarkorange: beta=0.0104, p=0.0292, FDR=0.469
- Ventral Striatum MElavenderblush3: beta=-0.12, p=0.00744, FDR=0.328
- Ventral Striatum MElightcyan1: beta=-0.0968, p=0.0301, FDR=0.438
- Ventral Striatum MElightgreen: beta=-0.112, p=0.0131, FDR=0.328

## Tier 2 and ML-prioritized enrichment

- Frontal Cortex MEsalmon / all_Tier2_genes: overlap=25, p=7.23e-08, FDR=1.88e-06
- Hippocampal Formation MEgreen / same_region_Tier2_genes: overlap=20, p=2.65e-09, FDR=1.83e-07
- Insula MEdarkgrey / all_Tier2_genes: overlap=34, p=0.00239, FDR=0.0692
- Ventral Striatum MElightcyan1 / same_region_Tier2_genes: overlap=13, p=6.36e-05, FDR=0.00159
- Ventral Striatum MEdarkred / high_plus_moderate_ML_genes: overlap=8, p=0.000128, FDR=0.00638
- Ventral Striatum MElavenderblush3 / same_region_Tier2_genes: overlap=4, p=6.26e-05, FDR=0.00159

## Published AHBA module overlap

The AHBA overlap step tests each MDD-derived module against published AHBA modules and highlights M26, M28, M2, M6, and M7.

- Frontal Cortex MEdarkgrey vs M2: overlap=10, p=0.0015, FDR=0.0171
- Frontal Cortex MEdarkgrey vs M26: overlap=45, p=1.74e-19, FDR=1.31e-17
- Frontal Cortex MEdarkorange vs M26: overlap=52, p=2.74e-10, FDR=1.14e-08
- Frontal Cortex MEdarkorange vs M28: overlap=37, p=1.45e-05, FDR=0.000303
- Frontal Cortex MEyellow vs M26: overlap=58, p=3.83e-09, FDR=1.33e-07
- Frontal Cortex MEbrown vs M6: overlap=95, p=4.75e-22, FDR=3.95e-20
- Frontal Cortex MEbrown vs M7: overlap=20, p=0.00224, FDR=0.0242
- Frontal Cortex MEgreenyellow vs M28: overlap=28, p=5.18e-05, FDR=0.000938
- Frontal Cortex MEred vs M26: overlap=27, p=5.98e-05, FDR=0.00104
- Frontal Cortex MEblack vs M28: overlap=14, p=0.00077, FDR=0.0103
- Frontal Cortex MEblack vs M6: overlap=18, p=0.000379, FDR=0.00584
- Frontal Cortex MEblue vs M26: overlap=18, p=0.00421, FDR=0.0427
- Frontal Cortex MEblue vs M28: overlap=20, p=0.000103, FDR=0.00175
- Frontal Cortex MEtan vs M6: overlap=53, p=1.59e-13, FDR=8.81e-12
- Frontal Cortex MEsaddlebrown vs M6: overlap=10, p=3.2e-06, FDR=7.19e-05
- Frontal Cortex MEdarkgreen vs M26: overlap=11, p=2.17e-05, FDR=0.000421
- Frontal Cortex MEsienna3 vs M6: overlap=12, p=1.04e-07, FDR=2.99e-06
- Frontal Cortex MElightsteelblue1 vs M7: overlap=3, p=0.00391, FDR=0.0402
- Frontal Cortex MEfloralwhite vs M6: overlap=6, p=0.000673, FDR=0.00933
- Hippocampal Formation MEyellowgreen vs M26: overlap=12, p=2.27e-06, FDR=9.14e-05
- Hippocampal Formation MEfirebrick4 vs M2: overlap=16, p=2.88e-10, FDR=1.7e-08
- Hippocampal Formation MEmediumpurple2 vs M26: overlap=17, p=7.23e-08, FDR=3.42e-06
- Hippocampal Formation MEdarkorange2 vs M6: overlap=5, p=0.0053, FDR=0.0946
- Hippocampal Formation MEdarkgreen vs M2: overlap=41, p=2.58e-31, FDR=8.01e-29
- Hippocampal Formation MEgreenyellow vs M26: overlap=29, p=7.29e-07, FDR=3.17e-05
- Hippocampal Formation MEblack vs M28: overlap=95, p=4.57e-10, FDR=2.62e-08
- Hippocampal Formation MElightcyan vs M6: overlap=29, p=5.49e-13, FDR=4.26e-11
- Hippocampal Formation MElightcyan vs M7: overlap=12, p=5.33e-07, FDR=2.37e-05
- Hippocampal Formation MEbrown vs M26: overlap=15, p=0.00508, FDR=0.0913
- Hippocampal Formation MEbrown vs M28: overlap=19, p=3.9e-05, FDR=0.00139
- Hippocampal Formation MEdarkred vs M6: overlap=28, p=2.87e-19, FDR=3.91e-17
- Hippocampal Formation MEdarkred vs M7: overlap=5, p=0.00545, FDR=0.0963
- Hippocampal Formation MElightcyan1 vs M26: overlap=23, p=7.85e-26, FDR=1.55e-23
- Hippocampal Formation MEskyblue vs M6: overlap=21, p=7.34e-16, FDR=8.87e-14
- Hippocampal Formation MEgrey60 vs M6: overlap=5, p=0.000575, FDR=0.0145
- Hippocampal Formation MEwhite vs M6: overlap=21, p=1.91e-12, FDR=1.34e-10
- Hippocampal Formation MEfloralwhite vs M6: overlap=13, p=6.98e-11, FDR=4.22e-09
- Hippocampal Formation MEplum3 vs M6: overlap=7, p=8.97e-06, FDR=0.000337
- Hippocampal Formation MEdarkseagreen3 vs M26: overlap=5, p=0.000346, FDR=0.00942
- Hippocampal Formation MEdarkslateblue vs M2: overlap=19, p=2.14e-13, FDR=1.79e-11
- Hippocampal Formation MEsalmon4 vs M2: overlap=5, p=6.91e-06, FDR=0.000269
- Hippocampal Formation MEmagenta vs M6: overlap=23, p=0.0019, FDR=0.0414
- Hippocampal Formation MEmagenta vs M7: overlap=13, p=0.000136, FDR=0.00406
- Hippocampal Formation MEdarkviolet vs M26: overlap=10, p=0.00192, FDR=0.0414
- Hippocampal Formation MElightyellow vs M26: overlap=14, p=0.0045, FDR=0.083
- Hippocampal Formation MEdarkseagreen4 vs M7: overlap=4, p=0.00222, FDR=0.0459
- Hippocampal Formation MEthistle3 vs M7: overlap=6, p=1.12e-06, FDR=4.69e-05
- Insula MEblue vs M26: overlap=34, p=0.000111, FDR=0.00234
- Insula MEblue vs M28: overlap=35, p=1.11e-05, FDR=0.000258
- Insula MEblack vs M26: overlap=29, p=0.000148, FDR=0.00298
- Insula MEbrown vs M26: overlap=69, p=1.66e-10, FDR=6.43e-09
- Insula MEbrown vs M28: overlap=71, p=6.7e-13, FDR=3.89e-11
- Insula MEsalmon vs M28: overlap=76, p=8.27e-06, FDR=0.000197
- Insula MEmagenta vs M6: overlap=38, p=1.73e-10, FDR=6.43e-09
- Insula MEmagenta vs M7: overlap=10, p=0.00427, FDR=0.0574
- Insula MEcyan vs M6: overlap=57, p=1.26e-11, FDR=5.83e-10
- Insula MElightyellow vs M2: overlap=5, p=0.000727, FDR=0.0127
- Insula MEdarkgrey vs M26: overlap=10, p=0.00495, FDR=0.0657
- Insula MEorange vs M2: overlap=6, p=4.32e-06, FDR=0.000108
- Insula MEsaddlebrown vs M6: overlap=5, p=0.00567, FDR=0.0725
- Insula MEsaddlebrown vs M7: overlap=3, p=0.00571, FDR=0.0725
- Ventral Striatum MEantiquewhite4 vs M6: overlap=11, p=0.000401, FDR=0.0101
- Ventral Striatum MEmediumpurple2 vs M28: overlap=56, p=6.32e-06, FDR=0.000236
- Ventral Striatum MEsalmon vs M26: overlap=69, p=1.17e-67, FDR=4.59e-65
- Ventral Striatum MEdarkgreen vs M6: overlap=19, p=0.000582, FDR=0.0145
- Ventral Striatum MEdarkgreen vs M7: overlap=19, p=5.5e-10, FDR=3.59e-08
- Ventral Striatum MEbrown vs M26: overlap=19, p=0.00492, FDR=0.0908
- Ventral Striatum MEbrown vs M28: overlap=29, p=3.2e-09, FDR=2.01e-07
- Ventral Striatum MEgreen vs M6: overlap=46, p=3.73e-07, FDR=1.89e-05
- Ventral Striatum MEgreen vs M7: overlap=20, p=0.000294, FDR=0.0077
- Ventral Striatum MEdarkolivegreen4 vs M26: overlap=7, p=6.75e-07, FDR=3.11e-05
- Ventral Striatum MEdarkturquoise vs M26: overlap=15, p=0.00178, FDR=0.0383
- Ventral Striatum MElightgreen vs M6: overlap=20, p=0.000604, FDR=0.0148
- Ventral Striatum MEthistle1 vs M26: overlap=14, p=2.53e-15, FDR=2.84e-13
- Ventral Striatum MEskyblue1 vs M7: overlap=12, p=1.07e-06, FDR=4.66e-05
- Ventral Striatum MEdarkolivegreen vs M7: overlap=7, p=0.000678, FDR=0.0164
- Ventral Striatum MEskyblue2 vs M6: overlap=5, p=0.00128, FDR=0.03
- Ventral Striatum MEpalevioletred3 vs M7: overlap=6, p=1.79e-05, FDR=0.000625
- Ventral Striatum MElavenderblush3 vs M6: overlap=8, p=2.68e-05, FDR=0.000908

## Important Module Summary

- Frontal Cortex MEdarkgreen: Exploratory diagnosis signal with candidate-gene or AHBA support. Top hubs: DLG2-AS2;GPR18;LINC01004;LOC107985903;LOC107985904;LOC107986276;LOC107987081;LOC644285;LSAMP-AS1;TNFSF4
- Frontal Cortex MEdarkmagenta: Exploratory raw-p Diagnosis association only. Top hubs: RPL10;RPL10A;RPL17;RPL18;RPL26;RPL35;RPLP0;RPS10;RPS17;RPS4X
- Frontal Cortex MEblack: Exploratory diagnosis signal with candidate-gene or AHBA support. Top hubs: CELF2-AS2;EML6-AS1;LOC105370588;LOC105371068;LOC105375724;LOC107984615;LOC107985137;LOC112268125;LOC112268151;LOC112268323
- Frontal Cortex MEblue: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: LINC01000;LINC01567;LOC100190986;LOC105371800;LOC107987258;POGZ;SCUBE3;SPDYE1;TNRC17;ZNF814
- Frontal Cortex MElightsteelblue1: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: DAB1;EPB41L4B;HPCAL1;LOC105372244;MARCHF4;PCDH19;RIMBP2;RNF165;SHISA9;SORBS2
- Frontal Cortex MEgreenyellow: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ARID1A;ATN1;BAZ2A;BCL9L;CIC;MINK1;MTA1;NCOR2;SKI;TNRC18
- Frontal Cortex MEdarkgrey: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ARHGAP6;DPYSL3;GABRQ;HTR2C;LINC01833;SCN9A;SIX3;SIX3-AS1;SOX11;TMEM255A
- Frontal Cortex MEdarkorange: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: AEBP1;ATF7-NPFF;ECHDC2;FLNA;KIAA0754;MOV10;NINL;PLEKHA7;SLC12A4;ZNF462
- Frontal Cortex MEsienna3: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: B4GALNT1;CRACDL;DLGAP2;HRH2;KCNC1;MGAT5B;PKNOX2;PPFIA4;R3HDM1;RBFOX3
- Frontal Cortex MEsaddlebrown: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: HSD3BP4;LINC00622;LINC02296;LINC02669;LOC107984120;LOC107984192;OSBPL6;SLC36A1;STARD4-AS1;ZNF385B
- Frontal Cortex MEbrown: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ADGRL2;GRIN2A;KCNA2;KCNQ5;KCNS2;MIR600HG;NLK;SCN8A;STRBP;ZNF365
- Frontal Cortex MEfloralwhite: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: BEGAIN;CSNK1G2;DEAF1;HSD11B1L;KIFC2;LRP3;MAZ;MEIS3;PLCH2;TSPOAP1
- Frontal Cortex MEred: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: AP1G2;DNAH1;L3MBTL1;LENG8;NPIPA5;PKD1P4-NPIPA8;PKD1P5-LOC105376752;STAG3L3;TRIM74;ZNF692
- Frontal Cortex MEyellow: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ADGRA3;BMPR1B;IQCK;MSI2;NPAS3;PDLIM5;PPARA;STON2;YAP1;ZBTB20
- Frontal Cortex MEsalmon: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: RNU4-1;RNU4-2;RNY1;SCARNA4;SCARNA8;SIRT4;SNORA11;SNORA44;SNORA6;SNORD8
- Frontal Cortex MEtan: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ADD2;ATCAY;CYFIP2;KIF3C;MADD;PHF24;PI4KA;PRKAR1B;SLC6A17;UNC13A
- Frontal Cortex MElightcyan1: Nominal candidate-gene or AHBA context only. Top hubs: C2CD4C;LOC105371730;MED19;MIEN1;MIF;MIF-AS1;MZT2B;RPS6KA4;TGFBR3L;TTC9B
- Frontal Cortex MEivory: Nominal candidate-gene or AHBA context only. Top hubs: ATP2B4;CACNA1E;EFNB3;ENTREP2;FXYD6-FXYD2;GABRA5;RAP2B;RASGEF1C;STUM;TENM3
- Frontal Cortex MEcyan: Nominal candidate-gene or AHBA context only. Top hubs: C3;CSF1R;DOCK2;FYB1;ITGAM;LAIR1;NCKAP1L;SLCO2B1;SYK;TBXAS1
- Frontal Cortex MEyellowgreen: Nominal candidate-gene or AHBA context only. Top hubs: ANP32B;C22orf46P;HNRNPA1;HNRNPA1P10;HNRNPH3;PNN;PRKCSH;RPLP1;SMC1A;SRRT
- Frontal Cortex MEwhite: Nominal candidate-gene or AHBA context only. Top hubs: CAMTA2;DLG4;DMTN;KCNQ2;MAPK8IP3;MAPT;MARK2;PACSIN1;RUSC2;SYNPO
- Hippocampal Formation MEgreen: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: C10orf90;ENOSF1;FNBP1;KCNH8;LOC728158;PDE8A;PHLDB1;PLEKHG3;PLEKHH1;USP54
- Hippocampal Formation MEthistle3: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ABCC12;AKAP5;CHRM4;PRKCB;PTGIR;RAB15;RFC5;SNX32;SYT5;TMEM271
- Hippocampal Formation MEdarkslateblue: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ANKRD33B;B4GALNT1;CACNG8;CAMK2B;NOS1AP;PGBD5;PLPPR2;RAB3D;RTN4RL1;SMPD3
- Hippocampal Formation MEmagenta: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ADGRB2;AGAP2;AGAP2-AS1;DLG4;IQSEC2;JPH4;KCNQ2;PIK3R2;SEMA6B;SYNGAP1
- Hippocampal Formation MEdarkgreen: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: CYFIP2;GABRB3;GRIA1;HNRNPLL;KCNG3;PTPRO;RALGPS2;ST6GALNAC5;SUSD4;TOMM34
- Hippocampal Formation MEplum3: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ARL4A;CA11;FRMPD2B;IL1RAPL2;LINC00943;LOC105369302;PKIG;RB1CC1;SEC1P;STYXL2
- Hippocampal Formation MEgreenyellow: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ATP13A4;GNG12;GOLIM4;GRAMD2B;ID4;METTL7A;PLSCR4;PON2;TLR4;YAP1
- Hippocampal Formation MEyellowgreen: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: C7;FBLN1;FRZB;GPX8;MARVELD2;MGP;PRLR;SLC13A4;SYNE2;TENT5C
- Hippocampal Formation MElightcyan1: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ADGB;CFAP157;CFAP43;CFAP45;DNAH11;FAM216B;FOXJ1;MAP3K19;ODAD2;SPAG17
- Hippocampal Formation MEsalmon4: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: BCL11B;CTXND1;FIBCD1;FNDC1;FNDC1-AS1;LOC101929667;LOC105375153;LOC105379315;LOC107986697;WIPF3
- Hippocampal Formation MEmediumpurple2: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: CD44;GFAP;HSPB8;LOC105370803;MASP1;PIRT;RASL12;SPARC;STEAP3;WDR49
- Hippocampal Formation MEdarkseagreen4: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: EIF4EBP2;EMX2;EMX2OS;ETNPPL;FGD6;GPT2;MRO;NR2E1;RAB31;S1PR1
- Hippocampal Formation MEdarkorange2: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: BTBD9;CAPN5;CELF4;GPRIN1;LOC105378751;LOC112268263;NHSL2;PEA15;PPL;STIM1
- Hippocampal Formation MEdarkviolet: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: CACNA1B-AS1;CELF2-AS1;DANT2;LOC105371077;LOC105372546;LOC105377775;LOC107984317;LOC107985137;LOC112268151;NPTN-IT1
- Hippocampal Formation MElightyellow: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: EYA2;LAMB2;LOC105372850;MAPKAPK3;MLC1;PBXIP1;RAB34;RHOC;RREB1;UTRN
- Hippocampal Formation MEfirebrick4: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ACVR1;CCNJL;GHSR;IGSF10;LATS2;LINC01322;LINC02338;LOC101928196;OSR1;SHISA6
- Hippocampal Formation MEblack: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: BRD4;CC2D1A;FBRSL1;HIVEP3;KMT2D;MLXIP;NCOR2;PRR12;RAI1;TNRC18
- Hippocampal Formation MEbrown: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: DANT1;LINC01567;LOC105371271;LOC105371800;LOC105376328;LOC107984984;MAK;RGL4;SPDYE1;TFAP2E
- Hippocampal Formation MElightcyan: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ANO3;AP3S1;CHL1;DCLK1;DIRAS2;ELAVL4;HCN1;MAN1A1;RTN1;SYT16
- Hippocampal Formation MEgrey60: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: FGF12-AS1;HECW1-IT1;LOC101927182;LOC105372624;LOC105372669;LOC112268076;MAGI1-IT1;MIAT;PTPRT-AS1;SLC36A1
- Hippocampal Formation MEwhite: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: CACNB4;HLF;HTR5A-AS1;MRPS30-DT;NUAK1;PPM1L;PTPRR;RGS4;STRBP;VSNL1
- Hippocampal Formation MEdarkseagreen3: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: BMP7;CFLAR;FNDC3B;IL17RD;LOC101929555;NHSL1;NKAIN3;NKAIN3-IT1;RXRG;TRPS1
- Hippocampal Formation MEdarkred: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ARPP21;FAM163B;HIVEP1;JCAD;LOC101929256;MLIP;NGEF;R3HDM1;ST6GAL2;TESPA1
- Hippocampal Formation MEfloralwhite: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: LOC100288175;LOC107984207;RPLP0P2;SEMA4A;SLC22A10;SLC22A24;SLC22A25;SLC22A9;SYTL2;VSTM2B-DT
- Hippocampal Formation MEskyblue: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: LINC01938;LINC02143;LINC03000;LOC105377702;LRRC4C;MIR600HG;OSBPL3;PTPRT;SORCS1;TOX
- Hippocampal Formation MEhoneydew: Nominal candidate-gene or AHBA context only. Top hubs: CABP1;CAMKK2;DRP2;F12;FARP1;KALRN;KCNJ4;NAGPA;POU6F1;TMEM151B
- Hippocampal Formation MEcoral3: Nominal candidate-gene or AHBA context only. Top hubs: ABHD8;CRIP2;DVL1;GUK1;LOC112268170;PDXP;PPP1R37;TLE5;UBXN6;WSB2
- Hippocampal Formation MEthistle2: Nominal candidate-gene or AHBA context only. Top hubs: ABI3BP;C2orf88;CHST6;GALNT15;LPAR4;RFX3-DT;SPATA13;TJP2;TPPP3;USH1C
- Hippocampal Formation MEdarkolivegreen: Nominal candidate-gene or AHBA context only. Top hubs: C3orf52;GCSAM;GFRA1;GPR78;GRIK1;LPO;NDST4;NPSR1-AS1;SLC9C1;TBILA
- Hippocampal Formation MEthistle: Nominal candidate-gene or AHBA context only. Top hubs: AK1;CA10;H1-1;KLHDC8B;NR2F1-AS1;PHACTR2;SPATS2L;TINCR;UNC13C;UNC5D
- Hippocampal Formation MEblue: Nominal candidate-gene or AHBA context only. Top hubs: ATAD1;ATP5PB;CAND1;CCT6A;CCT8;DNAJA2;OPA1;RLIM;RUFY3;TM2D2
- Hippocampal Formation MElightpink4: Nominal candidate-gene or AHBA context only. Top hubs: AGBL4;LINC01773;LOC101927383;LOC101927871;LOC105378823;LOC105379170;LOC107986526;MIR3659HG;POU3F1;ST8SIA2
- Hippocampal Formation MEcoral1: Nominal candidate-gene or AHBA context only. Top hubs: BRINP3;FXYD6;FXYD6-FXYD2;KIAA1217;KLF8;LDAF1;LOC105371659;LOC105371737;NTSR1;SEC24D
- Hippocampal Formation MEantiquewhite2: Nominal candidate-gene or AHBA context only. Top hubs: BLM;CHRDL1;CSGALNACT1;CSPG5;DIO2;FAM184A;GLDC;LOC107985165;RLBP1;THSD1
- Hippocampal Formation MEplum: Nominal candidate-gene or AHBA context only. Top hubs: FGFR3;LOC100129034;NEK6;PLXNB2;PTCH1;SNTA1;TNFSF12-TNFSF13;TNS3;VEGFB;ZNF423
- Insula MEdarkorange: Exploratory raw-p Diagnosis association only. Top hubs: DDX3Y;EIF1AY;KDM5D;NLGN4Y;TTTY10;TTTY14;TTTY15;TXLNGY;USP9Y;UTY
- Insula MEorange: Exploratory diagnosis signal with candidate-gene or AHBA support. Top hubs: CCN4;GALNT14;GRIK3;HS3ST4;NR4A2;OLFML2B;PTPRU;SHB;SMYD1;TMEM233
- Insula MEbrown: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: AP1G2;CCDC180;CROCCP3;DNAH1;LAS1L;LENG8;LINC01000;MICAL1;PHYKPL;SNED1
- Insula MElightyellow: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: DGKB;EPHA6;GRIA4;HS6ST3;NPR3;PIP5K1B;RASGRF2;RELN;THSD7A;VSTM2A
- Insula MEmagenta: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: CD200;ELMOD1;GABRG2;GPR158;HCN1;NCEH1;NPTN;OXR1;PPP3CA;RTN3
- Insula MEsaddlebrown: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: CAMTA1;CD101;KCNA2;KCNQ5;LRATD1;PREP;PTGFRN;RTKN2;SCN8A;STS
- Insula MEsalmon: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: AHDC1;ATN1;KDM4B;KMT2D;MINK1;PRR12;PTPN23;SETD1B;SRCAP;TBCD
- Insula MEdarkgrey: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: KCNIP4-IT1;LOC105370588;LOC105374056;LOC107986348;LOC107986350;LOC112268277;MAPK10-AS1;ST7-AS2;TUSC7;UBE2E2
- Insula MEblue: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: CACNA1C-AS4;LGALSL-DT;LOC105371077;LOC107984131;LOC107984445;LOC107984615;LOC107985137;LOC107985285;LOC107986977;LOC107987247
- Insula MEblack: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: BMPR1B;IQCK;P2RY1;PLPP3;PON2;SLC1A2;SLC39A12;SLC7A11;TLR4;YAP1
- Insula MEcyan: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ABCA3;ATP1A3;ATP2B3;CLSTN1;CLSTN3;CNTNAP1;FAIM2;PTPRN;SEZ6L2;SLC6A17
- Insula MEdarkmagenta: Nominal candidate-gene or AHBA context only. Top hubs: ADAMTS9-AS1;C6orf141;DIAPH2;DOK6;LINC00488;LOC101929473;NPBWR1;PIK3C2G;PPP4R4;TMEM215
- Insula MEsienna3: Nominal candidate-gene or AHBA context only. Top hubs: ATP5IF1;BAD;H2AC25;LOC105374114;MZT2B;NDUFS7;NME3;RABAC1;STMN3;TMED9
- Insula MEplum1: Nominal candidate-gene or AHBA context only. Top hubs: B4GALT2;DVL1;FBXL19;GP1BB;HDAC5;LRFN4;LRP3;MAPK8IP1;NELFB;PARD6A
- Insula MEdarkolivegreen: Nominal candidate-gene or AHBA context only. Top hubs: AJAP1;C4orf51;CELF2;COPG2;COPG2IT1;GRIN2B;KCNJ6;MEST;MRTFB;TANC2
- Insula MEturquoise: Nominal candidate-gene or AHBA context only. Top hubs: AASDHPPT;DNAJA2;MATR3;MORF4L2;MRPL3;NUDT21;PRKACB;RANBP6;RNF6;SEC23A
- Insula MEmidnightblue: Nominal candidate-gene or AHBA context only. Top hubs: ATP6V0E1;BNIP2;GATM;GPM6B;HADHB;MTMR10;QKI;RFTN2;SGK3;YES1
- Ventral Striatum MElavenderblush3: Exploratory diagnosis signal with candidate-gene or AHBA support. Top hubs: CACNA2D1;CAP2;CLCN4;CLTC;DMXL2;GLRB;ITFG1;ITFG1-AS1;PREPL;SCN2A
- Ventral Striatum MElightgreen: Exploratory diagnosis signal with candidate-gene or AHBA support. Top hubs: GDAP1;GPRASP2;GRIA2;KLHL7;LOC112268471;PNMA2;PRKACB;RLIM;SLC2A13;SLC4A10
- Ventral Striatum MElightcyan1: Exploratory diagnosis signal with candidate-gene or AHBA support. Top hubs: FAM169A;FAM20A;HACD3;HSPE1-MOB4;KLHL11;MOB4;NUDT21;PRKAR1A;RAB18;VMA21
- Ventral Striatum MEskyblue2: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ARPP19;B3GNT2;CDC40;CYRIB;MCTP1;PDP1;SACS;SEMA3A;SLMAP;TOMM70
- Ventral Striatum MEsalmon: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ADGB;CAPSL;CFAP107;CFAP45;CFAP52;LINC01765;MAP3K19;MYLK3;ODAD2;TEKT1
- Ventral Striatum MEantiquewhite4: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: CAMK2D;CDH10;GNG2;GRK3;KCNH5;KCNJ6;MYO16;SEZ6L;TAFA1;TAFA2
- Ventral Striatum MEbrown: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: BTBD9-AS1;DANT2;IGFBP-AS1;LOC105371068;LOC105373873;LOC105373874;LOC107984615;LOC107984805;LOC107985363;XKR6
- Ventral Striatum MEmediumpurple2: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: CREBBP;GRIK5;HIC2;KMT2D;POM121;PRR12;PRRC2A;SETD1B;SRCAP;WNK2
- Ventral Striatum MEdarkolivegreen: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: CALN1;CAMK4;GRM5-AS1;KCNQ5;LMO7;LRRC7;MPPE1;OSBPL8;RNF150;SGIP1
- Ventral Striatum MEpalevioletred3: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: BAP1;C20orf27;CCM2;CGREF1;FOXH1;KIFC2;MAP3K10;RUNDC3A;ST6GALNAC6;TMEM179
- Ventral Striatum MEdarkturquoise: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: BMPR1B;CPE;DDAH1;GABRG1;GLUD1;HSDL2;PDLIM5;SLC15A2;SLC1A3;TTPA
- Ventral Striatum MEthistle1: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: AK7;CCDC146;CFAP47;CFAP54;CFAP58;DCDC1;DNAI7;TEX9;ZBBX;ZMYND10
- Ventral Striatum MEdarkolivegreen4: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ARHGAP18;COL21A1;FABP5;FAM181A-AS1;FAM47E-STBD1;FRMD3;KLHDC8B;PLAAT5;SPATS2L;TC2N
- Ventral Striatum MEdarkred: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ADD3;ANO6;ERLIN2;GOLIM4;MDFIC;PTTG1IP;RASSF8;REST;SPPL2A;TMX1
- Ventral Striatum MEdarkgreen: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: CNR1;DLC1;GPR88;ITPR1;KCNAB1;MME;PDE10A;PRKCB;RGS7BP;TUNAR
- Ventral Striatum MEskyblue1: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: ATP2A2;ATP2B2;ATP2B3;GNB5;KCNB1;MPRIP;PPP2R2C;SLC12A5;SORBS2;TANC2
- Ventral Striatum MEgreen: Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association. Top hubs: C2CD2L;EPB41L1;FBXL16;FBXO41;IQSEC2;MINK1;NLGN2;PIK3R2;STX1B;TMEM8B
- Ventral Striatum MEsteelblue: Nominal candidate-gene or AHBA context only. Top hubs: CACNA2D1-AS1;DLG2-AS2;LOC101927647;LOC105374152;LOC105374162;LOC105374310;LOC107984425;LOC112268009;TAS2R31;TMEM14EP
- Ventral Striatum MEbrown2: Nominal candidate-gene or AHBA context only. Top hubs: ADGRA1;CHD5;GOLGA7B;KIF21B;MAST1;PCDHGC4;PLXNA1;SHANK1;SLIT1;SRGAP3
- Ventral Striatum MEorangered3: Nominal candidate-gene or AHBA context only. Top hubs: ALG1L9P;HDAC10;LOC105372249;MAPK12;MEG9;MOV10;NINL;PLEKHG5;TSNAXIP1;VILL
- Ventral Striatum MEskyblue3: Nominal candidate-gene or AHBA context only. Top hubs: AFAP1;C4orf51;CADM3;CADM3-AS1;FAM78B;LINC01123;SLC24A3;SLITRK5;SPOCK2;USP11
- Ventral Striatum MEcoral1: Nominal candidate-gene or AHBA context only. Top hubs: DAB2IP;EPB41;LRTM2;NOB1;NRGN;PPP6R1;PTPN7;RHOBTB2;RYR3;TMOD1
- Ventral Striatum MEpurple: Nominal candidate-gene or AHBA context only. Top hubs: CD53;CYBB;DOCK2;FYB1;IKZF1;NCKAP1L;SYK;TFEC;TLR1;TLR6
- Ventral Striatum MEpaleturquoise: Nominal candidate-gene or AHBA context only. Top hubs: CNTNAP3;ERICH3-AS1;LOC101929617;LOC105369672;LOC105369673;LOC105374681;LOC107985664;LOC285762;RALGPS2-AS1;SYBU-AS1
- Ventral Striatum MElightpink4: Nominal candidate-gene or AHBA context only. Top hubs: ITGA7;LOC105379524;MTSS2;RNU1-1;RNU1-2;RNU1-3;RNU1-4;RNVU1-7;SEC22B4P;ZNF423

## Limitations

- This is exploratory and should be interpreted as module discovery, not validation.
- Some regions have modest sample sizes, so module stability may be limited.
- Residualized analyses are sensitivity analyses and may use Age/Sex if those covariates are available.
- No external validation cohort has been used yet.

## Main output files

- `mdd_wgcna_parameter_sensitivity.csv`
- `best_solution_gene_to_module_assignments.csv`
- `best_solution_module_size_table.csv`
- `best_solution_module_eigengenes.csv`
- `module_trait_association_results.csv`
- `module_enrichment_for_Tier2_and_ML_genes.csv`
- `mdd_modules_vs_published_AHBA_modules_enrichment.csv`
- `module_hub_genes_top50.csv`
- `important_modules_summary.csv`
- `mdd_wgcna_manifest.json`
