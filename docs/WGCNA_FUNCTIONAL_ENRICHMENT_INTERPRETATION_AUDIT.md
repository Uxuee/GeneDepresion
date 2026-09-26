# WGCNA Functional-Enrichment and Biological-Interpretation Audit

## Executive summary

Phase 3C audited the fixed Phase-3B enrichment output without rerunning WGCNA or enrichment. The ten prespecified modules contained 200 terms with BH q<0.05 and 51 additional terms with 0.05≤q<0.10. Functional annotation remains separate from Diagnosis association: zero modules have Diagnosis BH q<0.10.

## Audit scope

Git baseline: `dd59b1dc2b64d681ec9e7b8f0092b41b005471ef`. Databases were restricted to GO BP/MF/CC, Reactome, and MSigDB Hallmark. The audit added gene tracing, redundancy grouping, family-concentration flags, and publication wording guardrails; it did not calculate new enrichment p-values.

## Prioritized modules

| region | module | module_genes_total | module_genes_mapped | mapping_fraction | Diagnosis_beta | Diagnosis_raw_p | Diagnosis_BH_q | structural_core_gene_count | Tier2_gene_count | corrected_ML_high_gene_count | corrected_ML_high_plus_moderate_gene_count |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Ventral Striatum | MElavenderblush3 | 35 | 34 | 0.9714285714285714 | -0.1202869351796773 | 0.0074412142761905 | 0.3278102226356968 | 35 | 4 | 0 | 0 |
| Ventral Striatum | MEdarkred | 227 | 217 | 0.9559471365638766 | 0.0218305133174333 | 0.6318333095761178 | 0.942208467476964 | 50 | 13 | 0 | 6 |
| Ventral Striatum | MElightgreen | 299 | 296 | 0.9899665551839464 | -0.1122054088716714 | 0.0131124089054278 | 0.3278102226356968 | 50 | 8 | 0 | 0 |
| Ventral Striatum | MElightcyan1 | 586 | 576 | 0.9829351535836176 | -0.0967591970949214 | 0.0300868059667943 | 0.4383577572829973 | 50 | 27 | 1 | 3 |
| Insula | MEorange | 53 | 53 | 1.0 | -0.093854516508098 | 0.0323771955035717 | 0.469469334801791 | 44 | 2 | 0 | 1 |
| Hippocampal Formation | MEgreen | 857 | 796 | 0.9288214702450408 | 0.0829309270459142 | 0.0736470402599981 | 0.9327758355003674 | 50 | 35 | 0 | 3 |
| dlPFC BA8/9 | MEblue | 703 | 611 | 0.8691322901849218 | -0.0945655364898348 | 0.0289793135594273 | 0.2723941494602633 | 50 | 37 | 0 | 2 |
| dlPFC BA8/9 | MEorange | 23 | 23 | 1.0 | -0.1068614219237258 | 0.015407352197724 | 0.2723941494602633 | 22 | 5 | 0 | 0 |
| Cingulate gyrus 25 | MEwhite | 50 | 50 | 1.0 | 0.1486780438879465 | 0.0467459884076952 | 0.8967327691588728 | 47 | 5 | 0 | 1 |
| OFC BA11 | MEtan | 41 | 39 | 0.951219512195122 | 0.0540101494727599 | 0.2355236279078876 | 0.5494934282476818 | 11 | 5 | 2 | 3 |

## Mapping quality

Phase-3B reported 2,611/2,783 unique mapped identifiers (93.8%) across the prioritized collection. Per-module accounting is retained below; unmapped identifiers are listed rather than discarded.

| region | module | module_genes_total | module_genes_mapped | module_genes_unmapped | mapping_percent | mapping_bias_caution |
|---|---|---|---|---|---|---|
| Ventral Striatum | MElavenderblush3 | 35 | 34 | 1 | 97.14285714285714 | At least 90% mapped; residual unmapped identifiers remain listed. |
| Ventral Striatum | MEdarkred | 227 | 217 | 10 | 95.59471365638767 | At least 90% mapped; residual unmapped identifiers remain listed. |
| Ventral Striatum | MElightgreen | 299 | 296 | 3 | 98.99665551839463 | At least 90% mapped; residual unmapped identifiers remain listed. |
| Ventral Striatum | MElightcyan1 | 586 | 576 | 10 | 98.29351535836176 | At least 90% mapped; residual unmapped identifiers remain listed. |
| Insula | MEorange | 53 | 53 | 0 | 100.0 | At least 90% mapped; residual unmapped identifiers remain listed. |
| Hippocampal Formation | MEgreen | 857 | 796 | 61 | 92.88214702450408 | At least 90% mapped; residual unmapped identifiers remain listed. |
| dlPFC BA8/9 | MEblue | 703 | 611 | 92 | 86.91322901849217 | Below 90% mapping; interpretation may underrepresent unmapped identifiers. |
| dlPFC BA8/9 | MEorange | 23 | 23 | 0 | 100.0 | At least 90% mapped; residual unmapped identifiers remain listed. |
| Cingulate gyrus 25 | MEwhite | 50 | 50 | 0 | 100.0 | At least 90% mapped; residual unmapped identifiers remain listed. |
| OFC BA11 | MEtan | 41 | 39 | 2 | 95.1219512195122 | At least 90% mapped; residual unmapped identifiers remain listed. |

## FDR-supported functional terms

There are 200 BH-q<0.05 terms. A separate file contains the 51 suggestive terms; the categories are never merged. Every row retains exact overlap genes, module and background sizes, raw p, BH q, and ratios.

## Term redundancy

The 251 q<0.10 terms collapse into 71 overlap/hierarchy clusters using exact overlap, overlap-gene Jaccard ≥0.80, or an immediate GO parent-child relation. Cluster counts are descriptive; related terms are not treated as independent findings.

## Gene drivers

| region | module | gene_symbol | fdr05_term_count | q_lt_0_10_term_count | gene_family |
|---|---|---|---|---|---|
| dlPFC BA8/9 | MEorange | FOS | 73 | 83 | other |
| dlPFC BA8/9 | MEorange | EGR1 | 59 | 71 | other |
| dlPFC BA8/9 | MEorange | EGR2 | 54 | 58 | other |
| dlPFC BA8/9 | MEorange | JUNB | 52 | 57 | other |
| dlPFC BA8/9 | MEorange | EGR3 | 52 | 54 | other |
| dlPFC BA8/9 | MEorange | NR4A3 | 51 | 58 | other |
| dlPFC BA8/9 | MEorange | DUSP4 | 46 | 50 | other |
| dlPFC BA8/9 | MEorange | DUSP6 | 44 | 49 | other |
| dlPFC BA8/9 | MEorange | PER2 | 40 | 46 | other |
| dlPFC BA8/9 | MEorange | BHLHE40 | 38 | 40 | other |
| dlPFC BA8/9 | MEorange | ETV5 | 38 | 38 | other |
| dlPFC BA8/9 | MEorange | DUSP1 | 37 | 43 | other |
| dlPFC BA8/9 | MEorange | CSRNP1 | 33 | 33 | other |
| dlPFC BA8/9 | MEorange | CCN1 | 29 | 34 | other |
| dlPFC BA8/9 | MEorange | DUSP5 | 28 | 32 | other |
| dlPFC BA8/9 | MEblue | TAS2R43 | 27 | 27 | taste_receptor |
| dlPFC BA8/9 | MEblue | TAS2R31 | 27 | 27 | taste_receptor |
| dlPFC BA8/9 | MEblue | TAS2R46 | 27 | 27 | taste_receptor |
| dlPFC BA8/9 | MEblue | TAS2R50 | 27 | 27 | taste_receptor |
| dlPFC BA8/9 | MEblue | TAS2R3 | 27 | 27 | taste_receptor |

## Gene-family effects

Recognized family concentration was flagged for taste/olfactory receptors, protocadherins, histones, ribosomal genes, immunoglobulins, and small nuclear/nucleolar RNA genes. A family flag does not declare an artifact; it limits how literally broad ontology labels should be interpreted.

| region | old_module | term | BH_q | dominant_gene_family | dominant_family_fraction | overlap_gene_symbols |
|---|---|---|---|---|---|---|
| dlPFC BA8/9 | MEblue | sensory perception of taste | 3.1029260742217435e-10 | taste_receptor | 1.0 | TAS2R10;TAS2R13;TAS2R14;TAS2R19;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7;TAS2R9 |
| dlPFC BA8/9 | MEblue | detection of chemical stimulus involved in sensory perception of bitter taste | 1.3124629813900428e-09 | taste_receptor | 1.0 | TAS2R10;TAS2R13;TAS2R14;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7;TAS2R9 |
| dlPFC BA8/9 | MEblue | detection of chemical stimulus involved in sensory perception of taste | 1.3124629813900428e-09 | taste_receptor | 1.0 | TAS2R10;TAS2R13;TAS2R14;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7;TAS2R9 |
| dlPFC BA8/9 | MEblue | sensory perception of bitter taste | 1.3124629813900428e-09 | taste_receptor | 1.0 | TAS2R10;TAS2R13;TAS2R14;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7;TAS2R9 |
| dlPFC BA8/9 | MEblue | sensory perception of chemical stimulus | 4.447891204412151e-08 | taste_receptor | 0.75 | OR2D2;OR2D3;OR6A2;OR6X1;TAS2R10;TAS2R13;TAS2R14;TAS2R19;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7;TAS2R9 |
| dlPFC BA8/9 | MEblue | detection of chemical stimulus involved in sensory perception | 1.07478107918881e-07 | taste_receptor | 0.7333333333333333 | OR2D2;OR2D3;OR6A2;OR6X1;TAS2R10;TAS2R13;TAS2R14;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7;TAS2R9 |
| dlPFC BA8/9 | MEblue | detection of chemical stimulus | 8.597608054260924e-07 | taste_receptor | 0.7333333333333333 | OR2D2;OR2D3;OR6A2;OR6X1;TAS2R10;TAS2R13;TAS2R14;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7;TAS2R9 |
| dlPFC BA8/9 | MEblue | detection of stimulus involved in sensory perception | 4.346225564394925e-06 | taste_receptor | 0.7333333333333333 | OR2D2;OR2D3;OR6A2;OR6X1;TAS2R10;TAS2R13;TAS2R14;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7;TAS2R9 |
| dlPFC BA8/9 | MEblue | taste receptor activity | 2.8590144217881757e-12 | taste_receptor | 1.0 | TAS2R10;TAS2R13;TAS2R14;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7;TAS2R9 |
| dlPFC BA8/9 | MEblue | bitter taste receptor activity | 2.8590144217881757e-12 | taste_receptor | 1.0 | TAS2R10;TAS2R13;TAS2R14;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7;TAS2R9 |
| dlPFC BA8/9 | MEblue | Class C/3 (Metabotropic glutamate/pheromone receptors) | 1.2575869494340457e-13 | taste_receptor | 1.0 | TAS2R10;TAS2R13;TAS2R14;TAS2R19;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7;TAS2R9 |
| dlPFC BA8/9 | MEblue | Sensory perception of sweet, bitter, and umami (glutamate) taste | 1.5762669231543378e-10 | taste_receptor | 1.0 | TAS2R10;TAS2R13;TAS2R14;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7 |
| dlPFC BA8/9 | MEblue | Sensory perception of taste | 3.560161284924573e-10 | taste_receptor | 1.0 | TAS2R10;TAS2R13;TAS2R14;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7 |
| dlPFC BA8/9 | MEblue | G alpha (i) signalling events | 1.2498240349688793e-06 | taste_receptor | 0.8 | GPR18;P2RY14;RRH;TAS2R10;TAS2R13;TAS2R14;TAS2R19;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7;TAS2R9 |
| dlPFC BA8/9 | MEblue | Sensory Perception | 8.393005843255522e-06 | taste_receptor | 0.7142857142857143 | OR2D2;OR2D3;OR6A2;OR6X1;TAS2R10;TAS2R13;TAS2R14;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7 |
| Cingulate gyrus 25 | MEwhite | RNA processing | 2.68259748164069e-46 | small_nuclear_nucleolar_RNA | 0.975609756097561 | MIR1843;RNU4-1;RNU5A-1;RNU5B-1;RNU6-7;RNU6-9;RNU6ATAC;SCARNA3;SCARNA4;SCARNA8;SNORA1;SNORA10;SNORA10B;SNORA11;SNORA17A;SNORA21;SNORA33;SNORA38B;SNORA44;SNORA54;SNORA58;SNORA6;SNORA63;SNORA63E;SNORA65;SNORA67;SNORA68;SNORA68B;SNORA71B;SNORA71E;SNORA73A;SNORA75;SNORA77B;SNORA80E;SNORA9;SNORD10;SNORD17;SNORD22;SNORD67;SNORD89;SNORD90 |
| Cingulate gyrus 25 | MEwhite | spliceosomal tri-snRNP complex assembly | 8.93018121508832e-08 | small_nuclear_nucleolar_RNA | 1.0 | RNU4-1;RNU5A-1;RNU5B-1;RNU6-7;RNU6-9;RNU6ATAC |
| Cingulate gyrus 25 | MEwhite | spliceosomal snRNP assembly | 8.93018121508832e-08 | small_nuclear_nucleolar_RNA | 1.0 | RNU4-1;RNU5A-1;RNU5B-1;RNU6-7;RNU6-9;RNU6ATAC |
| Cingulate gyrus 25 | MEwhite | protein-RNA complex assembly | 4.754013049358493e-05 | small_nuclear_nucleolar_RNA | 1.0 | RNU4-1;RNU5A-1;RNU5B-1;RNU6-7;RNU6-9;RNU6ATAC |
| Cingulate gyrus 25 | MEwhite | protein-RNA complex organization | 4.887649948102264e-05 | small_nuclear_nucleolar_RNA | 1.0 | RNU4-1;RNU5A-1;RNU5B-1;RNU6-7;RNU6-9;RNU6ATAC |
| Cingulate gyrus 25 | MEwhite | ribonucleoprotein complex biogenesis | 0.0002047641349264609 | small_nuclear_nucleolar_RNA | 1.0 | RNU4-1;RNU5A-1;RNU5B-1;RNU6-7;RNU6-9;RNU6ATAC |
| Cingulate gyrus 25 | MEwhite | RNA splicing, via transesterification reactions | 0.0002047641349264609 | small_nuclear_nucleolar_RNA | 1.0 | RNU4-1;RNU5A-1;RNU5B-1;RNU6-7;RNU6-9;RNU6ATAC |
| Cingulate gyrus 25 | MEwhite | RNA splicing, via transesterification reactions with bulged adenosine as nucleophile | 0.0002047641349264609 | small_nuclear_nucleolar_RNA | 1.0 | RNU4-1;RNU5A-1;RNU5B-1;RNU6-7;RNU6-9;RNU6ATAC |
| Cingulate gyrus 25 | MEwhite | mRNA splicing, via spliceosome | 0.0002047641349264609 | small_nuclear_nucleolar_RNA | 1.0 | RNU4-1;RNU5A-1;RNU5B-1;RNU6-7;RNU6-9;RNU6ATAC |
| Cingulate gyrus 25 | MEwhite | mRNA processing | 0.0006948053611969317 | small_nuclear_nucleolar_RNA | 1.0 | RNU4-1;RNU5A-1;RNU5B-1;RNU6-7;RNU6-9;RNU6ATAC |
| Cingulate gyrus 25 | MEwhite | RNA splicing | 0.0011633705983651 | small_nuclear_nucleolar_RNA | 1.0 | RNU4-1;RNU5A-1;RNU5B-1;RNU6-7;RNU6-9;RNU6ATAC |
| Cingulate gyrus 25 | MEwhite | mRNA metabolic process | 0.0050712076267792 | small_nuclear_nucleolar_RNA | 1.0 | RNU4-1;RNU5A-1;RNU5B-1;RNU6-7;RNU6-9;RNU6ATAC |
| Cingulate gyrus 25 | MEwhite | RNA binding | 0.0117371626125935 | small_nuclear_nucleolar_RNA | 1.0 | RNU4-1;RNU6-7;RNU6-9;RNU6ATAC |
| Cingulate gyrus 25 | MEwhite | nucleolus | 1.050207635653775e-33 | small_nuclear_nucleolar_RNA | 0.9714285714285714 | MIR1843;SCARNA3;SCARNA4;SCARNA8;SNORA1;SNORA10;SNORA10B;SNORA11;SNORA17A;SNORA21;SNORA33;SNORA38B;SNORA44;SNORA54;SNORA58;SNORA6;SNORA63;SNORA63E;SNORA65;SNORA67;SNORA68;SNORA68B;SNORA71B;SNORA71E;SNORA73A;SNORA75;SNORA77B;SNORA80E;SNORA9;SNORD10;SNORD17;SNORD22;SNORD67;SNORD89;SNORD90 |
| Cingulate gyrus 25 | MEwhite | Cajal body | 1.432754430427883e-06 | small_nuclear_nucleolar_RNA | 0.8 | MIR1843;SCARNA3;SCARNA4;SCARNA8;SNORA38B |

## Cross-database convergence

| region | module | dominant_theme | databases_with_q_lt_0_10 | cross_database_class |
|---|---|---|---|---|
| Ventral Striatum | MElavenderblush3 | SYNAPTIC_PRESYNAPTIC | GO;MSigDB Hallmark;Reactome | MULTI_DATABASE_CONVERGENCE |
| Ventral Striatum | MEdarkred | EXTRACELLULAR_MATRIX_WOUND | GO | WEAK/REDUNDANT_SUPPORT |
| Ventral Striatum | MElightgreen | SYNAPTIC_PRESYNAPTIC | GO | SINGLE_DATABASE_SUPPORT |
| Ventral Striatum | MElightcyan1 | CHOLESTEROL_LIPID | MSigDB Hallmark;Reactome | MULTI_DATABASE_CONVERGENCE |
| Insula | MEorange | NO_CLEAR_THEME |  | NO_CLEAR_THEME |
| Hippocampal Formation | MEgreen | ACTIN_CELL_CORTEX | GO | SINGLE_DATABASE_SUPPORT |
| dlPFC BA8/9 | MEblue | RECEPTOR_GPCR_FAMILY | GO;Reactome | MULTI_DATABASE_CONVERGENCE |
| dlPFC BA8/9 | MEorange | SIGNALING_IMMEDIATE_EARLY | GO;MSigDB Hallmark;Reactome | MULTI_DATABASE_CONVERGENCE |
| Cingulate gyrus 25 | MEwhite | RNA_PROCESSING_NUCLEOLAR | GO;Reactome | MULTI_DATABASE_CONVERGENCE |
| OFC BA11 | MEtan | ANGIOGENESIS_VASCULAR | MSigDB Hallmark | WEAK/REDUNDANT_SUPPORT |

## ML-gene contribution to functional themes

The table `ml_gene_function_overlap.csv` reports both global-union and same-region corrected-ML status for unique functional-driver genes. The global columns match the gene universe used by the fixed Phase-3B module-over-representation analysis; same-region columns prevent cross-region priority from being mistaken for local evidence. ML overlap and pathway annotation are overlapping evidence layers, not independent replications.

## Hub/core contribution

`hub_function_overlap.csv` gives median abs(kME), top-50-hub fraction, and structural-core fraction for each q<0.10 term. Missing kME values remain missing; no inferential hub-enrichment test was introduced.

## Module-by-module interpretation

### Ventral Striatum MElavenderblush3

- Annotation: **SYNAPTIC_PRESYNAPTIC**; category: **MAIN_TEXT_CANDIDATE**.
- Diagnosis: BH q=0.328; not BH-FDR significant.
- Corrected ML overlap: high=0; high+moderate=0; high+moderate q=1.
- Exact corrected-ML genes in module: high [none]; high+moderate [none].
- Leading fixed terms: GO: low-density lipoprotein particle receptor binding (q=0.00138; n=3; genes=CLTC;HSP90B1;SYT1); GO: lipoprotein particle receptor binding (q=0.00138; n=3; genes=CLTC;HSP90B1;SYT1); GO: dense core granule (q=0.0281; n=3; genes=CACNA2D1;DMXL2;SYT1); MSigDB Hallmark: HALLMARK_PI3K_AKT_MTOR_SIGNALING (q=0.0339; n=3; genes=CLTC;HSP90B1;PRKAA2); Reactome: Transmission across Chemical Synapses (q=0.0433; n=5; genes=CACNA2D1;GLRB;NSF;PRKAA2;SYT1)
- Main caution: Module has no more than 50 mapped genes; term dependence is substantial.
- Safe wording: The Ventral Striatum MElavenderblush3 gene set was annotated to synaptic presynaptic terms in the fixed Phase-3B enrichment analysis.
- Wording to avoid: Avoid claiming that synaptic presynaptic is dysregulated, causal, replicated, or a validated MDD mechanism in Ventral Striatum.

### Ventral Striatum MEdarkred

- Annotation: **EXTRACELLULAR_MATRIX_WOUND**; category: **SUPPLEMENTARY_CONTEXT**.
- Diagnosis: BH q=0.942; not BH-FDR significant.
- Corrected ML overlap: high=0; high+moderate=6; high+moderate q=0.109.
- Exact corrected-ML genes in module: high [none]; high+moderate [ADH1B;ADIRF;ARHGAP9;DPP4;SMIM10;TCN2].
- Leading fixed terms: GO: regulation of extracellular matrix organization (q=0.0911; n=6; genes=DDR2;DPP4;HPSE;LAMB2;TGFB2;TNFRSF1B); GO: positive regulation of response to wounding (q=0.0911; n=6; genes=ANO6;DDR2;EMILIN2;F2R;HPSE;SCARF1); GO: regulation of extracellular matrix disassembly (q=0.0911; n=4; genes=DDR2;DPP4;HPSE;TGFB2)
- Main caution: Support is confined to one annotation database.
- Safe wording: The Ventral Striatum MEdarkred gene set was annotated to extracellular matrix wound terms in the fixed Phase-3B enrichment analysis.
- Wording to avoid: Avoid claiming that extracellular matrix wound is dysregulated, causal, replicated, or a validated MDD mechanism in Ventral Striatum.

### Ventral Striatum MElightgreen

- Annotation: **SYNAPTIC_PRESYNAPTIC**; category: **MAIN_TEXT_CANDIDATE**.
- Diagnosis: BH q=0.328; not BH-FDR significant.
- Corrected ML overlap: high=0; high+moderate=0; high+moderate q=1.
- Exact corrected-ML genes in module: high [none]; high+moderate [none].
- Leading fixed terms: GO: synaptic vesicle membrane (q=0.000664; n=13; genes=ATP6AP2;ATP6V1A;ATP6V1B2;ATP6V1C1;ATP8A1;OPRK1;SCAMP1;SLC35F1;SVOP;SYNGR3;SYNPR;SYP;SYT4); GO: exocytic vesicle membrane (q=0.000664; n=13; genes=ATP6AP2;ATP6V1A;ATP6V1B2;ATP6V1C1;ATP8A1;OPRK1;SCAMP1;SLC35F1;SVOP;SYNGR3;SYNPR;SYP;SYT4); GO: presynapse (q=0.0101; n=26; genes=AP3S1;ATP6AP2;ATP6V1A;ATP6V1B2;ATP6V1C1;ATP8A1;CADPS;CNTN1;FZD3;GAD1;GAD2;GAP43;NAPB;NECAP1;OPRK1;RAB3C;SCAMP1;SLC35F1;SLC4A10;SLC8A1;SVOP;SYNGR3;SYNPR;SYP;SYT4;YWHAH); GO: exocytic vesicle (q=0.0101; n=15; genes=ATP6AP2;ATP6V1A;ATP6V1B2;ATP6V1C1;ATP8A1;OPRK1;RAB3C;SCAMP1;SLC35F1;SVOP;SYNGR3;SYNPR;SYP;SYT4;SYTL5); GO: dendrite (q=0.0101; n=24; genes=AKAP5;ASIC2;ATP6AP2;CHL1;EPHA7;FZD3;GABRA1;GABRA4;GABRG2;GAP43;GRIA2;HTR2C;INPP5F;KCNIP1;KIFAP3;MAP2K1;NEGR1;NSG1;OPRK1;RTN1;SLC4A10;SLC8A1;SYT4;UHMK1)
- Main caution: Support is confined to one annotation database.
- Safe wording: The Ventral Striatum MElightgreen gene set was annotated to synaptic presynaptic terms in the fixed Phase-3B enrichment analysis.
- Wording to avoid: Avoid claiming that synaptic presynaptic is dysregulated, causal, replicated, or a validated MDD mechanism in Ventral Striatum.

### Ventral Striatum MElightcyan1

- Annotation: **CHOLESTEROL_LIPID**; category: **MAIN_TEXT_CANDIDATE**.
- Diagnosis: BH q=0.438; not BH-FDR significant.
- Corrected ML overlap: high=1; high+moderate=3; high+moderate q=1.
- Exact corrected-ML genes in module: high [LINC01361]; high+moderate [FSBP;LINC01361;PPARG].
- Leading fixed terms: MSigDB Hallmark: HALLMARK_CHOLESTEROL_HOMEOSTASIS (q=0.0219; n=9; genes=ACAT2;CYP51A1;GLDC;HMGCS1;IDI1;PNRC1;PPARG;SQLE;TMEM97); Reactome: Cholesterol biosynthesis (q=0.0568; n=6; genes=ACAT2;CYP51A1;HMGCS1;IDI1;MSMO1;SQLE); Reactome: Protein ubiquitination (q=0.0568; n=7; genes=BCL10;H2BC11;H2BC3;H2BC5;HLA-A;UBE2D3;UBE2W)
- Main caution: Interpret as module annotation only; it is not MDD-association evidence.
- Safe wording: The Ventral Striatum MElightcyan1 gene set was annotated to cholesterol lipid terms in the fixed Phase-3B enrichment analysis.
- Wording to avoid: Avoid claiming that cholesterol lipid is dysregulated, causal, replicated, or a validated MDD mechanism in Ventral Striatum.

### Insula MEorange

- Annotation: **NO_CLEAR_THEME**; category: **AVOID_INTERPRETATION**.
- Diagnosis: BH q=0.469; not BH-FDR significant.
- Corrected ML overlap: high=0; high+moderate=1; high+moderate q=1.
- Exact corrected-ML genes in module: high [none]; high+moderate [DPP4].
- Leading fixed terms: No term reached BH q<0.10.
- Main caution: No Phase-3B term reached BH q<0.10.
- Safe wording: The Insula MEorange gene set was annotated to no clear theme terms in the fixed Phase-3B enrichment analysis.
- Wording to avoid: Avoid claiming that no clear theme is dysregulated, causal, replicated, or a validated MDD mechanism in Insula.

### Hippocampal Formation MEgreen

- Annotation: **ACTIN_CELL_CORTEX**; category: **MAIN_TEXT_CANDIDATE**.
- Diagnosis: BH q=0.933; not BH-FDR significant.
- Corrected ML overlap: high=0; high+moderate=3; high+moderate q=1.
- Exact corrected-ML genes in module: high [none]; high+moderate [GAU1;LOC105376380;TRPV3].
- Leading fixed terms: GO: cell cortex (q=0.00362; n=23; genes=CDH1;CNKSR1;CTTNBP2;EPB41L2;FGF1;FNBP1;GIPC1;GPSM2;HIP1;HIP1R;IQGAP1;LLGL2;MYO1D;MYO6;MYO9B;PHLDB1;PIEZO2;PLA2G4C;RAPGEF3;RTKN;SEPTIN1;SEPTIN12;SEPTIN4); GO: actin binding (q=0.0131; n=27; genes=AIF1L;AMOTL2;BIN1;CCDC88A;CLMN;CTNNA3;DAAM2;EPB41L2;FMNL2;GIPC1;HIP1;HIP1R;IQGAP1;LDB3;LIMA1;MARCKSL1;MYO1D;MYO6;MYO9B;PIP;PLEKHG3;SHTN1;SSH3;SYNE3;TPRN;UACA;WASF2); GO: actin filament binding (q=0.0159; n=16; genes=AIF1L;AMOTL2;BIN1;CLMN;CTNNA3;FMNL2;HIP1;HIP1R;IQGAP1;LIMA1;MARCKSL1;MYO1D;MYO6;MYO9B;SHTN1;SYNE3); GO: semaphorin receptor binding (q=0.0242; n=6; genes=PLXNB1;SEMA3B;SEMA4C;SEMA4D;SEMA6A;SEMA6D); GO: carbohydrate binding (q=0.0242; n=17; genes=CD22;CNTN2;ENPP1;GALK1;GALNT12;GALNT8;GANC;KLRC2;KLRC3;KLRC4-KLRK1;KLRK1;MAN2A1;MAN2A2;PKD1L2;PLOD3;PTN;SFTPA2)
- Main caution: Support is confined to one annotation database.
- Safe wording: The Hippocampal Formation MEgreen gene set was annotated to actin cell cortex terms in the fixed Phase-3B enrichment analysis.
- Wording to avoid: Avoid claiming that actin cell cortex is dysregulated, causal, replicated, or a validated MDD mechanism in Hippocampal Formation.

### dlPFC BA8/9 MEblue

- Annotation: **RECEPTOR_GPCR_FAMILY**; category: **SUPPLEMENTARY_CONTEXT**.
- Diagnosis: BH q=0.272; not BH-FDR significant.
- Corrected ML overlap: high=0; high+moderate=2; high+moderate q=1.
- Exact corrected-ML genes in module: high [none]; high+moderate [LOC105377975;METTL14-DT].
- Leading fixed terms: Reactome: Class C/3 (Metabotropic glutamate/pheromone receptors) (q=1.26e-13; n=12; genes=TAS2R10;TAS2R13;TAS2R14;TAS2R19;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7;TAS2R9); GO: taste receptor activity (q=2.86e-12; n=11; genes=TAS2R10;TAS2R13;TAS2R14;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7;TAS2R9); GO: bitter taste receptor activity (q=2.86e-12; n=11; genes=TAS2R10;TAS2R13;TAS2R14;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7;TAS2R9); GO: G protein-coupled receptor activity (q=4.44e-11; n=26; genes=ACKR1;ACKR4;GNRHR;GPR171;GPR18;GPR87;HTR2B;OR2D2;OR2D3;OR6A2;OR6X1;P2RY14;RRH;TAS2R10;TAS2R13;TAS2R14;TAS2R19;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7;TAS2R9;TRHR); Reactome: Sensory perception of sweet, bitter, and umami (glutamate) taste (q=1.58e-10; n=10; genes=TAS2R10;TAS2R13;TAS2R14;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7)
- Main caution: Most supporting genes belong to one recognized gene family.
- Safe wording: The dlPFC BA8/9 MEblue gene set was annotated to receptor gpcr family terms in the fixed Phase-3B enrichment analysis.
- Wording to avoid: Avoid claiming that receptor gpcr family is dysregulated, causal, replicated, or a validated MDD mechanism in dlPFC BA8/9.

### dlPFC BA8/9 MEorange

- Annotation: **SIGNALING_IMMEDIATE_EARLY**; category: **SUPPLEMENTARY_CONTEXT**.
- Diagnosis: BH q=0.272; not BH-FDR significant.
- Corrected ML overlap: high=0; high+moderate=0; high+moderate q=1.
- Exact corrected-ML genes in module: high [none]; high+moderate [none].
- Leading fixed terms: Reactome: Nuclear Events (kinase and transcription factor activation) (q=5.34e-11; n=8; genes=ARC;DUSP4;DUSP6;EGR1;EGR2;EGR3;FOS;JUNB); Reactome: Signaling by NTRK1 (TRKA) (q=8.22e-11; n=8; genes=ARC;DUSP4;DUSP6;EGR1;EGR2;EGR3;FOS;JUNB); Reactome: Signaling by NTRKs (q=2.29e-10; n=8; genes=ARC;DUSP4;DUSP6;EGR1;EGR2;EGR3;FOS;JUNB); Reactome: NGF-stimulated transcription (q=3.75e-08; n=6; genes=ARC;EGR1;EGR2;EGR3;FOS;JUNB); MSigDB Hallmark: HALLMARK_TNFA_SIGNALING_VIA_NFKB (q=5.13e-08; n=11; genes=BHLHE40;CCN1;DUSP1;DUSP4;DUSP5;EGR1;EGR2;EGR3;FOS;JUNB;NR4A3)
- Main caution: Module has no more than 50 mapped genes; term dependence is substantial.
- Safe wording: The dlPFC BA8/9 MEorange gene set was annotated to signaling immediate early terms in the fixed Phase-3B enrichment analysis.
- Wording to avoid: Avoid claiming that signaling immediate early is dysregulated, causal, replicated, or a validated MDD mechanism in dlPFC BA8/9.

### Cingulate gyrus 25 MEwhite

- Annotation: **RNA_PROCESSING_NUCLEOLAR**; category: **SUPPLEMENTARY_CONTEXT**.
- Diagnosis: BH q=0.897; not BH-FDR significant.
- Corrected ML overlap: high=0; high+moderate=1; high+moderate q=1.
- Exact corrected-ML genes in module: high [none]; high+moderate [EXOC3-AS1].
- Leading fixed terms: GO: RNA processing (q=2.68e-46; n=41; genes=MIR1843;RNU4-1;RNU5A-1;RNU5B-1;RNU6-7;RNU6-9;RNU6ATAC;SCARNA3;SCARNA4;SCARNA8;SNORA1;SNORA10;SNORA10B;SNORA11;SNORA17A;SNORA21;SNORA33;SNORA38B;SNORA44;SNORA54;SNORA58;SNORA6;SNORA63;SNORA63E;SNORA65;SNORA67;SNORA68;SNORA68B;SNORA71B;SNORA71E;SNORA73A;SNORA75;SNORA77B;SNORA80E;SNORA9;SNORD10;SNORD17;SNORD22;SNORD67;SNORD89;SNORD90); GO: nucleolus (q=1.05e-33; n=35; genes=MIR1843;SCARNA3;SCARNA4;SCARNA8;SNORA1;SNORA10;SNORA10B;SNORA11;SNORA17A;SNORA21;SNORA33;SNORA38B;SNORA44;SNORA54;SNORA58;SNORA6;SNORA63;SNORA63E;SNORA65;SNORA67;SNORA68;SNORA68B;SNORA71B;SNORA71E;SNORA73A;SNORA75;SNORA77B;SNORA80E;SNORA9;SNORD10;SNORD17;SNORD22;SNORD67;SNORD89;SNORD90); GO: spliceosomal snRNP assembly (q=8.93e-08; n=6; genes=RNU4-1;RNU5A-1;RNU5B-1;RNU6-7;RNU6-9;RNU6ATAC); GO: spliceosomal tri-snRNP complex assembly (q=8.93e-08; n=6; genes=RNU4-1;RNU5A-1;RNU5B-1;RNU6-7;RNU6-9;RNU6ATAC); GO: Cajal body (q=1.43e-06; n=5; genes=MIR1843;SCARNA3;SCARNA4;SCARNA8;SNORA38B)
- Main caution: Most supporting genes belong to one recognized gene family.
- Safe wording: The Cingulate gyrus 25 MEwhite gene set was annotated to rna processing nucleolar terms in the fixed Phase-3B enrichment analysis.
- Wording to avoid: Avoid claiming that rna processing nucleolar is dysregulated, causal, replicated, or a validated MDD mechanism in Cingulate gyrus 25.

### OFC BA11 MEtan

- Annotation: **ANGIOGENESIS_VASCULAR**; category: **SUPPLEMENTARY_CONTEXT**.
- Diagnosis: BH q=0.549; not BH-FDR significant.
- Corrected ML overlap: high=2; high+moderate=3; high+moderate q=0.0686.
- Exact corrected-ML genes in module: high [LINC01361;LOC107985795]; high+moderate [GPR6;LINC01361;LOC107985795].
- Leading fixed terms: MSigDB Hallmark: HALLMARK_ANGIOGENESIS (q=0.024; n=2; genes=COL5A2;POSTN)
- Main caution: Theme is driven by no more than three distinct genes.
- Safe wording: The OFC BA11 MEtan gene set was annotated to angiogenesis vascular terms in the fixed Phase-3B enrichment analysis.
- Wording to avoid: Avoid claiming that angiogenesis vascular is dysregulated, causal, replicated, or a validated MDD mechanism in OFC BA11.

## Claims suitable for main text

| region | module | dominant_biological_theme | safe_wording | interpretation_caution |
|---|---|---|---|---|
| Ventral Striatum | MElavenderblush3 | SYNAPTIC_PRESYNAPTIC | The Ventral Striatum MElavenderblush3 gene set was annotated to synaptic presynaptic terms in the fixed Phase-3B enrichment analysis. | Module has no more than 50 mapped genes; term dependence is substantial. |
| Ventral Striatum | MElightgreen | SYNAPTIC_PRESYNAPTIC | The Ventral Striatum MElightgreen gene set was annotated to synaptic presynaptic terms in the fixed Phase-3B enrichment analysis. | Support is confined to one annotation database. |
| Ventral Striatum | MElightcyan1 | CHOLESTEROL_LIPID | The Ventral Striatum MElightcyan1 gene set was annotated to cholesterol lipid terms in the fixed Phase-3B enrichment analysis. | Interpret as module annotation only; it is not MDD-association evidence. |
| Hippocampal Formation | MEgreen | ACTIN_CELL_CORTEX | The Hippocampal Formation MEgreen gene set was annotated to actin cell cortex terms in the fixed Phase-3B enrichment analysis. | Support is confined to one annotation database. |

## Claims suitable only for supplement

| region | module | dominant_biological_theme | safe_wording | interpretation_caution |
|---|---|---|---|---|
| Ventral Striatum | MEdarkred | EXTRACELLULAR_MATRIX_WOUND | The Ventral Striatum MEdarkred gene set was annotated to extracellular matrix wound terms in the fixed Phase-3B enrichment analysis. | Support is confined to one annotation database. |
| dlPFC BA8/9 | MEblue | RECEPTOR_GPCR_FAMILY | The dlPFC BA8/9 MEblue gene set was annotated to receptor gpcr family terms in the fixed Phase-3B enrichment analysis. | Most supporting genes belong to one recognized gene family. |
| dlPFC BA8/9 | MEorange | SIGNALING_IMMEDIATE_EARLY | The dlPFC BA8/9 MEorange gene set was annotated to signaling immediate early terms in the fixed Phase-3B enrichment analysis. | Module has no more than 50 mapped genes; term dependence is substantial. |
| Cingulate gyrus 25 | MEwhite | RNA_PROCESSING_NUCLEOLAR | The Cingulate gyrus 25 MEwhite gene set was annotated to rna processing nucleolar terms in the fixed Phase-3B enrichment analysis. | Most supporting genes belong to one recognized gene family. |
| OFC BA11 | MEtan | ANGIOGENESIS_VASCULAR | The OFC BA11 MEtan gene set was annotated to angiogenesis vascular terms in the fixed Phase-3B enrichment analysis. | Theme is driven by no more than three distinct genes. |

## Claims to avoid

| region | module | dominant_biological_theme | wording_to_avoid |
|---|---|---|---|
| Insula | MEorange | NO_CLEAR_THEME | Avoid claiming that no clear theme is dysregulated, causal, replicated, or a validated MDD mechanism in Insula. |

Across all modules, avoid causal, biomarker, replicated-mechanism, or Diagnosis-significant wording. Annotation coherence does not override the absence of module-eigengene Diagnosis BH-FDR significance.

## Remaining biological uncertainties

- Functional databases are overlapping and hierarchically redundant.
- Small modules and small overlap sets can yield very low q-values while remaining gene-concentrated.
- Several gene-level Diagnosis statistics are unavailable outside the Tier2/corrected-ML table and were not reconstructed.
- Mapping improved substantially, but unmapped identifiers remain and may not be biologically random.
- No literature search or external biological validation was used to rescue or strengthen a theme.

## Final decision table

| region | module | dominant_biological_theme | supporting_databases | key_driver_genes | overlap_size | mapping_coverage | ML_high_plus_moderate_q | Diagnosis_FDR_status | interpretation_caution | publication_use_category |
|---|---|---|---|---|---|---|---|---|---|---|
| Ventral Striatum | MElavenderblush3 | SYNAPTIC_PRESYNAPTIC | GO;Reactome | ATP6V1H;CACNA2D1;CAP2;CLCN4;DMXL2;GLRB;NSF;NTNG1;PRKAA2;SYT1 | 10 | 0.9714285714285714 | 1.0 | NOT_BH_FDR_SIGNIFICANT | Module has no more than 50 mapped genes; term dependence is substantial. | MAIN_TEXT_CANDIDATE |
| Ventral Striatum | MEdarkred | EXTRACELLULAR_MATRIX_WOUND | GO | ANO6;DDR2;DPP4;EMILIN2;F2R;HPSE;LAMB2;SCARF1;TGFB2;TNFRSF1B | 10 | 0.9559471365638766 | 0.1085068144870678 | NOT_BH_FDR_SIGNIFICANT | Support is confined to one annotation database. | SUPPLEMENTARY_CONTEXT |
| Ventral Striatum | MElightgreen | SYNAPTIC_PRESYNAPTIC | GO | ADAM22;AKAP5;AP3S1;ASIC2;ATP6AP2;ATP6V1A;ATP6V1B2;ATP6V1C1;ATP8A1;CADPS;CNTN1;EPHA7;FZD3;GABRA1;GABRA4;GABRG2;GAD1;GAD2;GAP43;GRIA2;MAP2K1;NAPB;NECAP1;NEGR1;OPRK1;PAK3;PJA2;PSD3;RAB3C;RAPGEF4;RTN1;SCAMP1;SLC35F1;SLC4A10;SLC8A1;SVOP;SYNGR3;SYNPR;SYP;SYT4;SYTL5;YWHAH;YWHAZ | 43 | 0.9899665551839464 | 1.0 | NOT_BH_FDR_SIGNIFICANT | Support is confined to one annotation database. | MAIN_TEXT_CANDIDATE |
| Ventral Striatum | MElightcyan1 | CHOLESTEROL_LIPID | MSigDB Hallmark;Reactome | ACAT2;CYP51A1;GLDC;HMGCS1;IDI1;MSMO1;PNRC1;PPARG;SQLE;TMEM97 | 10 | 0.9829351535836176 | 1.0 | NOT_BH_FDR_SIGNIFICANT | Interpret as module annotation only; it is not MDD-association evidence. | MAIN_TEXT_CANDIDATE |
| Insula | MEorange | NO_CLEAR_THEME |  |  | 0 | 1.0 | 1.0 | NOT_BH_FDR_SIGNIFICANT | No Phase-3B term reached BH q<0.10. | AVOID_INTERPRETATION |
| Hippocampal Formation | MEgreen | ACTIN_CELL_CORTEX | GO | AIF1L;AMOTL2;BIN1;CCDC88A;CDH1;CLMN;CNKSR1;CTNNA3;CTTNBP2;DAAM2;EPB41L2;FGF1;FMNL2;FNBP1;GIPC1;GPSM2;HIP1;HIP1R;IQGAP1;LDB3;LIMA1;LLGL2;MARCKSL1;MYO1D;MYO6;MYO9B;PHLDB1;PIEZO2;PIP;PLA2G4C;PLEKHG3;RAPGEF3;RTKN;SEPTIN1;SEPTIN12;SEPTIN4;SHTN1;SSH3;SYNE3;TPRN;UACA;WASF2 | 42 | 0.9288214702450408 | 1.0 | NOT_BH_FDR_SIGNIFICANT | Support is confined to one annotation database. | MAIN_TEXT_CANDIDATE |
| dlPFC BA8/9 | MEblue | RECEPTOR_GPCR_FAMILY | GO;Reactome | ACKR1;ACKR4;CDON;GNRHR;GPR171;GPR18;GPR87;HTR2B;IL18R1;OR2D2;OR2D3;OR6A2;OR6X1;P2RY14;PDC;PDE10A;PDE7B;RRH;TAS2R10;TAS2R13;TAS2R14;TAS2R19;TAS2R3;TAS2R31;TAS2R43;TAS2R46;TAS2R5;TAS2R50;TAS2R7;TAS2R9;TRHR;VEGFA | 32 | 0.8691322901849218 | 1.0 | NOT_BH_FDR_SIGNIFICANT | Most supporting genes belong to one recognized gene family. | SUPPLEMENTARY_CONTEXT |
| dlPFC BA8/9 | MEorange | SIGNALING_IMMEDIATE_EARLY | GO;Reactome | ARC;BHLHE40;CCN1;CSRNP1;DUSP1;DUSP4;DUSP5;DUSP6;EGR1;EGR2;EGR3;ETV5;FOS;JUNB;NR4A3;PER2;SPRY4;TRIB2 | 18 | 1.0 | 1.0 | NOT_BH_FDR_SIGNIFICANT | Module has no more than 50 mapped genes; term dependence is substantial. | SUPPLEMENTARY_CONTEXT |
| Cingulate gyrus 25 | MEwhite | RNA_PROCESSING_NUCLEOLAR | GO | MIR1843;RNU4-1;RNU5A-1;RNU5B-1;RNU6-7;RNU6-9;RNU6ATAC;SCARNA3;SCARNA4;SCARNA8;SNORA1;SNORA10;SNORA10B;SNORA11;SNORA17A;SNORA21;SNORA33;SNORA38B;SNORA44;SNORA54;SNORA58;SNORA6;SNORA63;SNORA63E;SNORA65;SNORA67;SNORA68;SNORA68B;SNORA71B;SNORA71E;SNORA73A;SNORA75;SNORA77B;SNORA80E;SNORA9;SNORD10;SNORD17;SNORD22;SNORD67;SNORD89;SNORD90 | 41 | 1.0 | 1.0 | NOT_BH_FDR_SIGNIFICANT | Most supporting genes belong to one recognized gene family. | SUPPLEMENTARY_CONTEXT |
| OFC BA11 | MEtan | ANGIOGENESIS_VASCULAR | MSigDB Hallmark | COL5A2;POSTN | 2 | 0.951219512195122 | 0.0686067539236891 | NOT_BH_FDR_SIGNIFICANT | Theme is driven by no more than three distinct genes. | SUPPLEMENTARY_CONTEXT |
