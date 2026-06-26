# Tier 2 Enrichment Report

## Scope

This enrichment analysis used the Length > 100 Tier 2 exploratory DE candidate gene sets.

Tier 2 definition:

```text
p_value_diagnosis < 0.01 and abs_standardized_effect_size_diagnosis > 0.8
```

Tier 2 genes are exploratory DE candidate genes, not FDR-significant DE genes.

The enrichment background for each region/subregion was the tested-gene universe from that same analysis, not all genes in the genome.

## Input And Mapping Summary

    analysis_level   region_or_subregion   direction n_input_genes
        four_group  broad Frontal Cortex         all             1
        four_group  broad Frontal Cortex down_in_MDD             0
        four_group  broad Frontal Cortex   up_in_MDD             1
 frontal_subregion    Cingulate gyrus 25         all           186
 frontal_subregion    Cingulate gyrus 25 down_in_MDD            49
 frontal_subregion    Cingulate gyrus 25   up_in_MDD           137
 frontal_subregion           dlPFC BA8/9         all            72
 frontal_subregion           dlPFC BA8/9 down_in_MDD            48
 frontal_subregion           dlPFC BA8/9   up_in_MDD            24
        four_group Hippocampal Formation         all           127
        four_group Hippocampal Formation down_in_MDD            47
        four_group Hippocampal Formation   up_in_MDD            80
        four_group                Insula         all            70
        four_group                Insula down_in_MDD            33
        four_group                Insula   up_in_MDD            37
 frontal_subregion              OFC BA11         all            70
 frontal_subregion              OFC BA11 down_in_MDD            32
 frontal_subregion              OFC BA11   up_in_MDD            38
        four_group      Ventral Striatum         all           146
        four_group      Ventral Striatum down_in_MDD            70
        four_group      Ventral Striatum   up_in_MDD            76
 n_mapped_input_genes n_background_genes n_mapped_background_genes
                    1              23941                     23941
                    0              23941                     23941
                    1              23941                     23941
                  186              23941                     23941
                   49              23941                     23941
                  137              23941                     23941
                   72              23941                     23941
                   48              23941                     23941
                   24              23941                     23941
                  127              24371                     24371
                   47              24371                     24371
                   80              24371                     24371
                   70              23090                     23090
                   33              23090                     23090
                   37              23090                     23090
                   70              23941                     23941
                   32              23941                     23941
                   38              23941                     23941
                  146              23917                     23917
                   70              23917                     23917
                   76              23917                     23917

Gene sets with fewer than 5 mapped input genes were skipped. Broad Frontal Cortex was retained for documentation only and is not useful for enrichment because it has only 1 Tier 2 gene.

## Main Enrichment Summary

    analysis_level   region_or_subregion   direction            database
 frontal_subregion              OFC BA11 down_in_MDD        MSigDB_C5_GO
 frontal_subregion           dlPFC BA8/9   up_in_MDD               GO_MF
 frontal_subregion              OFC BA11 down_in_MDD               GO_CC
 frontal_subregion           dlPFC BA8/9   up_in_MDD               GO_CC
        four_group Hippocampal Formation   up_in_MDD                KEGG
 frontal_subregion           dlPFC BA8/9         all               GO_CC
        four_group                Insula down_in_MDD   MSigDB_C2_curated
 frontal_subregion              OFC BA11   up_in_MDD MSigDB_C8_cell_type
 frontal_subregion              OFC BA11         all               GO_BP
 frontal_subregion    Cingulate gyrus 25 down_in_MDD               GO_MF
        four_group                Insula down_in_MDD               GO_MF
 frontal_subregion              OFC BA11 down_in_MDD                KEGG
 frontal_subregion           dlPFC BA8/9   up_in_MDD                KEGG
 frontal_subregion    Cingulate gyrus 25         all               GO_CC
 frontal_subregion    Cingulate gyrus 25 down_in_MDD MSigDB_C8_cell_type
 frontal_subregion    Cingulate gyrus 25   up_in_MDD               GO_MF
 frontal_subregion              OFC BA11 down_in_MDD               GO_BP
 frontal_subregion              OFC BA11 down_in_MDD     MSigDB_Hallmark
        four_group      Ventral Striatum   up_in_MDD        MSigDB_C5_GO
 frontal_subregion    Cingulate gyrus 25 down_in_MDD               GO_BP
        four_group Hippocampal Formation down_in_MDD               GO_MF
        four_group                Insula down_in_MDD     MSigDB_Hallmark
 frontal_subregion              OFC BA11         all               GO_CC
        four_group      Ventral Striatum   up_in_MDD MSigDB_C8_cell_type
 frontal_subregion           dlPFC BA8/9 down_in_MDD     MSigDB_Hallmark
        four_group                Insula         all MSigDB_C8_cell_type
        four_group                Insula down_in_MDD               GO_CC
        four_group                Insula down_in_MDD MSigDB_C8_cell_type
 frontal_subregion              OFC BA11         all        MSigDB_C5_GO
 frontal_subregion              OFC BA11 down_in_MDD MSigDB_C8_cell_type
 n_terms_padj005 n_terms_padj010
               4              42
               0              26
              12              18
              13              17
               6              17
               0              16
               5              12
               1               9
               0               8
               5               7
               0               6
               0               6
               1               5
               1               4
               3               4
               0               4
               3               4
               0               4
               4               4
               0               3
               3               3
               1               3
               3               3
               2               3
               0               2
               2               2
               1               2
               1               2
               2               2
               0               2
                                               top_term_padj010 status
                                GOBP_RENAL_SODIUM_ION_TRANSPORT     ok
                          single-stranded telomeric DNA binding     ok
                    nuclear DNA-directed RNA polymerase complex     ok
                                           spliceosomal complex     ok
                        Progesterone-mediated oocyte maturation     ok
                                           spliceosomal complex     ok
                                MIKKELSEN_MEF_HCP_WITH_H3K27ME3     ok
                         DESCARTES_MAIN_FETAL_ENDOCARDIAL_CELLS     ok
                            cellular response to osmotic stress     ok
                              protein kinase inhibitor activity     ok
 transmembrane receptor protein serine/threonine kinase binding     ok
                                                 RNA polymerase     ok
                                                    Spliceosome     ok
                                 Sm-like protein family complex     ok
                               MANNO_MIDBRAIN_NEUROTYPES_HRGL2C     ok
                           aldo-keto reductase (NADPH) activity     ok
                                    renal sodium ion absorption     ok
                       HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY     ok
                                         HP_DECREASED_FERTILITY     ok
                            negative regulation of MAPK cascade     ok
      MAP kinase tyrosine/serine/threonine phosphatase activity     ok
                                    HALLMARK_TGF_BETA_SIGNALING     ok
                    nuclear DNA-directed RNA polymerase complex     ok
                 DESCARTES_FETAL_LUNG_CILIATED_EPITHELIAL_CELLS     ok
                                       HALLMARK_HEME_METABOLISM     ok
                 DESCARTES_FETAL_THYMUS_THYMIC_EPITHELIAL_CELLS     ok
                                           axon initial segment     ok
                 DESCARTES_FETAL_THYMUS_THYMIC_EPITHELIAL_CELLS     ok
                       GOBP_CELLULAR_RESPONSE_TO_OSMOTIC_STRESS     ok
                              HE_LIM_SUN_FETAL_LUNG_C2_DC1_CELL     ok

## Region/Subregion Interpretation

### Ventral Striatum

   direction            database n_terms_padj005 n_terms_padj010
   up_in_MDD        MSigDB_C5_GO               4               4
   up_in_MDD MSigDB_C8_cell_type               2               3
 down_in_MDD   MSigDB_C2_curated               1               2
   up_in_MDD               GO_MF               0               2
   up_in_MDD               GO_CC               0               1
                                top_term_padj010 status
                          HP_DECREASED_FERTILITY     ok
  DESCARTES_FETAL_LUNG_CILIATED_EPITHELIAL_CELLS     ok
                    GRATIAS_RETINOBLASTOMA_16Q24     ok
 all-trans-retinol dehydrogenase (NAD+) activity     ok
                                   motile cilium     ok

### Insula

   direction            database n_terms_padj005 n_terms_padj010
 down_in_MDD   MSigDB_C2_curated               5              12
 down_in_MDD               GO_MF               0               6
 down_in_MDD     MSigDB_Hallmark               1               3
         all MSigDB_C8_cell_type               2               2
 down_in_MDD               GO_CC               1               2
 down_in_MDD MSigDB_C8_cell_type               1               2
         all   MSigDB_C2_curated               1               1
         all     MSigDB_Hallmark               0               1
   up_in_MDD                KEGG               1               1
                                               top_term_padj010 status
                                MIKKELSEN_MEF_HCP_WITH_H3K27ME3     ok
 transmembrane receptor protein serine/threonine kinase binding     ok
                                    HALLMARK_TGF_BETA_SIGNALING     ok
                 DESCARTES_FETAL_THYMUS_THYMIC_EPITHELIAL_CELLS     ok
                                           axon initial segment     ok
                 DESCARTES_FETAL_THYMUS_THYMIC_EPITHELIAL_CELLS     ok
                                MIKKELSEN_MEF_HCP_WITH_H3K27ME3     ok
                                    HALLMARK_TGF_BETA_SIGNALING     ok
                                     PI3K-Akt signaling pathway     ok

### Hippocampal Formation

   direction database n_terms_padj005 n_terms_padj010
   up_in_MDD     KEGG               6              17
 down_in_MDD    GO_MF               3               3
         all    GO_CC               0               1
                                          top_term_padj010 status
                   Progesterone-mediated oocyte maturation     ok
 MAP kinase tyrosine/serine/threonine phosphatase activity     ok
                                              X chromosome     ok

### Cingulate gyrus 25

   direction            database n_terms_padj005 n_terms_padj010
 down_in_MDD               GO_MF               5               7
         all               GO_CC               1               4
 down_in_MDD MSigDB_C8_cell_type               3               4
   up_in_MDD               GO_MF               0               4
 down_in_MDD               GO_BP               0               3
         all                KEGG               0               1
 down_in_MDD     MSigDB_Hallmark               0               1
   up_in_MDD               GO_BP               0               1
   up_in_MDD               GO_CC               0               1
   up_in_MDD                KEGG               1               1
                                        top_term_padj010 status
                       protein kinase inhibitor activity     ok
                          Sm-like protein family complex     ok
                        MANNO_MIDBRAIN_NEUROTYPES_HRGL2C     ok
                    aldo-keto reductase (NADPH) activity     ok
                     negative regulation of MAPK cascade     ok
                                    Porphyrin metabolism     ok
                        HALLMARK_TNFA_SIGNALING_VIA_NFKB     ok
 regulation of release of cytochrome c from mitochondria     ok
                          Sm-like protein family complex     ok
                                    Porphyrin metabolism     ok

### dlPFC BA8/9

   direction        database n_terms_padj005 n_terms_padj010
   up_in_MDD           GO_MF               0              26
   up_in_MDD           GO_CC              13              17
         all           GO_CC               0              16
   up_in_MDD            KEGG               1               5
 down_in_MDD MSigDB_Hallmark               0               2
         all            KEGG               0               1
   up_in_MDD MSigDB_Hallmark               1               1
                      top_term_padj010 status
 single-stranded telomeric DNA binding     ok
                  spliceosomal complex     ok
                  spliceosomal complex     ok
                           Spliceosome     ok
              HALLMARK_HEME_METABOLISM     ok
                           Spliceosome     ok
               HALLMARK_MYC_TARGETS_V1     ok

### OFC BA11

   direction            database n_terms_padj005 n_terms_padj010
 down_in_MDD        MSigDB_C5_GO               4              42
 down_in_MDD               GO_CC              12              18
   up_in_MDD MSigDB_C8_cell_type               1               9
         all               GO_BP               0               8
 down_in_MDD                KEGG               0               6
 down_in_MDD               GO_BP               3               4
 down_in_MDD     MSigDB_Hallmark               0               4
         all               GO_CC               3               3
         all        MSigDB_C5_GO               2               2
 down_in_MDD MSigDB_C8_cell_type               0               2
                            top_term_padj010 status
             GOBP_RENAL_SODIUM_ION_TRANSPORT     ok
 nuclear DNA-directed RNA polymerase complex     ok
      DESCARTES_MAIN_FETAL_ENDOCARDIAL_CELLS     ok
         cellular response to osmotic stress     ok
                              RNA polymerase     ok
                 renal sodium ion absorption     ok
    HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY     ok
 nuclear DNA-directed RNA polymerase complex     ok
    GOBP_CELLULAR_RESPONSE_TO_OSMOTIC_STRESS     ok
           HE_LIM_SUN_FETAL_LUNG_C2_DC1_CELL     ok

## Up-In-MDD Versus Down-In-MDD Patterns

   region_or_subregion   direction            database n_terms_padj010
    Cingulate gyrus 25 down_in_MDD               GO_MF               7
    Cingulate gyrus 25 down_in_MDD MSigDB_C8_cell_type               4
    Cingulate gyrus 25 down_in_MDD               GO_BP               3
    Cingulate gyrus 25 down_in_MDD     MSigDB_Hallmark               1
    Cingulate gyrus 25   up_in_MDD               GO_MF               4
    Cingulate gyrus 25   up_in_MDD               GO_BP               1
    Cingulate gyrus 25   up_in_MDD               GO_CC               1
    Cingulate gyrus 25   up_in_MDD                KEGG               1
           dlPFC BA8/9 down_in_MDD     MSigDB_Hallmark               2
           dlPFC BA8/9   up_in_MDD               GO_MF              26
           dlPFC BA8/9   up_in_MDD               GO_CC              17
           dlPFC BA8/9   up_in_MDD                KEGG               5
           dlPFC BA8/9   up_in_MDD     MSigDB_Hallmark               1
 Hippocampal Formation down_in_MDD               GO_MF               3
 Hippocampal Formation   up_in_MDD                KEGG              17
                Insula down_in_MDD   MSigDB_C2_curated              12
                Insula down_in_MDD               GO_MF               6
                Insula down_in_MDD     MSigDB_Hallmark               3
                Insula down_in_MDD               GO_CC               2
                Insula down_in_MDD MSigDB_C8_cell_type               2
                Insula   up_in_MDD                KEGG               1
              OFC BA11 down_in_MDD        MSigDB_C5_GO              42
              OFC BA11 down_in_MDD               GO_CC              18
              OFC BA11 down_in_MDD                KEGG               6
              OFC BA11 down_in_MDD               GO_BP               4
              OFC BA11 down_in_MDD     MSigDB_Hallmark               4
              OFC BA11 down_in_MDD MSigDB_C8_cell_type               2
              OFC BA11 down_in_MDD   MSigDB_C2_curated               1
              OFC BA11   up_in_MDD MSigDB_C8_cell_type               9
      Ventral Striatum down_in_MDD   MSigDB_C2_curated               2
                                               top_term_padj010
                              protein kinase inhibitor activity
                               MANNO_MIDBRAIN_NEUROTYPES_HRGL2C
                            negative regulation of MAPK cascade
                               HALLMARK_TNFA_SIGNALING_VIA_NFKB
                           aldo-keto reductase (NADPH) activity
        regulation of release of cytochrome c from mitochondria
                                 Sm-like protein family complex
                                           Porphyrin metabolism
                                       HALLMARK_HEME_METABOLISM
                          single-stranded telomeric DNA binding
                                           spliceosomal complex
                                                    Spliceosome
                                        HALLMARK_MYC_TARGETS_V1
      MAP kinase tyrosine/serine/threonine phosphatase activity
                        Progesterone-mediated oocyte maturation
                                MIKKELSEN_MEF_HCP_WITH_H3K27ME3
 transmembrane receptor protein serine/threonine kinase binding
                                    HALLMARK_TGF_BETA_SIGNALING
                                           axon initial segment
                 DESCARTES_FETAL_THYMUS_THYMIC_EPITHELIAL_CELLS
                                     PI3K-Akt signaling pathway
                                GOBP_RENAL_SODIUM_ION_TRANSPORT
                    nuclear DNA-directed RNA polymerase complex
                                                 RNA polymerase
                                    renal sodium ion absorption
                       HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY
                              HE_LIM_SUN_FETAL_LUNG_C2_DC1_CELL
                                   HAMAI_APOPTOSIS_VIA_TRAIL_DN
                         DESCARTES_MAIN_FETAL_ENDOCARDIAL_CELLS
                                   GRATIAS_RETINOBLASTOMA_16Q24

## Cautious Interpretation

These enrichment results are exploratory and depend on the Tier 2 candidate definition. They should be used to prioritize biological hypotheses and downstream analyses, not as final proof of pathway-level disease mechanisms.

The p-value correction sensitivity audit showed that larger single-gene counts mainly appear under effect-size-gated or direction-specific correction scopes. Those remain sensitivity results and do not replace Tier 2.

## Recommended Next Step

After this ORA enrichment, run ranked GSEA using the Tier 4 signed standardized log-p signatures.

## Outputs

- `tables/`: full, adjusted p < 0.05, and adjusted p < 0.10 enrichment tables.
- `plots/`: dotplots, barplots, and enrichment maps when enough terms were available.
- `unmapped_genes/`: unmapped input and background genes for each gene set.
- `enrichment_summary_by_gene_set.csv`: summary of mapping and term counts.
- `top_enrichment_terms_all_sets.csv`: top 20 terms from each successful enrichment analysis.
