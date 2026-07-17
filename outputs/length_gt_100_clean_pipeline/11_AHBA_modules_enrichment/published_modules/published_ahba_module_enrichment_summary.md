# Published AHBA Canonical Module Enrichment Summary

## Source

Hawrylycz et al. 2015 Supplementary Table 6 was downloaded from the official Springer/Nature supplementary materials and standardized into:

```text
D:\genes\GeneDepresion\outputs\length_gt_100_clean_pipeline\11_AHBA_modules_enrichment\published_modules\ahba_published_gene_to_module_assignments.csv
```

Official source URL:

```text
https://static-content.springer.com/esm/art%3A10.1038%2Fnn.4171/MediaObjects/41593_2015_BFnn4171_MOESM97_ESM.csv
```

## Published Module Table

- Published AHBA module-assignment rows: 15,627
- Unique published AHBA module genes: 15,627
- Published modules: 32
- MDD-tested ? published AHBA module background genes: 12,042

## All-Region Candidate Overlap

| candidate_set      |   original_gene_region_pair_count |   unique_gene_count |   unique_genes_overlapping_published_AHBA_background |   modules_tested |
|:-------------------|----------------------------------:|--------------------:|-----------------------------------------------------:|-----------------:|
| all_Tier2          |                               671 |                 640 |                                                  276 |               32 |
| high_confidence    |                                35 |                  35 |                                                   12 |               32 |
| high_plus_moderate |                               165 |                 154 |                                                   56 |               32 |

## Module Sizes In The MDD-Tested Background

| published_module   |   background_gene_count |
|:-------------------|------------------------:|
| M30                |                    1438 |
| M24                |                     982 |
| M19                |                     899 |
| M23                |                     838 |
| M27                |                     713 |
| M20                |                     685 |
| M26                |                     554 |
| M28                |                     523 |
| M17                |                     485 |
| M31                |                     466 |
| M6                 |                     455 |
| M29                |                     446 |
| M4                 |                     418 |
| M1                 |                     358 |
| M11                |                     245 |
| M21                |                     224 |
| M3                 |                     218 |
| M18                |                     193 |
| M8                 |                     180 |
| M15                |                     172 |
| M9                 |                     163 |
| M25                |                     160 |
| M2                 |                     159 |
| M5                 |                     149 |
| M32                |                     145 |
| M22                |                     145 |
| M7                 |                     139 |
| M10                |                     136 |
| M14                |                     128 |
| M16                |                     100 |
| M13                |                      81 |
| M12                |                      45 |

## Enrichment Result Counts

- Nominal enrichment rows with p < 0.05: 41
- Rows with BH-FDR < 0.10: 6
- Rows with BH-FDR < 0.05: 3

## Top Enrichment Rows

| candidate_set      | scope                                              | published_module   |   candidate_genes_in_background |   module_genes_in_background |   overlap_genes |   p_value_hypergeometric |      fdr_BH | overlap_gene_symbols                                                 |
|:-------------------|:---------------------------------------------------|:-------------------|--------------------------------:|-----------------------------:|----------------:|-------------------------:|------------:|:---------------------------------------------------------------------|
| all_Tier2          | region:Ventral Striatum|direction:up_in_MDD        | M26                |                              30 |                          554 |               9 |              5.19607e-06 | 0.000166274 | C4ORF47;CDNF;FBXL13;MORN5;MYCBPAP;RPGR;TNFRSF10B;TSGA10;TTC23L       |
| all_Tier2          | region:Cingulate gyrus 25|direction:down_in_MDD    | M2                 |                              27 |                          159 |               4 |              0.000405302 | 0.0129697   | DUSP4;NR4A3;PTPRE;TG                                                 |
| all_Tier2          | region:Ventral Striatum|direction:down_in_MDD      | M6                 |                              51 |                          455 |               8 |              0.000595408 | 0.0190531   | DIRAS2;GLDC;KCNV1;KIAA0513;KLHL18;MCHR2;PRKAA2;SHC3                  |
| all_Tier2          | region:Hippocampal Formation|direction:down_in_MDD | M7                 |                              24 |                          139 |               3 |              0.00255056  | 0.051894    | DUSP14;PRMT7;PRR7                                                    |
| all_Tier2          | region:Hippocampal Formation|direction:down_in_MDD | M28                |                              24 |                          523 |               5 |              0.00324338  | 0.051894    | C9ORF163;GMEB2;NOL6;PNKP;SPATA2L                                     |
| high_plus_moderate | region:Ventral Striatum|direction:up_in_MDD        | M26                |                               7 |                          554 |               3 |              0.00294979  | 0.0943933   | MORN5;RPGR;TSGA10                                                    |
| all_Tier2          | region:Ventral Striatum                            | M26                |                              81 |                          554 |              10 |              0.00388688  | 0.12438     | C4ORF47;CDNF;FBXL13;MORN5;MYCBPAP;RPGR;TFDP1;TNFRSF10B;TSGA10;TTC23L |
| all_Tier2          | region:OFC BA11|direction:down_in_MDD              | M28                |                               9 |                          523 |               3 |              0.005623    | 0.179936    | ETV3;FAM43B;NOL6                                                     |
| all_Tier2          | region:Ventral Striatum                            | M6                 |                              81 |                          455 |               8 |              0.0112679   | 0.180287    | DIRAS2;GLDC;KCNV1;KIAA0513;KLHL18;MCHR2;PRKAA2;SHC3                  |
| high_plus_moderate | region:Ventral Striatum                            | M26                |                               9 |                          554 |               3 |              0.00660676  | 0.211416    | MORN5;RPGR;TSGA10                                                    |
| high_plus_moderate | region:Cingulate gyrus 25|direction:down_in_MDD    | M24                |                               2 |                          982 |               2 |              0.00664384  | 0.212603    | BACH2;VPS54                                                          |
| high_plus_moderate | region:Hippocampal Formation|direction:up_in_MDD   | M24                |                               2 |                          982 |               2 |              0.00664384  | 0.212603    | TOR2A;ZNF738                                                         |

## Interpretation

Published AHBA canonical modules provide healthy-brain coexpression context for ML-prioritized MDD candidate genes. They do not validate these genes as disease biomarkers.

The enrichment uses unique gene symbols within each candidate set, not duplicated gene-region pairs. The background is the intersection of the MDD clean-pipeline tested-gene universe and genes present in the published AHBA module assignment table.
