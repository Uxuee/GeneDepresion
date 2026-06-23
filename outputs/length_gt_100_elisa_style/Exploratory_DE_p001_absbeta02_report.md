# Exploratory DE Filter: p < 0.01 and |beta| > 0.2

## Definition

Exploratory DE candidates are genes with:

`p_value_diagnosis < 0.01` and `abs(beta_diagnosis) > 0.2`

The script used `beta_diagnosis` as the diagnosis beta/effect-size column and `p_value_diagnosis` as the diagnosis p-value column. `fdr_BH` was retained as a reference column.

## Important Interpretation

These are nominal/exploratory DE candidates, not final FDR-significant DE genes.

FDR is not the effect size. `beta_diagnosis` is the effect size, `p_value_diagnosis` is the nominal statistical evidence, and FDR is the multiple-testing-corrected p-value.

## Summary By Region

| brain_region          |   genes_tested |   exploratory_DE_genes |   up_in_MDD |   down_in_MDD |   median_abs_beta |   max_abs_beta |   min_p_value |   number_with_FDR_lt_005 |   number_with_FDR_lt_010 |
|:----------------------|---------------:|-----------------------:|------------:|--------------:|------------------:|---------------:|--------------:|-------------------------:|-------------------------:|
| Frontal_Cortex        |          23941 |                     99 |          16 |            83 |          0.241918 |       0.528167 |   7.83571e-05 |                        0 |                        0 |
| Insula                |          23090 |                     92 |          54 |            38 |          0.263585 |       0.564171 |   0.00010202  |                        0 |                        0 |
| Ventral_Striatum      |          23917 |                    167 |         116 |            51 |          0.26459  |       0.539819 |   7.96917e-06 |                        0 |                        0 |
| Hippocampal_Formation |          24371 |                     70 |          51 |            19 |          0.273618 |       0.625618 |   6.03267e-05 |                        0 |                        0 |

## Top 20 By P-Value

### Frontal_Cortex

| rank | gene | beta_diagnosis | p_value_diagnosis | abs_beta_diagnosis | direction |
|---:|---|---:|---:|---:|---|
| 1 | C19orf73 | 0.2008 | 7.836e-05 | 0.2008 | up_in_MDD |
| 2 | QPRT | 0.2255 | 0.000126 | 0.2255 | up_in_MDD |
| 3 | DDT | 0.2107 | 0.0002076 | 0.2107 | up_in_MDD |
| 4 | LOC101929347 | -0.2644 | 0.0002567 | 0.2644 | down_in_MDD |
| 5 | LOC107984015 | -0.252 | 0.0003227 | 0.252 | down_in_MDD |
| 6 | SURF1 | 0.2523 | 0.0003518 | 0.2523 | up_in_MDD |
| 7 | LINC00710 | -0.3327 | 0.0005088 | 0.3327 | down_in_MDD |
| 8 | LINC00473 | -0.3995 | 0.0005223 | 0.3995 | down_in_MDD |
| 9 | LOC105371630 | -0.2305 | 0.0005508 | 0.2305 | down_in_MDD |
| 10 | LOC107985130 | -0.2099 | 0.0005727 | 0.2099 | down_in_MDD |
| 11 | PKIA-AS1 | -0.2247 | 0.0005896 | 0.2247 | down_in_MDD |
| 12 | LOC112268044 | -0.2702 | 0.0006431 | 0.2702 | down_in_MDD |
| 13 | LOC105373894 | -0.3061 | 0.0008923 | 0.3061 | down_in_MDD |
| 14 | LOC101928929 | -0.2046 | 0.0008928 | 0.2046 | down_in_MDD |
| 15 | LOC107986526 | -0.3199 | 0.0009152 | 0.3199 | down_in_MDD |
| 16 | LOC105375097 | -0.2589 | 0.001303 | 0.2589 | down_in_MDD |
| 17 | PPP1R14B-AS1 | 0.2267 | 0.001346 | 0.2267 | up_in_MDD |
| 18 | LOC100506474 | -0.2144 | 0.001466 | 0.2144 | down_in_MDD |
| 19 | LOC105370972 | -0.3161 | 0.001689 | 0.3161 | down_in_MDD |
| 20 | OR10A4 | -0.2526 | 0.001781 | 0.2526 | down_in_MDD |

### Hippocampal_Formation

| rank | gene | beta_diagnosis | p_value_diagnosis | abs_beta_diagnosis | direction |
|---:|---|---:|---:|---:|---|
| 1 | NOL6 | -0.2022 | 6.033e-05 | 0.2022 | down_in_MDD |
| 2 | HLA-F-AS1 | 0.2518 | 0.0001173 | 0.2518 | up_in_MDD |
| 3 | ZNF663P | 0.3435 | 0.0001823 | 0.3435 | up_in_MDD |
| 4 | LINC02137 | -0.2516 | 0.0003982 | 0.2516 | down_in_MDD |
| 5 | LOC105377663 | 0.227 | 0.0008446 | 0.227 | up_in_MDD |
| 6 | ARR3 | 0.2263 | 0.00127 | 0.2263 | up_in_MDD |
| 7 | MAGIX | -0.215 | 0.00133 | 0.215 | down_in_MDD |
| 8 | SNORA25 | 0.4375 | 0.001469 | 0.4375 | up_in_MDD |
| 9 | LINC01010 | 0.2805 | 0.001693 | 0.2805 | up_in_MDD |
| 10 | SCRG1 | 0.2694 | 0.001803 | 0.2694 | up_in_MDD |
| 11 | MRTFA-AS1 | -0.2557 | 0.001932 | 0.2557 | down_in_MDD |
| 12 | SAMMSON | 0.3937 | 0.002052 | 0.3937 | up_in_MDD |
| 13 | MGST1 | 0.2844 | 0.002335 | 0.2844 | up_in_MDD |
| 14 | LYPLA2 | -0.2239 | 0.002495 | 0.2239 | down_in_MDD |
| 15 | LOC105377623 | 0.302 | 0.002673 | 0.302 | up_in_MDD |
| 16 | DUSP5 | -0.5602 | 0.002771 | 0.5602 | down_in_MDD |
| 17 | TMEM161B-DT | 0.2923 | 0.002885 | 0.2923 | up_in_MDD |
| 18 | IFI35 | 0.2311 | 0.002973 | 0.2311 | up_in_MDD |
| 19 | APOBEC2 | 0.272 | 0.003344 | 0.272 | up_in_MDD |
| 20 | LOC107987057 | 0.2051 | 0.003554 | 0.2051 | up_in_MDD |

### Insula

| rank | gene | beta_diagnosis | p_value_diagnosis | abs_beta_diagnosis | direction |
|---:|---|---:|---:|---:|---|
| 1 | UST | -0.3162 | 0.000102 | 0.3162 | down_in_MDD |
| 2 | POU6F2-AS2 | -0.3481 | 0.0001189 | 0.3481 | down_in_MDD |
| 3 | DPP4 | -0.3623 | 0.0001719 | 0.3623 | down_in_MDD |
| 4 | BMS1P2-AGAP9 | 0.3437 | 0.0002088 | 0.3437 | up_in_MDD |
| 5 | BMS1P2 | 0.3115 | 0.0002102 | 0.3115 | up_in_MDD |
| 6 | TMEM161B-DT | 0.3442 | 0.0003056 | 0.3442 | up_in_MDD |
| 7 | LAMA3 | 0.2101 | 0.0003634 | 0.2101 | up_in_MDD |
| 8 | LGI4 | 0.2198 | 0.0003639 | 0.2198 | up_in_MDD |
| 9 | TEK | 0.2659 | 0.0004748 | 0.2659 | up_in_MDD |
| 10 | GABRG3 | -0.2781 | 0.0005329 | 0.2781 | down_in_MDD |
| 11 | LOC101928268 | 0.2823 | 0.0006427 | 0.2823 | up_in_MDD |
| 12 | KCNA1 | -0.3179 | 0.0008301 | 0.3179 | down_in_MDD |
| 13 | LOC107986284 | -0.3022 | 0.001471 | 0.3022 | down_in_MDD |
| 14 | VWC2L | -0.2681 | 0.001504 | 0.2681 | down_in_MDD |
| 15 | SLC2A14 | -0.2382 | 0.001563 | 0.2382 | down_in_MDD |
| 16 | LOC112268239 | 0.368 | 0.001583 | 0.368 | up_in_MDD |
| 17 | FAM43B | -0.2628 | 0.001615 | 0.2628 | down_in_MDD |
| 18 | ARL14EP-DT | 0.5642 | 0.001901 | 0.5642 | up_in_MDD |
| 19 | NBPF11 | 0.207 | 0.002123 | 0.207 | up_in_MDD |
| 20 | DCLK3 | 0.3075 | 0.002237 | 0.3075 | up_in_MDD |

### Ventral_Striatum

| rank | gene | beta_diagnosis | p_value_diagnosis | abs_beta_diagnosis | direction |
|---:|---|---:|---:|---:|---|
| 1 | TSGA10 | 0.234 | 7.969e-06 | 0.234 | up_in_MDD |
| 2 | BMS1P1 | 0.2235 | 1.214e-05 | 0.2235 | up_in_MDD |
| 3 | LOC107987147 | 0.3302 | 2.385e-05 | 0.3302 | up_in_MDD |
| 4 | MLLT11 | -0.2553 | 8.944e-05 | 0.2553 | down_in_MDD |
| 5 | NPTX2 | -0.5389 | 0.0001219 | 0.5389 | down_in_MDD |
| 6 | BMS1P2 | 0.3096 | 0.0001977 | 0.3096 | up_in_MDD |
| 7 | MORN5 | 0.5003 | 0.0002906 | 0.5003 | up_in_MDD |
| 8 | PSMG3 | -0.2371 | 0.0002958 | 0.2371 | down_in_MDD |
| 9 | COQ7-DT | 0.3262 | 0.0004843 | 0.3262 | up_in_MDD |
| 10 | LOC112268145 | 0.2147 | 0.0004878 | 0.2147 | up_in_MDD |
| 11 | NIPAL1 | -0.2139 | 0.0005264 | 0.2139 | down_in_MDD |
| 12 | MOSPD3 | 0.2126 | 0.0005745 | 0.2126 | up_in_MDD |
| 13 | POTEF | 0.2138 | 0.0005896 | 0.2138 | up_in_MDD |
| 14 | SNAP25-AS1 | -0.216 | 0.0006242 | 0.216 | down_in_MDD |
| 15 | GPR22 | -0.3081 | 0.0006509 | 0.3081 | down_in_MDD |
| 16 | SHC3 | -0.2144 | 0.0006604 | 0.2144 | down_in_MDD |
| 17 | LOC105371067 | 0.2505 | 0.0006834 | 0.2505 | up_in_MDD |
| 18 | SPACDR | 0.334 | 0.000699 | 0.334 | up_in_MDD |
| 19 | ANKRD34C | -0.2107 | 0.0007945 | 0.2107 | down_in_MDD |
| 20 | AGAP9 | 0.2275 | 0.0007995 | 0.2275 | up_in_MDD |


## Top 20 By Absolute Beta

### Frontal_Cortex

| rank | gene | beta_diagnosis | p_value_diagnosis | abs_beta_diagnosis | direction |
|---:|---|---:|---:|---:|---|
| 1 | LOC107986426 | -0.5282 | 0.009071 | 0.5282 | down_in_MDD |
| 2 | LINC00499 | 0.4408 | 0.00842 | 0.4408 | up_in_MDD |
| 3 | RNVU1-1 | 0.4382 | 0.002645 | 0.4382 | up_in_MDD |
| 4 | NR4A3 | -0.408 | 0.006901 | 0.408 | down_in_MDD |
| 5 | LINC00473 | -0.3995 | 0.0005223 | 0.3995 | down_in_MDD |
| 6 | LOC101928937 | -0.3877 | 0.001984 | 0.3877 | down_in_MDD |
| 7 | HTR2A-AS1 | -0.3734 | 0.003666 | 0.3734 | down_in_MDD |
| 8 | ARL14EP-DT | 0.3702 | 0.003839 | 0.3702 | up_in_MDD |
| 9 | LOC105369673 | -0.3409 | 0.005618 | 0.3409 | down_in_MDD |
| 10 | LINC00710 | -0.3327 | 0.0005088 | 0.3327 | down_in_MDD |
| 11 | LOC105374934 | -0.3283 | 0.004268 | 0.3283 | down_in_MDD |
| 12 | LOC107986526 | -0.3199 | 0.0009152 | 0.3199 | down_in_MDD |
| 13 | LOC105370972 | -0.3161 | 0.001689 | 0.3161 | down_in_MDD |
| 14 | PRKCE-AS1 | -0.3143 | 0.00627 | 0.3143 | down_in_MDD |
| 15 | HECW1-IT1 | -0.313 | 0.003772 | 0.313 | down_in_MDD |
| 16 | LOC102724906 | -0.3089 | 0.003878 | 0.3089 | down_in_MDD |
| 17 | LOC105373894 | -0.3061 | 0.0008923 | 0.3061 | down_in_MDD |
| 18 | LOC105375999 | 0.3031 | 0.009785 | 0.3031 | up_in_MDD |
| 19 | ARL17A | -0.2927 | 0.005323 | 0.2927 | down_in_MDD |
| 20 | LOC105378115 | -0.2927 | 0.00406 | 0.2927 | down_in_MDD |

### Hippocampal_Formation

| rank | gene | beta_diagnosis | p_value_diagnosis | abs_beta_diagnosis | direction |
|---:|---|---:|---:|---:|---|
| 1 | FAUP4 | -0.6256 | 0.008525 | 0.6256 | down_in_MDD |
| 2 | ARL14EP-DT | 0.5619 | 0.004751 | 0.5619 | up_in_MDD |
| 3 | DUSP5 | -0.5602 | 0.002771 | 0.5602 | down_in_MDD |
| 4 | NPTX2 | -0.5178 | 0.007327 | 0.5178 | down_in_MDD |
| 5 | LOC105376547 | 0.4829 | 0.005359 | 0.4829 | up_in_MDD |
| 6 | LINC01480 | 0.4783 | 0.006564 | 0.4783 | up_in_MDD |
| 7 | SNORA25 | 0.4375 | 0.001469 | 0.4375 | up_in_MDD |
| 8 | LINC01992 | 0.4241 | 0.004479 | 0.4241 | up_in_MDD |
| 9 | LOC107986295 | 0.4075 | 0.00471 | 0.4075 | up_in_MDD |
| 10 | GPC5-AS2 | 0.4011 | 0.009813 | 0.4011 | up_in_MDD |
| 11 | SAMMSON | 0.3937 | 0.002052 | 0.3937 | up_in_MDD |
| 12 | LINC02440 | -0.3933 | 0.005345 | 0.3933 | down_in_MDD |
| 13 | LOC105377432 | 0.3812 | 0.007249 | 0.3812 | up_in_MDD |
| 14 | DUSP5-DT | -0.3572 | 0.004263 | 0.3572 | down_in_MDD |
| 15 | LOC105370516 | 0.3572 | 0.004956 | 0.3572 | up_in_MDD |
| 16 | DDX43 | 0.3513 | 0.005086 | 0.3513 | up_in_MDD |
| 17 | LOC105378995 | 0.3506 | 0.00838 | 0.3506 | up_in_MDD |
| 18 | LOC112268213 | -0.3452 | 0.006705 | 0.3452 | down_in_MDD |
| 19 | ZNF663P | 0.3435 | 0.0001823 | 0.3435 | up_in_MDD |
| 20 | DNAL4 | -0.3326 | 0.005286 | 0.3326 | down_in_MDD |

### Insula

| rank | gene | beta_diagnosis | p_value_diagnosis | abs_beta_diagnosis | direction |
|---:|---|---:|---:|---:|---|
| 1 | ARL14EP-DT | 0.5642 | 0.001901 | 0.5642 | up_in_MDD |
| 2 | CCN2 | -0.5222 | 0.008736 | 0.5222 | down_in_MDD |
| 3 | LOC101929710 | 0.4937 | 0.00548 | 0.4937 | up_in_MDD |
| 4 | CNGB1 | 0.4455 | 0.003031 | 0.4455 | up_in_MDD |
| 5 | MAGI2-AS2 | -0.4137 | 0.007247 | 0.4137 | down_in_MDD |
| 6 | CLUL1 | -0.4013 | 0.008806 | 0.4013 | down_in_MDD |
| 7 | LOC105376380 | 0.3807 | 0.003981 | 0.3807 | up_in_MDD |
| 8 | LINC00343 | -0.3688 | 0.005767 | 0.3688 | down_in_MDD |
| 9 | LOC112268239 | 0.368 | 0.001583 | 0.368 | up_in_MDD |
| 10 | LOC101930048 | 0.3678 | 0.005865 | 0.3678 | up_in_MDD |
| 11 | DPP4 | -0.3623 | 0.0001719 | 0.3623 | down_in_MDD |
| 12 | POU6F2-AS2 | -0.3481 | 0.0001189 | 0.3481 | down_in_MDD |
| 13 | LOC101928092 | 0.347 | 0.00302 | 0.347 | up_in_MDD |
| 14 | TMEM161B-DT | 0.3442 | 0.0003056 | 0.3442 | up_in_MDD |
| 15 | BMS1P2-AGAP9 | 0.3437 | 0.0002088 | 0.3437 | up_in_MDD |
| 16 | TPM2 | 0.3425 | 0.007612 | 0.3425 | up_in_MDD |
| 17 | PNMA5 | 0.3348 | 0.009143 | 0.3348 | up_in_MDD |
| 18 | GSTM2P1 | 0.3271 | 0.003214 | 0.3271 | up_in_MDD |
| 19 | CALB2 | 0.3248 | 0.007897 | 0.3248 | up_in_MDD |
| 20 | GPR3 | -0.3244 | 0.009721 | 0.3244 | down_in_MDD |

### Ventral_Striatum

| rank | gene | beta_diagnosis | p_value_diagnosis | abs_beta_diagnosis | direction |
|---:|---|---:|---:|---:|---|
| 1 | ARL14EP-DT | 0.5398 | 0.007797 | 0.5398 | up_in_MDD |
| 2 | NPTX2 | -0.5389 | 0.0001219 | 0.5389 | down_in_MDD |
| 3 | MORN5 | 0.5003 | 0.0002906 | 0.5003 | up_in_MDD |
| 4 | LOC105377106 | 0.4747 | 0.001356 | 0.4747 | up_in_MDD |
| 5 | DUSP6 | -0.4675 | 0.005964 | 0.4675 | down_in_MDD |
| 6 | RGS2 | -0.4465 | 0.006186 | 0.4465 | down_in_MDD |
| 7 | LOC105376520 | 0.4301 | 0.001034 | 0.4301 | up_in_MDD |
| 8 | LOC112268307 | 0.4243 | 0.001943 | 0.4243 | up_in_MDD |
| 9 | STIM1-AS1 | 0.4183 | 0.002228 | 0.4183 | up_in_MDD |
| 10 | LINC03002 | 0.4122 | 0.002677 | 0.4122 | up_in_MDD |
| 11 | ARMC3 | 0.4074 | 0.008335 | 0.4074 | up_in_MDD |
| 12 | LOC105376547 | 0.3943 | 0.00913 | 0.3943 | up_in_MDD |
| 13 | RARB-AS1 | 0.3902 | 0.003565 | 0.3902 | up_in_MDD |
| 14 | SOX10 | 0.3896 | 0.003392 | 0.3896 | up_in_MDD |
| 15 | RAB27B | -0.3692 | 0.005206 | 0.3692 | down_in_MDD |
| 16 | FSHB | 0.3676 | 0.001312 | 0.3676 | up_in_MDD |
| 17 | LOC105370763 | 0.3665 | 0.009578 | 0.3665 | up_in_MDD |
| 18 | SLC47A2 | 0.3625 | 0.007524 | 0.3625 | up_in_MDD |
| 19 | CDC42EP1 | 0.3593 | 0.007569 | 0.3593 | up_in_MDD |
| 20 | PI16 | 0.3576 | 0.003328 | 0.3576 | up_in_MDD |


## Recommended Next Step

Use these exploratory DE candidates and ranked signatures for enrichment analysis, ML feature selection, and AHBA module mapping.
