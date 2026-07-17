# Stable Core Genes for Module-Level ML #2

This folder defines stable gene sets for future module-level ML #2 features/eigengenes. No classification was run in this step.

## Selection Principle

Stable genes were selected from WGCNA/module information only: non-grey best-solution module membership, presence in the normalized expression matrix, and module membership/hubness measured by kME or hub rank. Diagnosis p-values were not used for selection. Tier2 and ML-prioritized flags were added only after selection as annotation columns.

## Stability Method

Run-level gene-to-module assignment files found: 0.
Run-level gene-to-module assignment files were not saved for each parameter-grid run; assignment_stability_frequency is unavailable and stable cores use best-solution kME/hub rank.
Because run-level assignments were unavailable, stable_core_genes use |kME| >= 0.6; if a module had too few such genes, the top hub-ranked genes were added until up to 25 genes were available.

## Gene Set Types

- all_module_genes: every non-grey gene assigned to the prioritized module in the best solution.
- stable_core_genes: genes selected for future module score/eigengene features using kME/hubness stability rules.
- top_hub_genes: top 50 genes by abs(kME)/hub rank, or all genes if the module is smaller.

## Module Summary

| Region | Module | Module size | Stable core genes | Top hub genes | Median stable-core abs(kME) | ML #2 readiness |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Ventral Striatum | MElavenderblush3 | 35 | 35 | 35 | 0.796 | enough stable core genes |
| Ventral Striatum | MElightcyan1 | 586 | 50 | 50 | 0.900 | enough stable core genes |
| Ventral Striatum | MEdarkred | 227 | 50 | 50 | 0.835 | enough stable core genes |
| Ventral Striatum | MElightgreen | 299 | 50 | 50 | 0.900 | enough stable core genes |
| Insula | MEorange | 53 | 44 | 50 | 0.705 | enough stable core genes |
| Hippocampal Formation | MEgreen | 857 | 50 | 50 | 0.884 | enough stable core genes |
| Frontal Cortex | MEdarkgreen | 479 | 50 | 50 | 0.857 | enough stable core genes |

## How These Genes Can Be Used Later

For ML #2, each module can be represented by module scores computed from stable_core_genes, such as an eigengene/first principal component, mean scaled expression, or robust median scaled expression. The all_module_genes set remains available for sensitivity analyses, while top_hub_genes can be used for compact interpretable feature variants.

## Limitations

- Run-level assignment stability could not be estimated because parameter-grid gene assignment tables were not saved.
- kME values are available only for the saved top hub table, so genes outside the top hub table are retained in all_module_genes but are not selected into stable cores unless they appear in the hub/kME table.
- Stable core gene selection does not imply diagnosis association and should not be interpreted as biomarker validation.
- Future classification should be run in a separate step with strict cross-validation and no leakage from diagnosis labels into gene selection.

## Outputs

- stable_core_genes_by_module.csv
- stable_core_gene_summary_by_module.csv
- stable_core_gene_manifest.json
