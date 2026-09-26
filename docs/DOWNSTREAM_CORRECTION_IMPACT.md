# Downstream Impact of the ML-Priority Correction

This is a dependency audit, not a rerun.

| Artifact/workflow | Old list used? | Needs regeneration? | Reason |
|---|---|---|---|
| AHBA preprocessing/enrichment scripts 01/02/03/05 | Yes | Yes | Point to legacy Step 10 |
| AHBA report | Yes | Yes | Reports 35 high and 130 moderate |
| Broad/frontal WGCNA construction | ML lists used for annotation, not construction | Partial only | Preserve modules; refresh ML enrichment/labels |
| WGCNA biological annotation and figures | Yes | Yes | Manifests/scripts reference legacy ranking |
| Frontal top-15,000 sensitivity | Yes for ML overlay | Partial only | Preserve modules; refresh overlay if reported |
| Stable module-core extraction | Yes | Yes | Input is legacy final ranking |
| Module network/follow-up plots | Yes | Yes | Annotation scripts load legacy genes |
| Module-level ML2 | Potentially | Review, then rerun affected exports | Module scores may be independent; gene overlays may not be |
| Collected report images | Yes | Yes | Packaging includes legacy plots |
| Preliminary manuscript package | Yes | Yes | Inherits pre-correction figures/tables |
| Integrated/supervisor reports | Yes | Yes | Several name legacy directory/counts as current |
| Historical correction/code audits | Yes | No | Retain as clearly historical provenance |

## Direct code references queued for later update

- `code/Step6_MDD_WGCNA_modules.R`
- `code/Step6_Frontal_Subregion_WGCNA_modules.R`
- `code/Step7_Annotate_MDD_WGCNA_Modules.R`
- `code/Step8_Define_WGCNA_Module_ML2_Stable_Core_Genes.py`
- `code/ahba_modules/01_preprocess_ahba_microarray.py`
- `code/ahba_modules/02_build_ahba_modules.R`
- `code/ahba_modules/03_test_ahba_module_enrichment.py`
- `code/ahba_modules/05_test_published_ahba_module_enrichment.py`
- `code/create_frontal_subregion_wgcna_followups.py`
- `code/create_mdd_wgcna_full_module_networks.py`
- `code/copy_original_manuscript_figures.py`
- `code/length_gt_100_clean_pipeline/pipeline_lib.py`

They were not changed here because this task stops before downstream regeneration.
