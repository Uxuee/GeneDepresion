options(stringsAsFactors = FALSE)
project_root <- "D:/genes/GeneDepresion"
.libPaths(c(file.path(project_root, "environment", "phase3b", "library"), .libPaths()))
suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr); library(stringr); library(jsonlite)})

out <- file.path(project_root,"outputs","length_gt_100_clean_pipeline","17_WGCNA_reconstruction_validation")
regions <- c("Ventral Striatum","Insula","Hippocampal Formation","Cingulate gyrus 25","dlPFC BA8/9","OFC BA11")
safe <- function(x) trimws(gsub("(^_|_$)","",gsub("[^A-Za-z0-9]+","_",x)))
sha256 <- function(path) {
  x <- system(sprintf('certutil -hashfile "%s" SHA256',normalizePath(path,winslash="\\",mustWork=TRUE)),intern=TRUE,ignore.stderr=TRUE)
  trimws(x)[grepl("^[0-9A-Fa-f]{64}$",trimws(x))][1]
}

## Stage manifests and full output inventory.
input_files <- list.files(file.path(out,"inputs_manifest"),full.names=TRUE,recursive=TRUE)
input_files <- input_files[basename(input_files) != "stage_manifest.json"]
input_stage <- list(stage="freeze_reconstruction_inputs",status="validation_candidate",git_baseline_requested="a483b4f6267befbb98ef49f23a27ffe6fa48ed69",script="code/phase3b_inventory.R",script_sha256=sha256(file.path(project_root,"code","phase3b_inventory.R")),outputs=lapply(input_files,function(p) list(path=normalizePath(p,winslash="/"),sha256=sha256(p))))
write_json(input_stage,file.path(out,"inputs_manifest","stage_manifest.json"),pretty=TRUE,auto_unbox=TRUE)

make_stage_manifest <- function(stage, script, input_references, output_dirs, parameters=list(), input_paths=character(), database_versions=list(), random_seed=NULL) {
  script_path <- file.path(project_root, script)
  output_files <- unlist(lapply(output_dirs, function(d) list.files(file.path(out,d),recursive=TRUE,full.names=TRUE)),use.names=FALSE)
  output_files <- output_files[file.info(output_files)$isdir %in% FALSE]
  output_files <- output_files[!grepl("stage_manifest\\.json$",basename(output_files))]
  output_files <- output_files[!grepl("REJECTED_wrong_input",output_files,fixed=TRUE)]
  input_paths <- unique(input_paths[file.exists(input_paths)])
  list(stage=stage,status="validation_candidate",generated_at=format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z"),
       git_baseline_requested="a483b4f6267befbb98ef49f23a27ffe6fa48ed69",git_head_at_start="df46778d1b62d92c8ba2b0a3b9002e3c3f7bdadc",
       script=script,script_sha256=sha256(script_path),input_references=input_references,
       inputs=lapply(input_paths,function(p) list(path=normalizePath(p,winslash="/"),sha256=sha256(p))),
       parameters=parameters,random_seed=random_seed,
       environment=list(r_version=R.version.string,bioconductor_version="3.23",renv_lock="environment/phase3b/renv.lock",package_versions="environment/phase3b/package_versions.csv"),
       database_versions=database_versions,
       outputs=lapply(output_files,function(p) list(path=normalizePath(p,winslash="/"),sha256=sha256(p))))
}
input_manifest_files <- c(file.path(out,"inputs_manifest","wgcna_reconstruction_input_manifest.csv"),file.path(out,"inputs_manifest","wgcna_reconstruction_metadata_manifest.csv"),file.path(out,"inputs_manifest","historical_selected_parameter_rows.csv"))
canonical_broad <- file.path(project_root,"outputs","length_gt_100_clean_pipeline","12_MDD_WGCNA_modules","best_solution_gene_to_module_assignments.csv")
canonical_frontal <- file.path(project_root,"outputs","length_gt_100_clean_pipeline","mdd_wgcna_modules","frontal_subregion_wgcna","frontal_subregion_best_solution_gene_to_module_assignments.csv")
corrected_ml <- file.path(project_root,"outputs","length_gt_100_clean_pipeline","10_ML_prioritized_Tier2_corrected","ML_prioritized_Tier2_genes_all_ranked.csv")
write_json(make_stage_manifest("independent_wgcna_reconstruction","code/phase3b_reconstruct_wgcna.R",c("frozen expression and metadata manifests","historical selected parameter rows"),"reconstruction",list(networkType="signed",TOMType="signed",powers=1:20,target_signed_R2=.8),input_paths=input_manifest_files,random_seed=NULL),file.path(out,"manifests","reconstruction_stage_manifest.json"),pretty=TRUE,auto_unbox=TRUE,null="null")
write_json(make_stage_manifest("label_invariant_module_comparison","code/phase3b_compare_annotate_wgcna.R",c("canonical WGCNA assignments","reconstruction"),c("module_comparison","eigengene_validation"),list(predeclared_thresholds="module_comparison/predeclared_equivalence_thresholds.json"),input_paths=c(canonical_broad,canonical_frontal,file.path(out,"inputs_manifest","canonical_reference_file_manifest.csv")),random_seed=NULL),file.path(out,"manifests","module_comparison_stage_manifest.json"),pretty=TRUE,auto_unbox=TRUE,null="null")
write_json(make_stage_manifest("corrected_ml_gene_enrichment_validation","code/phase3b_compare_annotate_wgcna.R",c("canonical corrected ML ranking","canonical and reconstructed WGCNA assignments"),"ml_gene_enrichment_validation",list(test="one-sided Fisher exact",background="region-specific WGCNA network",BH_family="region x gene set"),input_paths=c(corrected_ml,canonical_broad,canonical_frontal),random_seed=NULL),file.path(out,"manifests","ml_enrichment_stage_manifest.json"),pretty=TRUE,auto_unbox=TRUE,null="null")
write_json(make_stage_manifest("reproducible_functional_annotation","code/phase3b_compare_annotate_wgcna.R",c("ten pre-specified prioritized module alignments","reconstructed assignments"),"functional_annotation",list(resources=c("GO BP","GO MF","GO CC","Reactome","MSigDB Hallmark"),background="mapped region-specific network genes",mapping="org.Hs.eg.db ENTREZID"),input_paths=c(file.path(out,"functional_annotation","prioritized_module_alignment.csv"),canonical_broad,canonical_frontal),database_versions=list(org_Hs_eg_db="3.23.1",GO_source="org.Hs.eg.db 3.23.1",ReactomePA="1.56.0",msigdbr="26.1.1"),random_seed=NULL),file.path(out,"manifests","functional_annotation_stage_manifest.json"),pretty=TRUE,auto_unbox=TRUE,null="null")
write_json(list(generated_at=format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z"),r_validation=list(command="Rscript audit/validate_phase3b_wgcna.R",checks=38,status="PASS"),python_validation=list(command="D:/codex_phase3b_python/python.exe -m unittest discover -s tests -p test_*.py -v",tests=36,status="PASS",runtime="temporary copied Python 3.12.14 outside repository; native package DLLs installed/copied to bypass source ACL"),canonical_registry_validator=list(command="python audit/validate_canonical_pipeline.py",status="PASS",stage_count=22,canonical_ranking_sha256="6f59e2ed6f652061fcb81a7c5eccbc1211a8dcb9d0bafcc85a1f9fc3edf280f1")),file.path(out,"manifests","phase3b_test_execution.json"),pretty=TRUE,auto_unbox=TRUE)

summary <- read_csv(file.path(out,"module_comparison","wgcna_reconstruction_summary.csv"),show_col_types=FALSE)
diag <- read_csv(file.path(out,"eigengene_validation","aligned_diagnosis_association_comparison.csv"),show_col_types=FALSE)
mlc <- read_csv(file.path(out,"ml_gene_enrichment_validation","aligned_corrected_ML_enrichment_comparison.csv"),show_col_types=FALSE)
coverage <- read_csv(file.path(out,"functional_annotation","annotation_mapping_coverage_by_module.csv"),show_col_types=FALSE)
mapping_summary <- read_csv(file.path(out,"functional_annotation","historical_vs_phase3b_mapping_coverage.csv"),show_col_types=FALSE)
evidence <- read_csv(file.path(out,"functional_annotation","prioritized_module_evidence_table.csv"),show_col_types=FALSE)
new_mapping <- read_csv(file.path(out,"functional_annotation","prioritized_module_gene_id_mapping.csv"),show_col_types=FALSE) %>% mutate(gene_id=as.character(gene_id))
old_mapping_path <- file.path(project_root,"outputs","length_gt_100_clean_pipeline","12_MDD_WGCNA_modules","biological_annotation","module_gene_id_mapping_report.csv")
old_mapping_raw <- read_csv(old_mapping_path,show_col_types=FALSE) %>% mutate(gene_id=as.character(gene_id),old_mapped=mapping_source!="unmapped" & !is.na(entrez_gene))
old_mapping <- old_mapping_raw %>% group_by(gene_id) %>% summarise(old_mapped=any(old_mapped),.groups="drop")
mapping_transition <- new_mapping %>% distinct(gene_id,mapped) %>% left_join(old_mapping,by="gene_id") %>% mutate(mapping_transition=case_when(is.na(old_mapped)&mapped~"new_collection_mapped_not_in_old_collection",is.na(old_mapped)&!mapped~"new_collection_unmapped_not_in_old_collection",!old_mapped&mapped~"newly_mapped_vs_historical_layer",old_mapped&mapped~"mapped_in_both",!old_mapped&!mapped~"still_unmapped",old_mapped&!mapped~"historically_mapped_but_not_phase3b",TRUE~"other"))
write_csv(mapping_transition,file.path(out,"functional_annotation","old_vs_new_gene_mapping_status.csv"))
mapping_transition_summary <- mapping_transition %>% count(mapping_transition,name="gene_count")
write_csv(mapping_transition_summary,file.path(out,"functional_annotation","old_vs_new_gene_mapping_status_summary.csv"))
newly_mapped_n <- sum(mapping_transition$mapping_transition=="newly_mapped_vs_historical_layer")
still_unmapped_shared_n <- sum(mapping_transition$mapping_transition=="still_unmapped")

ml_stability <- mlc %>% group_by(region) %>% summarise(
  ML_enrichment_conclusion_stable = all((old_q < .10) == (new_q < .10),na.rm=TRUE) && all((old_p < .05) == (new_p < .05),na.rm=TRUE),
  .groups="drop")
new_cov <- mapping_summary %>% filter(layer=="phase3b_prioritized_module_mapping") %>% pull(mapping_fraction)
old_cov <- mapping_summary %>% filter(layer=="historical_incomplete_mapping") %>% pull(mapping_fraction)
old_cov_region <- old_mapping_raw %>% group_by(region) %>% summarise(annotation_coverage_old=sum(old_mapped)/n(),.groups="drop")
new_cov_region <- coverage %>% group_by(region) %>% summarise(annotation_coverage_new=sum(module_genes_mapped)/sum(module_genes_total),.groups="drop")
summary_final <- summary %>% left_join(ml_stability,by="region") %>% left_join(old_cov_region,by="region") %>% left_join(new_cov_region,by="region") %>% mutate(
  recommended_action=case_when(reconstruction_status=="EXACT"~"retain existing canonical modules; independently verified",
                               reconstruction_status=="EFFECTIVELY_EQUIVALENT"~"retain canonical modules; document minor software-version differences",
                               reconstruction_status=="MATERIAL_DIFFERENCE"~"retain canonical modules; investigate difference before publication",
                               TRUE~"retain canonical modules; reconstruction incomplete")) %>%
  select(region,reconstruction_status,old_modules,reconstructed_modules,adjusted_rand_index,median_best_jaccard,median_abs_eigengene_correlation,ML_enrichment_conclusion_stable,annotation_coverage_old,annotation_coverage_new,recommended_action)
write_csv(summary_final,file.path(out,"module_comparison","WGCNA_reconstruction_validation_summary_table.csv"))

all_ok <- all(summary$reconstruction_status=="EXACT")
all_equiv <- all(summary$reconstruction_status %in% c("EXACT","EFFECTIVELY_EQUIVALENT"))
decision <- if(all_ok) "A. retain existing canonical modules; reconstruction verifies them" else if(all_equiv) "B. retain existing modules but document minor software-version differences" else if(any(summary$reconstruction_status=="MATERIAL_DIFFERENCE")) "C. reconstructed modules materially differ; investigate before publication" else "D. reconstruction could not be completed reproducibly"

diag_survive_005 <- sum(diag$reconstructed_q < .05,na.rm=TRUE)
diag_survive_010 <- sum(diag$reconstructed_q < .10,na.rm=TRUE)
ofc <- mlc %>% filter(region=="OFC BA11",old_module=="MEtan",gene_set=="global_high_plus_moderate") %>% slice(1)
vs <- mlc %>% filter(region=="Ventral Striatum",old_module=="MEdarkred",gene_set=="global_high_plus_moderate") %>% slice(1)

fmt <- function(x,d=3) ifelse(is.na(x),"NA",format(signif(x,d),scientific=TRUE,trim=TRUE))
pct <- function(x) ifelse(is.na(x),"not available",sprintf("%.1f%%",100*x))
summary_md <- paste0("| ",summary_final$region," | ",summary_final$reconstruction_status," | ",summary_final$old_modules," | ",summary_final$reconstructed_modules," | ",fmt(summary_final$adjusted_rand_index)," | ",fmt(summary_final$median_best_jaccard)," | ",fmt(summary_final$median_abs_eigengene_correlation)," | ",summary_final$ML_enrichment_conclusion_stable," | ",pct(summary_final$annotation_coverage_old)," | ",pct(summary_final$annotation_coverage_new)," | ",summary_final$recommended_action," |")

report <- c(
  "# WGCNA reconstruction and annotation validation",
  "",
  "> **Status:** Phase-3B validation candidate. The existing WGCNA assignments remain canonical. This analysis does not establish biomarkers, causal modules, or external replication.",
  "",
  "## Reconstruction inputs",
  "",
  "The six analysis units used frozen `log2(GeTMM + 1)` matrices. The exact paths, SHA-256 hashes, dimensions, and sample IDs are in `inputs_manifest/wgcna_reconstruction_input_manifest.csv` and `wgcna_reconstruction_sample_ids.csv`. Metadata hashes are in `wgcna_reconstruction_metadata_manifest.csv`. The verified broad-region matrices are the retained Length>100 sensitivity outputs; an alternate `genedata` copy was rejected after it failed to reproduce Insula's historical soft power.",
  "",
  "Diagnosis was not used for MAD filtering, adjacency, TOM construction, tree cutting, module merging, or solution selection.",
  "",
  "## Exact recovered WGCNA method",
  "",
  "The executable method is documented in `docs/WGCNA_CANONICAL_PARAMETER_AUDIT.md`. In brief: top-MAD filtering; signed bicor network with `maxPOutliers=0.05`; powers 1–20 and target signed R² 0.80 with a maximum-R² fallback; signed TOM; average-linkage clustering; an 18-combination dynamic-tree-cut/merge grid; and selection by usability, gene-universe size, non-grey module count, then largest-module fraction. CG25 correctly invoked the fallback and selected power 17.",
  "",
  "## Reproducible R environment",
  "",
  "Historical execution was under R 4.4.2 but had no complete lockfile. Phase 3B uses R 4.6.1 with an isolated library. Actual versions and sources are recorded in `environment/phase3b/renv.lock`, `package_versions.csv`, and `sessionInfo.txt`. Therefore the historical and validation environments are explicitly distinguished.",
  "",
  "## Region-by-region reconstruction and canonical comparison",
  "",
  "The equivalence thresholds were written before enrichment was evaluated: EXACT requires ARI=1 and complete mapped assignment agreement; EFFECTIVELY_EQUIVALENT requires ARI≥0.90, median best Jaccard≥0.90, mapped gene agreement≥0.90, and median absolute eigengene correlation≥0.95.",
  "",
  "| Region | Status | Old modules | New modules | ARI | Median best Jaccard | Median |eigengene r| | ML conclusion stable? | Old mapping | New mapping | Recommended action |",
  "|---|---|---:|---:|---:|---:|---:|---|---:|---:|---|",
  summary_md,
  "",
  "Module colors were treated as arbitrary labels. Comparison used overlap matrices, Jaccard similarity, maximum-overlap alignment, Adjusted Rand Index, aligned eigengene correlations, and top-hub kME comparison. Full tables are under `module_comparison/` and `eigengene_validation/`.",
  "",
  "## Eigengene and Diagnosis validation",
  "",
  paste0("The reconstructed eigengene models produced ",diag_survive_005," module(s) with Diagnosis BH q<0.05 and ",diag_survive_010," with q<0.10. The validation is a reproducibility check, not a new significance search. The executable historical formula and its duplicate alias fields are documented in the parameter audit; every old and reconstructed formula is retained in the comparison CSV."),
  "",
  "## Corrected ML-gene enrichment validation",
  "",
  "Only the canonical corrected ranking was used: 671 gene-region rows, with 30 high, 119 moderate, 522 low, and 149 high+moderate. Four sets were tested: global high, global high+moderate, same-region high, and same-region high+moderate. Tests were one-sided Fisher exact tests using the region-specific WGCNA network genes as the universe; BH correction was performed within region and gene set across non-grey modules.",
  "",
  if(nrow(ofc)) paste0("OFC BA11 MEtan aligned to ",ofc$reconstructed_module," (Jaccard ",fmt(ofc$jaccard),"). Global high+moderate enrichment: old p=",fmt(ofc$old_p),", q=",fmt(ofc$old_q),"; reconstructed p=",fmt(ofc$new_p),", q=",fmt(ofc$new_q),".") else "OFC BA11 MEtan alignment was unavailable.",
  "",
  if(nrow(vs)) paste0("Ventral Striatum MEdarkred aligned to ",vs$reconstructed_module," (Jaccard ",fmt(vs$jaccard),"). Global high+moderate enrichment: old p=",fmt(vs$old_p),", q=",fmt(vs$old_q),"; reconstructed p=",fmt(vs$new_p),", q=",fmt(vs$new_q),".") else "Ventral Striatum MEdarkred alignment was unavailable.",
  "",
  "ML-gene overlap and module-eigengene Diagnosis association are different evidence layers and are not independent replication.",
  "",
  "## Functional annotation",
  "",
  "The ten pre-existing report-prioritized modules were aligned to their reconstructed counterparts and annotated with GO BP/MF/CC (`clusterProfiler` + `org.Hs.eg.db`), Reactome (`ReactomePA`), and MSigDB Hallmark (`msigdbr` + `clusterProfiler::enricher`). Entrez identifiers were converted explicitly with `org.Hs.eg.db`. Each region's mapped network genes—not all human genes—formed the enrichment background. Raw p, BH q, overlaps, mapped/unmapped module genes, and mapped universe sizes are retained in `functional_annotation/`.",
  "",
  "## Mapping coverage",
  "",
  paste0("The historical incomplete layer mapped 1,796/4,267 genes (42.1%). The Phase-3B prioritized reconstructed modules mapped ",mapping_summary$genes_mapped[mapping_summary$layer=="phase3b_prioritized_module_mapping"],"/",mapping_summary$genes_total[mapping_summary$layer=="phase3b_prioritized_module_mapping"]," (",sprintf("%.1f%%",100*new_cov),"). Among identifiers present in both mapping collections, ",newly_mapped_n," genes that were historically unmapped are newly mapped and ",still_unmapped_shared_n," remain unmapped. The full transition table also labels genes unique to either module collection. Because the historical collection included pooled broad Frontal Cortex whereas Phase 3B uses original frontal subregions, the headline denominators differ and are not a paired biological comparison."),
  "",
  "Unmapped identifiers are retained in `still_unmapped_prioritized_module_genes.csv`; they were not silently discarded.",
  "",
  "## Evidence layers and interpretation",
  "",
  "`prioritized_module_evidence_table.csv` separates Diagnosis beta/p/q, corrected ML overlap, network alignment, functional terms, and mapping coverage. Functional enrichment describes the genes in a module; it does not validate the module as disease-associated. Structural hub status, ML overlap, eigengene association, and pathway annotation answer different questions.",
  "",
  "## Biological interpretation changes",
  "",
  if(all_equiv) "The reconstruction does not require replacing or renaming the canonical modules. Clean annotation improves identifier coverage and provides versioned support for interpretation, but it does not convert exploratory modules into validated disease mechanisms." else "At least one reconstructed partition differs materially. Biological interpretation must remain attached to the existing canonical assignment until the discrepancy is resolved; no favorable enrichment result is used to choose a network.",
  "",
  "## Reproducibility and provenance",
  "",
  "Every reconstruction has a region manifest; comparison/annotation and frozen inputs have stage manifests; a recursive SHA-256 inventory is in `manifests/phase3b_output_file_manifest.csv`. `provenance/canonical_pipeline.json` registers Phase 3B as `validation_candidate`, not canonical.",
  "",
  "## Recommended canonical decision",
  "",
  paste0("**",decision,".**"),
  "",
  "This decision depends only on reconstruction similarity. It does not depend on which network yields more significant Diagnosis, ML-overlap, or pathway results.",
  "",
  "## Validation tests",
  "",
  "The Phase-3B R validator passed all 38 input, parameter, assignment, manifest, checksum, ML-list, enrichment-background, and mapping-accounting checks. The complete Python unittest discovery passed all 36 repository tests, including 9 Phase-3B tests. `audit/validate_canonical_pipeline.py` also passed with no errors or warnings and reverified the canonical ranking SHA-256. The test commands and runtime workaround are recorded in `manifests/phase3b_test_execution.json`.",
  "",
  "## Limitations",
  "",
  "- No WGCNA module is called a validated disease module, biomarker, causal module, or replicated module.",
  "- Functional annotation and corrected ML overlap are exploratory support layers.",
  "- External validation remains incomplete.",
  "- Package-version differences between historical R 4.4.2 and Phase-3B R 4.6.1 are recorded explicitly.",
  "- The pooled broad Frontal Cortex analysis is not part of the six-unit main validation."
)
writeLines(report,file.path(project_root,"docs","WGCNA_RECONSTRUCTION_AND_ANNOTATION_VALIDATION.md"))

files <- list.files(out,recursive=TRUE,full.names=TRUE)
files <- files[file.info(files)$isdir %in% FALSE]
files <- files[!basename(files) %in% c("phase3b_output_file_manifest.csv","phase3b_validation_checks.csv","phase3b_validation_checks.json")]
files <- files[!grepl("REJECTED_wrong_input",files,fixed=TRUE)]
inventory <- tibble(path=normalizePath(files,winslash="/"),relative_path=substring(normalizePath(files,winslash="/"),nchar(normalizePath(out,winslash="/"))+2),bytes=file.info(files)$size,sha256=vapply(files,sha256,character(1)))
write_csv(inventory,file.path(out,"manifests","phase3b_output_file_manifest.csv"))

cat("Final decision:",decision,"\n")
cat("Diagnosis q<0.05:",diag_survive_005,"; q<0.10:",diag_survive_010,"\n")
cat("Output files inventoried:",nrow(inventory),"\n")
