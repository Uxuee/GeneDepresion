options(stringsAsFactors = FALSE)
root <- "D:/genes/GeneDepresion"
.libPaths(c(file.path(root,"environment","phase3b","library"),.libPaths()))
suppressPackageStartupMessages({library(readr);library(dplyr);library(jsonlite);library(stringr)})
out <- file.path(root,"outputs","length_gt_100_clean_pipeline","17_WGCNA_reconstruction_validation")
regions <- c("Ventral Striatum","Insula","Hippocampal Formation","Cingulate gyrus 25","dlPFC BA8/9","OFC BA11")
safe <- function(x) trimws(gsub("(^_|_$)","",gsub("[^A-Za-z0-9]+","_",x)))
sha <- function(path){x<-system(sprintf('certutil -hashfile "%s" SHA256',normalizePath(path,winslash="\\",mustWork=TRUE)),intern=TRUE);trimws(x)[grepl("^[0-9A-Fa-f]{64}$",trimws(x))][1]}
checks <- list(); add <- function(name,ok,detail="") checks[[length(checks)+1]] <<- tibble(check=name,passed=isTRUE(ok),detail=detail)

im <- read_csv(file.path(out,"inputs_manifest","wgcna_reconstruction_input_manifest.csv"),show_col_types=FALSE)
add("six frozen expression inputs",nrow(im)==6,paste(im$region,collapse="; "))
add("frozen expression SHA-256",all(vapply(seq_len(nrow(im)),function(i) sha(im$path[i])==im$sha256[i],logical(1))))

src <- readLines(file.path(root,"code","phase3b_reconstruct_wgcna.R"),warn=FALSE)
add("construction manifest declares no Diagnosis",any(grepl("construction_used_diagnosis = FALSE",src,fixed=TRUE)))
add("no Tier2 or corrected ML construction input",!any(grepl("Tier2|ML_prioritized",src)))

hist <- read_csv(file.path(out,"inputs_manifest","historical_selected_parameter_rows.csv"),show_col_types=FALSE)
for(region in regions){
  rd<-file.path(out,"reconstruction",safe(region)); add(paste(region,"complete"),file.exists(file.path(rd,"COMPLETE.txt")))
  p<-read_csv(file.path(rd,"selected_solution_parameters.csv"),show_col_types=FALSE); h<-hist%>%filter(region==!!region)%>%slice(1)
  add(paste(region,"selected parameters"),p$power==h$power && p$deepSplit==h$deepSplit && p$minModuleSize==h$minModuleSize && p$mergeCutHeight==h$mergeCutHeight && p$corType==h$corType)
  a<-read_csv(file.path(rd,"reconstructed_gene_to_module_assignments.csv"),show_col_types=FALSE)
  add(paste(region,"unique assignment per gene"),nrow(a)==n_distinct(a$gene_id) && all(startsWith(a$reconstructed_module,"ME")))
  m<-read_json(file.path(rd,"reconstruction_manifest.json"),simplifyVector=TRUE)
  add(paste(region,"noncanonical manifest"),m$status=="validation_candidate" && identical(m$construction_used_diagnosis,FALSE))
}

cmp<-read_csv(file.path(out,"module_comparison","wgcna_reconstruction_summary.csv"),show_col_types=FALSE)
add("six label-invariant comparisons",nrow(cmp)==6)
add("ARI bounded",all(cmp$adjusted_rand_index>=-1 & cmp$adjusted_rand_index<=1))
add("Jaccard bounded",all(cmp$median_best_jaccard>=0 & cmp$median_best_jaccard<=1))

ml<-read_csv(file.path(root,"outputs","length_gt_100_clean_pipeline","10_ML_prioritized_Tier2_corrected","ML_prioritized_Tier2_genes_all_ranked.csv"),show_col_types=FALSE)
ml_high <- sum(str_detect(ml$ML_priority_class,"High-confidence"))
ml_moderate <- sum(str_detect(ml$ML_priority_class,"Moderate-confidence"))
counts<-c(total=nrow(ml),high=ml_high,moderate=ml_moderate,low=nrow(ml)-ml_high-ml_moderate)
add("corrected ML 671/30/119/522",identical(as.integer(counts),c(671L,30L,119L,522L)),paste(counts,collapse="/"))
enr<-read_csv(file.path(out,"ml_gene_enrichment_validation","corrected_ML_gene_enrichment_old_and_reconstructed.csv"),show_col_types=FALSE)
add("enrichment accounting",all(enr$overlap_count<=pmin(enr$module_gene_count,enr$gene_set_count_in_background)) && all(enr$module_gene_count<=enr$network_universe_size))

cov<-read_csv(file.path(out,"functional_annotation","annotation_mapping_coverage_by_module.csv"),show_col_types=FALSE)
add("ten prioritized modules annotated",nrow(cov)==10)
add("mapping accounting",all(cov$module_genes_total==cov$module_genes_mapped+cov$module_genes_unmapped))

output_manifest_path <- file.path(out,"manifests","phase3b_output_file_manifest.csv")
output_manifest <- read_csv(output_manifest_path,show_col_types=FALSE)
output_files_exist <- file.exists(output_manifest$path)
output_hashes_match <- output_files_exist & vapply(seq_len(nrow(output_manifest)),function(i) {
  if (!output_files_exist[i]) return(FALSE)
  identical(sha(output_manifest$path[i]),output_manifest$sha256[i])
},logical(1))
add("output manifest has records",nrow(output_manifest)>0,paste(nrow(output_manifest),"files"))
add("output manifest files exist",all(output_files_exist),paste(output_manifest$relative_path[!output_files_exist],collapse="; "))
add("output manifest SHA-256",all(output_hashes_match),paste(output_manifest$relative_path[!output_hashes_match],collapse="; "))

result<-bind_rows(checks)
dir.create(file.path(out,"manifests"),recursive=TRUE,showWarnings=FALSE)
write_csv(result,file.path(out,"manifests","phase3b_validation_checks.csv"))
write_json(list(generated_at=format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z"),passed=all(result$passed),checks=result),file.path(out,"manifests","phase3b_validation_checks.json"),pretty=TRUE,auto_unbox=TRUE)
print(result,n=Inf)
if(!all(result$passed)) stop("Phase-3B validation failed: ",paste(result$check[!result$passed],collapse="; "))
cat("All",nrow(result),"Phase-3B validation checks passed.\n")
