## Phase 3B post-reconstruction validation and versioned annotation.
## Equivalence thresholds are declared before any enrichment is computed.

options(stringsAsFactors = FALSE)
project_root <- "D:/genes/GeneDepresion"
.libPaths(c(file.path(project_root, "environment", "phase3b", "library"), .libPaths()))
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(tibble)
  library(mclust); library(clue); library(ggplot2); library(jsonlite)
  library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)
  library(ReactomePA); library(msigdbr)
})
select <- dplyr::select
slice <- dplyr::slice
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise

out_root <- file.path(project_root, "outputs", "length_gt_100_clean_pipeline", "17_WGCNA_reconstruction_validation")
dirs <- c("module_comparison", "eigengene_validation", "ml_gene_enrichment_validation", "functional_annotation", "figures", "manifests")
for (d in dirs) dir.create(file.path(out_root, d), recursive = TRUE, showWarnings = FALSE)

## Prospective, biology-blind equivalence criteria.
thresholds <- list(
  exact_ari = 1.0,
  effectively_equivalent_ari_min = 0.90,
  effectively_equivalent_median_best_jaccard_min = 0.90,
  effectively_equivalent_gene_agreement_min = 0.90,
  effectively_equivalent_median_abs_eigengene_correlation_min = 0.95
)
write_json(thresholds, file.path(out_root, "module_comparison", "predeclared_equivalence_thresholds.json"), pretty = TRUE, auto_unbox = TRUE)

regions <- c("Ventral Striatum", "Insula", "Hippocampal Formation", "Cingulate gyrus 25", "dlPFC BA8/9", "OFC BA11")
safe <- function(x) trimws(gsub("(^_|_$)", "", gsub("[^A-Za-z0-9]+", "_", x)))
std <- function(x) {
  y <- trimws(gsub(" +", " ", gsub("[^a-z0-9]+", " ", tolower(x))))
  case_when(y %in% c("cg25", "cingulate gyrus 25") ~ "cingulate gyrus 25",
            y == "dlpfc ba8 9" ~ "dlpfc ba8 9", y == "ofc ba11" ~ "ofc ba11",
            y == "ventral striatum" ~ "ventral striatum", y == "insula" ~ "insula",
            y == "hippocampal formation" ~ "hippocampal formation", TRUE ~ y)
}

broad_root <- file.path(project_root, "outputs", "length_gt_100_clean_pipeline", "12_MDD_WGCNA_modules")
front_root <- file.path(project_root, "outputs", "length_gt_100_clean_pipeline", "mdd_wgcna_modules", "frontal_subregion_wgcna")
old_assign_b <- read_csv(file.path(broad_root, "best_solution_gene_to_module_assignments.csv"), show_col_types = FALSE)
old_assign_f <- read_csv(file.path(front_root, "frontal_subregion_best_solution_gene_to_module_assignments.csv"), show_col_types = FALSE)
old_assign <- bind_rows(old_assign_b, old_assign_f) %>% mutate(gene_id = as.character(gene_id)) %>% filter(region %in% regions)

old_mes_b <- read_csv(file.path(broad_root, "best_solution_module_eigengenes.csv"), show_col_types = FALSE)
old_mes_f <- read_csv(file.path(front_root, "frontal_subregion_module_eigengenes.csv"), show_col_types = FALSE)
old_traits <- bind_rows(
  read_csv(file.path(broad_root, "module_trait_association_results.csv"), show_col_types = FALSE),
  read_csv(file.path(front_root, "frontal_subregion_module_trait_association_results.csv"), show_col_types = FALSE)
) %>% filter(region %in% regions)
old_hubs <- bind_rows(
  read_csv(file.path(broad_root, "module_hub_genes_top50.csv"), show_col_types = FALSE),
  read_csv(file.path(front_root, "frontal_subregion_hub_genes_top50.csv"), show_col_types = FALSE)
) %>% mutate(gene_id = as.character(gene_id)) %>% filter(region %in% regions)

all_overlap <- list(); all_mapping <- list(); all_summary <- list(); all_eig <- list(); all_diag <- list(); all_hub <- list()
for (region in regions) {
  rd <- file.path(out_root, "reconstruction", safe(region))
  if (!file.exists(file.path(rd, "COMPLETE.txt"))) stop("Incomplete reconstruction: ", region)
  old <- old_assign %>% filter(region == !!region) %>% select(gene_id, old_module = module)
  new <- read_csv(file.path(rd, "reconstructed_gene_to_module_assignments.csv"), show_col_types = FALSE) %>% mutate(gene_id = as.character(gene_id)) %>% select(gene_id, reconstructed_module)
  both <- inner_join(old, new, by = "gene_id")
  tab <- table(both$old_module, both$reconstructed_module)
  ov <- as.data.frame(tab, stringsAsFactors = FALSE)
  names(ov) <- c("old_module", "reconstructed_module", "overlap")
  ov <- ov %>% filter(overlap > 0)
  os <- table(both$old_module); ns <- table(both$reconstructed_module)
  ov <- ov %>% mutate(region = region, old_size = as.integer(os[old_module]), reconstructed_size = as.integer(ns[reconstructed_module]),
                      union_size = old_size + reconstructed_size - overlap, jaccard = overlap / union_size)
  all_overlap[[region]] <- ov

  ## Maximum-overlap one-to-one assignment (padded square cost matrix).
  nr <- max(nrow(tab), ncol(tab)); score <- matrix(0, nr, nr); score[seq_len(nrow(tab)), seq_len(ncol(tab))] <- tab
  sol <- solve_LSAP(score, maximum = TRUE)
  pairs <- tibble(old_index = seq_len(nrow(tab)), new_index = as.integer(sol[seq_len(nrow(tab))])) %>%
    filter(new_index <= ncol(tab)) %>% transmute(region = region, old_module = rownames(tab)[old_index], reconstructed_module = colnames(tab)[new_index], overlap = score[cbind(old_index, new_index)]) %>%
    left_join(ov %>% select(region, old_module, reconstructed_module, old_size, reconstructed_size, union_size, jaccard), by = c("region", "old_module", "reconstructed_module"))
  ## Also retain each old module's best overlap to support complete alignment.
  best_each <- ov %>% group_by(old_module) %>% arrange(desc(overlap), desc(jaccard), reconstructed_module) %>% slice(1) %>% ungroup() %>% mutate(region = region)
  mapping <- best_each %>% select(region, old_module, reconstructed_module, overlap, old_size, reconstructed_size, union_size, jaccard) %>%
    mutate(in_one_to_one_mapping = paste(old_module, reconstructed_module) %in% paste(pairs$old_module, pairs$reconstructed_module))
  all_mapping[[region]] <- mapping
  ari <- mclust::adjustedRandIndex(as.character(both$old_module), as.character(both$reconstructed_module))
  gene_agree <- sum(mapping$overlap) / nrow(both)

  old_me <- if (region %in% unique(old_mes_b$region)) old_mes_b %>% filter(region == !!region) else old_mes_f %>% filter(region == !!region)
  sid_col <- intersect(c("SampleID", "sample_id", ".sample_id"), names(old_me))[1]
  new_me <- read_csv(file.path(rd, "reconstructed_module_eigengenes.csv"), show_col_types = FALSE)
  common_samples <- intersect(old_me[[sid_col]], new_me$sample_id)
  old_me <- old_me[match(common_samples, old_me[[sid_col]]), , drop = FALSE]; new_me <- new_me[match(common_samples, new_me$sample_id), , drop = FALSE]
  eig <- mapping %>% rowwise() %>% mutate(eigengene_correlation = if (old_module %in% names(old_me) && reconstructed_module %in% names(new_me)) cor(old_me[[old_module]], new_me[[reconstructed_module]], use = "pairwise.complete.obs") else NA_real_, abs_eigengene_correlation = abs(eigengene_correlation)) %>% ungroup()
  all_eig[[region]] <- eig
  med_eig <- median(eig$abs_eigengene_correlation[eig$old_module != "MEgrey"], na.rm = TRUE)
  exact <- abs(ari - 1) < 1e-12 && gene_agree == 1
  eq <- ari >= thresholds$effectively_equivalent_ari_min && median(mapping$jaccard[mapping$old_module != "MEgrey"], na.rm = TRUE) >= thresholds$effectively_equivalent_median_best_jaccard_min && gene_agree >= thresholds$effectively_equivalent_gene_agreement_min && med_eig >= thresholds$effectively_equivalent_median_abs_eigengene_correlation_min
  status <- if (exact) "EXACT" else if (eq) "EFFECTIVELY_EQUIVALENT" else "MATERIAL_DIFFERENCE"

  pars <- read_csv(file.path(rd, "selected_solution_parameters.csv"), show_col_types = FALSE)
  old_count <- n_distinct(old$old_module[old$old_module != "MEgrey"]); new_count <- n_distinct(new$reconstructed_module[new$reconstructed_module != "MEgrey"])
  all_summary[[region]] <- tibble(region = region, reconstruction_status = status, n_common_genes = nrow(both), old_modules = old_count, reconstructed_modules = new_count,
    old_grey_genes = sum(old$old_module == "MEgrey"), reconstructed_grey_genes = sum(new$reconstructed_module == "MEgrey"),
    adjusted_rand_index = ari, median_best_jaccard = median(mapping$jaccard[mapping$old_module != "MEgrey"], na.rm = TRUE),
    gene_assignment_agreement_from_best_mapping = gene_agree, median_abs_eigengene_correlation = med_eig,
    reconstructed_power = pars$power, reconstructed_deepSplit = pars$deepSplit, reconstructed_minModuleSize = pars$minModuleSize, reconstructed_mergeCutHeight = pars$mergeCutHeight)

  nd <- read_csv(file.path(rd, "reconstructed_module_trait_associations.csv"), show_col_types = FALSE)
  all_diag[[region]] <- mapping %>% select(region, old_module, reconstructed_module, jaccard) %>%
    left_join(old_traits %>% filter(region == !!region) %>% select(region, old_module = module, old_beta = beta_diagnosis, old_p = p_diagnosis, old_q = fdr_diagnosis_BH, old_formula = model_formula), by = c("region", "old_module")) %>%
    left_join(nd %>% select(region, reconstructed_module, reconstructed_beta = beta_diagnosis, reconstructed_p = p_diagnosis, reconstructed_q = fdr_diagnosis_BH, reconstructed_formula = model_formula), by = c("region", "reconstructed_module"))
  nh <- read_csv(file.path(rd, "reconstructed_module_hub_genes_top50.csv"), show_col_types = FALSE) %>% mutate(gene_id = as.character(gene_id))
  all_hub[[region]] <- mapping %>% filter(old_module != "MEgrey", reconstructed_module != "MEgrey") %>% rowwise() %>% do({
    z <- .; a <- old_hubs %>% filter(region == z$region, module == z$old_module) %>% select(gene_id, old_kME = kME)
    b <- nh %>% filter(reconstructed_module == z$reconstructed_module) %>% select(gene_id, reconstructed_kME = kME)
    ab <- inner_join(a,b,by="gene_id")
    tibble(region=z$region, old_module=z$old_module, reconstructed_module=z$reconstructed_module, common_top50_hubs=nrow(ab), kME_spearman=if(nrow(ab)>=3) cor(ab$old_kME,ab$reconstructed_kME,method="spearman") else NA_real_)
  }) %>% ungroup()
}

overlap <- bind_rows(all_overlap); mapping <- bind_rows(all_mapping); summary <- bind_rows(all_summary)
eig <- bind_rows(all_eig); diag <- bind_rows(all_diag); hub <- bind_rows(all_hub)
write_csv(overlap, file.path(out_root, "module_comparison", "old_by_reconstructed_module_overlap_and_jaccard.csv"))
write_csv(mapping, file.path(out_root, "module_comparison", "old_to_reconstructed_module_mapping.csv"))
write_csv(summary, file.path(out_root, "module_comparison", "wgcna_reconstruction_summary.csv"))
write_csv(eig, file.path(out_root, "eigengene_validation", "aligned_module_eigengene_correlations.csv"))
write_csv(diag, file.path(out_root, "eigengene_validation", "aligned_diagnosis_association_comparison.csv"))
write_csv(hub, file.path(out_root, "module_comparison", "aligned_top50_hub_kME_comparison.csv"))

## Corrected ML sets: 30 high, 119 moderate, 149 high+moderate.
ml_path <- file.path(project_root, "outputs", "length_gt_100_clean_pipeline", "10_ML_prioritized_Tier2_corrected", "ML_prioritized_Tier2_genes_all_ranked.csv")
ml <- read_csv(ml_path, show_col_types = FALSE) %>% mutate(gene_id = as.character(gene_id), region_key = std(region_or_subregion))
stopifnot(sum(str_detect(ml$ML_priority_class, "High-confidence")) == 30L)
stopifnot(sum(str_detect(ml$ML_priority_class, "Moderate-confidence")) == 119L)
stopifnot(sum(str_detect(ml$ML_priority_class, "High|Moderate")) == 149L)

enrich_one <- function(assign, region, version) {
  assign <- assign %>% transmute(gene_id = as.character(gene_id), module = if ("module" %in% names(assign)) module else reconstructed_module)
  bg <- unique(assign$gene_id); rk <- std(region)
  sets <- list(
    global_high = ml$gene_id[str_detect(ml$ML_priority_class, "High-confidence")],
    global_high_plus_moderate = ml$gene_id[str_detect(ml$ML_priority_class, "High|Moderate")],
    same_region_high = ml$gene_id[str_detect(ml$ML_priority_class, "High-confidence") & ml$region_key == rk],
    same_region_high_plus_moderate = ml$gene_id[str_detect(ml$ML_priority_class, "High|Moderate") & ml$region_key == rk]
  )
  sets <- lapply(sets, function(x) intersect(unique(x), bg))
  rows <- lapply(unique(assign$module), function(mod) lapply(names(sets), function(sn) {
    mg <- unique(assign$gene_id[assign$module == mod]); sg <- sets[[sn]]; ov <- intersect(mg, sg)
    mat <- matrix(c(length(ov), length(mg)-length(ov), length(sg)-length(ov), length(bg)-length(mg)-length(sg)+length(ov)), nrow=2)
    tibble(region=region, assignment_version=version, module=mod, gene_set=sn, network_universe_size=length(bg), module_gene_count=length(mg), gene_set_count_in_background=length(sg), overlap_count=length(ov), overlap_genes=paste(sort(ov),collapse=";"), p_value=fisher.test(mat,alternative="greater")$p.value)
  })) %>% bind_rows()
  rows %>% group_by(region, assignment_version, gene_set) %>% mutate(fdr_BH=if_else(module=="MEgrey",NA_real_,p.adjust(if_else(module=="MEgrey",NA_real_,p_value),"BH"))) %>% ungroup()
}

ml_enrich <- list()
for (region in regions) {
  old <- old_assign %>% filter(region == !!region) %>% select(gene_id, module)
  new <- read_csv(file.path(out_root, "reconstruction", safe(region), "reconstructed_gene_to_module_assignments.csv"), show_col_types=FALSE)
  ml_enrich[[paste0(region,"_old")]] <- enrich_one(old, region, "canonical_old")
  ml_enrich[[paste0(region,"_new")]] <- enrich_one(new, region, "reconstructed")
}
ml_enrich <- bind_rows(ml_enrich)
write_csv(ml_enrich, file.path(out_root, "ml_gene_enrichment_validation", "corrected_ML_gene_enrichment_old_and_reconstructed.csv"))
ml_compare <- mapping %>% select(region, old_module, reconstructed_module, jaccard) %>%
  left_join(ml_enrich %>% filter(assignment_version=="canonical_old") %>% select(region,old_module=module,gene_set,old_overlap=overlap_count,old_p=p_value,old_q=fdr_BH),by=c("region","old_module")) %>%
  left_join(ml_enrich %>% filter(assignment_version=="reconstructed") %>% select(region,reconstructed_module=module,gene_set,new_overlap=overlap_count,new_p=p_value,new_q=fdr_BH),by=c("region","reconstructed_module","gene_set"))
write_csv(ml_compare, file.path(out_root, "ml_gene_enrichment_validation", "aligned_corrected_ML_enrichment_comparison.csv"))

## Detailed annotation is limited prospectively to the ten report-prioritized modules.
priority <- tribble(
  ~region, ~old_module,
  "Ventral Striatum","MElavenderblush3", "Ventral Striatum","MEdarkred", "Ventral Striatum","MElightgreen", "Ventral Striatum","MElightcyan1",
  "Insula","MEorange", "Hippocampal Formation","MEgreen", "dlPFC BA8/9","MEblue", "dlPFC BA8/9","MEorange", "Cingulate gyrus 25","MEwhite", "OFC BA11","MEtan"
) %>% left_join(mapping %>% select(region,old_module,reconstructed_module,jaccard),by=c("region","old_module"))
write_csv(priority, file.path(out_root, "functional_annotation", "prioritized_module_alignment.csv"))

hallmark <- msigdbr(species="Homo sapiens", collection="H")
entrez_col <- intersect(c("ncbi_gene","entrez_gene"), names(hallmark))[1]
hallmark <- hallmark %>% transmute(term=gs_name, gene=as.character(.data[[entrez_col]])) %>% filter(!is.na(gene),gene!="") %>% distinct()
annotation_rows <- list(); map_rows <- list(); coverage_rows <- list()
for (i in seq_len(nrow(priority))) {
  z <- priority[i,]; rd <- file.path(out_root,"reconstruction",safe(z$region))
  a <- read_csv(file.path(rd,"reconstructed_gene_to_module_assignments.csv"),show_col_types=FALSE) %>% mutate(gene_id=as.character(gene_id))
  bg <- unique(a$gene_id); genes <- unique(a$gene_id[a$reconstructed_module==z$reconstructed_module])
  ann <- AnnotationDbi::select(org.Hs.eg.db, keys=unique(c(bg,genes)), keytype="ENTREZID", columns=c("SYMBOL","ENTREZID")) %>% as_tibble() %>% distinct(ENTREZID,.keep_all=TRUE)
  bgm <- intersect(bg,ann$ENTREZID[!is.na(ann$SYMBOL)]); gm <- intersect(genes,bgm)
  map_rows[[i]] <- tibble(region=z$region,old_module=z$old_module,reconstructed_module=z$reconstructed_module,gene_id=genes,mapped=genes%in%gm,symbol=ann$SYMBOL[match(genes,ann$ENTREZID)])
  coverage_rows[[i]] <- tibble(region=z$region,old_module=z$old_module,reconstructed_module=z$reconstructed_module,module_genes_total=length(genes),module_genes_mapped=length(gm),module_genes_unmapped=length(setdiff(genes,gm)),mapping_fraction=length(gm)/length(genes),universe_genes=length(bg),universe_mapped=length(bgm),universe_mapping_fraction=length(bgm)/length(bg),mapping_method="org.Hs.eg.db ENTREZID->SYMBOL")
  add_res <- function(obj, db, ont=NA_character_) {
    if (is.null(obj) || nrow(as.data.frame(obj))==0) return(NULL)
    as_tibble(as.data.frame(obj)) %>% transmute(region=z$region,old_module=z$old_module,reconstructed_module=z$reconstructed_module,database=db,ontology=ont,term_id=as.character(ID),term=as.character(Description),gene_ratio=as.character(GeneRatio),background_ratio=as.character(BgRatio),raw_p=as.numeric(pvalue),BH_q=as.numeric(p.adjust),overlap_genes=as.character(geneID),overlap_count=as.integer(Count),module_genes_total=length(genes),module_genes_mapped=length(gm),universe_genes=length(bg),universe_mapped=length(bgm))
  }
  for (ont in c("BP","MF","CC")) {
    eg <- suppressMessages(enrichGO(gene=gm,universe=bgm,OrgDb=org.Hs.eg.db,keyType="ENTREZID",ont=ont,pAdjustMethod="BH",pvalueCutoff=1,qvalueCutoff=1,readable=FALSE))
    annotation_rows[[length(annotation_rows)+1]] <- add_res(eg,"GO",ont)
  }
  er <- suppressMessages(enrichPathway(gene=gm,universe=bgm,organism="human",pAdjustMethod="BH",pvalueCutoff=1,qvalueCutoff=1,readable=FALSE))
  annotation_rows[[length(annotation_rows)+1]] <- add_res(er,"Reactome","Pathway")
  eh <- suppressMessages(enricher(gene=gm,universe=bgm,TERM2GENE=hallmark,pAdjustMethod="BH",pvalueCutoff=1,qvalueCutoff=1))
  annotation_rows[[length(annotation_rows)+1]] <- add_res(eh,"MSigDB Hallmark","Hallmark")
}
mapping_quality <- bind_rows(map_rows); coverage <- bind_rows(coverage_rows); functional <- bind_rows(annotation_rows)
write_csv(mapping_quality,file.path(out_root,"functional_annotation","prioritized_module_gene_id_mapping.csv"))
write_csv(coverage,file.path(out_root,"functional_annotation","annotation_mapping_coverage_by_module.csv"))
write_csv(functional,file.path(out_root,"functional_annotation","functional_enrichment_GO_Reactome_Hallmark.csv"))
write_csv(functional %>% group_by(region,old_module,database,ontology) %>% arrange(BH_q,raw_p) %>% slice_head(n=10) %>% ungroup(),file.path(out_root,"functional_annotation","top_functional_terms_by_prioritized_module.csv"))

new_total <- n_distinct(mapping_quality$gene_id); new_mapped <- n_distinct(mapping_quality$gene_id[mapping_quality$mapped])
mapping_summary <- tibble(layer=c("historical_incomplete_mapping","phase3b_prioritized_module_mapping"),genes_total=c(4267L,new_total),genes_mapped=c(1796L,new_mapped),mapping_fraction=c(1796/4267,new_mapped/new_total),method=c("historical annotation layer without org.Hs.eg.db","org.Hs.eg.db ENTREZID->SYMBOL; Phase-3B reconstructed prioritized modules"))
write_csv(mapping_summary,file.path(out_root,"functional_annotation","historical_vs_phase3b_mapping_coverage.csv"))
write_csv(mapping_quality %>% filter(!mapped),file.path(out_root,"functional_annotation","still_unmapped_prioritized_module_genes.csv"))

## Review figures only.
p1 <- ggplot(summary,aes(x=region,y=adjusted_rand_index,fill=reconstruction_status))+geom_col()+geom_hline(yintercept=.9,lty=2)+coord_cartesian(ylim=c(0,1))+coord_flip()+theme_bw()+labs(title="Canonical vs reconstructed WGCNA partition agreement",x=NULL,y="Adjusted Rand Index")
ggsave(file.path(out_root,"figures","module_partition_ARI_by_region.png"),p1,width=8,height=5,dpi=180)
p2 <- ggplot(eig %>% filter(old_module!="MEgrey"),aes(x=region,y=abs_eigengene_correlation))+geom_boxplot()+coord_flip()+theme_bw()+labs(title="Aligned module eigengene agreement",x=NULL,y="Absolute correlation")
ggsave(file.path(out_root,"figures","aligned_eigengene_correlation_by_region.png"),p2,width=8,height=5,dpi=180)
p3 <- ggplot(mapping_summary,aes(x=layer,y=mapping_fraction,fill=layer))+geom_col()+geom_text(aes(label=scales::percent(mapping_fraction,accuracy=.1)),vjust=-.2)+coord_cartesian(ylim=c(0,1))+theme_bw()+theme(legend.position="none",axis.text.x=element_text(angle=20,hjust=1))+labs(title="Annotation mapping coverage",x=NULL,y="Mapped fraction")
ggsave(file.path(out_root,"figures","annotation_mapping_coverage_comparison.png"),p3,width=7,height=5,dpi=180)

evidence <- priority %>%
  left_join(diag %>% select(region,old_module,reconstructed_module,reconstructed_beta,reconstructed_p,reconstructed_q),by=c("region","old_module","reconstructed_module")) %>%
  left_join(ml_enrich %>% filter(assignment_version=="reconstructed",gene_set=="global_high_plus_moderate") %>% select(region,reconstructed_module=module,ml_high_moderate_overlap=overlap_count,ml_overlap_p=p_value,ml_overlap_q=fdr_BH),by=c("region","reconstructed_module")) %>%
  left_join(coverage %>% select(region,old_module,reconstructed_module,module_genes_total,module_genes_mapped,mapping_fraction),by=c("region","old_module","reconstructed_module")) %>%
  left_join(functional %>% group_by(region,old_module,reconstructed_module) %>% arrange(BH_q,raw_p) %>% summarise(top_functional_annotation=paste(head(paste0(database,": ",term," (q=",signif(BH_q,3),")"),3),collapse="; "),.groups="drop"),by=c("region","old_module","reconstructed_module"))
write_csv(evidence,file.path(out_root,"functional_annotation","prioritized_module_evidence_table.csv"))

manifest <- list(status="validation_candidate",git_baseline_requested="a483b4f6267befbb98ef49f23a27ffe6fa48ed69",git_head_at_start="df46778d1b62d92c8ba2b0a3b9002e3c3f7bdadc",equivalence_thresholds=thresholds,
  corrected_ml_counts=list(total=671,high=30,moderate=119,low=522,high_plus_moderate=149),
  enrichment=list(test="one-sided Fisher exact",background="region-specific WGCNA network gene universe",BH_family="within region and gene-set across non-grey modules"),
  annotation=list(background="mapped genes in each region-specific tested WGCNA network",id_conversion="org.Hs.eg.db ENTREZID to SYMBOL",resources=c("GO BP","GO MF","GO CC","Reactome","MSigDB Hallmark")),
  sessionInfo=capture.output(sessionInfo()))
write_json(manifest,file.path(out_root,"manifests","phase3b_comparison_annotation_manifest.json"),pretty=TRUE,auto_unbox=TRUE,null="null")
cat("Phase-3B comparisons and annotations complete.\n")
print(summary)
cat("Maximum reconstructed Diagnosis BH q:", min(diag$reconstructed_q,na.rm=TRUE),"(minimum q; no-survival if >=0.05)\n")
