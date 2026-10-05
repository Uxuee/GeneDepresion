#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("RobustRankAggreg", quietly = TRUE)) {
    stop("Required canonical package RobustRankAggreg is unavailable; refusing fallback aggregation.")
  }
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else file.path("outputs", "length_gt_100_clean_pipeline", "18_Tier2_ML_RRA_prioritization")
rank_path <- file.path(input_dir, "Tier2_fold_gene_rankings.csv")
performance_path <- file.path(input_dir, "Tier2_outer_performance.csv")
unit_path <- file.path(input_dir, "Tier2_analysis_unit_manifest.csv")
tier2_path <- file.path("outputs", "length_gt_100_clean_pipeline", "07_Tier2_candidate_gene_sets", "Tier2_candidate_genes_all_regions.csv")
for (path in c(rank_path, performance_path, unit_path, tier2_path)) {
  if (!file.exists(path)) stop("Missing required input: ", path)
}

rankings <- read.csv(rank_path, stringsAsFactors = FALSE, check.names = FALSE)
performance <- read.csv(performance_path, stringsAsFactors = FALSE, check.names = FALSE)
units <- read.csv(unit_path, stringsAsFactors = FALSE, check.names = FALSE)
tier2 <- read.csv(tier2_path, stringsAsFactors = FALSE, check.names = FALSE)
ranking_flag <- tolower(trimws(as.character(performance$contributes_gene_ranking)))
if (any(!ranking_flag %in% c("true", "false"))) {
  stop("contributes_gene_ranking contains values other than True/False")
}
performance$contributes_gene_ranking <- ranking_flag == "true"
rankings$gene_id <- as.character(rankings$gene_id)
tier2$gene_id <- as.character(tier2$gene_id)
tier2$region <- tier2$region_or_subregion
tier2$region[tier2$region == "Cingulate gyrus 25"] <- "CG25"

required <- c("fit_id", "region", "model", "gene_id", "within_fit_rank")
if (!all(required %in% colnames(rankings))) stop("Fold-ranking schema is incomplete")
if (any(duplicated(rankings[c("fit_id", "gene_id")]))) stop("Duplicate gene within a fit ranking")
models_expected <- c("logistic_l1", "elastic_net", "logistic_l2", "linear_svm", "random_forest")
if (any(performance$model == "dummy_most_frequent" & performance$contributes_gene_ranking)) {
  stop("Dummy classifier is incorrectly marked as ranking evidence")
}

level1_all <- list()
level2_all <- list()
accounting_all <- list()

for (region in units$region) {
  candidates <- tier2[tier2$region == region, c("gene_id", "gene_symbol"), drop = FALSE]
  n_universe <- nrow(candidates)
  expected_n <- units$n_tier2_candidates[units$region == region]
  if (length(expected_n) != 1 || n_universe != expected_n) stop("Tier2 universe mismatch for ", region)
  region_rankings <- rankings[rankings$region == region, , drop = FALSE]
  model_tables <- list()

  for (model in models_expected) {
    model_data <- region_rankings[region_rankings$model == model, , drop = FALSE]
    fit_ids <- unique(as.character(performance$fit_id[performance$region == region & performance$model == model]))
    if (length(fit_ids) != 25) stop(region, " / ", model, ": expected 25 fits, found ", length(fit_ids))
    gene_lists <- lapply(fit_ids, function(fit_id) {
      part <- model_data[model_data$fit_id == fit_id, , drop = FALSE]
      part <- part[order(part$within_fit_rank, part$gene_id), , drop = FALSE]
      as.character(part$gene_id)
    })
    if (any(vapply(gene_lists, anyDuplicated, integer(1)) > 0)) stop("Non-unique RRA input list")
    if (any(!unlist(gene_lists, use.names = FALSE) %in% candidates$gene_id)) stop("Non-Tier2 gene in RRA input")
    names(gene_lists) <- fit_ids
    aggregated <- RobustRankAggreg::aggregateRanks(
      glist = gene_lists, N = n_universe, method = "RRA", full = FALSE, exact = FALSE
    )
    colnames(aggregated)[1:2] <- c("gene_id", "rra_score")
    aggregated$gene_id <- as.character(aggregated$gene_id)
    aggregated <- aggregated[order(aggregated$rra_score, aggregated$gene_id), , drop = FALSE]
    aggregated$model_rank <- seq_len(nrow(aggregated))
    appearances <- table(model_data$gene_id)
    aggregated$outer_list_appearance_count <- as.integer(appearances[aggregated$gene_id])
    aggregated$outer_list_appearance_proportion <- aggregated$outer_list_appearance_count / 25
    aggregated$region <- region
    aggregated$model <- model
    aggregated$n_outer_lists <- 25
    aggregated$tier2_universe_size <- n_universe
    aggregated$gene_symbol <- candidates$gene_symbol[match(aggregated$gene_id, candidates$gene_id)]
    level1_all[[paste(region, model, sep = "::")]] <- aggregated
    model_tables[[model]] <- aggregated
  }

  model_lists <- lapply(models_expected, function(model) as.character(model_tables[[model]]$gene_id))
  names(model_lists) <- models_expected
  consensus <- RobustRankAggreg::aggregateRanks(
    glist = model_lists, N = n_universe, method = "RRA", full = FALSE, exact = FALSE
  )
  colnames(consensus)[1:2] <- c("gene_id", "consensus_rra_score")
  consensus$gene_id <- as.character(consensus$gene_id)
  consensus <- consensus[order(consensus$consensus_rra_score, consensus$gene_id), , drop = FALSE]
  consensus$consensus_rank <- seq_len(nrow(consensus))
  consensus$region <- region
  consensus$gene_symbol <- candidates$gene_symbol[match(consensus$gene_id, candidates$gene_id)]
  consensus$model_appearance_count <- rowSums(sapply(model_tables, function(tab) consensus$gene_id %in% tab$gene_id))
  consensus$tier2_universe_size <- n_universe
  for (model in models_expected) {
    lookup_rank <- setNames(model_tables[[model]]$model_rank, model_tables[[model]]$gene_id)
    lookup_score <- setNames(model_tables[[model]]$rra_score, model_tables[[model]]$gene_id)
    consensus[[paste0(model, "_rank")]] <- as.integer(lookup_rank[consensus$gene_id])
    consensus[[paste0(model, "_rra_score")]] <- as.numeric(lookup_score[consensus$gene_id])
  }
  level2_all[[region]] <- consensus
  any_evidence <- candidates$gene_id %in% unique(region_rankings$gene_id)
  accounting_all[[region]] <- data.frame(
    region = region,
    total_tier2_candidates = n_universe,
    candidates_with_any_fold_ranking_evidence = sum(any_evidence),
    candidates_absent_from_every_fold_ranking = sum(!any_evidence),
    candidates_in_consensus_rra = nrow(consensus),
    stringsAsFactors = FALSE
  )
}

level1 <- do.call(rbind, level1_all)
level2 <- do.call(rbind, level2_all)
accounting <- do.call(rbind, accounting_all)
rownames(level1) <- NULL
rownames(level2) <- NULL
rownames(accounting) <- NULL
if (sum(accounting$total_tier2_candidates) != 671) stop("Tier2 accounting does not total 671")
write.csv(level1, file.path(input_dir, "Tier2_model_level_rra_rankings.csv"), row.names = FALSE, na = "")
write.csv(level2, file.path(input_dir, "Tier2_consensus_rra_rankings.csv"), row.names = FALSE, na = "")
write.csv(accounting, file.path(input_dir, "Tier2_candidate_accounting.csv"), row.names = FALSE, na = "")
metadata_lines <- c(
  "package=RobustRankAggreg",
  paste0("package_version=", as.character(utils::packageVersion("RobustRankAggreg"))),
  "function=RobustRankAggreg::aggregateRanks",
  "method=RRA",
  "exact=FALSE",
  "level1=25 outer-fit Tier2 rankings per region/model",
  "level2=five model-level Tier2 rankings per region",
  "dummy_contributes_gene_evidence=FALSE",
  "predictive_metric_weighting=FALSE",
  paste0("r_version=", R.version.string)
)
writeLines(metadata_lines, file.path(input_dir, "Tier2_rra_run_metadata.txt"))
cat("Tier2 RRA complete:", nrow(level1), "model-level rows and", nrow(level2), "consensus rows\n")
