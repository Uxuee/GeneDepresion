#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("RobustRankAggreg", quietly = TRUE)) {
    stop("Required canonical package RobustRankAggreg is unavailable; refusing fallback aggregation.")
  }
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else file.path("outputs", "length_gt_100_clean_pipeline", "17_ML_RRA_consensus")
rank_path <- file.path(input_dir, "fold_gene_rankings.csv")
performance_path <- file.path(input_dir, "outer_fold_performance.csv")
manifest_path <- file.path(input_dir, "run_manifest.json")
unit_manifest_path <- file.path(input_dir, "analysis_unit_manifest.csv")
if (!file.exists(rank_path)) stop("Missing fold rankings: ", rank_path)
if (!file.exists(performance_path)) stop("Missing performance ledger: ", performance_path)
if (!file.exists(manifest_path)) stop("Missing run manifest: ", manifest_path)
if (!file.exists(unit_manifest_path)) stop("Missing analysis-unit manifest: ", unit_manifest_path)

rankings <- read.csv(rank_path, stringsAsFactors = FALSE, check.names = FALSE)
performance <- read.csv(performance_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("fit_id", "region", "model", "gene_id", "within_fit_rank")
if (!all(required %in% colnames(rankings))) stop("Fold-ranking schema is incomplete")
if (any(duplicated(rankings[c("fit_id", "gene_id")]))) stop("Duplicate gene within a fit ranking")
rankings$gene_id <- as.character(rankings$gene_id)
unit_manifest <- read.csv(unit_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)

models_expected <- c("logistic_l1", "elastic_net", "logistic_l2", "linear_svm", "random_forest")
regions <- unique(rankings$region)
level1_all <- list()
level2_all <- list()

for (region in regions) {
  region_data <- rankings[rankings$region == region, , drop = FALSE]
  # N is the actual expression universe, never the smaller union of selected lists.
  unit_row <- unit_manifest[unit_manifest$region == region, , drop = FALSE]
  if (nrow(unit_row) != 1) stop("Could not recover unique n_genes for ", region)
  n_universe <- as.integer(unit_row$n_genes[[1]])

  model_tables <- list()
  for (model in models_expected) {
    model_data <- region_data[region_data$model == model, , drop = FALSE]
    fit_ids <- performance$fit_id[performance$region == region & performance$model == model]
    fit_ids <- unique(as.character(fit_ids))
    if (length(fit_ids) != 25) stop(region, " / ", model, ": expected 25 outer lists, found ", length(fit_ids))
    gene_lists <- lapply(fit_ids, function(fit_id) {
      part <- model_data[model_data$fit_id == fit_id, , drop = FALSE]
      part <- part[order(part$within_fit_rank, part$gene_id), , drop = FALSE]
      as.character(part$gene_id)
    })
    if (any(vapply(gene_lists, anyDuplicated, integer(1)) > 0)) stop("Non-unique RRA list")
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
    aggregated$outer_list_appearance_proportion <- aggregated$outer_list_appearance_count / length(fit_ids)
    aggregated$region <- region
    aggregated$model <- model
    aggregated$n_outer_lists <- length(fit_ids)
    aggregated$gene_universe_size <- n_universe
    level1_all[[paste(region, model, sep = "::")]] <- aggregated
    model_tables[[model]] <- aggregated
  }

  if (!identical(sort(names(model_tables)), sort(models_expected))) stop("All five model families are required")
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
  consensus$model_appearance_count <- rowSums(sapply(model_tables, function(tab) consensus$gene_id %in% tab$gene_id))
  consensus$gene_universe_size <- n_universe
  for (model in models_expected) {
    lookup_rank <- setNames(model_tables[[model]]$model_rank, model_tables[[model]]$gene_id)
    lookup_score <- setNames(model_tables[[model]]$rra_score, model_tables[[model]]$gene_id)
    consensus[[paste0(model, "_rank")]] <- as.integer(lookup_rank[consensus$gene_id])
    consensus[[paste0(model, "_rra_score")]] <- as.numeric(lookup_score[consensus$gene_id])
  }
  level2_all[[region]] <- consensus
}

level1 <- do.call(rbind, level1_all)
level2 <- do.call(rbind, level2_all)
rownames(level1) <- NULL
rownames(level2) <- NULL
write.csv(level1, file.path(input_dir, "model_level_rra_rankings.csv"), row.names = FALSE, na = "")
write.csv(level2, file.path(input_dir, "consensus_rra_rankings.csv"), row.names = FALSE, na = "")

metadata <- list(
  package = "RobustRankAggreg",
  package_version = as.character(utils::packageVersion("RobustRankAggreg")),
  function_name = "RobustRankAggreg::aggregateRanks",
  method = "RRA",
  exact = FALSE,
  level1 = "25 outer-fold gene lists per region/model",
  level2 = "five complete model-level RRA lists per region",
  predictive_metric_weighting = FALSE,
  r_version = R.version.string
)
writeLines(capture.output(dput(metadata)), file.path(input_dir, "rra_run_metadata.dput"))
cat("RRA complete:", nrow(level1), "model-level rows and", nrow(level2), "consensus rows\n")
