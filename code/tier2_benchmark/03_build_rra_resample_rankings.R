#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("RobustRankAggreg", quietly = TRUE)) {
    stop("RobustRankAggreg is required; no fallback aggregation is permitted")
  }
})

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args) >= 1) args[[1]] else file.path("outputs", "length_gt_100_clean_pipeline", "19_Tier2_ML_benchmark")
primary <- file.path("outputs", "length_gt_100_clean_pipeline", "18_Tier2_ML_RRA_prioritization")
rankings <- read.csv(file.path(primary, "Tier2_fold_gene_rankings.csv"), stringsAsFactors = FALSE, check.names = FALSE)
performance <- read.csv(file.path(primary, "Tier2_outer_performance.csv"), stringsAsFactors = FALSE, check.names = FALSE)
tier2 <- read.csv(file.path("outputs", "length_gt_100_clean_pipeline", "07_Tier2_candidate_gene_sets", "Tier2_candidate_genes_all_regions.csv"), stringsAsFactors = FALSE)
tier2$region <- tier2$region_or_subregion
tier2$region[tier2$region == "Cingulate gyrus 25"] <- "CG25"
tier2$gene_id <- as.character(tier2$gene_id)
rankings$gene_id <- as.character(rankings$gene_id)
models <- c("logistic_l1", "elastic_net", "logistic_l2", "linear_svm", "random_forest")
all_results <- list()

for (region in unique(tier2$region)) {
  candidates <- tier2[tier2$region == region, c("gene_id", "gene_symbol")]
  n_universe <- nrow(candidates)
  region_perf <- performance[performance$region == region & performance$model %in% models, ]
  splits <- unique(region_perf[c("repeat", "fold")])
  if (nrow(splits) != 25) stop(region, ": expected 25 outer splits")
  for (row_index in seq_len(nrow(splits))) {
    repeat_id <- splits[["repeat"]][[row_index]]
    fold_id <- splits$fold[[row_index]]
    lists <- lapply(models, function(model) {
      fit <- region_perf[region_perf[["repeat"]] == repeat_id & region_perf$fold == fold_id & region_perf$model == model, ]
      if (nrow(fit) != 1) stop("Missing unique fit for resample RRA")
      part <- rankings[rankings$fit_id == fit$fit_id[[1]], ]
      part <- part[order(part$within_fit_rank, part$gene_id), ]
      as.character(part$gene_id)
    })
    names(lists) <- models
    aggregated <- RobustRankAggreg::aggregateRanks(
      glist = lists, N = n_universe, method = "RRA", full = FALSE, exact = FALSE
    )
    colnames(aggregated)[1:2] <- c("gene_id", "resample_rra_score")
    aggregated$gene_id <- as.character(aggregated$gene_id)
    aggregated <- aggregated[order(aggregated$resample_rra_score, aggregated$gene_id), ]
    aggregated$resample_rank <- seq_len(nrow(aggregated))
    aggregated$region <- region
    aggregated[["repeat"]] <- repeat_id
    aggregated$fold <- fold_id
    aggregated$outer_seed <- unique(region_perf$outer_seed[region_perf[["repeat"]] == repeat_id & region_perf$fold == fold_id])[[1]]
    aggregated$gene_symbol <- candidates$gene_symbol[match(aggregated$gene_id, candidates$gene_id)]
    all_results[[paste(region, repeat_id, fold_id, sep = "::")]] <- aggregated
  }
}

result <- do.call(rbind, all_results)
rownames(result) <- NULL
if (nrow(result) != 671 * 25) stop("Unexpected split-level RRA row count")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(result, file.path(output, "RRA_resample_consensus_rankings.csv"), row.names = FALSE, na = "")
cat("Built", nrow(result), "split-level RRA rows across 150 region-resamples\n")
