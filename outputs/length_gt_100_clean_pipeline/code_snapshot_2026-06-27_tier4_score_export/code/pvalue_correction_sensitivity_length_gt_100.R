## pvalue_correction_sensitivity_length_gt_100.R
##
## Purpose:
##   Compare plausible p-value correction methods and correction scopes for the
##   Length > 100 standardized-effect-size MDD analysis.
##
## Important:
##   This is a sensitivity/audit branch only. It does not rerun LMMs and it
##   does not overwrite the current Tier 1/Tier 2/Tier 3 outputs.

rm(list = ls())
options(stringsAsFactors = FALSE)

has_qvalue <- requireNamespace("qvalue", quietly = TRUE)
has_ihw <- requireNamespace("IHW", quietly = TRUE)

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(script_path)) {
  script_path <- normalizePath("code/pvalue_correction_sensitivity_length_gt_100.R", mustWork = FALSE)
}
repo_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

root <- file.path(repo_dir, "outputs", "length_gt_100_elisa_style")
out_dir <- file.path(root, "pvalue_correction_sensitivity")
candidate_dir <- file.path(out_dir, "candidate_gene_lists")
plot_dir <- file.path(out_dir, "plots")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

four_group_file <- file.path(root, "standardized_effect_size_4groups", "LMM_all_genes_with_standardized_effect_size_4groups.csv")
frontal_subregion_file <- file.path(root, "frontal_subregion_sensitivity", "frontal_subregion_LMM_all_genes_with_standardized_effect_size.csv")

for (path in c(four_group_file, frontal_subregion_file)) {
  if (!file.exists(path)) {
    stop("Required input file was not found: ", path)
  }
}

message("Reading four-group standardized-effect-size table:")
message(four_group_file)
four_group <- read.csv(four_group_file, check.names = FALSE)

message("Reading frontal-subregion standardized-effect-size table:")
message(frontal_subregion_file)
frontal_subregion <- read.csv(frontal_subregion_file, check.names = FALSE)

as_true_flag <- function(x) {
  if (is.logical(x)) {
    return(!is.na(x) & x)
  }
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

short_scope <- function(x) {
  switch(
    x,
    per_region_all_genes = "perRegion",
    all_four_grouped_regions_together = "all4",
    effect_size_gated_correction = "esGate",
    direction_specific_correction = "dir",
    effect_size_gated_plus_direction_specific_correction = "esGateDir",
    safe_filename(x)
  )
}

short_method <- function(x) {
  switch(
    x,
    Bonferroni = "bonf",
    Holm = "holm",
    Hochberg = "hoch",
    BH_FDR = "BH",
    BY_FDR = "BY",
    Storey_qvalue = "Storey",
    IHW = "IHW",
    Independent_filter_mean_expression_20pct_BH = "IF20BH",
    safe_filename(x)
  )
}

standard_direction <- function(x) {
  ifelse(is.na(x), NA_character_, ifelse(x > 0, "up_in_MDD", ifelse(x < 0, "down_in_MDD", "no_change")))
}

prepare_input <- function(df, analysis_level) {
  if (analysis_level == "four_group") {
    region_col <- if ("grouped_region" %in% colnames(df)) "grouped_region" else "brain_region"
  } else {
    region_col <- "original_frontal_subregion"
  }

  needed <- c(
    region_col, "gene_id", "gene_symbol", "beta_diagnosis",
    "standardized_effect_size_diagnosis",
    "abs_standardized_effect_size_diagnosis", "p_value_diagnosis"
  )
  missing_cols <- setdiff(needed, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Input table is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  out <- data.frame(
    analysis_level = analysis_level,
    region_or_subregion = df[[region_col]],
    gene_id = as.character(df$gene_id),
    gene_symbol = df$gene_symbol,
    beta_diagnosis = suppressWarnings(as.numeric(df$beta_diagnosis)),
    standardized_effect_size_diagnosis = suppressWarnings(as.numeric(df$standardized_effect_size_diagnosis)),
    abs_standardized_effect_size_diagnosis = suppressWarnings(as.numeric(df$abs_standardized_effect_size_diagnosis)),
    p_value_diagnosis = suppressWarnings(as.numeric(df$p_value_diagnosis)),
    direction = if ("direction" %in% colnames(df)) df$direction else standard_direction(df$standardized_effect_size_diagnosis),
    mean_expression_gene = if ("mean_expression_gene" %in% colnames(df)) suppressWarnings(as.numeric(df$mean_expression_gene)) else NA_real_,
    stringsAsFactors = FALSE
  )

  if (analysis_level == "four_group" && "nominal_p001_absStdES08" %in% colnames(df)) {
    out$current_Tier2_flag <- as_true_flag(df$nominal_p001_absStdES08)
  } else if (analysis_level == "frontal_subregion" && "Tier2_nominal_p001_absStdES08" %in% colnames(df)) {
    out$current_Tier2_flag <- as_true_flag(df$Tier2_nominal_p001_absStdES08)
  } else {
  out$current_Tier2_flag <- out$p_value_diagnosis < 0.01 & out$abs_standardized_effect_size_diagnosis > 0.8
  }

  out$direction <- standard_direction(out$standardized_effect_size_diagnosis)
  out$absStdES08_flag <- out$abs_standardized_effect_size_diagnosis > 0.8
  out$current_Tier2_flag[is.na(out$current_Tier2_flag)] <- FALSE
  out$absStdES08_flag[is.na(out$absStdES08_flag)] <- FALSE
  out
}

four_group_prepared <- prepare_input(four_group, "four_group")
frontal_prepared <- prepare_input(frontal_subregion, "frontal_subregion")

adjust_p_values <- function(p_values, method, covariate = NULL) {
  adjusted <- rep(NA_real_, length(p_values))
  valid <- is.finite(p_values) & !is.na(p_values) & p_values >= 0 & p_values <= 1

  if (!any(valid)) {
    return(adjusted)
  }

  if (method %in% c("Bonferroni", "Holm", "Hochberg", "BH_FDR", "BY_FDR")) {
    p_method <- switch(
      method,
      Bonferroni = "bonferroni",
      Holm = "holm",
      Hochberg = "hochberg",
      BH_FDR = "BH",
      BY_FDR = "BY"
    )
    adjusted[valid] <- p.adjust(p_values[valid], method = p_method)
    return(adjusted)
  }

  if (method == "Storey_qvalue") {
    if (!has_qvalue || sum(valid) < 2) {
      return(adjusted)
    }
    adjusted[valid] <- tryCatch(
      qvalue::qvalue(p_values[valid])$qvalues,
      error = function(e) rep(NA_real_, sum(valid))
    )
    return(adjusted)
  }

  if (method == "IHW") {
    cov_valid <- valid & is.finite(covariate) & !is.na(covariate)
    if (!has_ihw || sum(cov_valid) < 50 || length(unique(covariate[cov_valid])) < 2) {
      return(adjusted)
    }
    adjusted[cov_valid] <- tryCatch(
      IHW::adj_pvalues(IHW::ihw(p_values[cov_valid] ~ covariate[cov_valid], alpha = 0.1)),
      error = function(e) rep(NA_real_, sum(cov_valid))
    )
    return(adjusted)
  }

  if (method == "Independent_filter_mean_expression_20pct_BH") {
    cov_valid <- valid & is.finite(covariate) & !is.na(covariate)
    if (sum(cov_valid) < 2 || length(unique(covariate[cov_valid])) < 2) {
      return(adjusted)
    }
    cutoff <- stats::quantile(covariate[cov_valid], probs = 0.20, na.rm = TRUE, names = FALSE)
    filtered <- cov_valid & covariate >= cutoff
    if (sum(filtered) < 2) {
      return(adjusted)
    }
    adjusted[filtered] <- p.adjust(p_values[filtered], method = "BH")
    return(adjusted)
  }

  stop("Unknown correction method: ", method)
}

correction_methods <- c(
  "Bonferroni",
  "Holm",
  "Hochberg",
  "BH_FDR",
  "BY_FDR",
  if (has_qvalue) "Storey_qvalue" else character(),
  if (has_ihw) "IHW" else character(),
  "Independent_filter_mean_expression_20pct_BH"
)

make_scope_groups <- function(df, scope_name, include_all_regions_scope = FALSE) {
  groups <- list()
  regions <- sort(unique(df$region_or_subregion))

  if (scope_name == "per_region_all_genes") {
    for (region in regions) {
      idx <- which(df$region_or_subregion == region)
      groups[[length(groups) + 1]] <- list(region = region, idx = idx)
    }
  }

  if (scope_name == "all_four_grouped_regions_together") {
    if (include_all_regions_scope) {
      for (region in regions) {
        groups[[length(groups) + 1]] <- list(region = region, idx = which(seq_len(nrow(df)) > 0))
      }
    }
  }

  if (scope_name == "effect_size_gated_correction") {
    for (region in regions) {
      idx <- which(df$region_or_subregion == region & df$absStdES08_flag)
      groups[[length(groups) + 1]] <- list(region = region, idx = idx)
    }
  }

  if (scope_name == "direction_specific_correction") {
    for (region in regions) {
      for (direction in c("up_in_MDD", "down_in_MDD")) {
        idx <- which(df$region_or_subregion == region & df$direction == direction)
        groups[[length(groups) + 1]] <- list(region = region, idx = idx)
      }
    }
  }

  if (scope_name == "effect_size_gated_plus_direction_specific_correction") {
    for (region in regions) {
      for (direction in c("up_in_MDD", "down_in_MDD")) {
        idx <- which(df$region_or_subregion == region & df$direction == direction & df$absStdES08_flag)
        groups[[length(groups) + 1]] <- list(region = region, idx = idx)
      }
    }
  }

  groups
}

run_correction_sensitivity <- function(df, analysis_level) {
  message("Running correction sensitivity for ", analysis_level)

  scopes <- c(
    "per_region_all_genes",
    if (analysis_level == "four_group") "all_four_grouped_regions_together" else character(),
    "effect_size_gated_correction",
    "direction_specific_correction",
    "effect_size_gated_plus_direction_specific_correction"
  )

  corrected_rows <- list()
  summary_rows <- list()
  overlap_rows <- list()
  gene_list_counter <- 0L

  region_stats <- lapply(split(df, df$region_or_subregion), function(x) {
    data.frame(
      region_or_subregion = x$region_or_subregion[1],
      n_absStdES08 = sum(x$absStdES08_flag, na.rm = TRUE),
      n_current_Tier2 = sum(x$current_Tier2_flag, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  region_stats <- do.call(rbind, region_stats)

  for (scope in scopes) {
    include_all_regions_scope <- scope == "all_four_grouped_regions_together"
    scope_groups <- make_scope_groups(df, scope, include_all_regions_scope = include_all_regions_scope)

    for (method in correction_methods) {
      adjusted <- rep(NA_real_, nrow(df))
      genes_corrected_over <- setNames(rep(0L, length(unique(df$region_or_subregion))), unique(df$region_or_subregion))

      for (group in scope_groups) {
        idx <- group$idx
        if (length(idx) == 0) {
          next
        }
        group_adj <- adjust_p_values(df$p_value_diagnosis[idx], method, covariate = df$mean_expression_gene[idx])
        adjusted[idx] <- group_adj

        corrected_idx <- idx[is.finite(group_adj) & !is.na(group_adj)]
        corrected_regions <- unique(df$region_or_subregion[corrected_idx])
        for (region in corrected_regions) {
          genes_corrected_over[region] <- genes_corrected_over[region] + sum(df$region_or_subregion[corrected_idx] == region)
        }
      }

      in_table <- is.finite(adjusted) & !is.na(adjusted)
      if (any(in_table)) {
        corrected_part <- data.frame(
          analysis_level = df$analysis_level[in_table],
          region_or_subregion = df$region_or_subregion[in_table],
          gene_id = df$gene_id[in_table],
          gene_symbol = df$gene_symbol[in_table],
          beta_diagnosis = df$beta_diagnosis[in_table],
          standardized_effect_size_diagnosis = df$standardized_effect_size_diagnosis[in_table],
          abs_standardized_effect_size_diagnosis = df$abs_standardized_effect_size_diagnosis[in_table],
          p_value_diagnosis = df$p_value_diagnosis[in_table],
          current_Tier2_flag = df$current_Tier2_flag[in_table],
          correction_scope = scope,
          correction_method = method,
          adjusted_p_value = adjusted[in_table],
          direction = df$direction[in_table],
          stringsAsFactors = FALSE
        )
        corrected_rows[[length(corrected_rows) + 1]] <- corrected_part
      }

      for (region in sort(unique(df$region_or_subregion))) {
        region_idx <- df$region_or_subregion == region
        adjusted_set <- region_idx & is.finite(adjusted) & !is.na(adjusted) &
          adjusted < 0.10 & df$absStdES08_flag
        current_set <- region_idx & df$current_Tier2_flag

        n_union <- length(union(df$gene_id[current_set], df$gene_id[adjusted_set]))
        n_overlap <- length(intersect(df$gene_id[current_set], df$gene_id[adjusted_set]))
        overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
          analysis_level = analysis_level,
          region_or_subregion = region,
          correction_scope = scope,
          correction_method = method,
          adjusted_threshold = "adjusted_p_lt_0.10_absStdES08",
          current_Tier2_genes = sum(current_set, na.rm = TRUE),
          adjusted_set_genes = sum(adjusted_set, na.rm = TRUE),
          overlap_with_current_Tier2 = n_overlap,
          genes_unique_to_current_Tier2 = length(setdiff(df$gene_id[current_set], df$gene_id[adjusted_set])),
          genes_unique_to_adjusted_set = length(setdiff(df$gene_id[adjusted_set], df$gene_id[current_set])),
          jaccard_index = ifelse(n_union > 0, n_overlap / n_union, NA_real_),
          stringsAsFactors = FALSE
        )

        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          analysis_level = analysis_level,
          region_or_subregion = region,
          correction_scope = scope,
          correction_method = method,
          n_genes_corrected_over = unname(genes_corrected_over[region]),
          n_absStdES08 = region_stats$n_absStdES08[match(region, region_stats$region_or_subregion)],
          n_current_Tier2 = region_stats$n_current_Tier2[match(region, region_stats$region_or_subregion)],
          n_adj005_absStdES08 = sum(region_idx & is.finite(adjusted) & adjusted < 0.05 & df$absStdES08_flag, na.rm = TRUE),
          n_adj010_absStdES08 = sum(region_idx & is.finite(adjusted) & adjusted < 0.10 & df$absStdES08_flag, na.rm = TRUE),
          n_adj020_absStdES08 = sum(region_idx & is.finite(adjusted) & adjusted < 0.20 & df$absStdES08_flag, na.rm = TRUE),
          n_up_adj010_absStdES08 = sum(region_idx & is.finite(adjusted) & adjusted < 0.10 & df$absStdES08_flag & df$direction == "up_in_MDD", na.rm = TRUE),
          n_down_adj010_absStdES08 = sum(region_idx & is.finite(adjusted) & adjusted < 0.10 & df$absStdES08_flag & df$direction == "down_in_MDD", na.rm = TRUE),
          stringsAsFactors = FALSE
        )

        for (threshold in c(0.05, 0.10, 0.20)) {
        gene_mask <- region_idx & is.finite(adjusted) & !is.na(adjusted) &
            adjusted < threshold & df$absStdES08_flag
          gene_mask[is.na(gene_mask)] <- FALSE
          n_gene_list <- sum(gene_mask)
          gene_list <- data.frame(
            gene_id = df$gene_id[gene_mask],
            gene_symbol = df$gene_symbol[gene_mask],
            region_or_subregion = df$region_or_subregion[gene_mask],
            correction_scope = rep(scope, n_gene_list),
            correction_method = rep(method, n_gene_list),
            beta_diagnosis = df$beta_diagnosis[gene_mask],
            standardized_effect_size_diagnosis = df$standardized_effect_size_diagnosis[gene_mask],
            abs_standardized_effect_size_diagnosis = df$abs_standardized_effect_size_diagnosis[gene_mask],
            p_value_diagnosis = df$p_value_diagnosis[gene_mask],
            adjusted_p_value = adjusted[gene_mask],
            direction = df$direction[gene_mask],
            current_Tier2_flag = df$current_Tier2_flag[gene_mask],
            stringsAsFactors = FALSE
          )
          if (nrow(gene_list) > 0) {
            gene_list <- gene_list[order(gene_list$adjusted_p_value, gene_list$p_value_diagnosis, na.last = TRUE), , drop = FALSE]
          }
          file_name <- paste(
            safe_filename(ifelse(analysis_level == "frontal_subregion", "frontSub", "fourGrp")),
            substr(safe_filename(region), 1, 28),
            short_scope(scope),
            short_method(method),
            paste0("adjp_lt_", gsub("\\.", "", sprintf("%.2f", threshold)), "_absStdES08.csv"),
            sep = "__"
          )
          write.csv(gene_list, file.path(candidate_dir, file_name), row.names = FALSE)
          gene_list_counter <- gene_list_counter + 1L
        }
      }
    }
  }

  list(
    corrected = do.call(rbind, corrected_rows),
    summary = do.call(rbind, summary_rows),
    overlap = do.call(rbind, overlap_rows),
    n_gene_lists = gene_list_counter
  )
}

four_results <- run_correction_sensitivity(four_group_prepared, "four_group")
frontal_results <- run_correction_sensitivity(frontal_prepared, "frontal_subregion")

write.csv(
  four_results$corrected,
  file.path(out_dir, "four_group_pvalue_correction_sensitivity_full.csv"),
  row.names = FALSE
)
write.csv(
  frontal_results$corrected,
  file.path(out_dir, "frontal_subregion_pvalue_correction_sensitivity_full.csv"),
  row.names = FALSE
)
write.csv(
  four_results$summary,
  file.path(out_dir, "four_group_correction_sensitivity_summary.csv"),
  row.names = FALSE
)
write.csv(
  frontal_results$summary,
  file.path(out_dir, "frontal_subregion_correction_sensitivity_summary.csv"),
  row.names = FALSE
)

combined_overlap <- rbind(four_results$overlap, frontal_results$overlap)
write.csv(
  combined_overlap,
  file.path(out_dir, "correction_vs_current_Tier2_overlap.csv"),
  row.names = FALSE
)

combined_summary <- rbind(four_results$summary, frontal_results$summary)

plot_png <- function(path, width = 1600, height = 1000, code) {
  png(path, width = width, height = height)
  tryCatch(code(), finally = dev.off())
}

plot_png(file.path(plot_dir, "barplot_adj010_counts_by_method_scope.png"), height = 1200, code = function() {
  y <- combined_summary$n_adj010_absStdES08
  labels <- paste(combined_summary$region_or_subregion, combined_summary$correction_method, combined_summary$correction_scope, sep = " | ")
  top <- order(y, decreasing = TRUE)[seq_len(min(60, length(y)))]
  par(mar = c(14, 5, 4, 1))
  barplot(
    y[top],
    names.arg = labels[top],
    las = 2,
    cex.names = 0.55,
    col = "#4C78A8",
    ylab = "Adjusted p < 0.10 and absStdES > 0.8 gene count",
    main = "Top correction sensitivity counts"
  )
})

plot_png(file.path(plot_dir, "heatmap_table_adj010_counts_by_method_scope.png"), width = 1800, height = 1200, code = function() {
  keys <- paste(combined_summary$correction_method, combined_summary$correction_scope, sep = "\n")
  regions <- unique(combined_summary$region_or_subregion)
  key_levels <- unique(keys)
  mat <- matrix(0, nrow = length(regions), ncol = length(key_levels), dimnames = list(regions, key_levels))
  for (i in seq_len(nrow(combined_summary))) {
    mat[combined_summary$region_or_subregion[i], keys[i]] <- combined_summary$n_adj010_absStdES08[i]
  }
  par(mar = c(12, 12, 4, 2))
  image(
    x = seq_len(ncol(mat)),
    y = seq_len(nrow(mat)),
    z = t(mat[nrow(mat):1, , drop = FALSE]),
    axes = FALSE,
    col = colorRampPalette(c("white", "#F2C14E", "#D95F02"))(100),
    xlab = "",
    ylab = "",
    main = "n_adj010_absStdES08 by correction method and scope"
  )
  axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.55)
  axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.75)
  for (x in seq_len(ncol(mat))) {
    for (y_pos in seq_len(nrow(mat))) {
      text(x, y_pos, labels = mat[nrow(mat) - y_pos + 1, x], cex = 0.65)
    }
  }
})

plot_sample_scatter <- function(corrected_df, path, title_text) {
  sample_df <- corrected_df[corrected_df$correction_scope %in% c("per_region_all_genes", "effect_size_gated_correction") &
                              corrected_df$correction_method %in% c("Bonferroni", "BH_FDR", "Storey_qvalue", "IHW"), , drop = FALSE]
  if (nrow(sample_df) > 60000) {
    set.seed(100)
    sample_df <- sample_df[sample(seq_len(nrow(sample_df)), 60000), , drop = FALSE]
  }
  plot_png(path, width = 1300, height = 1000, code = function() {
    cols <- as.integer(factor(paste(sample_df$correction_method, sample_df$correction_scope)))
    plot(
      sample_df$p_value_diagnosis,
      sample_df$adjusted_p_value,
      pch = 16,
      cex = 0.35,
      col = adjustcolor(cols, alpha.f = 0.35),
      xlab = "p_value_diagnosis",
      ylab = "adjusted_p_value",
      main = title_text
    )
    abline(0, 1, col = "gray40", lty = 2)
    legend("bottomright", legend = levels(factor(paste(sample_df$correction_method, sample_df$correction_scope))), pch = 16, col = seq_along(levels(factor(paste(sample_df$correction_method, sample_df$correction_scope)))), cex = 0.65)
  })
}
plot_sample_scatter(four_results$corrected, file.path(plot_dir, "four_group_pvalue_vs_adjusted_pvalue_scatter.png"), "Four-group p-value vs adjusted p-value")
plot_sample_scatter(frontal_results$corrected, file.path(plot_dir, "frontal_subregion_pvalue_vs_adjusted_pvalue_scatter.png"), "Frontal-subregion p-value vs adjusted p-value")

plot_png(file.path(plot_dir, "current_Tier2_vs_adjusted010_absStdES08_counts.png"), width = 1600, height = 1000, code = function() {
  best_per_region <- aggregate(
    n_adj010_absStdES08 ~ analysis_level + region_or_subregion,
    data = combined_summary,
    FUN = max
  )
  tier2 <- aggregate(
    n_current_Tier2 ~ analysis_level + region_or_subregion,
    data = combined_summary,
    FUN = max
  )
  merged <- merge(best_per_region, tier2, by = c("analysis_level", "region_or_subregion"))
  par(mar = c(10, 5, 4, 1))
  mat <- rbind(merged$n_current_Tier2, merged$n_adj010_absStdES08)
  colnames(mat) <- merged$region_or_subregion
  barplot(
    mat,
    beside = TRUE,
    las = 2,
    col = c("#4C78A8", "#F58518"),
    ylab = "Gene count",
    main = "Current Tier 2 count vs maximum adjusted p < 0.10 sensitivity count"
  )
  legend("topright", legend = c("Current Tier 2", "Max adjusted p < 0.10"), fill = c("#4C78A8", "#F58518"))
})

summarize_extremes <- function(summary_df) {
  positive <- summary_df[summary_df$n_adj010_absStdES08 > 0, , drop = FALSE]
  largest <- summary_df[order(-summary_df$n_adj010_absStdES08), , drop = FALSE][1:min(10, nrow(summary_df)), ]
  smallest <- if (nrow(positive) > 0) positive[order(positive$n_adj010_absStdES08), , drop = FALSE][1:min(10, nrow(positive)), ] else positive
  list(largest = largest, smallest_positive = smallest)
}

four_extremes <- summarize_extremes(four_results$summary)
frontal_extremes <- summarize_extremes(frontal_results$summary)

top_overlap <- combined_overlap[order(-combined_overlap$jaccard_index, -combined_overlap$overlap_with_current_Tier2, na.last = TRUE), , drop = FALSE]
top_overlap <- head(top_overlap, 20)

report_path <- file.path(out_dir, "Pvalue_Correction_Sensitivity_Report.md")
report <- c(
  "# P-Value Correction Sensitivity Audit",
  "",
  "## Purpose",
  "",
  "This audit tests whether the number of MDD candidate genes is sensitive to the p-value correction method and correction scope. It uses existing LMM p-values and standardized effect sizes only; no LMMs were rerun.",
  "",
  "## Important Warning",
  "",
  "The current Tier 2 definition is preserved and remains the working exploratory DE candidate definition unless we explicitly decide otherwise.",
  "",
  "Current Tier 2: `p_value_diagnosis < 0.01` and `abs_standardized_effect_size_diagnosis > 0.8`.",
  "",
  "## Methods Tested",
  "",
  paste0("Correction methods included: ", paste(correction_methods, collapse = ", "), "."),
  "",
  "Correction scopes included:",
  "",
  "- Per-region/per-subregion correction across all genes.",
  "- Four grouped regions corrected together, for the four-group analysis only.",
  "- Effect-size-gated correction: first filter to `abs_standardized_effect_size_diagnosis > 0.8`, then correct p-values within that subset. This is not globally controlled FDR across all tested genes.",
  "- Direction-specific correction: split up-in-MDD and down-in-MDD genes, then correct separately. This is a sensitivity result only.",
  "- Effect-size-gated plus direction-specific correction: first filter by effect size, then split by direction and correct separately. This is a sensitivity result only.",
  "",
  "Independent filtering used `mean_expression_gene` where available and applied BH after retaining genes at or above the 20th percentile of mean expression within the tested correction set.",
  "",
  "## Main Four-Group Results",
  "",
  "Largest four-group adjusted p < 0.10 plus absStdES > 0.8 counts:",
  "",
  paste(capture.output(print(four_extremes$largest[, c("region_or_subregion", "correction_scope", "correction_method", "n_adj010_absStdES08", "n_current_Tier2")], row.names = FALSE)), collapse = "\n"),
  "",
  "Smallest positive four-group adjusted p < 0.10 plus absStdES > 0.8 counts:",
  "",
  paste(capture.output(print(four_extremes$smallest_positive[, c("region_or_subregion", "correction_scope", "correction_method", "n_adj010_absStdES08", "n_current_Tier2")], row.names = FALSE)), collapse = "\n"),
  "",
  "## Frontal Subregion Results",
  "",
  "Largest frontal-subregion adjusted p < 0.10 plus absStdES > 0.8 counts:",
  "",
  paste(capture.output(print(frontal_extremes$largest[, c("region_or_subregion", "correction_scope", "correction_method", "n_adj010_absStdES08", "n_current_Tier2")], row.names = FALSE)), collapse = "\n"),
  "",
  "These results show whether Cingulate gyrus 25, dlPFC BA8/9, or OFC BA11 gains corrected genes under any alternative correction. Any gains under effect-size-gated or direction-specific scopes should be interpreted as sensitivity findings, not replacement primary results.",
  "",
  "## Comparison With Current Tier 2",
  "",
  "Highest overlaps with current Tier 2 at adjusted p < 0.10 plus absStdES > 0.8:",
  "",
  paste(capture.output(print(top_overlap[, c("analysis_level", "region_or_subregion", "correction_scope", "correction_method", "current_Tier2_genes", "adjusted_set_genes", "overlap_with_current_Tier2", "genes_unique_to_current_Tier2", "genes_unique_to_adjusted_set", "jaccard_index")], row.names = FALSE)), collapse = "\n"),
  "",
  "## Recommendation",
  "",
  "Do not choose a correction only because it gives more genes. Keep Tier 2 as the exploratory working DE candidate set unless a different correction is methodologically justified and documented. Alternative corrected sets should be used as sensitivity results, while Tier 2 remains the practical exploratory set for enrichment, ML preparation, and AHBA mapping.",
  "",
  "## Suggested Wording",
  "",
  "\"Because the exact correction implementation in the reference workflow may differ from the written description, we performed a p-value correction sensitivity audit. The current Tier 2 exploratory DE candidate definition, p < 0.01 and |standardized effect size| > 0.8, was preserved. Alternative correction scopes and methods were evaluated only as sensitivity analyses and were not used to overwrite the main working gene sets.\""
)
writeLines(report, report_path)

message("")
message("Saved p-value correction sensitivity audit to: ", out_dir)
message("Four-group summary rows: ", nrow(four_results$summary))
message("Frontal-subregion summary rows: ", nrow(frontal_results$summary))
message("Candidate gene-list files created: ", four_results$n_gene_lists + frontal_results$n_gene_lists)
message("qvalue available: ", has_qvalue)
message("IHW available: ", has_ihw)
