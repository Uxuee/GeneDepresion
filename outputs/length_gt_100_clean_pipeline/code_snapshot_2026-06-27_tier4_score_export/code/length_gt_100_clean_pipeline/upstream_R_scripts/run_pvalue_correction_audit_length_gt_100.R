## Fixed p-value correction audit for Length > 100 Elisa-style LMM results.
##
## This script does not rerun normalization and does not rerun the LMM.
## It reads the existing LMM table plus expression covariates computed from the
## normalized log2(GeTMM + 1) matrices, then writes FIXED audit outputs.

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)
} else {
  script_path <- normalizePath("code/run_pvalue_correction_audit_length_gt_100.R", mustWork = FALSE)
}
repo_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
base_dir <- file.path(repo_dir, "outputs", "length_gt_100_elisa_style")
audit_dir <- file.path(base_dir, "pvalue_correction_audit")

lmm_file <- file.path(base_dir, "LMM_all_genes_length_gt_100.csv")
covariate_file <- file.path(audit_dir, "expression_filter_covariates_length_gt_100.csv")

corrected_file <- file.path(audit_dir, "LMM_all_genes_with_pvalue_corrections_FIXED.csv")
merged_file <- file.path(audit_dir, "LMM_all_genes_with_expression_covariates.csv")
summary_file <- file.path(audit_dir, "pvalue_correction_summary_by_region_FIXED.csv")
report_file <- file.path(audit_dir, "Pvalue_Correction_Audit_Report_FIXED.md")
plots_dir <- file.path(audit_dir, "plots_FIXED")
candidate_dir <- file.path(audit_dir, "candidate_gene_lists_FIXED")

dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)

baseline_col <- "exploratory_p001_absbeta02"

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]+", "_", as.character(x))
  gsub("^_|_$", "", x)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

bh_adjust <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- !is.na(p)
  out[ok] <- p.adjust(p[ok], method = "BH")
  out
}

by_adjust <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- !is.na(p)
  out[ok] <- p.adjust(p[ok], method = "BY")
  out
}

storey_adjust <- function(p) {
  out <- list(qvalues = rep(NA_real_, length(p)), pi0 = NA_real_, error = NA_character_)
  ok <- !is.na(p) & p >= 0 & p <= 1
  if (sum(ok) < 2 || !requireNamespace("qvalue", quietly = TRUE)) {
    return(out)
  }
  result <- tryCatch(
    qvalue::qvalue(p[ok]),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    out$error <- conditionMessage(result)
    return(out)
  }
  out$qvalues[ok] <- result$qvalues
  out$pi0 <- result$pi0
  out
}

run_ihw <- function(region_data, covariate_column) {
  out <- rep(NA_real_, nrow(region_data))
  if (!requireNamespace("IHW", quietly = TRUE) || is.na(covariate_column)) {
    return(out)
  }

  covariate <- safe_numeric(region_data[[covariate_column]])
  p_values <- safe_numeric(region_data$p_value_diagnosis)
  ok <- !is.na(p_values) & p_values >= 0 & p_values <= 1 & !is.na(covariate)
  if (sum(ok) < 50 || length(unique(covariate[ok])) < 2) {
    return(out)
  }

  ihw_data <- data.frame(
    pvalue = p_values[ok],
    covariate = covariate[ok],
    stringsAsFactors = FALSE
  )

  fit <- tryCatch(
    IHW::ihw(pvalue ~ covariate, data = ihw_data, alpha = 0.10),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(out)
  }
  out[ok] <- IHW::adj_pvalues(fit)
  out
}

build_filter_specs <- function(region_data) {
  specs <- list()

  add_quantile_filters <- function(column, label) {
    values <- safe_numeric(region_data[[column]])
    if (all(is.na(values))) return(NULL)
    for (q in c(0.10, 0.20, 0.30)) {
      threshold <- as.numeric(stats::quantile(values, probs = q, na.rm = TRUE))
      specs[[length(specs) + 1]] <<- list(
        name = paste0(label, "_above_", q * 100, "th_percentile"),
        column = column,
        threshold = threshold,
        keep = values > threshold,
        aggressiveness = q
      )
    }
  }

  add_quantile_filters("mean_expression_overall", "mean_expression")
  add_quantile_filters("expression_variance_overall", "expression_variance")

  detection <- safe_numeric(region_data$detection_rate_overall)
  if (!all(is.na(detection))) {
    for (threshold in c(0.20, 0.50, 0.80)) {
      specs[[length(specs) + 1]] <- list(
        name = paste0("detection_rate_ge_", threshold),
        column = "detection_rate_overall",
        threshold = threshold,
        keep = detection >= threshold,
        aggressiveness = threshold
      )
    }
  }

  specs
}

write_candidate_file <- function(data, region, label, mask, sort_column = "p_value_diagnosis") {
  mask <- !is.na(mask) & mask
  out <- data[mask, , drop = FALSE]
  if (nrow(out) > 0 && sort_column %in% colnames(out)) {
    out <- out[order(out[[sort_column]], na.last = TRUE), , drop = FALSE]
  }
  write.csv(
    out,
    file.path(candidate_dir, paste0(safe_name(region), "_", label, ".csv")),
    row.names = FALSE
  )
}

message("Reading LMM results: ", lmm_file)
lmm <- read.csv(lmm_file, check.names = FALSE)
message("Reading expression covariates: ", covariate_file)
covariates <- read.csv(covariate_file, check.names = FALSE)

required_lmm <- c("brain_region", "gene_id", "beta_diagnosis", "p_value_diagnosis")
missing_lmm <- setdiff(required_lmm, colnames(lmm))
if (length(missing_lmm) > 0) {
  stop("LMM file is missing required columns: ", paste(missing_lmm, collapse = ", "))
}

required_covariates <- c(
  "brain_region",
  "gene_id",
  "mean_expression_overall",
  "expression_variance_overall",
  "detection_rate_overall"
)
missing_covariates <- setdiff(required_covariates, colnames(covariates))
if (length(missing_covariates) > 0) {
  stop("Covariate file is missing required columns: ", paste(missing_covariates, collapse = ", "))
}

lmm_rows_before_merge <- nrow(lmm)
lmm$gene_id <- as.character(lmm$gene_id)
covariates$gene_id <- as.character(covariates$gene_id)

covariate_cols <- c(
  "brain_region",
  "gene_id",
  "mean_expression_overall",
  "expression_variance_overall",
  "detection_rate_overall",
  "n_samples_for_covariate"
)
if ("gene_symbol" %in% colnames(covariates) && !("gene_symbol" %in% colnames(lmm))) {
  covariate_cols <- c(covariate_cols, "gene_symbol")
}

merged <- merge(
  lmm,
  covariates[, covariate_cols, drop = FALSE],
  by = c("brain_region", "gene_id"),
  all.x = TRUE,
  sort = FALSE
)

merge_summary <- data.frame(
  lmm_rows_before_merge = lmm_rows_before_merge,
  lmm_rows_after_merge = nrow(merged),
  nonmissing_mean_expression_overall = sum(!is.na(merged$mean_expression_overall)),
  nonmissing_expression_variance_overall = sum(!is.na(merged$expression_variance_overall)),
  nonmissing_detection_rate_overall = sum(!is.na(merged$detection_rate_overall)),
  stringsAsFactors = FALSE
)
write.csv(merge_summary, file.path(audit_dir, "merge_expression_covariates_summary.csv"), row.names = FALSE)
write.csv(merged, merged_file, row.names = FALSE)

message("Rows before merge: ", lmm_rows_before_merge)
message("Rows after merge:  ", nrow(merged))
message("Non-missing mean expression: ", merge_summary$nonmissing_mean_expression_overall)
message("Non-missing variance:        ", merge_summary$nonmissing_expression_variance_overall)
message("Non-missing detection rate:  ", merge_summary$nonmissing_detection_rate_overall)

merged$beta_diagnosis <- safe_numeric(merged$beta_diagnosis)
merged$p_value_diagnosis <- safe_numeric(merged$p_value_diagnosis)
merged$abs_beta_diagnosis <- abs(merged$beta_diagnosis)
merged$signed_logp_score <- sign(merged$beta_diagnosis) * -log10(pmax(merged$p_value_diagnosis, .Machine$double.xmin))
merged[[baseline_col]] <- merged$p_value_diagnosis < 0.01 & merged$abs_beta_diagnosis > 0.2

merged$fdr_BH <- NA_real_
merged$fdr_BY <- NA_real_
merged$qvalue_Storey <- NA_real_
merged$pi0_Storey <- NA_real_
merged$ihw_adj_pvalue <- NA_real_

qvalue_available <- requireNamespace("qvalue", quietly = TRUE)
ihw_available <- requireNamespace("IHW", quietly = TRUE)

filter_summary_rows <- list()
summary_rows <- list()
storey_errors <- list()

for (region in unique(merged$brain_region)) {
  message("Processing region: ", region)
  region_index <- which(merged$brain_region == region)
  region_data <- merged[region_index, , drop = FALSE]

  merged$fdr_BH[region_index] <- bh_adjust(region_data$p_value_diagnosis)
  merged$fdr_BY[region_index] <- by_adjust(region_data$p_value_diagnosis)

  storey_result <- storey_adjust(region_data$p_value_diagnosis)
  merged$qvalue_Storey[region_index] <- storey_result$qvalues
  merged$pi0_Storey[region_index] <- storey_result$pi0
  if (!is.na(storey_result$error)) {
    storey_errors[[length(storey_errors) + 1]] <- data.frame(
      brain_region = region,
      storey_error = storey_result$error,
      stringsAsFactors = FALSE
    )
  }

  ihw_covariate <- NA_character_
  if (sum(!is.na(region_data$mean_expression_overall)) > 0) {
    ihw_covariate <- "mean_expression_overall"
  } else if (sum(!is.na(region_data$detection_rate_overall)) > 0) {
    ihw_covariate <- "detection_rate_overall"
  }
  merged$ihw_adj_pvalue[region_index] <- run_ihw(region_data, ihw_covariate)

  filter_specs <- build_filter_specs(region_data)
  region_filter_rows <- list()

  for (spec in filter_specs) {
    keep <- !is.na(spec$keep) & spec$keep
    kept_p <- region_data$p_value_diagnosis[keep]

    bh <- rep(NA_real_, nrow(region_data))
    storey_q <- rep(NA_real_, nrow(region_data))
    if (sum(keep) > 0) {
      bh[keep] <- bh_adjust(kept_p)
      if (qvalue_available) {
        storey_filtered <- storey_adjust(kept_p)
        storey_q[keep] <- storey_filtered$qvalues
      }
    }

    bh_col <- paste0("BH_after_filter_", safe_name(spec$name))
    storey_col <- paste0("Storey_after_filter_", safe_name(spec$name))
    merged[[bh_col]] <- NA_real_
    merged[[storey_col]] <- NA_real_
    merged[[bh_col]][region_index] <- bh
    merged[[storey_col]][region_index] <- storey_q

    row <- data.frame(
      brain_region = region,
      filter_name = spec$name,
      filter_column = spec$column,
      filter_threshold = spec$threshold,
      filter_aggressiveness = spec$aggressiveness,
      genes_tested_after_filter = sum(keep),
      exploratory_p001_absbeta02_after_filter = sum(region_data[[baseline_col]][keep], na.rm = TRUE),
      BH_after_filter_lt_005 = sum(bh < 0.05, na.rm = TRUE),
      BH_after_filter_lt_010 = sum(bh < 0.10, na.rm = TRUE),
      Storey_after_filter_lt_005 = sum(storey_q < 0.05, na.rm = TRUE),
      Storey_after_filter_lt_010 = sum(storey_q < 0.10, na.rm = TRUE),
      bh_adjusted_column = bh_col,
      storey_adjusted_column = storey_col,
      stringsAsFactors = FALSE
    )
    region_filter_rows[[length(region_filter_rows) + 1]] <- row
    filter_summary_rows[[length(filter_summary_rows) + 1]] <- row
  }

  if (length(region_filter_rows) == 0) {
    best_filter <- NULL
  } else {
    region_filter_summary <- do.call(rbind, region_filter_rows)
    region_filter_summary <- region_filter_summary[order(
      -region_filter_summary$BH_after_filter_lt_010,
      -region_filter_summary$BH_after_filter_lt_005,
      region_filter_summary$filter_aggressiveness,
      -region_filter_summary$genes_tested_after_filter
    ), , drop = FALSE]
    best_filter <- region_filter_summary[1, , drop = FALSE]
  }

  group <- merged[region_index, , drop = FALSE]
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    brain_region = region,
    genes_tested = nrow(group),
    raw_p_lt_005 = sum(group$p_value_diagnosis < 0.05, na.rm = TRUE),
    raw_p_lt_001 = sum(group$p_value_diagnosis < 0.01, na.rm = TRUE),
    raw_p_lt_0001 = sum(group$p_value_diagnosis < 0.001, na.rm = TRUE),
    exploratory_p001_absbeta02 = sum(group[[baseline_col]], na.rm = TRUE),
    BH_FDR_lt_005 = sum(group$fdr_BH < 0.05, na.rm = TRUE),
    BH_FDR_lt_010 = sum(group$fdr_BH < 0.10, na.rm = TRUE),
    BY_FDR_lt_005 = sum(group$fdr_BY < 0.05, na.rm = TRUE),
    BY_FDR_lt_010 = sum(group$fdr_BY < 0.10, na.rm = TRUE),
    Storey_q_lt_005 = sum(group$qvalue_Storey < 0.05, na.rm = TRUE),
    Storey_q_lt_010 = sum(group$qvalue_Storey < 0.10, na.rm = TRUE),
    pi0_Storey = unique(group$pi0_Storey[!is.na(group$pi0_Storey)])[1],
    best_independent_filter_name = if (is.null(best_filter)) NA_character_ else best_filter$filter_name,
    genes_tested_after_best_filter = if (is.null(best_filter)) NA_integer_ else best_filter$genes_tested_after_filter,
    BH_after_best_filter_lt_005 = if (is.null(best_filter)) NA_integer_ else best_filter$BH_after_filter_lt_005,
    BH_after_best_filter_lt_010 = if (is.null(best_filter)) NA_integer_ else best_filter$BH_after_filter_lt_010,
    Storey_after_best_filter_lt_005 = if (is.null(best_filter)) NA_integer_ else best_filter$Storey_after_filter_lt_005,
    Storey_after_best_filter_lt_010 = if (is.null(best_filter)) NA_integer_ else best_filter$Storey_after_filter_lt_010,
    IHW_lt_005 = sum(group$ihw_adj_pvalue < 0.05, na.rm = TRUE),
    IHW_lt_010 = sum(group$ihw_adj_pvalue < 0.10, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  write_candidate_file(group, region, "exploratory_p001_absbeta02_FIXED", group[[baseline_col]])
  write_candidate_file(group, region, "BH_FDR_lt_005_FIXED", group$fdr_BH < 0.05, "fdr_BH")
  write_candidate_file(group, region, "BH_FDR_lt_010_FIXED", group$fdr_BH < 0.10, "fdr_BH")
  write_candidate_file(group, region, "BY_FDR_lt_005_FIXED", group$fdr_BY < 0.05, "fdr_BY")
  write_candidate_file(group, region, "BY_FDR_lt_010_FIXED", group$fdr_BY < 0.10, "fdr_BY")
  write_candidate_file(group, region, "Storey_q_lt_005_FIXED", group$qvalue_Storey < 0.05, "qvalue_Storey")
  write_candidate_file(group, region, "Storey_q_lt_010_FIXED", group$qvalue_Storey < 0.10, "qvalue_Storey")
  write_candidate_file(group, region, "IHW_lt_005_FIXED", group$ihw_adj_pvalue < 0.05, "ihw_adj_pvalue")
  write_candidate_file(group, region, "IHW_lt_010_FIXED", group$ihw_adj_pvalue < 0.10, "ihw_adj_pvalue")

  if (!is.null(best_filter)) {
    bh_col <- best_filter$bh_adjusted_column
    storey_col <- best_filter$storey_adjusted_column
    write_candidate_file(group, region, "best_independent_filter_BH_lt_005_FIXED", group[[bh_col]] < 0.05, bh_col)
    write_candidate_file(group, region, "best_independent_filter_BH_lt_010_FIXED", group[[bh_col]] < 0.10, bh_col)
    write_candidate_file(group, region, "best_independent_filter_Storey_lt_005_FIXED", group[[storey_col]] < 0.05, storey_col)
    write_candidate_file(group, region, "best_independent_filter_Storey_lt_010_FIXED", group[[storey_col]] < 0.10, storey_col)
  }
}

filter_summary <- if (length(filter_summary_rows) > 0) do.call(rbind, filter_summary_rows) else data.frame()
summary_table <- do.call(rbind, summary_rows)

write.csv(filter_summary, file.path(audit_dir, "independent_filtering_summary_FIXED.csv"), row.names = FALSE)
if (length(storey_errors) > 0) {
  write.csv(do.call(rbind, storey_errors), file.path(audit_dir, "storey_qvalue_errors_FIXED.csv"), row.names = FALSE)
}

write.csv(merged, corrected_file, row.names = FALSE)
write.csv(summary_table, summary_file, row.names = FALSE)

make_plots <- function(region, group) {
  prefix <- safe_name(region)
  p <- pmax(group$p_value_diagnosis, .Machine$double.xmin)
  minus_logp <- -log10(p)

  png(file.path(plots_dir, paste0(prefix, "_pvalue_histogram_FIXED.png")), width = 900, height = 650)
  hist(group$p_value_diagnosis, breaks = 50, col = "gray75", border = "white", main = paste("P-value histogram:", region), xlab = "p_value_diagnosis")
  dev.off()

  observed <- -log10(sort(p))
  expected <- -log10((seq_along(observed) - 0.5) / length(observed))
  png(file.path(plots_dir, paste0(prefix, "_qq_plot_FIXED.png")), width = 800, height = 800)
  plot(expected, observed, pch = 16, cex = 0.45, col = "#3b6ea8", xlab = "Expected -log10(p)", ylab = "Observed -log10(p)", main = paste("QQ plot:", region))
  abline(0, 1, lty = 2, col = "gray30")
  dev.off()

  png(file.path(plots_dir, paste0(prefix, "_beta_histogram_FIXED.png")), width = 900, height = 650)
  hist(group$beta_diagnosis, breaks = 50, col = "gray70", border = "white", main = paste("Beta histogram:", region), xlab = "beta_diagnosis")
  dev.off()

  png(file.path(plots_dir, paste0(prefix, "_abs_beta_histogram_FIXED.png")), width = 900, height = 650)
  hist(group$abs_beta_diagnosis, breaks = 50, col = "gray70", border = "white", main = paste("Absolute beta histogram:", region), xlab = "abs(beta_diagnosis)")
  dev.off()

  png(file.path(plots_dir, paste0(prefix, "_mean_expression_vs_minus_log10_p_FIXED.png")), width = 900, height = 650)
  plot(group$mean_expression_overall, minus_logp, pch = 16, cex = 0.45, col = "#3b6ea8", xlab = "mean_expression_overall", ylab = "-log10(p_value_diagnosis)", main = paste("Mean expression vs -log10(p):", region))
  dev.off()

  png(file.path(plots_dir, paste0(prefix, "_expression_variance_vs_minus_log10_p_FIXED.png")), width = 900, height = 650)
  plot(group$expression_variance_overall, minus_logp, pch = 16, cex = 0.45, col = "#3b6ea8", xlab = "expression_variance_overall", ylab = "-log10(p_value_diagnosis)", main = paste("Expression variance vs -log10(p):", region))
  dev.off()

  png(file.path(plots_dir, paste0(prefix, "_detection_rate_vs_minus_log10_p_FIXED.png")), width = 900, height = 650)
  plot(group$detection_rate_overall, minus_logp, pch = 16, cex = 0.45, col = "#3b6ea8", xlab = "detection_rate_overall", ylab = "-log10(p_value_diagnosis)", main = paste("Detection rate vs -log10(p):", region))
  dev.off()

  png(file.path(plots_dir, paste0(prefix, "_volcano_exploratory_highlighted_FIXED.png")), width = 900, height = 750)
  plot(group$beta_diagnosis, minus_logp, pch = 16, cex = 0.45, col = ifelse(group[[baseline_col]], "#b42318", "#747b85"), xlab = "beta_diagnosis", ylab = "-log10(p_value_diagnosis)", main = paste("Volcano:", region))
  abline(h = -log10(0.01), v = c(-0.2, 0.2), lty = 2, col = "gray30")
  dev.off()
}

for (region in unique(merged$brain_region)) {
  make_plots(region, merged[merged$brain_region == region, , drop = FALSE])
}

format_markdown_table <- function(df) {
  if (requireNamespace("knitr", quietly = TRUE)) {
    return(paste(capture.output(knitr::kable(df, format = "pipe")), collapse = "\n"))
  }
  paste(capture.output(print(df, row.names = FALSE)), collapse = "\n")
}

strict_recommendations <- paste(vapply(seq_len(nrow(summary_table)), function(i) {
  row <- summary_table[i, ]
  if (row$BH_FDR_lt_010 > 0) {
    paste0("- ", row$brain_region, ": use BH-FDR < 0.10 for strict reporting; keep exploratory list for enrichment and ML.")
  } else if (row$Storey_q_lt_010 > 0) {
    paste0("- ", row$brain_region, ": Storey q < 0.10 recovers candidates, but label this as q-value based and verify biologically.")
  } else {
    paste0("- ", row$brain_region, ": no corrected list at 0.10; use exploratory/ranked lists only for exploratory enrichment, ML, and AHBA mapping.")
  }
}, character(1)), collapse = "\n")

report <- paste0(
  "# Fixed P-Value Correction Audit: Length > 100\n\n",
  "## Previous Issue\n\n",
  "- `qvalue` was missing, so Storey q-value was skipped.\n",
  "- `IHW` was missing, so independent hypothesis weighting was skipped.\n",
  "- Expression covariates were missing from the LMM table, so independent filtering could not run.\n\n",
  "## What Was Fixed\n\n",
  "- Added `code/setup_pvalue_audit_R_packages.R` to install/check `qvalue` and `IHW`.\n",
  "- Added expression covariates computed from normalized log2(GeTMM + 1) matrices: mean expression, expression variance, and detection rate.\n",
  "- Reran p-value corrections with BH, BY, Storey q-value, independent filtering + BH, independent filtering + Storey q-value, and IHW where available.\n\n",
  "## Preserved Baseline\n\n",
  "`p_value_diagnosis < 0.01 and abs(beta_diagnosis) > 0.2`\n\n",
  "The original exploratory DE report and tables were not overwritten.\n\n",
  "## Definitions\n\n",
  "- `beta_diagnosis` = effect size.\n",
  "- `p_value_diagnosis` = nominal p-value.\n",
  "- FDR / q-value / IHW adjusted p-value = multiple-testing correction.\n\n",
  "## Merge Check\n\n",
  format_markdown_table(merge_summary), "\n\n",
  "## Results By Region\n\n",
  format_markdown_table(summary_table), "\n\n",
  "## Method Availability\n\n",
  "- Storey q-value available: ", qvalue_available, "\n",
  "- IHW available: ", ihw_available, "\n\n",
  "## Interpretation\n\n",
  "This audit compares correction strategies; it does not replace the current exploratory DE baseline unless the comparison supports that decision. ",
  "Independent filtering now runs because the independent covariates are computed from expression matrices without using diagnosis, beta, direction, or p-values.\n\n",
  "## Recommendations\n\n",
  "1. Strict reporting list: use BH-FDR < 0.10 where available; otherwise state that no corrected DE genes were recovered.\n",
  "2. Exploratory DE list: keep using `p_value_diagnosis < 0.01 and abs(beta_diagnosis) > 0.2`.\n",
  "3. Enrichment-analysis input: use the exploratory DE list plus ranked full signatures.\n",
  "4. ML feature-selection input: use ranked `abs_beta_diagnosis` and `signed_logp_score`, with exploratory candidates as a focused candidate set.\n",
  "5. AHBA module-mapping input: use exploratory candidates and ranked signatures, clearly labeled as exploratory.\n\n",
  strict_recommendations, "\n"
)

writeLines(report, report_file)

message("Wrote fixed corrected table: ", corrected_file)
message("Wrote fixed summary: ", summary_file)
message("Wrote fixed report: ", report_file)
