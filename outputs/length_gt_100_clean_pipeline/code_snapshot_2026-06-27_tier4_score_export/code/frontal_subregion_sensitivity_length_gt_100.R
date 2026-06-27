## frontal_subregion_sensitivity_length_gt_100.R
##
## Purpose:
##   Split only the Length > 100 Frontal Cortex samples back into their
##   original GEO tissue labels and test whether the broad Frontal Cortex
##   grouping diluted the MDD differential-expression signal.
##
## Important:
##   This script does not redo normalization and does not overwrite the
##   four-group results. It subsets the existing Length > 100 log2(GeTMM + 1)
##   normalized Frontal Cortex matrix.

rm(list = ls())
options(stringsAsFactors = FALSE)

required_packages <- c("stats", "graphics", "utils")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

has_qvalue <- requireNamespace("qvalue", quietly = TRUE)
has_ihw <- requireNamespace("IHW", quietly = TRUE)

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(script_path)) {
  script_path <- normalizePath("code/frontal_subregion_sensitivity_length_gt_100.R", mustWork = FALSE)
}
repo_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

root <- file.path(repo_dir, "outputs", "length_gt_100_elisa_style")
step1_root <- file.path(root, "step1_from_raw", "Length_gt_100")
out_dir <- file.path(root, "frontal_subregion_sensitivity")
candidate_dir <- file.path(out_dir, "candidate_gene_lists")
matrix_dir <- file.path(out_dir, "subregion_normalized_expression_matrices")
plot_dir <- file.path(out_dir, "plots")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

metadata_file <- file.path(step1_root, "01_metadata", "keep_all_samples", "metadata_frontal_cortex.csv")
expression_file <- file.path(root, "normalized_expression_length_gt_100_by_region", "normalized_GeTMM_log2_frontal_cortex.csv")
grouped_lmm_file <- file.path(root, "LMM_all_genes_length_gt_100.csv")
grouped_standardized_file <- file.path(root, "standardized_effect_size_4groups", "LMM_all_genes_with_standardized_effect_size_4groups.csv")

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Required input file was not found: ", path)
  }
  TRUE
}
invisible(vapply(c(metadata_file, expression_file, grouped_lmm_file, grouped_standardized_file), stop_if_missing, logical(1)))

message("Reading metadata: ", metadata_file)
metadata <- read.csv(metadata_file, check.names = FALSE)

message("Reading normalized expression matrix: ", expression_file)
expression_matrix <- read.csv(expression_file, row.names = 1, check.names = FALSE)
expression_matrix <- as.matrix(expression_matrix)
storage.mode(expression_matrix) <- "numeric"

message("Reading grouped four-region LMM table: ", grouped_lmm_file)
grouped_lmm <- read.csv(grouped_lmm_file, check.names = FALSE)

message("Reading grouped standardized-effect table: ", grouped_standardized_file)
grouped_standardized <- read.csv(grouped_standardized_file, check.names = FALSE)

candidate_metadata_columns <- intersect(
  c(
    "tissue", "tissue_label", "geo_tissue", "original_tissue",
    "brain_region_original", "source_name_ch1", "characteristics",
    "exact_region", "OriginalRegion", "tissue:ch1"
  ),
  colnames(metadata)
)
message("Candidate original-tissue metadata columns:")
print(candidate_metadata_columns)

original_region_column <- if ("OriginalRegion" %in% colnames(metadata)) {
  "OriginalRegion"
} else if ("tissue" %in% colnames(metadata)) {
  "tissue"
} else if ("tissue:ch1" %in% colnames(metadata)) {
  "tissue:ch1"
} else {
  stop("No usable original tissue column found.")
}
message("Final original-tissue column used: ", original_region_column)

metadata_candidate_report <- data.frame(
  candidate_column = candidate_metadata_columns,
  used_for_original_frontal_subregion = candidate_metadata_columns == original_region_column,
  stringsAsFactors = FALSE
)
write.csv(
  metadata_candidate_report,
  file.path(out_dir, "metadata_original_region_candidate_columns.csv"),
  row.names = FALSE
)

required_metadata_columns <- c("SampleID", "Diagnosis", "Age", "Sex", "RIN", "BrainID", "BrainArea", original_region_column)
missing_metadata_columns <- setdiff(required_metadata_columns, colnames(metadata))
if (length(missing_metadata_columns) > 0) {
  stop("Metadata is missing required columns: ", paste(missing_metadata_columns, collapse = ", "))
}

frontal_metadata <- metadata[metadata$BrainArea == "Frontal_Cortex", , drop = FALSE]
frontal_metadata$original_frontal_subregion <- frontal_metadata[[original_region_column]]

frontal_grouping_map <- unique(data.frame(
  original_tissue_label = frontal_metadata$original_frontal_subregion,
  grouped_region = "Frontal Cortex",
  metadata_column_used = original_region_column,
  stringsAsFactors = FALSE
))
frontal_grouping_map <- frontal_grouping_map[order(frontal_grouping_map$original_tissue_label), , drop = FALSE]
write.csv(frontal_grouping_map, file.path(out_dir, "frontal_grouping_map.csv"), row.names = FALSE)

count_levels <- function(values) {
  tab <- table(values, useNA = "ifany")
  paste(paste(names(tab), as.integer(tab), sep = "="), collapse = "; ")
}

sample_count_rows <- lapply(split(frontal_metadata, frontal_metadata$original_frontal_subregion), function(x) {
  data.frame(
    original_frontal_subregion = x$original_frontal_subregion[1],
    total_samples = nrow(x),
    MDD_samples = sum(x$Diagnosis == "MDD", na.rm = TRUE),
    Control_samples = sum(x$Diagnosis %in% c("CTRL", "Control"), na.rm = TRUE),
    sex_counts = if ("Sex" %in% colnames(x)) count_levels(x$Sex) else NA_character_,
    unique_subjects = if ("BrainID" %in% colnames(x)) length(unique(x$BrainID)) else NA_integer_,
    batch_counts = if ("batch" %in% colnames(x)) count_levels(x$batch) else NA_character_,
    repeated_subjects = if ("BrainID" %in% colnames(x)) sum(table(x$BrainID) > 1) else NA_integer_,
    stringsAsFactors = FALSE
  )
})
sample_counts <- do.call(rbind, sample_count_rows)
sample_counts <- sample_counts[order(sample_counts$original_frontal_subregion), , drop = FALSE]
write.csv(sample_counts, file.path(out_dir, "frontal_subregion_sample_counts.csv"), row.names = FALSE)

primary_pass <- sample_counts$MDD_samples >= 5 & sample_counts$Control_samples >= 5
relaxed_pass <- sample_counts$MDD_samples >= 3 & sample_counts$Control_samples >= 3
sample_counts$passes_primary_5_per_group <- primary_pass
sample_counts$passes_relaxed_3_per_group <- relaxed_pass
write.csv(sample_counts, file.path(out_dir, "frontal_subregion_sample_counts.csv"), row.names = FALSE)

if (any(primary_pass)) {
  analyzable_subregions <- sample_counts$original_frontal_subregion[primary_pass]
  sample_rule_used <- "primary: at least 5 MDD and at least 5 Control samples"
} else if (any(relaxed_pass)) {
  analyzable_subregions <- sample_counts$original_frontal_subregion[relaxed_pass]
  sample_rule_used <- "relaxed: at least 3 MDD and at least 3 Control samples"
} else {
  analyzable_subregions <- character()
  sample_rule_used <- "none: no subregion had at least 3 MDD and at least 3 Control samples"
}
message("Sample rule used: ", sample_rule_used)
message("Analyzable subregions:")
print(analyzable_subregions)

if (length(analyzable_subregions) == 0) {
  stop("No frontal subregions are analyzable under the primary or relaxed sample rule.")
}

common_samples <- intersect(colnames(expression_matrix), frontal_metadata$SampleID)
if (length(common_samples) == 0) {
  stop("No sample IDs matched between frontal expression matrix and metadata.")
}
expression_matrix <- expression_matrix[, common_samples, drop = FALSE]
frontal_metadata <- frontal_metadata[match(common_samples, frontal_metadata$SampleID), , drop = FALSE]
if (!all(colnames(expression_matrix) == frontal_metadata$SampleID)) {
  stop("Expression matrix and metadata are not aligned.")
}

gene_info <- unique(grouped_lmm[, intersect(c("gene_id", "gene_symbol", "gene_length"), colnames(grouped_lmm)), drop = FALSE])
gene_info$gene_id <- as.character(gene_info$gene_id)

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

standard_direction <- function(beta) {
  ifelse(is.na(beta), NA_character_, ifelse(beta > 0, "up_in_MDD", ifelse(beta < 0, "down_in_MDD", "no_change")))
}

as_true_flag <- function(x) {
  if (is.logical(x)) {
    return(!is.na(x) & x)
  }
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}

normalize_region_name <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

adjust_storey <- function(p_values) {
  out <- rep(NA_real_, length(p_values))
  valid <- is.finite(p_values) & !is.na(p_values) & p_values >= 0 & p_values <= 1
  if (!has_qvalue || sum(valid) < 2) {
    return(out)
  }
  q <- tryCatch(qvalue::qvalue(p_values[valid])$qvalues, error = function(e) rep(NA_real_, sum(valid)))
  out[valid] <- q
  out
}

adjust_ihw <- function(p_values, covariate) {
  out <- rep(NA_real_, length(p_values))
  valid <- is.finite(p_values) & !is.na(p_values) & p_values >= 0 & p_values <= 1 &
    is.finite(covariate) & !is.na(covariate)
  if (!has_ihw || sum(valid) < 50 || length(unique(covariate[valid])) < 2) {
    return(out)
  }
  adj <- tryCatch(
    IHW::adj_pvalues(IHW::ihw(p_values[valid] ~ covariate[valid], alpha = 0.1)),
    error = function(e) rep(NA_real_, sum(valid))
  )
  out[valid] <- adj
  out
}

fit_gene_lm_matrix <- function(y_matrix, model_metadata) {
  model_metadata$Diagnosis <- factor(model_metadata$Diagnosis, levels = c("CTRL", "MDD"))
  model_metadata$Sex <- factor(model_metadata$Sex)
  model_metadata$Age <- suppressWarnings(as.numeric(model_metadata$Age))
  model_metadata$RIN <- suppressWarnings(as.numeric(model_metadata$RIN))

  complete_rows <- complete.cases(model_metadata[, c("Diagnosis", "Age", "Sex", "RIN")])
  model_metadata <- model_metadata[complete_rows, , drop = FALSE]
  y_matrix <- y_matrix[, model_metadata$SampleID, drop = FALSE]

  design <- model.matrix(~ Diagnosis + Age + Sex + RIN, data = model_metadata)
  diagnosis_column <- "DiagnosisMDD"
  if (!(diagnosis_column %in% colnames(design))) {
    stop("DiagnosisMDD coefficient was not present in the subregion design matrix.")
  }

  qr_design <- qr(design)
  rank_design <- qr_design$rank
  df_residual <- nrow(design) - rank_design
  if (df_residual <= 0) {
    stop("No residual degrees of freedom for subregion model.")
  }

  fit <- lm.fit(design, t(y_matrix))
  coefficients <- fit$coefficients
  fitted_values <- design %*% coefficients
  residuals <- t(y_matrix) - fitted_values
  rss <- colSums(residuals^2, na.rm = TRUE)
  sigma2 <- rss / df_residual
  xtx_inv <- chol2inv(qr.R(qr_design))
  diagnosis_index <- match(diagnosis_column, colnames(design))
  beta <- coefficients[diagnosis_index, ]
  standard_error <- sqrt(sigma2 * xtx_inv[diagnosis_index, diagnosis_index])
  statistic <- beta / standard_error
  p_value <- 2 * stats::pt(abs(statistic), df = df_residual, lower.tail = FALSE)

  data.frame(
    gene_id = rownames(y_matrix),
    beta_diagnosis = as.numeric(beta),
    standard_error_diagnosis = as.numeric(standard_error),
    statistic_diagnosis = as.numeric(statistic),
    p_value_diagnosis = as.numeric(p_value),
    model_formula_used = "expression_value ~ Diagnosis + Age + Sex + RIN",
    model_note = "BrainID random effect not estimable within subregion because each BrainID appears once; OriginalRegion fixed effect is constant after splitting.",
    residual_degrees_of_freedom = df_residual,
    stringsAsFactors = FALSE
  )
}

plot_png <- function(path, code) {
  png(path, width = 1200, height = 900)
  tryCatch(code(), finally = dev.off())
}

make_plots <- function(results, subregion_safe, subregion_label) {
  valid <- is.finite(results$p_value_diagnosis) & results$p_value_diagnosis > 0
  neg_logp <- rep(NA_real_, nrow(results))
  neg_logp[valid] <- -log10(results$p_value_diagnosis[valid])

  plot_png(file.path(plot_dir, paste0(subregion_safe, "_volcano_beta_highlight_Tier3.png")), function() {
    cols <- ifelse(results$Tier3_current_baseline_p001_absBeta02, "#D55E00", "gray70")
    plot(results$beta_diagnosis, neg_logp, pch = 16, cex = 0.45, col = cols,
         xlab = "beta_diagnosis", ylab = "-log10(p_value_diagnosis)",
         main = paste("Tier 3 beta volcano:", subregion_label))
    abline(h = -log10(0.01), lty = 2, col = "gray30")
    abline(v = c(-0.2, 0.2), lty = 2, col = "gray30")
  })

  plot_png(file.path(plot_dir, paste0(subregion_safe, "_volcano_standardized_highlight_Tier2.png")), function() {
    cols <- ifelse(results$Tier2_nominal_p001_absStdES08, "#0072B2", "gray70")
    plot(results$standardized_effect_size_diagnosis, neg_logp, pch = 16, cex = 0.45, col = cols,
         xlab = "standardized_effect_size_diagnosis", ylab = "-log10(p_value_diagnosis)",
         main = paste("Tier 2 standardized volcano:", subregion_label))
    abline(h = -log10(0.01), lty = 2, col = "gray30")
    abline(v = c(-0.8, 0.8), lty = 2, col = "gray30")
  })

  plot_png(file.path(plot_dir, paste0(subregion_safe, "_histogram_beta_diagnosis.png")), function() {
    hist(results$beta_diagnosis, breaks = 60, col = "gray80", border = "white",
         main = paste("beta_diagnosis:", subregion_label), xlab = "beta_diagnosis")
  })

  plot_png(file.path(plot_dir, paste0(subregion_safe, "_histogram_standardized_effect_size.png")), function() {
    hist(results$standardized_effect_size_diagnosis, breaks = 60, col = "gray80", border = "white",
         main = paste("standardized effect size:", subregion_label),
         xlab = "standardized_effect_size_diagnosis")
  })

  plot_png(file.path(plot_dir, paste0(subregion_safe, "_histogram_p_values.png")), function() {
    hist(results$p_value_diagnosis[valid], breaks = 50, col = "gray80", border = "white",
         main = paste("Diagnosis p-values:", subregion_label), xlab = "p_value_diagnosis")
  })

  plot_png(file.path(plot_dir, paste0(subregion_safe, "_qq_plot_p_values.png")), function() {
    p <- sort(results$p_value_diagnosis[valid])
    expected <- -log10(ppoints(length(p)))
    observed <- -log10(p)
    plot(expected, observed, pch = 16, cex = 0.5, col = "gray50",
         xlab = "Expected -log10(p)", ylab = "Observed -log10(p)",
         main = paste("QQ plot:", subregion_label))
    abline(0, 1, col = "#D55E00", lwd = 2)
  })

  plot_png(file.path(plot_dir, paste0(subregion_safe, "_sd_expression_vs_abs_beta.png")), function() {
    plot(results$sd_expression_gene, abs(results$beta_diagnosis), pch = 16, cex = 0.45, col = "gray60",
         xlab = "sd_expression_gene", ylab = "abs(beta_diagnosis)",
         main = paste("SD vs abs beta:", subregion_label))
  })

  plot_png(file.path(plot_dir, paste0(subregion_safe, "_sd_expression_vs_abs_standardized_effect_size.png")), function() {
    plot(results$sd_expression_gene, results$abs_standardized_effect_size_diagnosis,
         pch = 16, cex = 0.45, col = "gray60",
         xlab = "sd_expression_gene", ylab = "abs_standardized_effect_size_diagnosis",
         main = paste("SD vs abs standardized effect:", subregion_label))
  })
}

all_results <- list()
sign_check_rows <- list()
summary_rows <- list()
model_feasibility_rows <- list()

for (subregion in analyzable_subregions) {
  message("")
  message("============================================================")
  message("Analyzing frontal subregion: ", subregion)
  message("============================================================")

  sub_meta <- frontal_metadata[frontal_metadata$original_frontal_subregion == subregion, , drop = FALSE]
  sub_samples <- intersect(colnames(expression_matrix), sub_meta$SampleID)
  sub_meta <- sub_meta[match(sub_samples, sub_meta$SampleID), , drop = FALSE]
  sub_expr <- expression_matrix[, sub_samples, drop = FALSE]

  sub_safe <- safe_filename(subregion)
  write.csv(
    data.frame(gene_id = rownames(sub_expr), sub_expr, check.names = FALSE),
    file.path(matrix_dir, paste0(sub_safe, "_normalized_GeTMM_log2_length_gt_100.csv")),
    row.names = FALSE
  )

  repeated_brainids <- sum(table(sub_meta$BrainID) > 1)
  lmm_random_effect_estimable <- repeated_brainids > 0 && length(unique(sub_meta$BrainID)) < nrow(sub_meta)

  model_feasibility_rows[[subregion]] <- data.frame(
    original_frontal_subregion = subregion,
    n_samples = nrow(sub_meta),
    unique_BrainID = length(unique(sub_meta$BrainID)),
    repeated_BrainID = repeated_brainids,
    lmm_random_effect_estimable = lmm_random_effect_estimable,
    model_used = ifelse(
      lmm_random_effect_estimable,
      "LMM would be estimable, but this dataset did not trigger that case.",
      "fixed-effect LM on normalized expression because BrainID random effect is not estimable within this subregion"
    ),
    stringsAsFactors = FALSE
  )

  results <- fit_gene_lm_matrix(sub_expr, sub_meta)
  rownames(results) <- NULL
  results$original_frontal_subregion <- subregion
  results$n_samples <- nrow(sub_meta)
  results$n_MDD <- sum(sub_meta$Diagnosis == "MDD", na.rm = TRUE)
  results$n_Control <- sum(sub_meta$Diagnosis == "CTRL", na.rm = TRUE)
  results$direction <- standard_direction(results$beta_diagnosis)

  sub_expr_aligned <- sub_expr[results$gene_id, sub_meta$SampleID, drop = FALSE]
  results$sd_expression_gene <- apply(sub_expr_aligned, 1, stats::sd, na.rm = TRUE)
  results$mean_expression_gene <- rowMeans(sub_expr_aligned, na.rm = TRUE)
  results$n_samples_used_for_sd <- rowSums(!is.na(sub_expr_aligned))
  results$standardized_effect_size_diagnosis <- ifelse(
    is.finite(results$sd_expression_gene) & results$sd_expression_gene > 0,
    results$beta_diagnosis / results$sd_expression_gene,
    NA_real_
  )
  results$abs_standardized_effect_size_diagnosis <- abs(results$standardized_effect_size_diagnosis)
  results$abs_beta_diagnosis <- abs(results$beta_diagnosis)
  results$signed_logp_score <- sign(results$beta_diagnosis) * -log10(pmax(results$p_value_diagnosis, .Machine$double.xmin))
  results$signed_standardized_logp_score <- sign(results$standardized_effect_size_diagnosis) * -log10(pmax(results$p_value_diagnosis, .Machine$double.xmin))
  results$fdr_BH <- p.adjust(results$p_value_diagnosis, method = "BH")
  results$fdr_BY <- p.adjust(results$p_value_diagnosis, method = "BY")
  results$qvalue_Storey <- adjust_storey(results$p_value_diagnosis)
  results$ihw_adj_pvalue <- adjust_ihw(results$p_value_diagnosis, results$mean_expression_gene)

  results$Tier1_Storey_q010_absStdES08 <- results$qvalue_Storey < 0.10 & results$abs_standardized_effect_size_diagnosis > 0.8
  results$Tier2_nominal_p001_absStdES08 <- results$p_value_diagnosis < 0.01 & results$abs_standardized_effect_size_diagnosis > 0.8
  results$Tier3_current_baseline_p001_absBeta02 <- results$p_value_diagnosis < 0.01 & results$abs_beta_diagnosis > 0.2
  results$BH005_absStdES08 <- results$fdr_BH < 0.05 & results$abs_standardized_effect_size_diagnosis > 0.8
  results$BH010_absStdES08 <- results$fdr_BH < 0.10 & results$abs_standardized_effect_size_diagnosis > 0.8

  sign_valid <- is.finite(results$beta_diagnosis) & is.finite(results$standardized_effect_size_diagnosis) &
    results$beta_diagnosis != 0 & results$standardized_effect_size_diagnosis != 0
  sign_matches <- sign(results$beta_diagnosis[sign_valid]) == sign(results$standardized_effect_size_diagnosis[sign_valid])
  sign_check_rows[[subregion]] <- data.frame(
    original_frontal_subregion = subregion,
    rows_checked = sum(sign_valid),
    sign_mismatches = sum(!sign_matches),
    sign_check_passed = all(sign_matches),
    stringsAsFactors = FALSE
  )

  results <- merge(results, gene_info, by = "gene_id", all.x = TRUE, sort = FALSE)
  ordered_columns <- c(
    "gene_id", "gene_symbol", "gene_length", "original_frontal_subregion",
    "beta_diagnosis", "standard_error_diagnosis", "statistic_diagnosis",
    "p_value_diagnosis", "fdr_BH", "fdr_BY", "qvalue_Storey", "ihw_adj_pvalue",
    "standardized_effect_size_diagnosis", "abs_standardized_effect_size_diagnosis",
    "sd_expression_gene", "mean_expression_gene", "n_samples_used_for_sd",
    "abs_beta_diagnosis", "signed_logp_score", "signed_standardized_logp_score",
    "direction", "n_samples", "n_MDD", "n_Control",
    "Tier1_Storey_q010_absStdES08", "Tier2_nominal_p001_absStdES08",
    "Tier3_current_baseline_p001_absBeta02", "BH005_absStdES08", "BH010_absStdES08",
    "model_formula_used", "model_note", "residual_degrees_of_freedom"
  )
  results <- results[, c(ordered_columns, setdiff(colnames(results), ordered_columns)), drop = FALSE]

  tier_files <- list(
    Tier1_Storey_q010_absStdES08 = results[results$Tier1_Storey_q010_absStdES08 %in% TRUE, ],
    Tier2_nominal_p001_absStdES08 = results[results$Tier2_nominal_p001_absStdES08 %in% TRUE, ],
    Tier3_current_baseline_p001_absBeta02 = results[results$Tier3_current_baseline_p001_absBeta02 %in% TRUE, ],
    BH005_absStdES08 = results[results$BH005_absStdES08 %in% TRUE, ],
    BH010_absStdES08 = results[results$BH010_absStdES08 %in% TRUE, ]
  )
  for (tier_name in names(tier_files)) {
    tier_df <- tier_files[[tier_name]]
    if (nrow(tier_df) > 0) {
      tier_df <- tier_df[order(tier_df$p_value_diagnosis, -tier_df$abs_standardized_effect_size_diagnosis, na.last = TRUE), ]
    }
    write.csv(
      tier_df,
      file.path(candidate_dir, paste0(sub_safe, "_", tier_name, ".csv")),
      row.names = FALSE
    )
  }

  summary_rows[[subregion]] <- data.frame(
    original_frontal_subregion = subregion,
    n_samples = nrow(sub_meta),
    n_MDD = sum(sub_meta$Diagnosis == "MDD", na.rm = TRUE),
    n_Control = sum(sub_meta$Diagnosis == "CTRL", na.rm = TRUE),
    genes_tested = nrow(results),
    Tier1_Storey_q010_absStdES08 = sum(results$Tier1_Storey_q010_absStdES08, na.rm = TRUE),
    Tier2_nominal_p001_absStdES08 = sum(results$Tier2_nominal_p001_absStdES08, na.rm = TRUE),
    Tier3_current_baseline_p001_absBeta02 = sum(results$Tier3_current_baseline_p001_absBeta02, na.rm = TRUE),
    BH005_absStdES08 = sum(results$BH005_absStdES08, na.rm = TRUE),
    BH010_absStdES08 = sum(results$BH010_absStdES08, na.rm = TRUE),
    median_abs_beta = median(results$abs_beta_diagnosis, na.rm = TRUE),
    max_abs_beta = max(results$abs_beta_diagnosis, na.rm = TRUE),
    median_abs_standardized_effect_size = median(results$abs_standardized_effect_size_diagnosis, na.rm = TRUE),
    max_abs_standardized_effect_size = max(results$abs_standardized_effect_size_diagnosis, na.rm = TRUE),
    up_in_MDD_Tier2 = sum(results$Tier2_nominal_p001_absStdES08 & results$direction == "up_in_MDD", na.rm = TRUE),
    down_in_MDD_Tier2 = sum(results$Tier2_nominal_p001_absStdES08 & results$direction == "down_in_MDD", na.rm = TRUE),
    up_in_MDD_Tier3 = sum(results$Tier3_current_baseline_p001_absBeta02 & results$direction == "up_in_MDD", na.rm = TRUE),
    down_in_MDD_Tier3 = sum(results$Tier3_current_baseline_p001_absBeta02 & results$direction == "down_in_MDD", na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  make_plots(results, sub_safe, subregion)
  all_results[[subregion]] <- results
}

all_results_df <- do.call(rbind, all_results)
rownames(all_results_df) <- NULL
write.csv(
  all_results_df,
  file.path(out_dir, "frontal_subregion_LMM_all_genes_with_standardized_effect_size.csv"),
  row.names = FALSE
)

sign_check <- do.call(rbind, sign_check_rows)
write.csv(sign_check, file.path(out_dir, "frontal_subregion_standardized_effect_size_sign_check.csv"), row.names = FALSE)

model_feasibility <- do.call(rbind, model_feasibility_rows)
write.csv(model_feasibility, file.path(out_dir, "frontal_subregion_model_feasibility.csv"), row.names = FALSE)

summary_df <- do.call(rbind, summary_rows)
summary_df <- summary_df[order(summary_df$original_frontal_subregion), , drop = FALSE]
write.csv(summary_df, file.path(out_dir, "frontal_subregion_summary.csv"), row.names = FALSE)

grouped_region_norm <- normalize_region_name(grouped_standardized$grouped_region)
brain_region_norm <- normalize_region_name(grouped_standardized$brain_region)
grouped_frontal <- grouped_standardized[grouped_region_norm == "frontal_cortex" | brain_region_norm == "frontal_cortex", , drop = FALSE]
grouped_tier2 <- grouped_frontal$gene_id[as_true_flag(grouped_frontal$nominal_p001_absStdES08)]
grouped_tier3 <- grouped_frontal$gene_id[as_true_flag(grouped_frontal$current_baseline_p001_absBeta02)]
grouped_top100_p <- head(grouped_frontal$gene_id[order(grouped_frontal$p_value_diagnosis, na.last = TRUE)], 100)
grouped_top100_signed_std <- head(grouped_frontal$gene_id[order(-abs(grouped_frontal$signed_standardized_logp_score), na.last = TRUE)], 100)

overlap_rows <- list()
for (subregion in names(all_results)) {
  sub_df <- all_results[[subregion]]
  sub_tier2 <- sub_df$gene_id[sub_df$Tier2_nominal_p001_absStdES08 %in% TRUE]
  sub_tier3 <- sub_df$gene_id[sub_df$Tier3_current_baseline_p001_absBeta02 %in% TRUE]
  sub_top100_p <- head(sub_df$gene_id[order(sub_df$p_value_diagnosis, na.last = TRUE)], 100)
  sub_top100_signed_std <- head(sub_df$gene_id[order(-abs(sub_df$signed_standardized_logp_score), na.last = TRUE)], 100)

  overlap_rows[[subregion]] <- data.frame(
    original_frontal_subregion = subregion,
    grouped_Frontal_Cortex_Tier2_genes = length(grouped_tier2),
    subregion_Tier2_genes = length(sub_tier2),
    Tier2_overlap_genes = length(intersect(grouped_tier2, sub_tier2)),
    grouped_Frontal_Cortex_Tier3_genes = length(grouped_tier3),
    subregion_Tier3_genes = length(sub_tier3),
    Tier3_overlap_genes = length(intersect(grouped_tier3, sub_tier3)),
    top100_pvalue_overlap_genes = length(intersect(grouped_top100_p, sub_top100_p)),
    top100_signed_standardized_logp_overlap_genes = length(intersect(grouped_top100_signed_std, sub_top100_signed_std)),
    Tier2_overlap_gene_symbols = paste(na.omit(sub_df$gene_symbol[sub_df$gene_id %in% intersect(grouped_tier2, sub_tier2)]), collapse = "; "),
    Tier3_overlap_gene_symbols = paste(na.omit(sub_df$gene_symbol[sub_df$gene_id %in% intersect(grouped_tier3, sub_tier3)]), collapse = "; "),
    stringsAsFactors = FALSE
  )
}
overlap_df <- do.call(rbind, overlap_rows)
write.csv(overlap_df, file.path(out_dir, "frontal_subregion_vs_grouped_overlap.csv"), row.names = FALSE)

top20_p_rows <- list()
top20_std_rows <- list()
for (subregion in names(all_results)) {
  sub_df <- all_results[[subregion]]
  p_candidates <- sub_df[sub_df$abs_standardized_effect_size_diagnosis > 0.8, , drop = FALSE]
  p_candidates <- p_candidates[order(p_candidates$p_value_diagnosis, na.last = TRUE), ]
  top20_p_rows[[subregion]] <- head(p_candidates, 20)
  std_candidates <- sub_df[order(-sub_df$abs_standardized_effect_size_diagnosis, sub_df$p_value_diagnosis, na.last = TRUE), ]
  top20_std_rows[[subregion]] <- head(std_candidates, 20)
}
top20_p_df <- do.call(rbind, top20_p_rows)
top20_std_df <- do.call(rbind, top20_std_rows)
write.csv(top20_p_df, file.path(out_dir, "top20_by_pvalue_absStdES08_by_subregion.csv"), row.names = FALSE)
write.csv(top20_std_df, file.path(out_dir, "top20_by_abs_standardized_effect_size_by_subregion.csv"), row.names = FALSE)

stronger_than_grouped <- summary_df$Tier2_nominal_p001_absStdES08 > 1
any_stronger <- any(stronger_than_grouped)

report_path <- file.path(out_dir, "Frontal_Subregion_Sensitivity_Report.md")
report <- c(
  "# Frontal Cortex Subregion Sensitivity Report",
  "",
  "## Scope",
  "",
  "The original four-group Length > 100 result was not overwritten. This analysis writes only to `outputs/length_gt_100_elisa_style/frontal_subregion_sensitivity/` plus the standalone script in `code/`.",
  "",
  "This sensitivity analysis was run because the grouped Frontal Cortex result had 99 current-baseline raw-beta genes but only 1 nominal standardized-effect-size gene and 0 corrected standardized-effect-size genes.",
  "",
  "## Original Tissue Labels Grouped Into Frontal Cortex",
  "",
  paste(capture.output(print(frontal_grouping_map[, c("original_tissue_label", "grouped_region")], row.names = FALSE)), collapse = "\n"),
  "",
  "## Metadata Column Used",
  "",
  paste0("Candidate columns printed by the script: `", paste(candidate_metadata_columns, collapse = "`, `"), "`."),
  paste0("Final column used: `", original_region_column, "`."),
  "",
  "## Sample Counts",
  "",
  paste(capture.output(print(sample_counts, row.names = FALSE)), collapse = "\n"),
  "",
  "## Model Feasibility",
  "",
  paste0("Sample rule used: ", sample_rule_used, "."),
  "",
  "Inside each original frontal subregion, `OriginalRegion` is constant. In this dataset every `BrainID` appears once per original frontal subregion, so the `(1 | BrainID)` random effect from the broad Frontal Cortex LMM is not estimable after splitting. The script therefore records the LMM feasibility check and uses the closest available normalized-expression diagnosis model: `expression_value ~ Diagnosis + Age + Sex + RIN`. This avoids silently running impossible mixed models.",
  "",
  paste(capture.output(print(model_feasibility, row.names = FALSE)), collapse = "\n"),
  "",
  "## Sign Check",
  "",
  paste(capture.output(print(sign_check, row.names = FALSE)), collapse = "\n"),
  "",
  "The sign check confirms whether `sign(beta_diagnosis)` matches `sign(standardized_effect_size_diagnosis)` for every row with valid gene-level SD.",
  "",
  "## Evidence Tier Summary",
  "",
  paste(capture.output(print(summary_df, row.names = FALSE)), collapse = "\n"),
  "",
  "Tier definitions:",
  "",
  "- Tier 1: Storey q < 0.10 and abs_standardized_effect_size_diagnosis > 0.8.",
  "- Tier 2: p_value_diagnosis < 0.01 and abs_standardized_effect_size_diagnosis > 0.8.",
  "- Tier 3: p_value_diagnosis < 0.01 and abs(beta_diagnosis) > 0.2.",
  "",
  "Tier 2 and Tier 3 are exploratory, nominal candidate signatures. They are not final significant DE genes.",
  "",
  "## Comparison With Grouped Frontal Cortex",
  "",
  paste(capture.output(print(overlap_df, row.names = FALSE)), collapse = "\n"),
  "",
  if (any_stronger) {
    "At least one original frontal subregion recovers a stronger nominal standardized-effect-size signal than the grouped Frontal Cortex result. This supports the interpretation that broad grouping can dilute a frontal signal."
  } else {
    "No original frontal subregion recovered a stronger nominal standardized-effect-size signal than the grouped Frontal Cortex result. This does not support grouping dilution as the main explanation."
  },
  "",
  "## Top 20 Genes By P-Value Among abs_standardized_effect_size_diagnosis > 0.8",
  "",
  paste(capture.output(print(top20_p_df[, intersect(c("original_frontal_subregion", "gene_id", "gene_symbol", "beta_diagnosis", "standardized_effect_size_diagnosis", "p_value_diagnosis", "qvalue_Storey", "direction"), colnames(top20_p_df))], row.names = FALSE)), collapse = "\n"),
  "",
  "## Top 20 Genes By Absolute Standardized Effect Size",
  "",
  paste(capture.output(print(top20_std_df[, intersect(c("original_frontal_subregion", "gene_id", "gene_symbol", "beta_diagnosis", "standardized_effect_size_diagnosis", "p_value_diagnosis", "qvalue_Storey", "direction"), colnames(top20_std_df))], row.names = FALSE)), collapse = "\n"),
  "",
  "## Recommendation",
  "",
  "1. Keep the broad Frontal Cortex result as the main comparable four-group result.",
  "2. Report this frontal subregion sensitivity result as exploratory because the split removes the estimable random-effect structure used in the broad Frontal Cortex LMM.",
  "3. Use any Tier 2 subregion candidates and ranked signatures for exploratory enrichment, ML feature selection, and AHBA/module mapping, not as final significant DE genes.",
  "4. Ventral Striatum should remain the main robust region from the four-group analysis because it retained high exploratory counts and the only Storey q < 0.10 standardized-effect-size corrected genes.",
  "",
  "Interpretation terms to use: corrected genes for strict reporting, nominal standardized-effect-size genes for exploratory biology, and ranked signatures for ML and AHBA module mapping."
)
writeLines(report, report_path)

message("")
message("Saved frontal subregion sensitivity outputs to: ", out_dir)
message("Summary:")
print(summary_df)
message("Sign check:")
print(sign_check)
