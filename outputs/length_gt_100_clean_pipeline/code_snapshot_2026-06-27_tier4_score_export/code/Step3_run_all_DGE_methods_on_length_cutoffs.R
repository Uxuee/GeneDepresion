## Step3_run_all_DGE_methods_on_length_cutoffs.R
##
## Purpose:
##   Run differential expression methods on every Step 1 length-cutoff output.
##
## This script compares:
##
##   A. Transcriptomics-appropriate methods
##      - limma-voom for Insula, Ventral_Striatum, Hippocampal_Formation
##      - dream for Frontal_Cortex
##
##   B. Normalized-expression diagnostic comparator methods
##      - ordinary LM on log2(GeTMM + 1) for small groups
##      - LMM on log2(GeTMM + 1) for Frontal_Cortex
##      - fixed-effect LM on log2(GeTMM + 1) for Frontal_Cortex
##
## Important:
##   1. This script does not run PCA.
##   2. This script does not run the combined all-region dream model.
##   3. This script expects Step1_length_cutoff_sensitivity_GeTMM.R to have
##      already created the length-cutoff input folders.
##   4. This script is not optimized for speed. It is written for readability.

rm(list = ls())
options(stringsAsFactors = FALSE)

## -------------------------------------------------------------------------
## 0. Packages
## -------------------------------------------------------------------------

required_packages <- c(
  "edgeR",
  "limma",
  "variancePartition",
  "BiocParallel",
  "lme4",
  "lmerTest"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

library(edgeR)
library(limma)
library(variancePartition)
library(BiocParallel)

## -------------------------------------------------------------------------
## 1. User settings
## -------------------------------------------------------------------------

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(script_path)) {
  script_path <- normalizePath("code/Step3_run_all_DGE_methods_on_length_cutoffs.R", mustWork = FALSE)
}

repo_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

step1_sensitivity_folder <- Sys.getenv(
  "STEP3_LENGTH_INPUT_ROOT",
  unset = file.path(repo_dir, "outputs", "step1_length_sensitivity")
)
output_root <- Sys.getenv(
  "STEP3_LENGTH_OUTPUT_ROOT",
  unset = file.path(repo_dir, "outputs", "step3_length_sensitivity")
)

length_cutoffs_env <- Sys.getenv("STEP3_LENGTH_CUTOFFS", unset = "")
if (nzchar(length_cutoffs_env)) {
  length_cutoffs <- suppressWarnings(as.numeric(strsplit(length_cutoffs_env, ",")[[1]]))
  length_cutoffs <- length_cutoffs[!is.na(length_cutoffs)]
  if (length(length_cutoffs) == 0) {
    stop("STEP3_LENGTH_CUTOFFS was set, but no numeric cutoffs were found.")
  }
} else {
  length_cutoffs <- c(20, 50, 100, 200, 500)
}

brain_areas <- data.frame(
  brain_area = c("Frontal_Cortex", "Insula", "Ventral_Striatum", "Hippocampal_Formation"),
  file_stub = c("frontal_cortex", "insula", "ventral_striatum", "hippocampal_formation"),
  is_frontal = c(TRUE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

diagnosis_coefficient <- "DiagnosisMDD"

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------------
## 2. Small helper functions
## -------------------------------------------------------------------------

make_length_label <- function(length_cutoff) {
  paste0("Length_gt_", length_cutoff)
}

read_csv_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("Matrix file was not found: ", path)
  }

  matrix_data <- read.csv(
    path,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  matrix_data <- as.matrix(matrix_data)
  storage.mode(matrix_data) <- "numeric"

  if (any(is.na(matrix_data))) {
    stop("Matrix contains NA values after numeric conversion: ", path)
  }

  matrix_data
}

read_metadata <- function(path, required_columns) {
  if (!file.exists(path)) {
    stop("Metadata file was not found: ", path)
  }

  metadata <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  missing_columns <- setdiff(required_columns, colnames(metadata))
  if (length(missing_columns) > 0) {
    stop("Metadata is missing required columns: ", paste(missing_columns, collapse = ", "))
  }

  if (any(duplicated(metadata$SampleID))) {
    stop("Metadata contains duplicated SampleID values: ", path)
  }

  metadata
}

align_matrix_and_metadata <- function(expression_or_count_matrix, metadata) {
  if (any(duplicated(colnames(expression_or_count_matrix)))) {
    stop("Expression/count matrix contains duplicated sample columns.")
  }

  common_samples <- intersect(colnames(expression_or_count_matrix), metadata$SampleID)
  message("Samples found in both matrix and metadata: ", length(common_samples))

  if (length(common_samples) == 0) {
    stop("No sample IDs matched between matrix columns and metadata rows.")
  }

  expression_or_count_matrix <- expression_or_count_matrix[, common_samples, drop = FALSE]
  metadata <- metadata[match(common_samples, metadata$SampleID), , drop = FALSE]

  message("First matrix sample IDs after matching:")
  print(head(colnames(expression_or_count_matrix)))

  message("First metadata SampleID values after matching:")
  print(head(metadata$SampleID))

  if (!all(colnames(expression_or_count_matrix) == metadata$SampleID)) {
    stop("Matrix columns and metadata rows are not aligned after matching.")
  }

  list(
    matrix = expression_or_count_matrix,
    metadata = metadata
  )
}

prepare_small_group_metadata <- function(metadata) {
  diagnosis_values <- sort(unique(metadata$Diagnosis))
  message("Diagnosis values found:")
  print(diagnosis_values)

  if (!identical(diagnosis_values, c("CTRL", "MDD"))) {
    stop("Expected Diagnosis values to be exactly CTRL and MDD.")
  }

  metadata$Diagnosis <- factor(metadata$Diagnosis, levels = c("CTRL", "MDD"))
  metadata$Sex <- factor(metadata$Sex)
  metadata$Age <- suppressWarnings(as.numeric(metadata$Age))
  metadata$RIN <- suppressWarnings(as.numeric(metadata$RIN))

  message("Diagnosis counts:")
  print(table(metadata$Diagnosis, useNA = "ifany"))

  message("Missing values in model covariates:")
  print(colSums(is.na(metadata[, c("Diagnosis", "Age", "Sex", "RIN")])))

  complete_rows <- complete.cases(metadata[, c("Diagnosis", "Age", "Sex", "RIN")])
  message("Samples before removing missing covariates: ", nrow(metadata))
  message("Samples after removing missing covariates:  ", sum(complete_rows))

  metadata[complete_rows, , drop = FALSE]
}

prepare_frontal_metadata <- function(metadata) {
  diagnosis_values <- sort(unique(metadata$Diagnosis))
  message("Diagnosis values found:")
  print(diagnosis_values)

  if (!identical(diagnosis_values, c("CTRL", "MDD"))) {
    stop("Expected Diagnosis values to be exactly CTRL and MDD.")
  }

  metadata$Diagnosis <- factor(metadata$Diagnosis, levels = c("CTRL", "MDD"))
  metadata$Sex <- factor(metadata$Sex)
  metadata$OriginalRegion <- factor(metadata$OriginalRegion)
  metadata$BrainID <- factor(metadata$BrainID)
  metadata$Age <- suppressWarnings(as.numeric(metadata$Age))
  metadata$RIN <- suppressWarnings(as.numeric(metadata$RIN))

  message("Diagnosis counts:")
  print(table(metadata$Diagnosis, useNA = "ifany"))

  message("OriginalRegion by Diagnosis:")
  print(table(metadata$OriginalRegion, metadata$Diagnosis, useNA = "ifany"))

  message("Number of samples: ", nrow(metadata))
  message("Number of unique BrainID donors: ", length(unique(metadata$BrainID)))
  message("Number of repeated BrainID donors: ", sum(table(metadata$BrainID) > 1))

  message("Missing values in model covariates:")
  print(colSums(is.na(metadata[, c("Diagnosis", "Age", "Sex", "RIN", "OriginalRegion", "BrainID")])))

  complete_rows <- complete.cases(metadata[, c("Diagnosis", "Age", "Sex", "RIN", "OriginalRegion", "BrainID")])
  message("Samples before removing missing covariates: ", nrow(metadata))
  message("Samples after removing missing covariates:  ", sum(complete_rows))

  metadata[complete_rows, , drop = FALSE]
}

save_histogram_and_top_tables <- function(results, output_folder, file_prefix, plot_title) {
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

  results_file <- file.path(output_folder, paste0(file_prefix, "_results.csv"))
  top_raw_file <- file.path(output_folder, paste0(file_prefix, "_top100_by_raw_p.csv"))
  top_fdr_file <- file.path(output_folder, paste0(file_prefix, "_top100_by_FDR.csv"))
  histogram_file <- file.path(output_folder, paste0(file_prefix, "_diagnosis_pvalue_histogram.png"))

  write.csv(results, results_file, row.names = FALSE)

  results_by_raw_p <- results[order(results$diagnosis_p_value, na.last = TRUE), , drop = FALSE]
  results_by_fdr <- results[order(results$diagnosis_FDR_BH, na.last = TRUE), , drop = FALSE]

  write.csv(head(results_by_raw_p, 100), top_raw_file, row.names = FALSE)
  write.csv(head(results_by_fdr, 100), top_fdr_file, row.names = FALSE)

  valid_p_values <- results$diagnosis_p_value[!is.na(results$diagnosis_p_value)]

  if (length(valid_p_values) == 0) {
    warning("No valid p-values were available for histogram: ", file_prefix)
  } else {
    png(filename = histogram_file, width = 900, height = 700)
    hist(
      valid_p_values,
      breaks = 50,
      main = plot_title,
      xlab = "Raw p-value for DiagnosisMDD",
      col = "gray80",
      border = "white"
    )
    dev.off()
  }

  message("Saved full results to: ", results_file)
  message("Saved top raw-p table to: ", top_raw_file)
  message("Saved top FDR table to: ", top_fdr_file)
}

make_summary_row <- function(length_cutoff, brain_area, method, results) {
  valid_p <- results$diagnosis_p_value[!is.na(results$diagnosis_p_value)]
  valid_fdr <- results$diagnosis_FDR_BH[!is.na(results$diagnosis_FDR_BH)]

  if (length(valid_p) == 0) {
    minimum_p_value <- NA_real_
    top_gene <- NA_character_
  } else {
    minimum_p_value <- min(valid_p)
    top_gene <- results$GeneID[which.min(results$diagnosis_p_value)]
  }

  if (length(valid_fdr) == 0) {
    minimum_fdr <- NA_real_
  } else {
    minimum_fdr <- min(valid_fdr)
  }

  convergence_count <- NA_integer_
  singular_count <- NA_integer_

  if ("convergence_warning" %in% colnames(results)) {
    convergence_count <- sum(results$convergence_warning, na.rm = TRUE)
  }

  if ("singular_fit_warning" %in% colnames(results)) {
    singular_count <- sum(results$singular_fit_warning, na.rm = TRUE)
  }

  data.frame(
    length_cutoff = length_cutoff,
    brain_area = brain_area,
    method = method,
    genes_tested = nrow(results),
    raw_p_lt_0.05 = sum(results$diagnosis_p_value < 0.05, na.rm = TRUE),
    raw_p_lt_0.01 = sum(results$diagnosis_p_value < 0.01, na.rm = TRUE),
    raw_p_lt_0.001 = sum(results$diagnosis_p_value < 0.001, na.rm = TRUE),
    FDR_lt_0.05 = sum(results$diagnosis_FDR_BH < 0.05, na.rm = TRUE),
    FDR_lt_0.10 = sum(results$diagnosis_FDR_BH < 0.10, na.rm = TRUE),
    minimum_p_value = minimum_p_value,
    minimum_FDR = minimum_fdr,
    top_gene = top_gene,
    convergence_warnings = convergence_count,
    singular_fit_warnings = singular_count,
    stringsAsFactors = FALSE
  )
}

make_top100_rows <- function(length_cutoff, brain_area, method, results) {
  ordered_results <- results[order(results$diagnosis_p_value, na.last = TRUE), , drop = FALSE]
  top_results <- head(ordered_results, 100)

  data.frame(
    length_cutoff = length_cutoff,
    brain_area = brain_area,
    method = method,
    rank_by_raw_p = seq_len(nrow(top_results)),
    GeneID = top_results$GeneID,
    stringsAsFactors = FALSE
  )
}

make_effect_rows <- function(length_cutoff, brain_area, method, results) {
  data.frame(
    length_cutoff = length_cutoff,
    brain_area = brain_area,
    method = method,
    GeneID = results$GeneID,
    diagnosis_effect = results$diagnosis_effect,
    stringsAsFactors = FALSE
  )
}

## -------------------------------------------------------------------------
## 3. Method A: limma-voom for small groups
## -------------------------------------------------------------------------

run_limma_voom <- function(length_cutoff, brain_area, file_stub, cutoff_input_folder, cutoff_output_folder) {
  method <- "limma_voom"
  message("")
  message("Running ", method, " | ", brain_area, " | Length > ", length_cutoff)

  counts_file <- file.path(
    cutoff_input_folder,
    "03_filtered_raw_counts",
    "keep_all_samples",
    paste0("filtered_raw_counts_", file_stub, ".csv")
  )

  metadata_file <- file.path(
    cutoff_input_folder,
    "01_metadata",
    "keep_all_samples",
    paste0("metadata_", file_stub, ".csv")
  )

  raw_counts <- read_csv_matrix(counts_file)
  metadata <- read_metadata(metadata_file, c("SampleID", "Diagnosis", "Age", "Sex", "RIN"))

  message("Raw count matrix dimensions before matching:")
  print(dim(raw_counts))
  message("Metadata rows before matching: ", nrow(metadata))

  aligned <- align_matrix_and_metadata(raw_counts, metadata)
  raw_counts <- aligned$matrix
  metadata <- prepare_small_group_metadata(aligned$metadata)
  raw_counts <- raw_counts[, metadata$SampleID, drop = FALSE]

  message("Final raw count matrix dimensions:")
  print(dim(raw_counts))

  model_formula <- ~ Diagnosis + Age + Sex + RIN
  design_matrix <- model.matrix(model_formula, data = metadata)

  message("Design matrix columns:")
  print(colnames(design_matrix))

  if (!(diagnosis_coefficient %in% colnames(design_matrix))) {
    stop("Could not find expected design column: ", diagnosis_coefficient)
  }

  dge_object <- edgeR::DGEList(counts = raw_counts)
  dge_object <- edgeR::calcNormFactors(dge_object, method = "TMM")

  voom_object <- limma::voom(dge_object, design = design_matrix, plot = FALSE)
  fit <- limma::lmFit(voom_object, design_matrix)
  fit <- limma::eBayes(fit)

  limma_results <- limma::topTable(
    fit,
    coef = diagnosis_coefficient,
    number = Inf,
    sort.by = "P"
  )

  results <- data.frame(
    GeneID = rownames(limma_results),
    length_cutoff = length_cutoff,
    brain_area = brain_area,
    method = method,
    diagnosis_term = diagnosis_coefficient,
    diagnosis_effect = limma_results$logFC,
    standard_error = NA_real_,
    statistic = limma_results$t,
    diagnosis_p_value = limma_results$P.Value,
    diagnosis_FDR_BH = limma_results$adj.P.Val,
    average_expression = limma_results$AveExpr,
    convergence_warning = 0,
    singular_fit_warning = 0,
    stringsAsFactors = FALSE
  )

  output_folder <- file.path(cutoff_output_folder, method)
  file_prefix <- paste0(method, "_", brain_area)
  save_histogram_and_top_tables(
    results = results,
    output_folder = output_folder,
    file_prefix = file_prefix,
    plot_title = paste("Diagnosis p-values:", brain_area, method, "Length >", length_cutoff)
  )

  results
}

## -------------------------------------------------------------------------
## 4. Method A: dream for Frontal Cortex
## -------------------------------------------------------------------------

run_dream <- function(length_cutoff, brain_area, file_stub, cutoff_input_folder, cutoff_output_folder) {
  method <- "dream"
  message("")
  message("Running ", method, " | ", brain_area, " | Length > ", length_cutoff)

  counts_file <- file.path(
    cutoff_input_folder,
    "03_filtered_raw_counts",
    "keep_all_samples",
    paste0("filtered_raw_counts_", file_stub, ".csv")
  )

  metadata_file <- file.path(
    cutoff_input_folder,
    "01_metadata",
    "keep_all_samples",
    paste0("metadata_", file_stub, ".csv")
  )

  raw_counts <- read_csv_matrix(counts_file)
  metadata <- read_metadata(metadata_file, c("SampleID", "Diagnosis", "Age", "Sex", "RIN", "OriginalRegion", "BrainID"))

  message("Raw count matrix dimensions before matching:")
  print(dim(raw_counts))
  message("Metadata rows before matching: ", nrow(metadata))

  aligned <- align_matrix_and_metadata(raw_counts, metadata)
  raw_counts <- aligned$matrix
  metadata <- prepare_frontal_metadata(aligned$metadata)
  raw_counts <- raw_counts[, metadata$SampleID, drop = FALSE]
  rownames(metadata) <- metadata$SampleID

  message("Final raw count matrix dimensions:")
  print(dim(raw_counts))

  model_formula <- ~ Diagnosis + Age + Sex + RIN + OriginalRegion + (1 | BrainID)
  message("dream model formula:")
  print(model_formula)

  dge_object <- edgeR::DGEList(counts = raw_counts)
  dge_object <- edgeR::calcNormFactors(dge_object, method = "TMM")

  serial_parameter <- BiocParallel::SerialParam()

  voom_object <- variancePartition::voomWithDreamWeights(
    dge_object,
    formula = model_formula,
    data = metadata,
    plot = FALSE,
    BPPARAM = serial_parameter
  )

  dream_fit <- variancePartition::dream(
    voom_object,
    formula = model_formula,
    data = metadata,
    BPPARAM = serial_parameter
  )

  if (!(diagnosis_coefficient %in% colnames(dream_fit$coefficients))) {
    stop("Could not find expected dream coefficient: ", diagnosis_coefficient)
  }

  dream_results <- variancePartition::topTable(
    dream_fit,
    coef = diagnosis_coefficient,
    number = Inf
  )

  standard_error <- rep(NA_real_, nrow(dream_results))
  if ("SE" %in% colnames(dream_results)) {
    standard_error <- dream_results$SE
  }

  results <- data.frame(
    GeneID = rownames(dream_results),
    length_cutoff = length_cutoff,
    brain_area = brain_area,
    method = method,
    diagnosis_term = diagnosis_coefficient,
    diagnosis_effect = dream_results$logFC,
    standard_error = standard_error,
    statistic = dream_results$t,
    diagnosis_p_value = dream_results$P.Value,
    diagnosis_FDR_BH = dream_results$adj.P.Val,
    average_expression = dream_results$AveExpr,
    convergence_warning = 0,
    singular_fit_warning = 0,
    stringsAsFactors = FALSE
  )

  output_folder <- file.path(cutoff_output_folder, method)
  file_prefix <- paste0(method, "_", brain_area)
  save_histogram_and_top_tables(
    results = results,
    output_folder = output_folder,
    file_prefix = file_prefix,
    plot_title = paste("Diagnosis p-values:", brain_area, method, "Length >", length_cutoff)
  )

  results
}

## -------------------------------------------------------------------------
## 5. Method B: ordinary LM on log2(GeTMM + 1)
## -------------------------------------------------------------------------

run_lm_getmm_log2 <- function(length_cutoff, brain_area, file_stub, cutoff_input_folder, cutoff_output_folder, frontal_fixed_effect) {
  if (frontal_fixed_effect) {
    method <- "lm_getmm_log2_frontal_fixed_effect"
  } else {
    method <- "lm_getmm_log2"
  }

  message("")
  message("Running ", method, " | ", brain_area, " | Length > ", length_cutoff)

  expression_file <- file.path(
    cutoff_input_folder,
    "02_normalized_expression",
    "keep_all_samples",
    paste0("normalized_GeTMM_log2_", file_stub, ".csv")
  )

  metadata_file <- file.path(
    cutoff_input_folder,
    "01_metadata",
    "keep_all_samples",
    paste0("metadata_", file_stub, ".csv")
  )

  expression_matrix <- read_csv_matrix(expression_file)

  if (frontal_fixed_effect) {
    metadata <- read_metadata(metadata_file, c("SampleID", "Diagnosis", "Age", "Sex", "RIN", "OriginalRegion", "BrainID"))
  } else {
    metadata <- read_metadata(metadata_file, c("SampleID", "Diagnosis", "Age", "Sex", "RIN"))
  }

  message("Expression matrix dimensions before matching:")
  print(dim(expression_matrix))
  message("Metadata rows before matching: ", nrow(metadata))

  aligned <- align_matrix_and_metadata(expression_matrix, metadata)
  expression_matrix <- aligned$matrix

  if (frontal_fixed_effect) {
    metadata <- prepare_frontal_metadata(aligned$metadata)
  } else {
    metadata <- prepare_small_group_metadata(aligned$metadata)
  }

  expression_matrix <- expression_matrix[, metadata$SampleID, drop = FALSE]

  message("Final expression matrix dimensions:")
  print(dim(expression_matrix))

  if (frontal_fixed_effect) {
    model_formula <- expression_value ~ Diagnosis + Age + Sex + RIN + OriginalRegion
  } else {
    model_formula <- expression_value ~ Diagnosis + Age + Sex + RIN
  }

  message("LM model formula:")
  print(model_formula)

  result_rows <- vector("list", nrow(expression_matrix))
  gene_ids <- rownames(expression_matrix)

  for (gene_index in seq_len(nrow(expression_matrix))) {
    gene_id <- gene_ids[gene_index]

    model_data <- metadata
    model_data$expression_value <- as.numeric(expression_matrix[gene_index, ])

    fit <- tryCatch(
      lm(model_formula, data = model_data),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      result_rows[[gene_index]] <- data.frame(
        GeneID = gene_id,
        diagnosis_effect = NA_real_,
        standard_error = NA_real_,
        statistic = NA_real_,
        diagnosis_p_value = NA_real_,
        fit_error = conditionMessage(fit),
        stringsAsFactors = FALSE
      )
      next
    }

    coefficient_table <- summary(fit)$coefficients

    if (!(diagnosis_coefficient %in% rownames(coefficient_table))) {
      result_rows[[gene_index]] <- data.frame(
        GeneID = gene_id,
        diagnosis_effect = NA_real_,
        standard_error = NA_real_,
        statistic = NA_real_,
        diagnosis_p_value = NA_real_,
        fit_error = paste("Missing coefficient:", diagnosis_coefficient),
        stringsAsFactors = FALSE
      )
      next
    }

    result_rows[[gene_index]] <- data.frame(
      GeneID = gene_id,
      diagnosis_effect = coefficient_table[diagnosis_coefficient, "Estimate"],
      standard_error = coefficient_table[diagnosis_coefficient, "Std. Error"],
      statistic = coefficient_table[diagnosis_coefficient, "t value"],
      diagnosis_p_value = coefficient_table[diagnosis_coefficient, "Pr(>|t|)"],
      fit_error = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  results <- do.call(rbind, result_rows)
  results$diagnosis_FDR_BH <- p.adjust(results$diagnosis_p_value, method = "BH")

  results <- data.frame(
    GeneID = results$GeneID,
    length_cutoff = length_cutoff,
    brain_area = brain_area,
    method = method,
    diagnosis_term = diagnosis_coefficient,
    diagnosis_effect = results$diagnosis_effect,
    standard_error = results$standard_error,
    statistic = results$statistic,
    diagnosis_p_value = results$diagnosis_p_value,
    diagnosis_FDR_BH = results$diagnosis_FDR_BH,
    average_expression = rowMeans(expression_matrix[results$GeneID, , drop = FALSE], na.rm = TRUE),
    convergence_warning = 0,
    singular_fit_warning = 0,
    fit_error = results$fit_error,
    stringsAsFactors = FALSE
  )

  output_folder <- file.path(cutoff_output_folder, method)
  file_prefix <- paste0(method, "_", brain_area)
  save_histogram_and_top_tables(
    results = results,
    output_folder = output_folder,
    file_prefix = file_prefix,
    plot_title = paste("Diagnosis p-values:", brain_area, method, "Length >", length_cutoff)
  )

  results
}

## -------------------------------------------------------------------------
## 6. Method B: LMM on log2(GeTMM + 1) for Frontal Cortex
## -------------------------------------------------------------------------

run_lmm_getmm_log2 <- function(length_cutoff, brain_area, file_stub, cutoff_input_folder, cutoff_output_folder) {
  method <- "lmm_getmm_log2"
  message("")
  message("Running ", method, " | ", brain_area, " | Length > ", length_cutoff)

  expression_file <- file.path(
    cutoff_input_folder,
    "02_normalized_expression",
    "keep_all_samples",
    paste0("normalized_GeTMM_log2_", file_stub, ".csv")
  )

  metadata_file <- file.path(
    cutoff_input_folder,
    "01_metadata",
    "keep_all_samples",
    paste0("metadata_", file_stub, ".csv")
  )

  expression_matrix <- read_csv_matrix(expression_file)
  metadata <- read_metadata(metadata_file, c("SampleID", "Diagnosis", "Age", "Sex", "RIN", "OriginalRegion", "BrainID"))

  message("Expression matrix dimensions before matching:")
  print(dim(expression_matrix))
  message("Metadata rows before matching: ", nrow(metadata))

  aligned <- align_matrix_and_metadata(expression_matrix, metadata)
  expression_matrix <- aligned$matrix
  metadata <- prepare_frontal_metadata(aligned$metadata)
  expression_matrix <- expression_matrix[, metadata$SampleID, drop = FALSE]

  message("Final expression matrix dimensions:")
  print(dim(expression_matrix))

  model_formula <- expression_value ~ Diagnosis + Age + Sex + RIN + OriginalRegion + (1 | BrainID)

  message("LMM model formula:")
  print(model_formula)

  result_rows <- vector("list", nrow(expression_matrix))
  gene_ids <- rownames(expression_matrix)

  for (gene_index in seq_len(nrow(expression_matrix))) {
    gene_id <- gene_ids[gene_index]

    model_data <- metadata
    model_data$expression_value <- as.numeric(expression_matrix[gene_index, ])

    warning_text <- character()

    fit <- tryCatch(
      withCallingHandlers(
        lmerTest::lmer(model_formula, data = model_data, REML = FALSE),
        warning = function(w) {
          warning_text <<- c(warning_text, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      result_rows[[gene_index]] <- data.frame(
        GeneID = gene_id,
        diagnosis_effect = NA_real_,
        standard_error = NA_real_,
        statistic = NA_real_,
        diagnosis_p_value = NA_real_,
        convergence_warning = length(warning_text) > 0,
        singular_fit_warning = NA,
        warning_text = paste(warning_text, collapse = " | "),
        fit_error = conditionMessage(fit),
        stringsAsFactors = FALSE
      )
      next
    }

    singular_fit <- lme4::isSingular(fit, tol = 1e-4)
    coefficient_table <- summary(fit)$coefficients
    p_value_column <- grep("^Pr\\(", colnames(coefficient_table), value = TRUE)

    if (length(p_value_column) == 0) {
      p_value_column <- NA_character_
    } else {
      p_value_column <- p_value_column[1]
    }

    if (!(diagnosis_coefficient %in% rownames(coefficient_table)) || is.na(p_value_column)) {
      result_rows[[gene_index]] <- data.frame(
        GeneID = gene_id,
        diagnosis_effect = NA_real_,
        standard_error = NA_real_,
        statistic = NA_real_,
        diagnosis_p_value = NA_real_,
        convergence_warning = length(warning_text) > 0,
        singular_fit_warning = singular_fit,
        warning_text = paste(warning_text, collapse = " | "),
        fit_error = paste("Missing coefficient or p-value:", diagnosis_coefficient),
        stringsAsFactors = FALSE
      )
      next
    }

    result_rows[[gene_index]] <- data.frame(
      GeneID = gene_id,
      diagnosis_effect = coefficient_table[diagnosis_coefficient, "Estimate"],
      standard_error = coefficient_table[diagnosis_coefficient, "Std. Error"],
      statistic = coefficient_table[diagnosis_coefficient, "t value"],
      diagnosis_p_value = coefficient_table[diagnosis_coefficient, p_value_column],
      convergence_warning = length(warning_text) > 0,
      singular_fit_warning = singular_fit,
      warning_text = paste(warning_text, collapse = " | "),
      fit_error = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  results <- do.call(rbind, result_rows)
  results$diagnosis_FDR_BH <- p.adjust(results$diagnosis_p_value, method = "BH")

  results <- data.frame(
    GeneID = results$GeneID,
    length_cutoff = length_cutoff,
    brain_area = brain_area,
    method = method,
    diagnosis_term = diagnosis_coefficient,
    diagnosis_effect = results$diagnosis_effect,
    standard_error = results$standard_error,
    statistic = results$statistic,
    diagnosis_p_value = results$diagnosis_p_value,
    diagnosis_FDR_BH = results$diagnosis_FDR_BH,
    average_expression = rowMeans(expression_matrix[results$GeneID, , drop = FALSE], na.rm = TRUE),
    convergence_warning = results$convergence_warning,
    singular_fit_warning = results$singular_fit_warning,
    warning_text = results$warning_text,
    fit_error = results$fit_error,
    stringsAsFactors = FALSE
  )

  output_folder <- file.path(cutoff_output_folder, method)
  file_prefix <- paste0(method, "_", brain_area)
  save_histogram_and_top_tables(
    results = results,
    output_folder = output_folder,
    file_prefix = file_prefix,
    plot_title = paste("Diagnosis p-values:", brain_area, method, "Length >", length_cutoff)
  )

  results
}

## -------------------------------------------------------------------------
## 7. Run every cutoff, brain area, and method
## -------------------------------------------------------------------------

all_summary_rows <- list()
all_top100_rows <- list()
all_effect_rows <- list()

for (length_cutoff in length_cutoffs) {
  length_label <- make_length_label(length_cutoff)

  message("")
  message("============================================================")
  message("Starting Step 3 methods for ", length_label)
  message("============================================================")

  cutoff_input_folder <- file.path(step1_sensitivity_folder, length_label)
  cutoff_output_folder <- file.path(output_root, length_label)

  if (!dir.exists(cutoff_input_folder)) {
    stop("Input folder does not exist. Run Step1_length_cutoff_sensitivity_GeTMM.R first: ",
         cutoff_input_folder)
  }

  dir.create(cutoff_output_folder, recursive = TRUE, showWarnings = FALSE)

  for (group_index in seq_len(nrow(brain_areas))) {
    brain_area <- brain_areas$brain_area[group_index]
    file_stub <- brain_areas$file_stub[group_index]
    is_frontal <- brain_areas$is_frontal[group_index]

    if (is_frontal) {
      dream_results <- run_dream(
        length_cutoff = length_cutoff,
        brain_area = brain_area,
        file_stub = file_stub,
        cutoff_input_folder = cutoff_input_folder,
        cutoff_output_folder = cutoff_output_folder
      )

      lmm_results <- run_lmm_getmm_log2(
        length_cutoff = length_cutoff,
        brain_area = brain_area,
        file_stub = file_stub,
        cutoff_input_folder = cutoff_input_folder,
        cutoff_output_folder = cutoff_output_folder
      )

      frontal_fixed_lm_results <- run_lm_getmm_log2(
        length_cutoff = length_cutoff,
        brain_area = brain_area,
        file_stub = file_stub,
        cutoff_input_folder = cutoff_input_folder,
        cutoff_output_folder = cutoff_output_folder,
        frontal_fixed_effect = TRUE
      )

      method_results <- list(
        dream = dream_results,
        lmm_getmm_log2 = lmm_results,
        lm_getmm_log2_frontal_fixed_effect = frontal_fixed_lm_results
      )
    } else {
      limma_voom_results <- run_limma_voom(
        length_cutoff = length_cutoff,
        brain_area = brain_area,
        file_stub = file_stub,
        cutoff_input_folder = cutoff_input_folder,
        cutoff_output_folder = cutoff_output_folder
      )

      lm_results <- run_lm_getmm_log2(
        length_cutoff = length_cutoff,
        brain_area = brain_area,
        file_stub = file_stub,
        cutoff_input_folder = cutoff_input_folder,
        cutoff_output_folder = cutoff_output_folder,
        frontal_fixed_effect = FALSE
      )

      method_results <- list(
        limma_voom = limma_voom_results,
        lm_getmm_log2 = lm_results
      )
    }

    for (method_name in names(method_results)) {
      results <- method_results[[method_name]]

      all_summary_rows[[length(all_summary_rows) + 1]] <- make_summary_row(
        length_cutoff = length_cutoff,
        brain_area = brain_area,
        method = method_name,
        results = results
      )

      all_top100_rows[[length(all_top100_rows) + 1]] <- make_top100_rows(
        length_cutoff = length_cutoff,
        brain_area = brain_area,
        method = method_name,
        results = results
      )

      all_effect_rows[[length(all_effect_rows) + 1]] <- make_effect_rows(
        length_cutoff = length_cutoff,
        brain_area = brain_area,
        method = method_name,
        results = results
      )
    }
  }
}

## -------------------------------------------------------------------------
## 8. Save the final all-method summary table
## -------------------------------------------------------------------------

summary_table <- do.call(rbind, all_summary_rows)

summary_file <- file.path(output_root, "Step3_length_cutoff_all_methods_summary.csv")
write.csv(summary_table, summary_file, row.names = FALSE)

message("")
message("Saved final Step 3 summary table:")
message(summary_file)

top100_table <- do.call(rbind, all_top100_rows)
effect_table <- do.call(rbind, all_effect_rows)

## -------------------------------------------------------------------------
## 9. Top gene stability versus Length > 20
## -------------------------------------------------------------------------

stability_rows <- list()

top100_groups <- unique(top100_table[, c("brain_area", "method", "length_cutoff"), drop = FALSE])

for (row_index in seq_len(nrow(top100_groups))) {
  brain_area <- top100_groups$brain_area[row_index]
  method <- top100_groups$method[row_index]
  length_cutoff <- top100_groups$length_cutoff[row_index]

  baseline_top100 <- top100_table[
    top100_table$brain_area == brain_area &
      top100_table$method == method &
      top100_table$length_cutoff == 20,
    "GeneID"
  ]

  current_top100 <- top100_table[
    top100_table$brain_area == brain_area &
      top100_table$method == method &
      top100_table$length_cutoff == length_cutoff,
    "GeneID"
  ]

  if (length(baseline_top100) == 0) {
    overlap_count <- NA_integer_
    percent_overlap <- NA_real_
  } else {
    overlap_count <- length(intersect(baseline_top100, current_top100))
    percent_overlap <- 100 * overlap_count / length(baseline_top100)
  }

  stability_rows[[length(stability_rows) + 1]] <- data.frame(
    brain_area = brain_area,
    method = method,
    length_cutoff = length_cutoff,
    overlap_with_baseline_top100 = overlap_count,
    percent_overlap_with_baseline_top100 = percent_overlap,
    stringsAsFactors = FALSE
  )
}

stability_table <- do.call(rbind, stability_rows)

stability_file <- file.path(output_root, "Step3_top_gene_stability_vs_Length_gt_20_all_methods.csv")
write.csv(stability_table, stability_file, row.names = FALSE)

message("Saved top-gene stability table:")
message(stability_file)

## -------------------------------------------------------------------------
## 10. Cross-method top 100 overlap
## -------------------------------------------------------------------------

cross_method_rows <- list()

length_area_pairs <- unique(top100_table[, c("length_cutoff", "brain_area"), drop = FALSE])

for (pair_index in seq_len(nrow(length_area_pairs))) {
  length_cutoff <- length_area_pairs$length_cutoff[pair_index]
  brain_area <- length_area_pairs$brain_area[pair_index]

  available_methods <- sort(unique(top100_table$method[
    top100_table$length_cutoff == length_cutoff &
      top100_table$brain_area == brain_area
  ]))

  if (length(available_methods) < 2) {
    next
  }

  method_pairs <- combn(available_methods, 2, simplify = FALSE)

  for (method_pair in method_pairs) {
    method_a <- method_pair[1]
    method_b <- method_pair[2]

    top_a <- top100_table$GeneID[
      top100_table$length_cutoff == length_cutoff &
        top100_table$brain_area == brain_area &
        top100_table$method == method_a
    ]

    top_b <- top100_table$GeneID[
      top100_table$length_cutoff == length_cutoff &
        top100_table$brain_area == brain_area &
        top100_table$method == method_b
    ]

    overlap_count <- length(intersect(top_a, top_b))
    percent_overlap <- 100 * overlap_count / min(length(top_a), length(top_b))

    cross_method_rows[[length(cross_method_rows) + 1]] <- data.frame(
      length_cutoff = length_cutoff,
      brain_area = brain_area,
      method_A = method_a,
      method_B = method_b,
      top100_overlap = overlap_count,
      percent_overlap = percent_overlap,
      stringsAsFactors = FALSE
    )
  }
}

cross_method_table <- do.call(rbind, cross_method_rows)

cross_method_file <- file.path(output_root, "Step3_cross_method_top100_overlap.csv")
write.csv(cross_method_table, cross_method_file, row.names = FALSE)

message("Saved cross-method top-100 overlap table:")
message(cross_method_file)

## -------------------------------------------------------------------------
## 11. Cross-method diagnosis-effect correlation
## -------------------------------------------------------------------------

effect_comparison_pairs <- data.frame(
  brain_area = c(
    "Insula",
    "Ventral_Striatum",
    "Hippocampal_Formation",
    "Frontal_Cortex",
    "Frontal_Cortex"
  ),
  method_A = c(
    "limma_voom",
    "limma_voom",
    "limma_voom",
    "dream",
    "lmm_getmm_log2"
  ),
  method_B = c(
    "lm_getmm_log2",
    "lm_getmm_log2",
    "lm_getmm_log2",
    "lmm_getmm_log2",
    "lm_getmm_log2_frontal_fixed_effect"
  ),
  stringsAsFactors = FALSE
)

effect_correlation_rows <- list()

for (length_cutoff in length_cutoffs) {
  for (pair_index in seq_len(nrow(effect_comparison_pairs))) {
    brain_area <- effect_comparison_pairs$brain_area[pair_index]
    method_a <- effect_comparison_pairs$method_A[pair_index]
    method_b <- effect_comparison_pairs$method_B[pair_index]

    effects_a <- effect_table[
      effect_table$length_cutoff == length_cutoff &
        effect_table$brain_area == brain_area &
        effect_table$method == method_a,
      c("GeneID", "diagnosis_effect"),
      drop = FALSE
    ]

    effects_b <- effect_table[
      effect_table$length_cutoff == length_cutoff &
        effect_table$brain_area == brain_area &
        effect_table$method == method_b,
      c("GeneID", "diagnosis_effect"),
      drop = FALSE
    ]

    merged_effects <- merge(effects_a, effects_b, by = "GeneID", suffixes = c("_A", "_B"))
    complete_effects <- complete.cases(merged_effects[, c("diagnosis_effect_A", "diagnosis_effect_B")])
    merged_effects <- merged_effects[complete_effects, , drop = FALSE]

    shared_genes <- nrow(merged_effects)

    if (shared_genes >= 2) {
      effect_correlation <- cor(
        merged_effects$diagnosis_effect_A,
        merged_effects$diagnosis_effect_B,
        use = "complete.obs"
      )
    } else {
      effect_correlation <- NA_real_
    }

    nonzero_effects <- merged_effects$diagnosis_effect_A != 0 &
      merged_effects$diagnosis_effect_B != 0

    if (sum(nonzero_effects) > 0) {
      same_direction_percent <- 100 * mean(
        sign(merged_effects$diagnosis_effect_A[nonzero_effects]) ==
          sign(merged_effects$diagnosis_effect_B[nonzero_effects])
      )
    } else {
      same_direction_percent <- NA_real_
    }

    effect_correlation_rows[[length(effect_correlation_rows) + 1]] <- data.frame(
      length_cutoff = length_cutoff,
      brain_area = brain_area,
      method_A = method_a,
      method_B = method_b,
      shared_genes = shared_genes,
      diagnosis_effect_correlation = effect_correlation,
      same_direction_percent = same_direction_percent,
      stringsAsFactors = FALSE
    )
  }
}

effect_correlation_table <- do.call(rbind, effect_correlation_rows)

effect_correlation_file <- file.path(output_root, "Step3_cross_method_effect_size_correlation.csv")
write.csv(effect_correlation_table, effect_correlation_file, row.names = FALSE)

message("Saved cross-method effect-size correlation table:")
message(effect_correlation_file)

## -------------------------------------------------------------------------
## 12. Finish
## -------------------------------------------------------------------------

message("")
message("Finished Step 3 length-cutoff sensitivity DGE analysis.")
message("Main output folder:")
message(normalizePath(output_root, mustWork = FALSE))
message("")
message("Important interpretation rule:")
message("Do not choose a cutoff automatically based only on the number of FDR-significant genes.")
message("Look for fewer tested genes, improved FDR, stable top genes, consistent effect directions,")
message("and no suspicious collapse of raw p-value signal.")
