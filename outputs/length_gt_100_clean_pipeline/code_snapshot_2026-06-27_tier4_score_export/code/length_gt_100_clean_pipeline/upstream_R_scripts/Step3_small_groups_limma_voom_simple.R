## Step3_small_groups_limma_voom_simple.R
##
## Purpose:
##   Differential expression for the smaller analysis groups using limma-voom.
##
## Groups:
##   - Insula
##   - Ventral_Striatum
##   - Hippocampal_Formation
##
## Important assumptions:
##   1. PCA/QC has already been done.
##   2. Step 1 preprocessing has already been done.
##   3. This script does not run PCA.
##   4. This script does not recreate the Step 1 GeTMM-log2 files.
##   5. limma-voom starts from filtered raw counts, not log2 GeTMM values.
##   6. Count matrix rows are genes and count matrix columns are samples.
##   7. Metadata has one row per sample and contains SampleID.

## -------------------------------------------------------------------------
## 0. Packages
## -------------------------------------------------------------------------

required_packages <- c("edgeR", "limma")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

library(edgeR)
library(limma)

## -------------------------------------------------------------------------
## 1. File locations
## -------------------------------------------------------------------------

step1_folder <- "D:/genes/genedata/Step1_GSE102556_GeTMM_brain_area_groups"
counts_folder <- file.path(step1_folder, "03_filtered_raw_counts", "keep_all_samples")
metadata_folder <- file.path(step1_folder, "01_metadata", "keep_all_samples")

output_folder <- file.path("outputs", "step3_limma_voom_dream")
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(output_folder)) {
  stop("Could not create output folder: ", output_folder)
}

analysis_groups <- data.frame(
  group_label = c("Insula", "Ventral_Striatum", "Hippocampal_Formation"),
  file_stub = c("insula", "ventral_striatum", "hippocampal_formation"),
  stringsAsFactors = FALSE
)

## This is the model used for each small group.
## Diagnosis is the term of interest.
## Age, Sex, and RIN are covariates chosen after QC/PCA.
model_formula <- ~ Diagnosis + Age + Sex + RIN

message("Model formula for small groups:")
print(model_formula)

## -------------------------------------------------------------------------
## 2. Run one explicit limma-voom analysis per group
## -------------------------------------------------------------------------

for (group_index in seq_len(nrow(analysis_groups))) {
  group_label <- analysis_groups$group_label[group_index]
  file_stub <- analysis_groups$file_stub[group_index]

  message("")
  message("============================================================")
  message("Starting limma-voom analysis for: ", group_label)
  message("============================================================")

  counts_file <- file.path(counts_folder, paste0("filtered_raw_counts_", file_stub, ".csv"))
  metadata_file <- file.path(metadata_folder, paste0("metadata_", file_stub, ".csv"))
  results_file <- file.path(output_folder, paste0("limma_voom_", group_label, "_results.csv"))
  histogram_file <- file.path(output_folder, paste0("limma_voom_", group_label, "_diagnosis_pvalue_histogram.png"))

  message("Counts file:    ", counts_file)
  message("Metadata file:  ", metadata_file)
  message("Results file:   ", results_file)
  message("Histogram file: ", histogram_file)

  if (!file.exists(counts_file)) {
    stop("Counts file was not found: ", counts_file)
  }

  if (!file.exists(metadata_file)) {
    stop("Metadata file was not found: ", metadata_file)
  }

  ## -----------------------------------------------------------------------
  ## 3. Load filtered raw counts and metadata
  ## -----------------------------------------------------------------------

  ## The first column contains gene IDs.
  ## Rows are genes. Columns are samples.
  raw_counts <- read.csv(
    counts_file,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  metadata <- read.csv(
    metadata_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  message("Raw count matrix rows (genes):   ", nrow(raw_counts))
  message("Raw count matrix columns(samples): ", ncol(raw_counts))
  message("Metadata rows:                   ", nrow(metadata))
  message("Metadata columns found:")
  print(names(metadata))

  message("First count sample IDs before matching:")
  print(head(colnames(raw_counts)))

  if ("SampleID" %in% names(metadata)) {
    message("First metadata SampleID values before matching:")
    print(head(metadata$SampleID))
  } else {
    stop("Metadata must contain a column named SampleID.")
  }

  ## -----------------------------------------------------------------------
  ## 4. Check required columns
  ## -----------------------------------------------------------------------

  required_columns <- c("SampleID", "Diagnosis", "Age", "Sex", "RIN")
  missing_columns <- setdiff(required_columns, names(metadata))

  if (length(missing_columns) > 0) {
    stop("Missing required metadata columns: ", paste(missing_columns, collapse = ", "))
  }

  if (any(duplicated(metadata$SampleID))) {
    stop("Metadata contains duplicated SampleID values.")
  }

  if (any(duplicated(colnames(raw_counts)))) {
    stop("Count matrix contains duplicated sample columns.")
  }

  ## -----------------------------------------------------------------------
  ## 5. Align count columns and metadata rows
  ## -----------------------------------------------------------------------

  common_samples <- intersect(colnames(raw_counts), metadata$SampleID)

  message("Samples found in both counts and metadata: ", length(common_samples))

  if (length(common_samples) == 0) {
    stop("No sample IDs matched between counts and metadata.")
  }

  raw_counts <- raw_counts[, common_samples, drop = FALSE]
  metadata <- metadata[match(common_samples, metadata$SampleID), , drop = FALSE]

  message("First count sample IDs after matching:")
  print(head(colnames(raw_counts)))

  message("First metadata SampleID values after matching:")
  print(head(metadata$SampleID))

  if (!all(colnames(raw_counts) == metadata$SampleID)) {
    stop("Counts and metadata are not aligned after matching.")
  }

  ## -----------------------------------------------------------------------
  ## 6. Prepare metadata for the design matrix
  ## -----------------------------------------------------------------------

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

  message("Sex counts:")
  print(table(metadata$Sex, useNA = "ifany"))

  message("Missing values in model covariates:")
  print(colSums(is.na(metadata[, c("Diagnosis", "Age", "Sex", "RIN")])))

  complete_rows <- complete.cases(metadata[, c("Diagnosis", "Age", "Sex", "RIN")])

  message("Samples before removing missing covariates: ", nrow(metadata))
  message("Samples after removing missing covariates:  ", sum(complete_rows))

  metadata <- metadata[complete_rows, , drop = FALSE]
  raw_counts <- raw_counts[, metadata$SampleID, drop = FALSE]

  if (nrow(metadata) < 3) {
    stop("Too few samples remain after removing missing covariates.")
  }

  ## -----------------------------------------------------------------------
  ## 7. Convert counts to a numeric matrix
  ## -----------------------------------------------------------------------

  raw_counts_matrix <- as.matrix(raw_counts)
  storage.mode(raw_counts_matrix) <- "numeric"

  if (any(is.na(raw_counts_matrix))) {
    stop("The count matrix contains NA values after numeric conversion.")
  }

  if (any(raw_counts_matrix < 0)) {
    stop("The count matrix contains negative values, which are invalid for raw counts.")
  }

  message("Final count matrix dimensions:")
  print(dim(raw_counts_matrix))

  ## -----------------------------------------------------------------------
  ## 8. Build the design matrix explicitly
  ## -----------------------------------------------------------------------

  design_matrix <- model.matrix(model_formula, data = metadata)

  message("Design matrix columns:")
  print(colnames(design_matrix))

  message("First rows of design matrix:")
  print(head(design_matrix))

  diagnosis_coefficient <- "DiagnosisMDD"

  if (!(diagnosis_coefficient %in% colnames(design_matrix))) {
    stop("Could not find expected design column: ", diagnosis_coefficient)
  }

  ## -----------------------------------------------------------------------
  ## 9. Run edgeR library-size normalization for voom
  ## -----------------------------------------------------------------------

  ## This does not recreate the Step 1 GeTMM-log2 output.
  ## This is the standard count-based scaling step used inside limma-voom.
  dge_object <- edgeR::DGEList(counts = raw_counts_matrix)
  dge_object <- edgeR::calcNormFactors(dge_object, method = "TMM")

  message("Library sizes and TMM factors:")
  print(head(dge_object$samples))

  ## -----------------------------------------------------------------------
  ## 10. Run voom and limma
  ## -----------------------------------------------------------------------

  voom_object <- limma::voom(
    dge_object,
    design = design_matrix,
    plot = FALSE
  )

  fit <- limma::lmFit(voom_object, design_matrix)
  fit <- limma::eBayes(fit)

  message("Model coefficients available:")
  print(colnames(fit$coefficients))

  ## -----------------------------------------------------------------------
  ## 11. Extract Diagnosis results
  ## -----------------------------------------------------------------------

  limma_results <- limma::topTable(
    fit,
    coef = diagnosis_coefficient,
    number = Inf,
    sort.by = "P"
  )

  final_results <- data.frame(
    GeneID = rownames(limma_results),
    diagnosis_term = diagnosis_coefficient,
    diagnosis_logFC = limma_results$logFC,
    average_expression = limma_results$AveExpr,
    moderated_t = limma_results$t,
    diagnosis_p_value = limma_results$P.Value,
    diagnosis_FDR_BH = limma_results$adj.P.Val,
    log_odds_DE = limma_results$B,
    stringsAsFactors = FALSE
  )

  message("Number of genes in result table: ", nrow(final_results))
  message("Genes with raw p < 0.05: ", sum(final_results$diagnosis_p_value < 0.05, na.rm = TRUE))
  message("Genes with FDR < 0.05:   ", sum(final_results$diagnosis_FDR_BH < 0.05, na.rm = TRUE))
  message("Genes with FDR < 0.10:   ", sum(final_results$diagnosis_FDR_BH < 0.10, na.rm = TRUE))

  ## -----------------------------------------------------------------------
  ## 12. Save CSV results
  ## -----------------------------------------------------------------------

  write.csv(final_results, results_file, row.names = FALSE)

  message("Saved results to:")
  message(results_file)

  ## -----------------------------------------------------------------------
  ## 13. Save a simple p-value histogram
  ## -----------------------------------------------------------------------

  png(filename = histogram_file, width = 900, height = 700)
  hist(
    final_results$diagnosis_p_value,
    breaks = 50,
    main = paste("Diagnosis p-values:", group_label),
    xlab = "Raw p-value for DiagnosisMDD",
    col = "gray80",
    border = "white"
  )
  dev.off()

  message("Saved p-value histogram to:")
  message(histogram_file)
}

message("")
message("Finished all small-group limma-voom analyses.")
