## LM_Small_Groups_simple.R
##
## Purpose:
##   Run one simple linear model per gene for the smaller analysis groups:
##   Insula, Ventral_Striatum, and Hippocampal_Formation.
##
## Important assumptions:
##   1. Each expression CSV has genes in rows and samples in columns.
##   2. The first column of each expression CSV contains gene IDs.
##   3. Each metadata CSV has one row per sample.
##   4. Each metadata CSV contains a SampleID column that matches expression column names.
##   5. PCA/QC has already been done. This script does not run PCA or make plots.

## -------------------------------------------------------------------------
## 1. File locations and group names
## -------------------------------------------------------------------------

## These folders point to the Step 1 outputs from the GeTMM-log2 preprocessing.
## These three groups do not appear to have the same within-group repeated donor
## issue seen in Frontal_Cortex, so this script uses simple linear models.
step1_folder <- "D:/genes/genedata/Step1_GSE102556_GeTMM_brain_area_groups"
expression_folder <- file.path(step1_folder, "02_normalized_expression", "keep_all_samples")
metadata_folder <- file.path(step1_folder, "01_metadata", "keep_all_samples")

## This script saves one differential-expression results table per group.
output_folder <- file.path("outputs", "differential_expression")
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

analysis_groups <- data.frame(
  group_label = c("Insula", "Ventral_Striatum", "Hippocampal_Formation"),
  file_stub = c("insula", "ventral_striatum", "hippocampal_formation"),
  stringsAsFactors = FALSE
)

## This is the model used for every gene in these groups.
## You can manually edit this formula later if your final covariate choice changes.
model_formula <- expression_value ~ Diagnosis + Age + Sex + RIN

message("Model formula:")
print(model_formula)

## -------------------------------------------------------------------------
## 2. Analyze each group
## -------------------------------------------------------------------------

for (group_index in seq_len(nrow(analysis_groups))) {
  group_name <- analysis_groups$group_label[group_index]
  file_stub <- analysis_groups$file_stub[group_index]

  message("")
  message("============================================================")
  message("Starting group: ", group_name)
  message("============================================================")

  expression_file <- file.path(expression_folder, paste0("normalized_GeTMM_log2_", file_stub, ".csv"))
  metadata_file <- file.path(metadata_folder, paste0("metadata_", file_stub, ".csv"))
  output_file <- file.path(output_folder, paste0("LM_", group_name, "_results.csv"))

  message("Expression file: ", expression_file)
  message("Metadata file:   ", metadata_file)
  message("Output file:     ", output_file)

  if (!file.exists(expression_file)) {
    stop("Expression file was not found: ", expression_file)
  }

  if (!file.exists(metadata_file)) {
    stop("Metadata file was not found: ", metadata_file)
  }

  ## -----------------------------------------------------------------------
  ## 3. Read expression and metadata
  ## -----------------------------------------------------------------------

  ## Expression is expected to be genes x samples.
  ## The first column is used as row names, so rownames(expression_matrix)
  ## should be gene IDs after reading the file.
  expression_matrix <- read.csv(
    expression_file,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  metadata <- read.csv(
    metadata_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  message("Expression rows (genes):     ", nrow(expression_matrix))
  message("Expression columns(samples): ", ncol(expression_matrix))
  message("Metadata rows:               ", nrow(metadata))
  message("Metadata columns found:")
  print(names(metadata))

  message("First expression sample IDs before matching:")
  print(head(colnames(expression_matrix)))

  if ("SampleID" %in% names(metadata)) {
    message("First metadata SampleID values before matching:")
    print(head(metadata$SampleID))
  } else {
    stop(metadata_file, " must contain a column named SampleID.")
  }

  ## -----------------------------------------------------------------------
  ## 4. Check required metadata columns
  ## -----------------------------------------------------------------------

  required_columns <- c("SampleID", "Diagnosis", "Age", "Sex", "RIN")

  missing_columns <- setdiff(required_columns, names(metadata))

  if (length(missing_columns) > 0) {
    stop("Missing required metadata columns in ", group_name, ": ", paste(missing_columns, collapse = ", "))
  }

  if (any(duplicated(metadata$SampleID))) {
    stop("Some SampleID values are duplicated in metadata for ", group_name, ". Please fix this before modeling.")
  }

  if (any(duplicated(colnames(expression_matrix)))) {
    stop("Some sample column names are duplicated in the expression matrix for ", group_name, ". Please fix this before modeling.")
  }

  ## -----------------------------------------------------------------------
  ## 5. Align metadata rows to expression columns
  ## -----------------------------------------------------------------------

  ## We keep only samples found in both files.
  ## Then we reorder metadata so its rows are in exactly the same order
  ## as the expression columns.
  common_samples <- intersect(colnames(expression_matrix), metadata$SampleID)

  message("Number of samples found in both expression and metadata: ", length(common_samples))

  if (length(common_samples) == 0) {
    stop("No matching sample IDs were found between expression columns and metadata$SampleID for ", group_name, ".")
  }

  expression_matrix <- expression_matrix[, common_samples, drop = FALSE]
  metadata <- metadata[match(common_samples, metadata$SampleID), , drop = FALSE]

  message("First expression sample IDs after matching:")
  print(head(colnames(expression_matrix)))

  message("First metadata SampleID values after matching:")
  print(head(metadata$SampleID))

  if (!all(colnames(expression_matrix) == metadata$SampleID)) {
    stop("Expression columns and metadata SampleID values are not aligned for ", group_name, ".")
  }

  ## -----------------------------------------------------------------------
  ## 6. Prepare covariates
  ## -----------------------------------------------------------------------

  ## Categorical variables are converted to factors.
  ## Numeric variables are converted to numeric.
  metadata$Diagnosis <- factor(metadata$Diagnosis)
  metadata$Sex <- factor(metadata$Sex)

  metadata$Age <- suppressWarnings(as.numeric(metadata$Age))
  metadata$RIN <- suppressWarnings(as.numeric(metadata$RIN))

  message("Diagnosis levels:")
  print(levels(metadata$Diagnosis))

  message("Sex levels:")
  print(levels(metadata$Sex))

  message("Missing values in model covariates:")
  print(colSums(is.na(metadata[, c("Diagnosis", "Age", "Sex", "RIN")])))

  ## Remove samples with missing covariates before fitting any gene models.
  ## This makes the sample set explicit instead of letting each model drop rows silently.
  complete_covariate_rows <- complete.cases(metadata[, c("Diagnosis", "Age", "Sex", "RIN")])

  message("Samples before removing missing covariates: ", nrow(metadata))
  message("Samples after removing missing covariates:  ", sum(complete_covariate_rows))

  metadata <- metadata[complete_covariate_rows, , drop = FALSE]
  expression_matrix <- expression_matrix[, metadata$SampleID, drop = FALSE]

  if (nrow(metadata) < 3) {
    stop("Too few samples remain after removing missing covariates for ", group_name, ".")
  }

  ## -----------------------------------------------------------------------
  ## 7. Choose the Diagnosis coefficient to extract
  ## -----------------------------------------------------------------------

  ## For a two-level Diagnosis factor, R usually names the diagnosis coefficient
  ## by combining the column name and the non-reference level.
  ##
  ## Example:
  ##   If levels(metadata$Diagnosis) are CTRL and MDD,
  ##   the coefficient is often named DiagnosisMDD.
  ##
  ## To keep this script easy to inspect, we find Diagnosis coefficients once
  ## from the first gene that fits. If more than one Diagnosis coefficient exists,
  ## the script stops and asks you to choose one manually.
  diagnosis_coefficient_name <- NA

  ## -----------------------------------------------------------------------
  ## 8. Fit one linear model per gene
  ## -----------------------------------------------------------------------

  number_of_genes <- nrow(expression_matrix)

  results <- data.frame(
    GeneID = rownames(expression_matrix),
    diagnosis_term = NA_character_,
    diagnosis_coefficient = NA_real_,
    diagnosis_standard_error = NA_real_,
    diagnosis_p_value = NA_real_,
    n_samples_used = NA_integer_,
    model_status = NA_character_,
    stringsAsFactors = FALSE
  )

  message("Starting per-gene linear models for ", group_name, ".")
  message("This may take a while because the script fits one model per gene.")

  for (gene_index in seq_len(number_of_genes)) {
    gene_id <- rownames(expression_matrix)[gene_index]

    if (gene_index %% 500 == 0) {
      message("Fitted genes so far: ", gene_index, " / ", number_of_genes)
    }

    ## Pull one gene's expression across all aligned samples.
    expression_value <- as.numeric(expression_matrix[gene_index, ])

    ## This model_data table is what lm sees for one gene.
    model_data <- data.frame(
      expression_value = expression_value,
      Diagnosis = metadata$Diagnosis,
      Age = metadata$Age,
      Sex = metadata$Sex,
      RIN = metadata$RIN
    )

    ## Remove samples with missing expression for this gene.
    model_data <- model_data[complete.cases(model_data), , drop = FALSE]
    results$n_samples_used[gene_index] <- nrow(model_data)

    if (nrow(model_data) < 3) {
      results$model_status[gene_index] <- "skipped_too_few_samples"
      next
    }

    fitted_model <- tryCatch(
      lm(model_formula, data = model_data),
      error = function(e) e
    )

    if (inherits(fitted_model, "error")) {
      results$model_status[gene_index] <- paste("model_error:", fitted_model$message)
      next
    }

    coefficient_table <- summary(fitted_model)$coefficients
    coefficient_names <- rownames(coefficient_table)

    if (is.na(diagnosis_coefficient_name)) {
      diagnosis_terms_found <- coefficient_names[grepl("^Diagnosis", coefficient_names)]

      if (length(diagnosis_terms_found) == 0) {
        stop("No Diagnosis coefficient was found in the fitted model for gene: ", gene_id)
      }

      if (length(diagnosis_terms_found) > 1) {
        stop(
          "More than one Diagnosis coefficient was found for ",
          group_name,
          ": ",
          paste(diagnosis_terms_found, collapse = ", "),
          ". Edit diagnosis_coefficient_name manually before running the script."
        )
      }

      diagnosis_coefficient_name <- diagnosis_terms_found
      message("Diagnosis coefficient selected for extraction: ", diagnosis_coefficient_name)
    }

    if (!(diagnosis_coefficient_name %in% coefficient_names)) {
      results$model_status[gene_index] <- "diagnosis_term_not_found"
      next
    }

    results$diagnosis_term[gene_index] <- diagnosis_coefficient_name
    results$diagnosis_coefficient[gene_index] <- coefficient_table[diagnosis_coefficient_name, "Estimate"]
    results$diagnosis_standard_error[gene_index] <- coefficient_table[diagnosis_coefficient_name, "Std. Error"]
    results$diagnosis_p_value[gene_index] <- coefficient_table[diagnosis_coefficient_name, "Pr(>|t|)"]
    results$model_status[gene_index] <- "ok"
  }

  ## -----------------------------------------------------------------------
  ## 9. Multiple-testing correction and save results
  ## -----------------------------------------------------------------------

  results$diagnosis_FDR_BH <- p.adjust(results$diagnosis_p_value, method = "BH")

  results <- results[order(results$diagnosis_FDR_BH, results$diagnosis_p_value), ]

  write.csv(results, output_file, row.names = FALSE)

  message("Finished linear-model analysis for ", group_name, ".")
  message("Results saved to: ", output_file)
  message("Model status counts:")
  print(table(results$model_status, useNA = "ifany"))
}

message("")
message("Finished all small-group linear-model analyses.")
