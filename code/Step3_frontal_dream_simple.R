## Step3_frontal_dream_simple.R
##
## Purpose:
##   Differential expression for Frontal_Cortex using dream.
##
## Why dream here?
##   Frontal_Cortex combines dlPFC, OFC, and Cg25.
##   That creates repeated BrainID values because the same donor can contribute
##   more than one frontal-region sample.
##
## Model:
##   expression ~ Diagnosis + Age + Sex + RIN + OriginalRegion + (1 | BrainID)
##
## Important assumptions:
##   1. PCA/QC has already been done.
##   2. Step 1 preprocessing has already been done.
##   3. This script does not run PCA.
##   4. This script does not recreate the Step 1 GeTMM-log2 files.
##   5. dream starts from filtered raw counts and estimates voom weights.
##   6. Count matrix rows are genes and count matrix columns are samples.
##   7. Metadata has one row per sample and contains SampleID.

## -------------------------------------------------------------------------
## 0. Packages
## -------------------------------------------------------------------------

required_packages <- c("edgeR", "limma", "variancePartition", "BiocParallel")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script.",
    "\nFor variancePartition/dream, use BiocManager::install('variancePartition')."
  )
}

library(edgeR)
library(limma)
library(variancePartition)
library(BiocParallel)

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

counts_file <- file.path(counts_folder, "filtered_raw_counts_frontal_cortex.csv")
metadata_file <- file.path(metadata_folder, "metadata_frontal_cortex.csv")
results_file <- file.path(output_folder, "dream_Frontal_Cortex_results.csv")
histogram_file <- file.path(output_folder, "dream_Frontal_Cortex_diagnosis_pvalue_histogram.png")

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

## -------------------------------------------------------------------------
## 2. Load filtered raw counts and metadata
## -------------------------------------------------------------------------

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

## -------------------------------------------------------------------------
## 3. Check required columns
## -------------------------------------------------------------------------

required_columns <- c("SampleID", "Diagnosis", "Age", "Sex", "RIN", "OriginalRegion", "BrainID")
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

## -------------------------------------------------------------------------
## 4. Align count columns and metadata rows
## -------------------------------------------------------------------------

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

## -------------------------------------------------------------------------
## 5. Prepare metadata for dream
## -------------------------------------------------------------------------

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

metadata <- metadata[complete_rows, , drop = FALSE]
raw_counts <- raw_counts[, metadata$SampleID, drop = FALSE]

if (nrow(metadata) < 3) {
  stop("Too few samples remain after removing missing covariates.")
}

## -------------------------------------------------------------------------
## 6. Convert counts to a numeric matrix
## -------------------------------------------------------------------------

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

## -------------------------------------------------------------------------
## 7. Write the dream model formula explicitly
## -------------------------------------------------------------------------

## DiagnosisMDD is the term of interest.
## OriginalRegion adjusts for dlPFC/OFC/Cg25 differences.
## (1 | BrainID) accounts for repeated donor samples.
model_formula <- ~ Diagnosis + Age + Sex + RIN + OriginalRegion + (1 | BrainID)

message("dream model formula:")
print(model_formula)

## -------------------------------------------------------------------------
## 8. Run edgeR library-size normalization for voom/dream
## -------------------------------------------------------------------------

## This does not recreate the Step 1 GeTMM-log2 output.
## This is the standard count-based scaling step used inside voom/dream.
dge_object <- edgeR::DGEList(counts = raw_counts_matrix)
dge_object <- edgeR::calcNormFactors(dge_object, method = "TMM")

message("Library sizes and TMM factors:")
print(head(dge_object$samples))

## -------------------------------------------------------------------------
## 9. Run voomWithDreamWeights and dream
## -------------------------------------------------------------------------

## SerialParam keeps this simple and avoids parallel execution.
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

dream_fit <- limma::eBayes(dream_fit)

message("Model coefficients available:")
print(colnames(dream_fit$coefficients))

diagnosis_coefficient <- "DiagnosisMDD"

if (!(diagnosis_coefficient %in% colnames(dream_fit$coefficients))) {
  stop("Could not find expected coefficient: ", diagnosis_coefficient)
}

## -------------------------------------------------------------------------
## 10. Extract Diagnosis results
## -------------------------------------------------------------------------

dream_results <- limma::topTable(
  dream_fit,
  coef = diagnosis_coefficient,
  number = Inf,
  sort.by = "P"
)

final_results <- data.frame(
  GeneID = rownames(dream_results),
  diagnosis_term = diagnosis_coefficient,
  diagnosis_coefficient = dream_results$logFC,
  average_expression = dream_results$AveExpr,
  moderated_t = dream_results$t,
  diagnosis_p_value = dream_results$P.Value,
  diagnosis_FDR_BH = dream_results$adj.P.Val,
  log_odds_DE = dream_results$B,
  stringsAsFactors = FALSE
)

message("Number of genes in result table: ", nrow(final_results))
message("Genes with raw p < 0.05: ", sum(final_results$diagnosis_p_value < 0.05, na.rm = TRUE))
message("Genes with FDR < 0.05:   ", sum(final_results$diagnosis_FDR_BH < 0.05, na.rm = TRUE))
message("Genes with FDR < 0.10:   ", sum(final_results$diagnosis_FDR_BH < 0.10, na.rm = TRUE))

## -------------------------------------------------------------------------
## 11. Save CSV results
## -------------------------------------------------------------------------

write.csv(final_results, results_file, row.names = FALSE)

message("Saved results to:")
message(results_file)

## -------------------------------------------------------------------------
## 12. Save a simple p-value histogram
## -------------------------------------------------------------------------

png(filename = histogram_file, width = 900, height = 700)
hist(
  final_results$diagnosis_p_value,
  breaks = 50,
  main = "Diagnosis p-values: Frontal_Cortex dream",
  xlab = "Raw p-value for DiagnosisMDD",
  col = "gray80",
  border = "white"
)
dev.off()

message("Saved p-value histogram to:")
message(histogram_file)

message("")
message("Finished Frontal_Cortex dream analysis.")
