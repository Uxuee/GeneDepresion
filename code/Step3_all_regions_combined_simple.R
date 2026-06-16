## Step3_all_regions_combined_simple.R
##
## Purpose:
##   Test for a shared depression effect across all available brain-area groups.
##
## Why combine regions?
##   Separate region analyses have small sample sizes.
##   A combined model can ask whether there is a diagnosis effect shared across
##   regions while adjusting for the brain area.
##
## Why dream here?
##   Donors can repeat across regions, so samples are not fully independent.
##   This script therefore uses dream with (1 | BrainID).
##
## Model:
##   expression ~ Diagnosis + BrainArea + Age + Sex + RIN + (1 | BrainID)
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

frontal_counts_file <- file.path(counts_folder, "filtered_raw_counts_frontal_cortex.csv")
insula_counts_file <- file.path(counts_folder, "filtered_raw_counts_insula.csv")
ventral_striatum_counts_file <- file.path(counts_folder, "filtered_raw_counts_ventral_striatum.csv")
hippocampal_counts_file <- file.path(counts_folder, "filtered_raw_counts_hippocampal_formation.csv")

frontal_metadata_file <- file.path(metadata_folder, "metadata_frontal_cortex.csv")
insula_metadata_file <- file.path(metadata_folder, "metadata_insula.csv")
ventral_striatum_metadata_file <- file.path(metadata_folder, "metadata_ventral_striatum.csv")
hippocampal_metadata_file <- file.path(metadata_folder, "metadata_hippocampal_formation.csv")

results_file <- file.path(output_folder, "dream_all_regions_combined_results.csv")
histogram_file <- file.path(output_folder, "dream_all_regions_combined_diagnosis_pvalue_histogram.png")

message("Output results file:   ", results_file)
message("Output histogram file: ", histogram_file)

input_files <- c(
  frontal_counts_file,
  insula_counts_file,
  ventral_striatum_counts_file,
  hippocampal_counts_file,
  frontal_metadata_file,
  insula_metadata_file,
  ventral_striatum_metadata_file,
  hippocampal_metadata_file
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop("Missing input files:\n", paste(missing_files, collapse = "\n"))
}

## -------------------------------------------------------------------------
## 2. Load each count matrix explicitly
## -------------------------------------------------------------------------

## The first column contains gene IDs.
## Rows are genes. Columns are samples.
frontal_counts <- read.csv(frontal_counts_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
insula_counts <- read.csv(insula_counts_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
ventral_striatum_counts <- read.csv(ventral_striatum_counts_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
hippocampal_counts <- read.csv(hippocampal_counts_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

message("Count matrix dimensions before combining:")
message("Frontal_Cortex:        ", paste(dim(frontal_counts), collapse = " x "))
message("Insula:                ", paste(dim(insula_counts), collapse = " x "))
message("Ventral_Striatum:      ", paste(dim(ventral_striatum_counts), collapse = " x "))
message("Hippocampal_Formation: ", paste(dim(hippocampal_counts), collapse = " x "))

## -------------------------------------------------------------------------
## 3. Load each metadata table explicitly
## -------------------------------------------------------------------------

frontal_metadata <- read.csv(frontal_metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
insula_metadata <- read.csv(insula_metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
ventral_striatum_metadata <- read.csv(ventral_striatum_metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
hippocampal_metadata <- read.csv(hippocampal_metadata_file, check.names = FALSE, stringsAsFactors = FALSE)

message("Metadata rows before combining:")
message("Frontal_Cortex:        ", nrow(frontal_metadata))
message("Insula:                ", nrow(insula_metadata))
message("Ventral_Striatum:      ", nrow(ventral_striatum_metadata))
message("Hippocampal_Formation: ", nrow(hippocampal_metadata))

## -------------------------------------------------------------------------
## 4. Keep only genes present in all four count matrices
## -------------------------------------------------------------------------

common_genes <- Reduce(
  intersect,
  list(
    rownames(frontal_counts),
    rownames(insula_counts),
    rownames(ventral_striatum_counts),
    rownames(hippocampal_counts)
  )
)

message("Genes present in all four groups: ", length(common_genes))

if (length(common_genes) == 0) {
  stop("No common genes were found across the four count matrices.")
}

frontal_counts <- frontal_counts[common_genes, , drop = FALSE]
insula_counts <- insula_counts[common_genes, , drop = FALSE]
ventral_striatum_counts <- ventral_striatum_counts[common_genes, , drop = FALSE]
hippocampal_counts <- hippocampal_counts[common_genes, , drop = FALSE]

## -------------------------------------------------------------------------
## 5. Combine count matrices and metadata
## -------------------------------------------------------------------------

combined_counts <- cbind(
  frontal_counts,
  insula_counts,
  ventral_striatum_counts,
  hippocampal_counts
)

combined_metadata <- rbind(
  frontal_metadata,
  insula_metadata,
  ventral_striatum_metadata,
  hippocampal_metadata
)

message("Combined count matrix dimensions:")
print(dim(combined_counts))

message("Combined metadata rows: ", nrow(combined_metadata))

if (any(duplicated(colnames(combined_counts)))) {
  stop("Combined count matrix contains duplicated sample IDs.")
}

if (any(duplicated(combined_metadata$SampleID))) {
  stop("Combined metadata contains duplicated SampleID values.")
}

## -------------------------------------------------------------------------
## 6. Check required metadata columns
## -------------------------------------------------------------------------

required_columns <- c("SampleID", "Diagnosis", "BrainArea", "Age", "Sex", "RIN", "BrainID")
missing_columns <- setdiff(required_columns, names(combined_metadata))

if (length(missing_columns) > 0) {
  stop("Missing required metadata columns: ", paste(missing_columns, collapse = ", "))
}

## -------------------------------------------------------------------------
## 7. Align count columns and metadata rows
## -------------------------------------------------------------------------

common_samples <- intersect(colnames(combined_counts), combined_metadata$SampleID)

message("Samples found in both combined counts and metadata: ", length(common_samples))

if (length(common_samples) == 0) {
  stop("No sample IDs matched between combined counts and combined metadata.")
}

combined_counts <- combined_counts[, common_samples, drop = FALSE]
combined_metadata <- combined_metadata[match(common_samples, combined_metadata$SampleID), , drop = FALSE]

message("First count sample IDs after matching:")
print(head(colnames(combined_counts)))

message("First metadata SampleID values after matching:")
print(head(combined_metadata$SampleID))

if (!all(colnames(combined_counts) == combined_metadata$SampleID)) {
  stop("Combined counts and metadata are not aligned after matching.")
}

## -------------------------------------------------------------------------
## 8. Prepare metadata for the combined dream model
## -------------------------------------------------------------------------

diagnosis_values <- sort(unique(combined_metadata$Diagnosis))
message("Diagnosis values found:")
print(diagnosis_values)

if (!identical(diagnosis_values, c("CTRL", "MDD"))) {
  stop("Expected Diagnosis values to be exactly CTRL and MDD.")
}

combined_metadata$Diagnosis <- factor(combined_metadata$Diagnosis, levels = c("CTRL", "MDD"))
combined_metadata$BrainArea <- factor(combined_metadata$BrainArea)
combined_metadata$Sex <- factor(combined_metadata$Sex)
combined_metadata$BrainID <- factor(combined_metadata$BrainID)
combined_metadata$Age <- suppressWarnings(as.numeric(combined_metadata$Age))
combined_metadata$RIN <- suppressWarnings(as.numeric(combined_metadata$RIN))

message("Diagnosis counts:")
print(table(combined_metadata$Diagnosis, useNA = "ifany"))

message("BrainArea counts:")
print(table(combined_metadata$BrainArea, useNA = "ifany"))

message("BrainArea by Diagnosis:")
print(table(combined_metadata$BrainArea, combined_metadata$Diagnosis, useNA = "ifany"))

message("Number of samples: ", nrow(combined_metadata))
message("Number of unique BrainID donors: ", length(unique(combined_metadata$BrainID)))
message("Number of repeated BrainID donors: ", sum(table(combined_metadata$BrainID) > 1))

message("Missing values in model covariates:")
print(colSums(is.na(combined_metadata[, c("Diagnosis", "BrainArea", "Age", "Sex", "RIN", "BrainID")])))

complete_rows <- complete.cases(combined_metadata[, c("Diagnosis", "BrainArea", "Age", "Sex", "RIN", "BrainID")])

message("Samples before removing missing covariates: ", nrow(combined_metadata))
message("Samples after removing missing covariates:  ", sum(complete_rows))

combined_metadata <- combined_metadata[complete_rows, , drop = FALSE]
combined_counts <- combined_counts[, combined_metadata$SampleID, drop = FALSE]

if (nrow(combined_metadata) < 3) {
  stop("Too few samples remain after removing missing covariates.")
}

## -------------------------------------------------------------------------
## 9. Convert counts to a numeric matrix
## -------------------------------------------------------------------------

combined_counts_matrix <- as.matrix(combined_counts)
storage.mode(combined_counts_matrix) <- "numeric"

if (any(is.na(combined_counts_matrix))) {
  stop("The combined count matrix contains NA values after numeric conversion.")
}

if (any(combined_counts_matrix < 0)) {
  stop("The combined count matrix contains negative values, which are invalid for raw counts.")
}

message("Final combined count matrix dimensions:")
print(dim(combined_counts_matrix))

## -------------------------------------------------------------------------
## 10. Write the combined dream model formula explicitly
## -------------------------------------------------------------------------

## We use BrainArea as a fixed covariate so the diagnosis effect is tested
## across regions while adjusting for broad region differences.
##
## We use (1 | BrainID) because some donors appear in more than one region.
## That means samples from the same donor are not fully independent.
model_formula <- ~ Diagnosis + BrainArea + Age + Sex + RIN + (1 | BrainID)

message("Combined dream model formula:")
print(model_formula)

## -------------------------------------------------------------------------
## 11. Run edgeR library-size normalization for voom/dream
## -------------------------------------------------------------------------

## This does not recreate the Step 1 GeTMM-log2 output.
## This is the standard count-based scaling step used inside voom/dream.
dge_object <- edgeR::DGEList(counts = combined_counts_matrix)
dge_object <- edgeR::calcNormFactors(dge_object, method = "TMM")

message("Library sizes and TMM factors:")
print(head(dge_object$samples))

## -------------------------------------------------------------------------
## 12. Run voomWithDreamWeights and dream
## -------------------------------------------------------------------------

## SerialParam keeps this simple and avoids parallel execution.
serial_parameter <- BiocParallel::SerialParam()

voom_object <- variancePartition::voomWithDreamWeights(
  dge_object,
  formula = model_formula,
  data = combined_metadata,
  plot = FALSE,
  BPPARAM = serial_parameter
)

dream_fit <- variancePartition::dream(
  voom_object,
  formula = model_formula,
  data = combined_metadata,
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
## 13. Extract Diagnosis results
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
## 14. Save CSV results
## -------------------------------------------------------------------------

write.csv(final_results, results_file, row.names = FALSE)

message("Saved results to:")
message(results_file)

## -------------------------------------------------------------------------
## 15. Save a simple p-value histogram
## -------------------------------------------------------------------------

png(filename = histogram_file, width = 900, height = 700)
hist(
  final_results$diagnosis_p_value,
  breaks = 50,
  main = "Diagnosis p-values: all regions combined dream",
  xlab = "Raw p-value for DiagnosisMDD",
  col = "gray80",
  border = "white"
)
dev.off()

message("Saved p-value histogram to:")
message(histogram_file)

message("")
message("Finished all-regions combined dream analysis.")
