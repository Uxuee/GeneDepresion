# ============================================================
# Step 2: Covariate exploration for GSE102556 Frontal_Cortex
# Input: log2(GeTMM + 1) expression matrix from Step 1
# Goal:
#   - Match expression samples to metadata
#   - Check diagnosis balance
#   - Check repeated donors
#   - Run PCA on top variable genes
#   - See which covariates explain PC structure
#   - Produce a model-decision summary for Step 3
# ============================================================

# ----------------------------
# 0. Packages
# ----------------------------
required_packages <- c("ggplot2", "dplyr", "readr", "tibble", "stringr", "tidyr")

missing_packages <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them first with:\ninstall.packages(c(",
    paste(sprintf('\"%s\"', missing_packages), collapse = ", "),
    "))"
  )
}

library(ggplot2)
library(dplyr)
library(readr)
library(tibble)
library(stringr)
library(tidyr)

# ----------------------------
# 1. Paths
# ----------------------------
base_dir <- "D:/genes/genedata"

expr_file <- file.path(
  base_dir,
  "Step1_GSE102556_GeTMM_brain_area_groups",
  "02_normalized_expression",
  "keep_all_samples",
  "normalized_GeTMM_log2_frontal_cortex.csv"
)

# IMPORTANT:
# If this auto-detection chooses the wrong metadata file,
# manually replace metadata_file with the correct path.
metadata_dir <- file.path(base_dir, "Step1_GSE102556_GeTMM_brain_area_groups", "01_metadata")

metadata_candidates <- list.files(
  metadata_dir,
  pattern = "frontal|Frontal|metadata|matched|keep_all",
  full.names = TRUE,
  recursive = TRUE
)

cat("\nMetadata candidates found:\n")
print(metadata_candidates)

if (length(metadata_candidates) == 0) {
  stop("No metadata candidates found. Please set metadata_file manually.")
}

metadata_file <- metadata_candidates[1]

cat("\nUsing expression file:\n", expr_file, "\n")
cat("\nUsing metadata file:\n", metadata_file, "\n")

out_dir <- file.path(
  base_dir,
  "Step2_covariate_exploration",
  "Frontal_Cortex_keep_all_samples"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 2. Read expression matrix
# ----------------------------
expr_df <- read.csv(expr_file, check.names = FALSE)

# First column is GeneID, even if it is called Unnamed: 0
colnames(expr_df)[1] <- "GeneID"

# Remove genes with duplicated IDs if any remain.
# Here we keep the first occurrence.
expr_df <- expr_df[!duplicated(expr_df$GeneID), ]

expr_mat <- as.matrix(expr_df[, -1])
rownames(expr_mat) <- as.character(expr_df$GeneID)
storage.mode(expr_mat) <- "numeric"

cat("\nExpression matrix dimensions:\n")
cat("Genes:", nrow(expr_mat), "\n")
cat("Samples:", ncol(expr_mat), "\n")

sample_ids_expr <- colnames(expr_mat)

# ----------------------------
# 3. Read metadata
# ----------------------------
metadata <- read.csv(metadata_file, check.names = FALSE)

# Try to standardize SampleID column if needed
possible_sample_cols <- c("SampleID", "geo_accession", "GSM", "sample", "Sample", "ID")

sample_col <- possible_sample_cols[possible_sample_cols %in% colnames(metadata)][1]

if (is.na(sample_col)) {
  stop(
    "Could not find a sample ID column in metadata. ",
    "Expected one of: ", paste(possible_sample_cols, collapse = ", ")
  )
}

metadata <- metadata %>%
  rename(SampleID = all_of(sample_col))

# Keep only metadata samples present in expression
metadata_fc <- metadata %>%
  filter(SampleID %in% sample_ids_expr)

# Reorder metadata to match expression columns
metadata_fc <- metadata_fc[match(sample_ids_expr, metadata_fc$SampleID), ]

if (!all(metadata_fc$SampleID == sample_ids_expr)) {
  stop("Metadata and expression sample order do not match.")
}

cat("\nMatched metadata rows:", nrow(metadata_fc), "\n")

# Save matched metadata
write.csv(
  metadata_fc,
  file.path(out_dir, "matched_metadata_frontal_cortex_keep_all_samples.csv"),
  row.names = FALSE
)

# ----------------------------
# 4. Basic metadata summary
# ----------------------------
metadata_summary <- data.frame(
  column = colnames(metadata_fc),
  n_missing = sapply(metadata_fc, function(x) sum(is.na(x) | x == "")),
  n_unique = sapply(metadata_fc, function(x) length(unique(x))),
  stringsAsFactors = FALSE
)

write.csv(
  metadata_summary,
  file.path(out_dir, "metadata_missingness_and_unique_counts.csv"),
  row.names = FALSE
)

# Diagnosis balance if available
if ("Diagnosis" %in% colnames(metadata_fc)) {
  diagnosis_balance <- metadata_fc %>%
    count(Diagnosis, name = "n_samples")

  write.csv(
    diagnosis_balance,
    file.path(out_dir, "diagnosis_balance.csv"),
    row.names = FALSE
  )

  print(diagnosis_balance)
}

# Diagnosis by original region
if (all(c("Diagnosis", "OriginalRegion") %in% colnames(metadata_fc))) {
  diagnosis_by_region <- metadata_fc %>%
    count(OriginalRegion, Diagnosis, name = "n_samples") %>%
    arrange(OriginalRegion, Diagnosis)

  write.csv(
    diagnosis_by_region,
    file.path(out_dir, "diagnosis_by_original_region.csv"),
    row.names = FALSE
  )

  p <- ggplot(diagnosis_by_region, aes(x = OriginalRegion, y = n_samples, fill = Diagnosis)) +
    geom_col(position = "dodge") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Diagnosis balance by original frontal region",
      x = "Original region",
      y = "Number of samples"
    )

  ggsave(
    file.path(out_dir, "diagnosis_by_original_region.png"),
    p,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# Repeated donor summary
if ("BrainID" %in% colnames(metadata_fc)) {
  donor_repetition <- metadata_fc %>%
    count(BrainID, name = "n_samples") %>%
    arrange(desc(n_samples))

  write.csv(
    donor_repetition,
    file.path(out_dir, "donor_repetition_from_matched_metadata.csv"),
    row.names = FALSE
  )

  cat("\nUnique donors:", nrow(donor_repetition), "\n")
  cat("Donors with repeated samples:", sum(donor_repetition$n_samples > 1), "\n")
}

# ----------------------------
# 5. Numeric covariate checks
# ----------------------------
numeric_covariates <- c("Age", "PMI", "pH", "RIN")
numeric_covariates <- numeric_covariates[numeric_covariates %in% colnames(metadata_fc)]

# Convert numeric covariates safely
for (cov in numeric_covariates) {
  metadata_fc[[cov]] <- suppressWarnings(as.numeric(metadata_fc[[cov]]))
}

if (length(numeric_covariates) > 0) {
  numeric_summary <- metadata_fc %>%
    summarise(
      across(
        all_of(numeric_covariates),
        list(
          n_missing = ~sum(is.na(.x)),
          mean = ~mean(.x, na.rm = TRUE),
          sd = ~sd(.x, na.rm = TRUE),
          min = ~min(.x, na.rm = TRUE),
          max = ~max(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )
    )

  write.csv(
    numeric_summary,
    file.path(out_dir, "numeric_covariate_summary.csv"),
    row.names = FALSE
  )
}

# Numeric covariates by diagnosis
if ("Diagnosis" %in% colnames(metadata_fc) && length(numeric_covariates) > 0) {
  for (cov in numeric_covariates) {
    p <- ggplot(metadata_fc, aes(x = Diagnosis, y = .data[[cov]])) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.7) +
      theme_bw() +
      labs(
        title = paste(cov, "by diagnosis"),
        x = "Diagnosis",
        y = cov
      )

    ggsave(
      file.path(out_dir, paste0("boxplot_", cov, "_by_diagnosis.png")),
      p,
      width = 5,
      height = 4,
      dpi = 300
    )
  }
}

# ----------------------------
# 6. PCA on top variable genes
# ----------------------------
# PCA should be done on samples, so we transpose expression matrix.
# To reduce noise, use top variable genes.
gene_variance <- apply(expr_mat, 1, var, na.rm = TRUE)

n_top_genes <- min(5000, length(gene_variance))
top_genes <- names(sort(gene_variance, decreasing = TRUE))[1:n_top_genes]

expr_top <- expr_mat[top_genes, ]

# Samples x genes
pca_input <- t(expr_top)

# Center genes, do not scale by default.
# Scaling can overweight low-variance genes, but here we already selected high-variance genes.
pca <- prcomp(pca_input, center = TRUE, scale. = FALSE)

pca_var <- (pca$sdev^2) / sum(pca$sdev^2)

pca_df <- as.data.frame(pca$x[, 1:10])
pca_df$SampleID <- rownames(pca_df)

pca_df <- pca_df %>%
  left_join(metadata_fc, by = "SampleID")

write.csv(
  pca_df,
  file.path(out_dir, "PCA_scores_with_metadata.csv"),
  row.names = FALSE
)

pca_variance_df <- data.frame(
  PC = paste0("PC", seq_along(pca_var)),
  variance_explained = pca_var
)

write.csv(
  pca_variance_df,
  file.path(out_dir, "PCA_variance_explained.csv"),
  row.names = FALSE
)

# Scree plot
scree <- pca_variance_df[1:20, ]

p_scree <- ggplot(scree, aes(x = factor(PC, levels = PC), y = variance_explained)) +
  geom_col() +
  theme_bw() +
  labs(
    title = "PCA variance explained",
    x = "Principal component",
    y = "Variance explained"
  )

ggsave(
  file.path(out_dir, "PCA_scree_plot_top5000_variable_genes.png"),
  p_scree,
  width = 7,
  height = 4,
  dpi = 300
)

# ----------------------------
# 7. PCA plots colored by covariates
# ----------------------------
categorical_covariates <- c("Diagnosis", "Sex", "OriginalRegion", "BrainArea", "BrainBank", "Batch")
categorical_covariates <- categorical_covariates[categorical_covariates %in% colnames(pca_df)]

for (cov in categorical_covariates) {
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[cov]])) +
    geom_point(size = 2.8, alpha = 0.85) +
    theme_bw() +
    labs(
      title = paste("PCA colored by", cov),
      x = paste0("PC1 (", round(100 * pca_var[1], 1), "%)"),
      y = paste0("PC2 (", round(100 * pca_var[2], 1), "%)"),
      color = cov
    )

  ggsave(
    file.path(out_dir, paste0("PCA_PC1_PC2_by_", cov, ".png")),
    p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

for (cov in numeric_covariates) {
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[cov]])) +
    geom_point(size = 2.8, alpha = 0.85) +
    theme_bw() +
    labs(
      title = paste("PCA colored by", cov),
      x = paste0("PC1 (", round(100 * pca_var[1], 1), "%)"),
      y = paste0("PC2 (", round(100 * pca_var[2], 1), "%)"),
      color = cov
    )

  ggsave(
    file.path(out_dir, paste0("PCA_PC1_PC2_by_", cov, ".png")),
    p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# ----------------------------
# 8. Test association between PCs and covariates
# ----------------------------
# This gives a rough guide:
# - For numeric covariates: linear model PC ~ covariate
# - For categorical covariates: ANOVA PC ~ covariate
# This is not the final DE model; it only helps decide covariates.

pc_cols <- paste0("PC", 1:10)

covariates_to_test <- c(
  "Diagnosis", "Sex", "Age", "PMI", "pH", "RIN",
  "OriginalRegion", "BrainArea", "BrainBank", "Batch"
)
covariates_to_test <- covariates_to_test[covariates_to_test %in% colnames(pca_df)]

pc_covariate_results <- list()

for (pc in pc_cols) {
  for (cov in covariates_to_test) {
    tmp <- pca_df[, c(pc, cov)]
    tmp <- tmp[complete.cases(tmp), ]

    if (nrow(tmp) < 5) next

    # If numeric, use lm. If categorical, use anova.
    if (is.numeric(tmp[[cov]])) {
      fit <- lm(tmp[[pc]] ~ tmp[[cov]])
      s <- summary(fit)
      pval <- coef(s)[2, 4]
      r2 <- s$r.squared
      test_type <- "linear_numeric"
    } else {
      tmp[[cov]] <- as.factor(tmp[[cov]])
      if (nlevels(tmp[[cov]]) < 2) next
      fit <- lm(tmp[[pc]] ~ tmp[[cov]])
      a <- anova(fit)
      pval <- a$`Pr(>F)`[1]
      r2 <- summary(fit)$r.squared
      test_type <- "anova_categorical"
    }

    pc_covariate_results[[length(pc_covariate_results) + 1]] <- data.frame(
      PC = pc,
      covariate = cov,
      test_type = test_type,
      p_value = pval,
      r_squared = r2
    )
  }
}

pc_covariate_results <- bind_rows(pc_covariate_results) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  arrange(FDR, desc(r_squared))

write.csv(
  pc_covariate_results,
  file.path(out_dir, "PC_covariate_association_results.csv"),
  row.names = FALSE
)

# Heatmap-like plot of R2 values
if (nrow(pc_covariate_results) > 0) {
  p_heat <- ggplot(pc_covariate_results, aes(x = PC, y = covariate, fill = r_squared)) +
    geom_tile() +
    theme_bw() +
    labs(
      title = "Covariate association with principal components",
      x = "Principal component",
      y = "Covariate",
      fill = "R²"
    )

  ggsave(
    file.path(out_dir, "PC_covariate_R2_heatmap.png"),
    p_heat,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# ----------------------------
# 9. Collinearity check among numeric covariates
# ----------------------------
if (length(numeric_covariates) >= 2) {
  numeric_meta <- metadata_fc[, numeric_covariates, drop = FALSE]
  numeric_cor <- cor(numeric_meta, use = "pairwise.complete.obs")

  write.csv(
    numeric_cor,
    file.path(out_dir, "numeric_covariate_correlation_matrix.csv")
  )
}

# Categorical covariate contingency tables
if (all(c("Diagnosis", "Sex") %in% colnames(metadata_fc))) {
  write.csv(
    as.data.frame.matrix(table(metadata_fc$Diagnosis, metadata_fc$Sex)),
    file.path(out_dir, "contingency_Diagnosis_by_Sex.csv")
  )
}

if (all(c("Diagnosis", "OriginalRegion") %in% colnames(metadata_fc))) {
  write.csv(
    as.data.frame.matrix(table(metadata_fc$Diagnosis, metadata_fc$OriginalRegion)),
    file.path(out_dir, "contingency_Diagnosis_by_OriginalRegion.csv")
  )
}

# ----------------------------
# 10. Draft model decision summary
# ----------------------------
model_notes <- c(
  "Recommended main analysis file: Frontal_Cortex keep_all_samples log2(GeTMM + 1).",
  "Because Frontal_Cortex combines OFC, dlPFC, and Cg25, repeated BrainID values are expected.",
  "If BrainID repeats, the main DE model should use a random intercept for BrainID.",
  "OriginalRegion should probably be included as a fixed effect because Frontal_Cortex contains multiple original frontal regions.",
  "Age, Sex, RIN, PMI, and pH should be considered if present and not too collinear.",
  "Diagnosis is the main biological variable of interest, not a nuisance covariate.",
  "PCA and PC-covariate association outputs should be inspected before finalizing the model."
)

writeLines(
  model_notes,
  con = file.path(out_dir, "draft_model_decision_notes.txt")
)

# ----------------------------
# 11. Session info
# ----------------------------
sink(file.path(out_dir, "session_info.txt"))
sessionInfo()
sink()

cat("\nDone. Step 2 covariate exploration outputs written to:\n")
cat(out_dir, "\n")
