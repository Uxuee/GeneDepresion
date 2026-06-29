#!/usr/bin/env Rscript

# Build AHBA gene coexpression modules from processed gene x sample matrix.
# Prefer WGCNA when installed. This script does not download AHBA data.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "code/ahba_modules/02_build_ahba_modules.R"
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)

clean_dir <- file.path(project_root, "outputs", "length_gt_100_clean_pipeline")
processed_dir <- file.path(clean_dir, "ahba_modules", "processed")
out_dir <- file.path(clean_dir, "ahba_modules", "modules")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr_path <- file.path(processed_dir, "ahba_gene_expression_matrix_gene_by_sample.csv")
if (!file.exists(expr_path)) {
  msg <- paste("Processed AHBA matrix not found:", expr_path, "\nRun python code/ahba_modules/01_preprocess_ahba_microarray.py first.")
  writeLines(msg, file.path(out_dir, "AHBA_MODULES_NOT_RUN.txt"))
  stop(msg)
}

expr <- read.csv(expr_path, row.names = 1, check.names = FALSE)
expr <- as.matrix(expr)
mode(expr) <- "numeric"

# WGCNA expects samples x genes.
datExpr <- t(expr)
datExpr <- datExpr[, colSums(is.na(datExpr)) == 0, drop = FALSE]
datExpr <- datExpr[, apply(datExpr, 2, sd, na.rm = TRUE) > 0, drop = FALSE]

params <- list(
  input_matrix = expr_path,
  n_samples = nrow(datExpr),
  n_genes = ncol(datExpr),
  method = NA,
  soft_power = NA,
  min_module_size = 30,
  merge_cut_height = 0.25
)

write_params <- function(params, path) {
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    writeLines(jsonlite::toJSON(params, pretty = TRUE, auto_unbox = TRUE), path)
  } else {
    capture.output(str(params), file = path)
  }
}

if (!requireNamespace("WGCNA", quietly = TRUE)) {
  msg <- paste(
    "WGCNA is not installed. Install it with BiocManager::install('WGCNA') and rerun.",
    "No fallback modules were generated to avoid changing the intended method."
  )
  writeLines(msg, file.path(out_dir, "WGCNA_NOT_INSTALLED.txt"))
  write_params(params, file.path(out_dir, "module_parameters.json"))
  stop(msg)
}

suppressPackageStartupMessages(library(WGCNA))
allowWGCNAThreads()

powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
fit <- sft$fitIndices
candidate <- fit$Power[which(fit$SFT.R.sq >= 0.80)[1]]
soft_power <- ifelse(is.na(candidate), 6, candidate)

net <- blockwiseModules(
  datExpr,
  power = soft_power,
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  verbose = 2
)

modules <- data.frame(
  gene_symbol = names(net$colors),
  ahba_module = as.character(net$colors),
  stringsAsFactors = FALSE
)
write.csv(modules, file.path(out_dir, "ahba_gene_to_module_assignments.csv"), row.names = FALSE)

module_sizes <- as.data.frame(table(modules$ahba_module), stringsAsFactors = FALSE)
names(module_sizes) <- c("ahba_module", "n_genes")
module_sizes <- module_sizes[order(module_sizes$n_genes, decreasing = TRUE), ]
write.csv(module_sizes, file.path(out_dir, "ahba_module_size_table.csv"), row.names = FALSE)

eigengenes <- as.data.frame(net$MEs)
eigengenes$ahba_sample_id <- rownames(eigengenes)
write.csv(eigengenes, file.path(out_dir, "ahba_module_eigengenes.csv"), row.names = FALSE)

params$method <- "WGCNA::blockwiseModules"
params$soft_power <- soft_power
params$selected_power_fit_table <- file.path(out_dir, "ahba_wgcna_soft_threshold_fit.csv")
write.csv(fit, params$selected_power_fit_table, row.names = FALSE)
write_params(params, file.path(out_dir, "module_parameters.json"))

report <- c(
  "# AHBA Module Construction Summary",
  "",
  paste("- Input matrix:", expr_path),
  paste("- Method:", params$method),
  paste("- Samples:", nrow(datExpr)),
  paste("- Genes:", ncol(datExpr)),
  paste("- Soft threshold power:", soft_power),
  "- Network type: signed",
  "- Minimum module size: 30",
  "- Merge cut height: 0.25",
  paste("- Modules:", length(unique(modules$ahba_module)))
)
writeLines(report, file.path(out_dir, "AHBA_module_construction_summary.md"))
