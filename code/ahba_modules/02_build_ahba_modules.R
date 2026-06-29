#!/usr/bin/env Rscript

# Build AHBA gene coexpression modules from processed gene x sample matrix.
# Prefer WGCNA when installed. This script does not download AHBA data.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "code/ahba_modules/02_build_ahba_modules.R"
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)

clean_dir <- file.path(project_root, "outputs", "length_gt_100_clean_pipeline")
source_length100_dir <- file.path(project_root, "outputs", "length_gt_100_elisa_style")
processed_dir <- file.path(clean_dir, "ahba_modules", "processed")
out_dir <- file.path(clean_dir, "ahba_modules", "modules")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

env_flag <- function(name, default = TRUE) {
  value <- Sys.getenv(name, unset = if (default) "true" else "false")
  tolower(value) %in% c("1", "true", "yes", "y")
}

env_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value)) default else value
}

restrict_to_mdd_tested <- env_flag("AHBA_RESTRICT_TO_MDD_TESTED", TRUE)
max_genes_for_wgcna <- env_int("AHBA_MAX_GENES_FOR_WGCNA", 10000)

expr_path <- file.path(processed_dir, "ahba_gene_expression_matrix_gene_by_sample.csv")
if (!file.exists(expr_path)) {
  msg <- paste("Processed AHBA matrix not found:", expr_path, "\nRun python code/ahba_modules/01_preprocess_ahba_microarray.py first.")
  writeLines(msg, file.path(out_dir, "AHBA_MODULES_NOT_RUN.txt"))
  stop(msg)
}

expr_df <- read.csv(expr_path, check.names = FALSE)
gene_symbols_raw <- expr_df[[1]]
keep_gene_symbols <- !is.na(gene_symbols_raw) & trimws(as.character(gene_symbols_raw)) != ""
removed_blank_gene_symbols <- sum(!keep_gene_symbols)
message("Removed blank/NA gene-symbol rows: ", removed_blank_gene_symbols)

expr_df <- expr_df[keep_gene_symbols, , drop = FALSE]
gene_symbols <- trimws(as.character(expr_df[[1]]))
if (anyDuplicated(gene_symbols) > 0) {
  stop("Duplicated gene symbols found after filtering blank rows. Please fix preprocessing before WGCNA.")
}

rownames(expr_df) <- gene_symbols
expr <- as.matrix(expr_df[, -1, drop = FALSE])
mode(expr) <- "numeric"

# WGCNA expects samples x genes.
datExpr <- t(expr)
datExpr <- datExpr[, colSums(is.na(datExpr)) == 0, drop = FALSE]
datExpr <- datExpr[, apply(datExpr, 2, sd, na.rm = TRUE) > 0, drop = FALSE]

load_mdd_tested_symbols <- function() {
  candidates <- c(
    file.path(clean_dir, "01_inputs", "tested_gene_universe_gene_symbols.csv"),
    file.path(clean_dir, "03_broad_region_LMM", "LMM_all_genes_length_gt_100.csv"),
    file.path(clean_dir, "04_standardized_effect_size", "LMM_all_genes_with_standardized_effect_size_4groups.csv"),
    file.path(source_length100_dir, "LMM_all_genes_length_gt_100.csv"),
    file.path(source_length100_dir, "standardized_effect_size_4groups", "LMM_all_genes_with_standardized_effect_size_4groups.csv")
  )
  for (path in candidates) {
    if (!file.exists(path)) next
    header <- names(read.csv(path, nrows = 1, check.names = FALSE))
    if (!"gene_symbol" %in% header) next
    tbl <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    genes <- unique(toupper(trimws(as.character(tbl$gene_symbol))))
    genes <- genes[!is.na(genes) & genes != ""]
    if (length(genes) > 0) {
      return(list(genes = genes, source = path))
    }
  }
  stop("Could not find MDD tested-gene table with gene_symbol for AHBA WGCNA restriction.")
}

gene_selection <- data.frame(
  step = "post_missingness_and_variance_filter",
  n_genes = ncol(datExpr),
  details = "AHBA genes with no missing values and nonzero variance",
  stringsAsFactors = FALSE
)

mdd_tested_source <- NA_character_
if (restrict_to_mdd_tested) {
  tested <- load_mdd_tested_symbols()
  mdd_tested_source <- tested$source
  before <- ncol(datExpr)
  keep <- toupper(colnames(datExpr)) %in% tested$genes
  datExpr <- datExpr[, keep, drop = FALSE]
  gene_selection <- rbind(
    gene_selection,
    data.frame(
      step = "restricted_to_mdd_tested_ahba_genes",
      n_genes = ncol(datExpr),
      details = paste("Kept", ncol(datExpr), "of", before, "AHBA genes using", mdd_tested_source),
      stringsAsFactors = FALSE
    )
  )
}

if (max_genes_for_wgcna > 0 && ncol(datExpr) > max_genes_for_wgcna) {
  before <- ncol(datExpr)
  gene_var <- apply(datExpr, 2, var, na.rm = TRUE)
  keep_names <- names(sort(gene_var, decreasing = TRUE))[seq_len(max_genes_for_wgcna)]
  datExpr <- datExpr[, keep_names, drop = FALSE]
  gene_selection <- rbind(
    gene_selection,
    data.frame(
      step = "restricted_to_top_variable_genes",
      n_genes = ncol(datExpr),
      details = paste("Kept top", max_genes_for_wgcna, "variable genes from", before, "MDD-tested AHBA genes"),
      stringsAsFactors = FALSE
    )
  )
}

if (ncol(datExpr) < 50) {
  stop("Too few genes remain for WGCNA after AHBA gene selection: ", ncol(datExpr))
}

write.csv(gene_selection, file.path(out_dir, "ahba_wgcna_gene_selection_audit.csv"), row.names = FALSE)

params <- list(
  input_matrix = expr_path,
  n_samples = nrow(datExpr),
  n_genes = ncol(datExpr),
  method = NA,
  soft_power = NA,
  min_module_size = 30,
  merge_cut_height = 0.25,
  restrict_to_mdd_tested_genes = restrict_to_mdd_tested,
  mdd_tested_gene_source = mdd_tested_source,
  max_genes_for_wgcna = max_genes_for_wgcna,
  gene_selection_audit = file.path(out_dir, "ahba_wgcna_gene_selection_audit.csv")
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
not_installed_marker <- file.path(out_dir, "WGCNA_NOT_INSTALLED.txt")
if (file.exists(not_installed_marker)) {
  file.remove(not_installed_marker)
}

powers <- c(1:10, seq(12, 20, by = 2))
soft_threshold_status <- "pickSoftThreshold_completed"
sft <- tryCatch(
  pickSoftThreshold(datExpr, powerVector = powers, verbose = 0),
  error = function(e) {
    soft_threshold_status <<- paste("pickSoftThreshold_failed_using_default_power_6:", conditionMessage(e))
    message(soft_threshold_status)
    NULL
  }
)
if (is.null(sft)) {
  fit <- data.frame()
  soft_power <- 6
} else {
  fit <- sft$fitIndices
  candidate <- fit$Power[which(fit$SFT.R.sq >= 0.80)[1]]
  soft_power <- ifelse(is.na(candidate), 6, candidate)
}

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
params$soft_threshold_status <- soft_threshold_status
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
  paste("- Restricted to MDD-tested genes:", restrict_to_mdd_tested),
  paste("- MDD tested-gene source:", mdd_tested_source),
  paste("- Maximum genes for WGCNA:", max_genes_for_wgcna),
  paste("- Gene-selection audit:", params$gene_selection_audit),
  paste("- Soft threshold power:", soft_power),
  "- Network type: signed",
  "- Minimum module size: 30",
  "- Merge cut height: 0.25",
  paste("- Modules:", length(unique(modules$ahba_module)))
)
writeLines(report, file.path(out_dir, "AHBA_module_construction_summary.md"))
