## Step6_MDD_WGCNA_modules.R
##
## Purpose:
##   Build MDD/control coexpression modules directly from the normalized
##   GeTMM-log2 expression matrices.
##
## Key method decision:
##   WGCNA is unsupervised. This script does NOT build modules only from
##   differentially expressed genes, Tier 2 genes, or ML-prioritized genes.
##   It builds modules from many expressed/informative genes, then tests
##   whether existing candidate genes are enriched inside the modules.

options(stringsAsFactors = FALSE)

## -------------------------------------------------------------------------
## 0. Packages
## -------------------------------------------------------------------------

required_packages <- c("WGCNA", "dynamicTreeCut", "readr", "dplyr", "tidyr", "stringr", "limma")
missing_required <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_required) > 0) {
  stop(
    "Missing required packages: ", paste(missing_required, collapse = ", "),
    "\nInstall them before running this script.",
    "\nExample: install.packages('WGCNA'); BiocManager::install('limma')"
  )
}

suppressPackageStartupMessages({
  library(WGCNA)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(limma)
})

has_broom <- requireNamespace("broom", quietly = TRUE)
if (!has_broom) {
  warning("Optional package 'broom' is not installed. Linear-model output will use base R.")
}

WGCNA::allowWGCNAThreads()

## -------------------------------------------------------------------------
## 1. File locations and run settings
## -------------------------------------------------------------------------

project_root <- "D:/genes/GeneDepresion"
pipeline_root <- file.path(project_root, "outputs", "length_gt_100_clean_pipeline")

metadata_file <- file.path(pipeline_root, "01_inputs", "sample_metadata_used.csv")
expression_folder <- file.path(pipeline_root, "02_normalized_expression")
tier2_file <- file.path(
  pipeline_root,
  "07_Tier2_candidate_gene_sets",
  "Tier2_candidate_genes_all_regions.csv"
)
ml_file <- file.path(
  pipeline_root,
  "10_ML_prioritized_Tier2_corrected",
  "ML_prioritized_Tier2_genes_all_ranked.csv"
)
ahba_file <- file.path(
  pipeline_root,
  "11_AHBA_modules_enrichment",
  "published_modules",
  "ahba_published_gene_to_module_assignments.csv"
)
gene_annotation_folder <- file.path(
  "D:/genes/genedata/Step1_GSE102556_GeTMM_brain_area_groups",
  "05_gene_annotation",
  "keep_all_samples"
)

output_folder <- file.path(pipeline_root, "12_MDD_WGCNA_modules")
plot_folder <- file.path(output_folder, "plots")
best_folder <- file.path(output_folder, "best_solutions")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(best_folder, recursive = TRUE, showWarnings = FALSE)

expression_files <- list.files(
  expression_folder,
  pattern = "^normalized_GeTMM_log2_.*\\.csv$",
  full.names = TRUE
)

gene_filter_targets <- c(5000, 10000, 15000)
gene_filter_targets_env <- Sys.getenv("STEP6_GENE_FILTER_TARGETS", unset = "")
if (nzchar(gene_filter_targets_env)) {
  parsed_filter_targets <- strsplit(gene_filter_targets_env, ",")[[1]] %>% trimws()
  parsed_filter_targets <- parsed_filter_targets[parsed_filter_targets != "all_usable_if_possible"]
  gene_filter_targets <- as.integer(parsed_filter_targets)
  gene_filter_targets <- gene_filter_targets[!is.na(gene_filter_targets)]
}
powers_to_test <- 1:20
target_signed_r2 <- 0.80
parameter_grid <- expand.grid(
  deepSplit = c(2, 3, 4),
  minModuleSize = c(20, 30, 50),
  mergeCutHeight = c(0.15, 0.25),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
parameter_grid_limit <- as.integer(Sys.getenv("STEP6_PARAM_GRID_LIMIT", unset = NA_character_))
if (!is.na(parameter_grid_limit) && parameter_grid_limit > 0) {
  parameter_grid <- parameter_grid[seq_len(min(parameter_grid_limit, nrow(parameter_grid))), , drop = FALSE]
}
max_regions <- as.integer(Sys.getenv("STEP6_MAX_REGIONS", unset = NA_character_))
if (!is.na(max_regions) && max_regions > 0) {
  expression_files <- expression_files[seq_len(min(max_regions, length(expression_files)))]
}
variant_mode <- Sys.getenv("STEP6_VARIANTS", unset = "all")
if (variant_mode == "unadjusted_and_residualized_if_covariates_available") {
  variant_mode <- "all"
}
highlight_ahba_modules <- c("M26", "M28", "M2", "M6", "M7")

run_control <- list(
  STEP6_MAX_REGIONS = Sys.getenv("STEP6_MAX_REGIONS", unset = ""),
  STEP6_GENE_FILTER_TARGETS = Sys.getenv("STEP6_GENE_FILTER_TARGETS", unset = ""),
  STEP6_PARAM_GRID_LIMIT = Sys.getenv("STEP6_PARAM_GRID_LIMIT", unset = ""),
  STEP6_VARIANTS = Sys.getenv("STEP6_VARIANTS", unset = "")
)
is_limited_run <- any(nzchar(unlist(run_control)))

message("Project root:    ", project_root)
message("Output folder:   ", output_folder)
message("Expression files found: ", length(expression_files))
message("Gene filter targets: ", paste(gene_filter_targets, collapse = ", "))
message("Parameter-grid rows: ", nrow(parameter_grid))
message("Variant mode: ", variant_mode)

required_files <- c(metadata_file, tier2_file, ml_file, ahba_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files:\n", paste(missing_files, collapse = "\n"))
}
if (length(expression_files) == 0) {
  stop("No normalized expression files were found in: ", expression_folder)
}

## -------------------------------------------------------------------------
## 2. Small helper functions
## -------------------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

clean_gene_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  x
}

standardize_region <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("_", " ") %>%
    str_replace_all("-", " ") %>%
    str_squish() %>%
    str_to_lower()
}

infer_region_from_file <- function(path) {
  basename(path) %>%
    str_remove("^normalized_GeTMM_log2_") %>%
    str_remove("\\.csv$") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
}

safe_name <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
}

find_first_column <- function(data, candidates) {
  matched <- candidates[candidates %in% names(data)]
  if (length(matched) > 0) matched[1] else NA_character_
}

read_region_annotation <- function(region) {
  file_stub <- region %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
  path <- file.path(gene_annotation_folder, paste0("gene_annotation_", file_stub, ".csv"))

  if (!file.exists(path)) {
    warning("Gene annotation file was not found for ", region, ": ", path)
    return(tibble(gene_id = character(), gene_symbol = character()))
  }

  readr::read_csv(path, show_col_types = FALSE) %>%
    transmute(
      gene_id = as.character(GeneID),
      gene_symbol = as.character(Symbol)
    ) %>%
    filter(!is.na(gene_id), gene_id != "") %>%
    distinct(gene_id, .keep_all = TRUE)
}

write_csv_safely <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(data, path)
  path
}

is_grey_module <- function(module_name) {
  !is.na(module_name) & module_name == "MEgrey"
}

collapse_semicolon <- function(x, max_items = Inf) {
  x <- unique(na.omit(as.character(x)))
  x <- x[x != ""]
  if (length(x) == 0) return("")
  if (is.finite(max_items)) x <- head(x, max_items)
  paste(sort(x), collapse = ";")
}

ensure_columns <- function(data, columns) {
  for (column in columns) {
    if (!(column %in% names(data))) {
      data[[column]] <- NA
    }
  }
  data
}

make_empty_row <- function(region, expression_file, analysis_variant, gene_filter_version,
                           n_samples, n_genes_input, notes) {
  tibble(
    region = region,
    expression_file = basename(expression_file),
    analysis_variant = analysis_variant,
    gene_filter_version = gene_filter_version,
    n_samples = n_samples,
    n_genes_input = n_genes_input,
    n_genes_after_filter = NA_integer_,
    power = NA_real_,
    power_selection_status = "not_run",
    deepSplit = NA_integer_,
    minModuleSize = NA_integer_,
    mergeCutHeight = NA_real_,
    corType = NA_character_,
    networkType = "signed",
    n_non_grey_modules = NA_integer_,
    largest_module_size = NA_integer_,
    largest_module_fraction = NA_real_,
    grey_gene_count = NA_integer_,
    grey_gene_fraction = NA_real_,
    median_module_size = NA_real_,
    usable_solution = FALSE,
    notes = notes
  )
}

read_expression_matrix <- function(path) {
  message("Reading expression file: ", basename(path))
  raw <- readr::read_csv(path, show_col_types = FALSE)
  if (ncol(raw) < 3) {
    stop("Expression file has too few columns: ", path)
  }

  gene_col <- find_first_column(raw, c("gene_symbol", "Gene", "gene", "symbol", "gene_id"))
  if (is.na(gene_col)) {
    gene_col <- names(raw)[1]
  }

  genes <- clean_gene_id(raw[[gene_col]])
  expr <- raw[, setdiff(names(raw), gene_col), drop = FALSE]
  expr <- as.data.frame(lapply(expr, as.numeric), check.names = FALSE)
  rownames(expr) <- genes
  expr <- expr[!is.na(rownames(expr)) & !duplicated(rownames(expr)), , drop = FALSE]
  expr <- as.matrix(expr)

  dat_expr <- t(expr)
  storage.mode(dat_expr) <- "numeric"
  dat_expr
}

prepare_metadata <- function(path) {
  metadata <- readr::read_csv(path, show_col_types = FALSE)

  sample_col <- find_first_column(
    metadata,
    c("sample_id", "SampleID", "sample", "geo_accession", "GSM", "Sample")
  )
  if (is.na(sample_col)) {
    stop("Could not find a sample ID column in metadata.")
  }

  metadata %>%
    mutate(
      .sample_id = as.character(.data[[sample_col]]),
      .Diagnosis_raw = as.character(Diagnosis %||% phenotype %||% `phenotype:ch1`),
      Diagnosis_clean = case_when(
        str_detect(str_to_lower(.Diagnosis_raw), "mdd|major|depress") ~ "MDD",
        str_detect(str_to_lower(.Diagnosis_raw), "ctrl|control") ~ "Control",
        TRUE ~ .Diagnosis_raw
      ),
      Diagnosis_numeric = case_when(
        Diagnosis_clean == "MDD" ~ 1,
        Diagnosis_clean == "Control" ~ 0,
        TRUE ~ NA_real_
      ),
      .region_standard = standardize_region(AnalysisGroup %||% BrainArea %||% tissue)
    )
}

numeric_or_factor_covariates <- function(metadata, candidate_covariates) {
  available <- candidate_covariates[candidate_covariates %in% names(metadata)]
  usable <- character(0)

  for (covar in available) {
    values <- metadata[[covar]]
    non_missing <- values[!is.na(values)]
    if (length(non_missing) < 3) next
    if (length(unique(non_missing)) < 2) next
    usable <- c(usable, covar)
  }

  usable
}

residualize_expression <- function(dat_expr, sample_metadata) {
  technical_covariates <- numeric_or_factor_covariates(
    sample_metadata,
    c("BrainBank", "SeqBatch", "RIN", "rin", "Age", "age", "Sex", "gender")
  )

  if (length(technical_covariates) == 0) {
    return(list(dat_expr = NULL, notes = "No usable technical/non-diagnosis covariates were found."))
  }

  complete_rows <- complete.cases(sample_metadata[, technical_covariates, drop = FALSE])
  if (sum(complete_rows) < 15) {
    return(list(dat_expr = NULL, notes = "Too few complete samples for residualization."))
  }

  covariate_data <- sample_metadata[complete_rows, technical_covariates, drop = FALSE]
  expression_complete <- dat_expr[complete_rows, , drop = FALSE]

  design <- model.matrix(
    as.formula(paste("~", paste(technical_covariates, collapse = " + "))),
    data = covariate_data
  )
  fit <- lm.fit(design, expression_complete)
  adjusted <- fit$residuals
  rownames(adjusted) <- rownames(expression_complete)
  colnames(adjusted) <- colnames(expression_complete)

  notes <- paste0(
    "Residualized sensitivity branch; Diagnosis was not regressed out; covariates used: ",
    paste(technical_covariates, collapse = ", ")
  )

  metadata_complete <- sample_metadata[complete_rows, , drop = FALSE]
  rownames(metadata_complete) <- rownames(expression_complete)

  list(dat_expr = adjusted, metadata = metadata_complete, notes = notes)
}

filter_genes_for_wgcna <- function(dat_expr, target_size) {
  usable <- dat_expr[, colSums(is.na(dat_expr)) == 0, drop = FALSE]
  gene_mad <- apply(usable, 2, mad, na.rm = TRUE)
  gene_var <- apply(usable, 2, var, na.rm = TRUE)
  keep <- is.finite(gene_mad) & is.finite(gene_var) & gene_var > 0
  usable <- usable[, keep, drop = FALSE]
  gene_mad <- gene_mad[keep]

  gsg <- WGCNA::goodSamplesGenes(usable, verbose = 0)
  usable <- usable[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  gene_mad <- gene_mad[colnames(usable)]

  if (is.na(target_size) || ncol(usable) <= target_size) {
    return(usable)
  }

  selected <- names(sort(gene_mad, decreasing = TRUE))[seq_len(target_size)]
  usable[, selected, drop = FALSE]
}

build_gene_filter_versions <- function(dat_expr) {
  usable_gene_count <- ncol(filter_genes_for_wgcna(dat_expr, NA_integer_))
  versions <- tibble(gene_filter_version = character(), target_size = integer())

  for (target in gene_filter_targets) {
    if (usable_gene_count >= target) {
      versions <- bind_rows(
        versions,
        tibble(gene_filter_version = paste0("top_", target, "_MAD"), target_size = target)
      )
    }
  }

  if (usable_gene_count <= 20000) {
    versions <- bind_rows(
      versions,
      tibble(gene_filter_version = "all_usable_genes", target_size = NA_integer_)
    )
  }

  versions
}

choose_soft_power <- function(dat_expr, cor_type = "bicor") {
  cor_function <- if (cor_type == "bicor") "bicor" else "cor"
  cor_options <- if (cor_type == "bicor") "maxPOutliers = 0.05, use = 'pairwise.complete.obs'" else "use = 'pairwise.complete.obs'"
  fit <- tryCatch(
    WGCNA::pickSoftThreshold(
      dat_expr,
      powerVector = powers_to_test,
      networkType = "signed",
      corFnc = cor_function,
      corOptions = cor_options,
      verbose = 0
    ),
    error = function(e) e
  )

  if (inherits(fit, "error") && cor_type == "bicor") {
    warning("bicor soft-threshold selection failed; falling back to Pearson. Error: ", fit$message)
    return(choose_soft_power(dat_expr, cor_type = "pearson"))
  }
  if (inherits(fit, "error")) {
    stop("Soft-threshold selection failed: ", fit$message)
  }

  indices <- as_tibble(fit$fitIndices)
  signed_col <- intersect(c("SFT.R.sq", "Rsquared.SFT"), names(indices))[1]
  if (is.na(signed_col)) {
    stop("Could not find signed R^2 column in pickSoftThreshold output.")
  }

  indices <- indices %>% mutate(.signed_r2 = .data[[signed_col]])
  reached <- indices %>% filter(.signed_r2 >= target_signed_r2)
  if (nrow(reached) > 0) {
    power <- reached$Power[1]
    status <- "first_power_reached_0.8"
  } else {
    best <- indices %>% arrange(desc(.signed_r2), Power) %>% slice(1)
    power <- best$Power
    status <- "no_power_reached_0.8"
  }

  list(power = power, status = status, fit = indices, corType = cor_type)
}

summarize_modules <- function(colors) {
  color_names <- if (is.character(colors)) colors else WGCNA::labels2colors(colors)
  sizes <- sort(table(color_names), decreasing = TRUE)
  non_grey_sizes <- sizes[names(sizes) != "grey"]
  grey_count <- as.integer(sizes["grey"] %||% 0)
  n_genes <- length(color_names)

  tibble(
    n_non_grey_modules = length(non_grey_sizes),
    largest_module_size = ifelse(length(non_grey_sizes) > 0, max(as.integer(non_grey_sizes)), 0L),
    largest_module_fraction = ifelse(length(non_grey_sizes) > 0, max(as.integer(non_grey_sizes)) / n_genes, 0),
    grey_gene_count = grey_count,
    grey_gene_fraction = grey_count / n_genes,
    median_module_size = ifelse(length(non_grey_sizes) > 0, median(as.integer(non_grey_sizes)), NA_real_)
  )
}

run_blockwise_modules <- function(dat_expr, power, grid_row, cor_type) {
  WGCNA::blockwiseModules(
    dat_expr,
    power = power,
    networkType = "signed",
    TOMType = "signed",
    corType = cor_type,
    maxPOutliers = ifelse(cor_type == "bicor", 0.05, 1),
    deepSplit = grid_row$deepSplit,
    minModuleSize = grid_row$minModuleSize,
    mergeCutHeight = grid_row$mergeCutHeight,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    maxBlockSize = max(5000, ncol(dat_expr) + 100),
    verbose = 3
  )
}

prepare_wgcna_network_once <- function(dat_expr, power, cor_type) {
  cor_function <- if (cor_type == "bicor") "bicor" else "cor"
  cor_options <- if (cor_type == "bicor") "maxPOutliers = 0.05, use = 'pairwise.complete.obs'" else "use = 'pairwise.complete.obs'"

  message("Preparing adjacency/TOM once for this region/variant/filter.")
  adjacency_matrix <- WGCNA::adjacency(
    dat_expr,
    power = power,
    type = "signed",
    corFnc = cor_function,
    corOptions = cor_options
  )
  tom <- WGCNA::TOMsimilarity(adjacency_matrix, TOMType = "signed", verbose = 1)
  rm(adjacency_matrix)
  gc()
  diss_tom <- 1 - tom
  gene_tree <- hclust(as.dist(diss_tom), method = "average")
  list(diss_tom = diss_tom, gene_tree = gene_tree)
}

run_grid_from_prepared_network <- function(dat_expr, prepared_network, grid_row) {
  dynamic_modules <- dynamicTreeCut::cutreeDynamic(
    dendro = prepared_network$gene_tree,
    distM = prepared_network$diss_tom,
    deepSplit = grid_row$deepSplit,
    pamRespectsDendro = FALSE,
    minClusterSize = grid_row$minModuleSize,
    verbose = 0
  )
  dynamic_colors <- WGCNA::labels2colors(dynamic_modules)
  merged <- WGCNA::mergeCloseModules(
    dat_expr,
    dynamic_colors,
    cutHeight = grid_row$mergeCutHeight,
    verbose = 0
  )
  list(
    colors = merged$colors,
    dendrograms = list(prepared_network$gene_tree),
    blockGenes = list(seq_len(ncol(dat_expr)))
  )
}

save_sample_dendrogram <- function(dat_expr, prefix) {
  path <- file.path(plot_folder, paste0(prefix, "_sample_clustering_dendrogram.png"))
  png(path, width = 1200, height = 800)
  sample_tree <- hclust(dist(dat_expr), method = "average")
  plot(sample_tree, main = "Sample clustering", xlab = "", sub = "", cex = 0.7)
  dev.off()
  path
}

save_soft_threshold_plot <- function(fit_indices, prefix) {
  path <- file.path(plot_folder, paste0(prefix, "_soft_threshold_diagnostic.png"))
  png(path, width = 1200, height = 600)
  par(mfrow = c(1, 2))
  plot(
    fit_indices$Power,
    fit_indices$.signed_r2,
    type = "b",
    xlab = "Soft-threshold power",
    ylab = "Signed scale-free topology R^2",
    main = "Scale independence"
  )
  abline(h = target_signed_r2, col = "red", lty = 2)
  mean_col <- intersect(c("mean.k.", "mean.k"), names(fit_indices))[1]
  if (!is.na(mean_col)) {
    plot(
      fit_indices$Power,
      fit_indices[[mean_col]],
      type = "b",
      xlab = "Soft-threshold power",
      ylab = "Mean connectivity",
      main = "Mean connectivity"
    )
  } else {
    plot.new()
    text(0.5, 0.5, "Mean connectivity column unavailable")
  }
  dev.off()
  path
}

save_module_dendrogram <- function(net, prefix) {
  path <- file.path(plot_folder, paste0(prefix, "_module_dendrogram_colors.png"))
  png(path, width = 1300, height = 800)
  colors <- WGCNA::labels2colors(net$colors)
  WGCNA::plotDendroAndColors(
    net$dendrograms[[1]],
    colors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = "Gene dendrogram and module colors"
  )
  dev.off()
  path
}

save_module_size_plot <- function(size_table, prefix) {
  path <- file.path(plot_folder, paste0(prefix, "_module_size_barplot.png"))
  png(path, width = 1200, height = 700)
  barplot(
    height = size_table$n_genes,
    names.arg = size_table$module,
    las = 2,
    col = size_table$module_color,
    main = "Module sizes",
    ylab = "Genes"
  )
  dev.off()
  path
}

save_trait_heatmap <- function(results, prefix) {
  path <- file.path(plot_folder, paste0(prefix, "_module_trait_heatmap.png"))
  png(path, width = 900, height = 1000)
  plot_data <- results %>% arrange(p_diagnosis)
  beta <- matrix(plot_data$beta_diagnosis, ncol = 1)
  rownames(beta) <- plot_data$module
  colnames(beta) <- "Diagnosis beta"
  if (nrow(beta) == 0 || !any(is.finite(beta))) {
    plot.new()
    text(0.5, 0.5, "No finite Diagnosis beta values available for module-trait heatmap.")
  } else {
    text_matrix <- matrix(
      ifelse(is.finite(plot_data$p_diagnosis), signif(plot_data$p_diagnosis, 2), "NA"),
      ncol = 1
    )
    rownames(text_matrix) <- rownames(beta)
    colnames(text_matrix) <- colnames(beta)
    WGCNA::labeledHeatmap(
      Matrix = beta,
      xLabels = colnames(beta),
      yLabels = rownames(beta),
      ySymbols = rownames(beta),
      colorLabels = FALSE,
      colors = WGCNA::blueWhiteRed(50),
      textMatrix = text_matrix,
      main = "Module association with Diagnosis"
    )
  }
  dev.off()
  path
}

save_dotplot <- function(data, x_col, y_col, size_col, color_col, title, path) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || nrow(data) == 0) {
    return(NA_character_)
  }
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])) +
    ggplot2::geom_point(ggplot2::aes(size = .data[[size_col]], color = .data[[color_col]]), alpha = 0.85) +
    ggplot2::scale_color_viridis_c(option = "magma", direction = -1) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = title, x = NULL, y = NULL, size = size_col, color = color_col)
  ggplot2::ggsave(path, p, width = 11, height = 7, dpi = 150)
  path
}

make_module_tables <- function(region, dat_expr, net, gene_annotation) {
  colors <- if (is.character(net$colors)) net$colors else WGCNA::labels2colors(net$colors)
  gene_to_module <- tibble(
    region = region,
    gene_id = colnames(dat_expr),
    module_label_numeric = match(colors, unique(colors)),
    module_color = colors,
    module = paste0("ME", colors)
  ) %>%
    left_join(gene_annotation, by = "gene_id") %>%
    mutate(gene = gene_id, .after = gene_symbol)

  size_table <- gene_to_module %>%
    count(region, module, module_color, name = "n_genes") %>%
    arrange(desc(n_genes))

  list(gene_to_module = gene_to_module, size_table = size_table)
}

fit_module_trait_models <- function(region, module_eigengenes, sample_metadata, size_table) {
  covariates <- numeric_or_factor_covariates(
    sample_metadata,
    c("Age", "age", "Sex", "gender", "BrainBank", "SeqBatch", "RIN", "rin", "PMI", "pmi")
  )

  results <- list()
  for (module_name in colnames(module_eigengenes)) {
    model_data <- sample_metadata %>%
      mutate(.ME = module_eigengenes[, module_name]) %>%
      filter(!is.na(Diagnosis_numeric), !is.na(.ME))

    full_covars <- covariates
    formula_text <- paste(".ME ~ Diagnosis_numeric", paste(full_covars, collapse = " + "), sep = ifelse(length(full_covars) > 0, " + ", ""))
    formula_obj <- as.formula(formula_text)

    model_notes <- "full_or_available_covariate_model"
    min_n_needed <- 5 + length(full_covars) * 2
    if (nrow(model_data) < min_n_needed) {
      formula_text <- ".ME ~ Diagnosis_numeric"
      formula_obj <- as.formula(formula_text)
      model_notes <- "reduced_model_due_to_sample_size"
    }

    fit <- tryCatch(lm(formula_obj, data = model_data), error = function(e) e)
    if (inherits(fit, "error")) {
      results[[module_name]] <- tibble(
        region = region,
        module = module_name,
        module_role = ifelse(is_grey_module(module_name), "unassigned_grey", "coexpression_module"),
        n_module_genes = size_table$n_genes[match(module_name, size_table$module)],
        beta_diagnosis = NA_real_,
        p_diagnosis = NA_real_,
        model_formula = formula_text,
        n_samples = nrow(model_data),
        notes = paste("model_failed:", fit$message)
      )
      next
    }

    coef_table <- summary(fit)$coefficients
    if ("Diagnosis_numeric" %in% rownames(coef_table)) {
      beta <- coef_table["Diagnosis_numeric", "Estimate"]
      pval <- coef_table["Diagnosis_numeric", "Pr(>|t|)"]
    } else {
      beta <- NA_real_
      pval <- NA_real_
      model_notes <- paste(model_notes, "Diagnosis term unavailable", sep = "; ")
    }

    results[[module_name]] <- tibble(
      region = region,
      module = module_name,
      module_role = ifelse(is_grey_module(module_name), "unassigned_grey", "coexpression_module"),
      n_module_genes = size_table$n_genes[match(module_name, size_table$module)],
      beta_diagnosis = beta,
      p_diagnosis = pval,
      model_formula = formula_text,
      n_samples = nrow(model_data),
      notes = model_notes
    )
  }

  bind_rows(results) %>%
    group_by(region) %>%
    mutate(
      fdr_diagnosis_BH = ifelse(
        is_grey_module(module),
        NA_real_,
        p.adjust(ifelse(is_grey_module(module), NA_real_, p_diagnosis), method = "BH")
      ),
      diagnosis_fdr_lt_0_10 = !is.na(fdr_diagnosis_BH) & fdr_diagnosis_BH < 0.10,
      diagnosis_raw_p_lt_0_05 = !is.na(p_diagnosis) & p_diagnosis < 0.05 & !is_grey_module(module)
    ) %>%
    ungroup() %>%
    select(region, module, module_role, n_module_genes, beta_diagnosis, p_diagnosis,
           fdr_diagnosis_BH, diagnosis_fdr_lt_0_10, diagnosis_raw_p_lt_0_05,
           model_formula, n_samples, notes)
}

test_enrichment <- function(region, gene_to_module, tier2_genes, ml_genes) {
  background <- unique(gene_to_module$gene_id)

  gene_sets <- list(
    all_Tier2_genes = unique(as.character(tier2_genes$gene_id)),
    high_confidence_ML_genes = unique(as.character((ml_genes %>% filter(str_detect(ML_priority_class, "High-confidence"))) $gene_id)),
    high_plus_moderate_ML_genes = unique(as.character((ml_genes %>% filter(str_detect(ML_priority_class, "High|Moderate"))) $gene_id))
  )

  region_key <- standardize_region(region)
  if ("region_or_subregion" %in% names(tier2_genes)) {
    same_region_tier2 <- tier2_genes %>%
      filter(standardize_region(region_or_subregion) == region_key) %>%
      pull(gene_id)
    gene_sets$same_region_Tier2_genes <- unique(same_region_tier2)
  }
  if ("region_or_subregion" %in% names(ml_genes)) {
    same_region_ml <- ml_genes %>%
      filter(
        standardize_region(region_or_subregion) == region_key,
        str_detect(ML_priority_class, "High|Moderate")
      ) %>%
      pull(gene_id)
    gene_sets$same_region_high_plus_moderate_ML_genes <- unique(same_region_ml)
  }

  gene_sets <- lapply(gene_sets, function(x) intersect(clean_gene_id(x), background))

  rows <- list()
  for (module_name in unique(gene_to_module$module)) {
    module_genes <- gene_to_module$gene_id[gene_to_module$module == module_name]
    for (set_name in names(gene_sets)) {
      set_genes <- gene_sets[[set_name]]
      overlap <- intersect(module_genes, set_genes)
      contingency <- matrix(
        c(
          length(overlap),
          length(module_genes) - length(overlap),
          length(set_genes) - length(overlap),
          length(background) - length(module_genes) - length(set_genes) + length(overlap)
        ),
        nrow = 2
      )
      pval <- tryCatch(fisher.test(contingency, alternative = "greater")$p.value, error = function(e) NA_real_)
      rows[[length(rows) + 1]] <- tibble(
        region = region,
        module = module_name,
        module_role = ifelse(is_grey_module(module_name), "unassigned_grey", "coexpression_module"),
        gene_set = set_name,
        gene_set_family = case_when(
          set_name == "all_Tier2_genes" ~ "all_Tier2_across_all_regions",
          set_name == "same_region_Tier2_genes" ~ "same_region_Tier2",
          set_name == "high_confidence_ML_genes" ~ "high_confidence_ML",
          set_name == "high_plus_moderate_ML_genes" ~ "high_plus_moderate_ML",
          set_name == "same_region_high_plus_moderate_ML_genes" ~ "same_region_high_plus_moderate_ML",
          TRUE ~ set_name
        ),
        module_gene_count = length(module_genes),
        gene_set_count_in_background = length(set_genes),
        overlap_count = length(overlap),
        overlap_genes = paste(sort(overlap), collapse = ";"),
        p_value = pval
      )
    }
  }

  bind_rows(rows) %>%
    group_by(region, gene_set) %>%
    mutate(
      fdr_BH = ifelse(
        is_grey_module(module),
        NA_real_,
        p.adjust(ifelse(is_grey_module(module), NA_real_, p_value), method = "BH")
      ),
      enrichment_fdr_lt_0_10 = !is.na(fdr_BH) & fdr_BH < 0.10,
      enrichment_raw_p_lt_0_05 = !is.na(p_value) & p_value < 0.05 & !is_grey_module(module)
    ) %>%
    ungroup()
}

test_ahba_overlap <- function(region, gene_to_module, ahba_modules) {
  mdd_with_symbols <- gene_to_module %>% filter(!is.na(gene_symbol), gene_symbol != "")
  background <- intersect(unique(mdd_with_symbols$gene_symbol), unique(ahba_modules$gene_symbol))
  ahba_modules <- ahba_modules %>% filter(gene_symbol %in% background)
  mdd_modules <- mdd_with_symbols %>% filter(gene_symbol %in% background)

  rows <- list()
  for (mdd_module in unique(mdd_modules$module)) {
    mdd_genes <- mdd_modules$gene_symbol[mdd_modules$module == mdd_module]
    for (ahba_module in sort(unique(ahba_modules$published_module))) {
      ahba_genes <- ahba_modules$gene_symbol[ahba_modules$published_module == ahba_module]
      overlap <- intersect(mdd_genes, ahba_genes)
      contingency <- matrix(
        c(
          length(overlap),
          length(mdd_genes) - length(overlap),
          length(ahba_genes) - length(overlap),
          length(background) - length(mdd_genes) - length(ahba_genes) + length(overlap)
        ),
        nrow = 2
      )
      pval <- tryCatch(fisher.test(contingency, alternative = "greater")$p.value, error = function(e) NA_real_)
      rows[[length(rows) + 1]] <- tibble(
        region = region,
        mdd_module = mdd_module,
        module_role = ifelse(is_grey_module(mdd_module), "unassigned_grey", "coexpression_module"),
        ahba_module = ahba_module,
        highlighted_ahba_module = ahba_module %in% highlight_ahba_modules,
        mdd_module_gene_count = length(mdd_genes),
        ahba_module_gene_count = length(ahba_genes),
        overlap_count = length(overlap),
        overlap_genes = paste(sort(overlap), collapse = ";"),
        p_value = pval
      )
    }
  }

  bind_rows(rows) %>%
    group_by(region) %>%
    mutate(
      fdr_BH = ifelse(
        is_grey_module(mdd_module),
        NA_real_,
        p.adjust(ifelse(is_grey_module(mdd_module), NA_real_, p_value), method = "BH")
      ),
      ahba_highlight_fdr_lt_0_10 = highlighted_ahba_module & !is.na(fdr_BH) & fdr_BH < 0.10,
      ahba_highlight_raw_p_lt_0_05 = highlighted_ahba_module & !is.na(p_value) & p_value < 0.05 & !is_grey_module(mdd_module)
    ) %>%
    ungroup()
}

calculate_hub_genes <- function(region, dat_expr, sample_metadata, gene_to_module, module_eigengenes) {
  diagnosis <- sample_metadata$Diagnosis_numeric
  rows <- list()

  for (module_name in unique(gene_to_module$module)) {
    if (module_name == "MEgrey") next
    module_genes <- gene_to_module$gene_id[gene_to_module$module == module_name]
    module_genes <- intersect(module_genes, colnames(dat_expr))
    if (length(module_genes) == 0 || !(module_name %in% colnames(module_eigengenes))) next

    kme <- cor(dat_expr[, module_genes, drop = FALSE], module_eigengenes[, module_name], use = "pairwise.complete.obs")
    gs <- rep(NA_real_, length(module_genes))
    if (sum(!is.na(diagnosis)) >= 10 && length(unique(na.omit(diagnosis))) == 2) {
      gs <- apply(dat_expr[, module_genes, drop = FALSE], 2, function(x) {
        suppressWarnings(cor(x, diagnosis, use = "pairwise.complete.obs"))
      })
    }

    rows[[module_name]] <- tibble(
      region = region,
      module = module_name,
      gene_id = module_genes,
      gene_symbol = gene_to_module$gene_symbol[match(module_genes, gene_to_module$gene_id)],
      kME = as.numeric(kme[, 1]),
      abs_kME = abs(as.numeric(kme[, 1])),
      gene_significance_diagnosis = as.numeric(gs),
      abs_gene_significance_diagnosis = abs(as.numeric(gs))
    ) %>%
      arrange(desc(abs_kME), desc(abs_gene_significance_diagnosis)) %>%
      slice_head(n = 50)
  }

  bind_rows(rows)
}

make_important_module_summary <- function(trait_results, enrichment_results, ahba_results, hub_genes) {
  if (nrow(trait_results) == 0) {
    return(tibble())
  }

  enrichment_flags <- enrichment_results %>%
    filter(!is_grey_module(module)) %>%
    group_by(region, module) %>%
    summarise(
      enriched_for_Tier2 = any(gene_set_family %in% c("all_Tier2_across_all_regions", "same_region_Tier2") & enrichment_fdr_lt_0_10, na.rm = TRUE),
      enriched_for_Tier2_exploratory = any(gene_set_family %in% c("all_Tier2_across_all_regions", "same_region_Tier2") & enrichment_raw_p_lt_0_05, na.rm = TRUE),
      enriched_for_high_moderate_ML = any(gene_set_family %in% c("high_plus_moderate_ML", "same_region_high_plus_moderate_ML") & enrichment_fdr_lt_0_10, na.rm = TRUE),
      enriched_for_high_moderate_ML_exploratory = any(gene_set_family %in% c("high_plus_moderate_ML", "same_region_high_plus_moderate_ML") & enrichment_raw_p_lt_0_05, na.rm = TRUE),
      enriched_for_high_confidence_ML = any(gene_set_family == "high_confidence_ML" & enrichment_fdr_lt_0_10, na.rm = TRUE),
      best_Tier2_or_ML_enrichment = collapse_semicolon(
        paste0(gene_set, " overlap=", overlap_count, " FDR=", signif(fdr_BH, 3))[
          !is.na(fdr_BH) & fdr_BH < 0.10
        ],
        max_items = 4
      ),
      exploratory_Tier2_or_ML_enrichment = collapse_semicolon(
        paste0(gene_set, " overlap=", overlap_count, " p=", signif(p_value, 3))[
          (is.na(fdr_BH) | fdr_BH >= 0.10) & enrichment_raw_p_lt_0_05
        ],
        max_items = 4
      ),
      .groups = "drop"
    )

  ahba_flags <- ahba_results %>%
    filter(!is_grey_module(mdd_module), highlighted_ahba_module) %>%
    group_by(region, module = mdd_module) %>%
    summarise(
      AHBA_overlap_modules = collapse_semicolon(
        paste0(ahba_module, " overlap=", overlap_count, " FDR=", signif(fdr_BH, 3))[
          ahba_highlight_fdr_lt_0_10
        ],
        max_items = 5
      ),
      AHBA_overlap_modules_exploratory = collapse_semicolon(
        paste0(ahba_module, " overlap=", overlap_count, " p=", signif(p_value, 3))[
          ahba_highlight_raw_p_lt_0_05 & !ahba_highlight_fdr_lt_0_10
        ],
        max_items = 5
      ),
      has_AHBA_highlight_overlap = any(ahba_highlight_fdr_lt_0_10, na.rm = TRUE),
      has_AHBA_highlight_overlap_exploratory = any(ahba_highlight_raw_p_lt_0_05, na.rm = TRUE),
      .groups = "drop"
    )

  hub_summary <- hub_genes %>%
    group_by(region, module) %>%
    arrange(desc(abs_kME), .by_group = TRUE) %>%
    summarise(
      top_hub_genes = collapse_semicolon(ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, gene_id), max_items = 10),
      .groups = "drop"
    )

  trait_results %>%
    left_join(enrichment_flags, by = c("region", "module")) %>%
    left_join(ahba_flags, by = c("region", "module")) %>%
    left_join(hub_summary, by = c("region", "module")) %>%
    mutate(
      across(
        c(enriched_for_Tier2, enriched_for_Tier2_exploratory,
          enriched_for_high_moderate_ML, enriched_for_high_moderate_ML_exploratory,
          enriched_for_high_confidence_ML, has_AHBA_highlight_overlap,
          has_AHBA_highlight_overlap_exploratory),
        ~ ifelse(is.na(.x), FALSE, .x)
      ),
      AHBA_overlap_modules = ifelse(is.na(AHBA_overlap_modules), "", AHBA_overlap_modules),
      AHBA_overlap_module_exploratory = ifelse(is.na(AHBA_overlap_modules_exploratory), "", AHBA_overlap_modules_exploratory),
      top_hub_genes = ifelse(is.na(top_hub_genes), "", top_hub_genes),
      interpretation = case_when(
        is_grey_module(module) ~ "Unassigned grey genes; retained in raw files but excluded from biological interpretation.",
        diagnosis_fdr_lt_0_10 & (enriched_for_Tier2 | enriched_for_high_moderate_ML | has_AHBA_highlight_overlap) ~ "Diagnosis-associated module with candidate-gene or AHBA support.",
        diagnosis_fdr_lt_0_10 ~ "Diagnosis-associated module eigengene.",
        diagnosis_raw_p_lt_0_05 & (enriched_for_Tier2 | enriched_for_high_moderate_ML | has_AHBA_highlight_overlap | enriched_for_Tier2_exploratory | enriched_for_high_moderate_ML_exploratory | has_AHBA_highlight_overlap_exploratory) ~ "Exploratory diagnosis signal with candidate-gene or AHBA support.",
        enriched_for_Tier2 | enriched_for_high_moderate_ML | enriched_for_high_confidence_ML | has_AHBA_highlight_overlap ~ "Candidate-gene or AHBA-supported coexpression module without corrected Diagnosis association.",
        diagnosis_raw_p_lt_0_05 ~ "Exploratory raw-p Diagnosis association only.",
        enriched_for_Tier2_exploratory | enriched_for_high_moderate_ML_exploratory | has_AHBA_highlight_overlap_exploratory ~ "Nominal candidate-gene or AHBA context only.",
        TRUE ~ "No highlighted Diagnosis, candidate-gene, or AHBA signal."
      ),
      caution_level = case_when(
        is_grey_module(module) ~ "unassigned_do_not_interpret",
        diagnosis_fdr_lt_0_10 & (enriched_for_Tier2 | enriched_for_high_moderate_ML | has_AHBA_highlight_overlap) ~ "strongest_exploratory_module_level_signal",
        diagnosis_fdr_lt_0_10 ~ "diagnosis_fdr_signal_no_enrichment_context",
        diagnosis_raw_p_lt_0_05 | enriched_for_Tier2 | enriched_for_high_moderate_ML | enriched_for_high_confidence_ML | has_AHBA_highlight_overlap ~ "exploratory_follow_up",
        enriched_for_Tier2_exploratory | enriched_for_high_moderate_ML_exploratory | has_AHBA_highlight_overlap_exploratory ~ "nominal_only_context",
        TRUE ~ "not_highlighted"
      ),
      important_module = !is_grey_module(module) & (
        diagnosis_fdr_lt_0_10 |
          diagnosis_raw_p_lt_0_05 |
          enriched_for_Tier2 |
          enriched_for_high_moderate_ML |
          enriched_for_high_confidence_ML |
          has_AHBA_highlight_overlap |
          enriched_for_Tier2_exploratory |
          enriched_for_high_moderate_ML_exploratory |
          has_AHBA_highlight_overlap_exploratory
      )
    ) %>%
    filter(important_module | is_grey_module(module)) %>%
    select(
      region,
      module,
      module_role,
      module_size = n_module_genes,
      diagnosis_beta = beta_diagnosis,
      diagnosis_p_value = p_diagnosis,
      diagnosis_BH_FDR = fdr_diagnosis_BH,
      diagnosis_fdr_lt_0_10,
      diagnosis_raw_p_lt_0_05,
      enriched_for_Tier2,
      enriched_for_high_moderate_ML,
      AHBA_overlap_modules,
      AHBA_overlap_module_exploratory,
      top_hub_genes,
      interpretation,
      caution_level,
      important_module
    ) %>%
    arrange(region, is_grey_module(module), caution_level, diagnosis_BH_FDR, diagnosis_p_value)
}

write_report <- function(summary_table, trait_results, enrichment_results, ahba_results, important_modules, output_path, full_report = FALSE) {
  trait_results <- ensure_columns(
    trait_results,
    c("region", "module", "module_role", "n_module_genes", "beta_diagnosis",
      "p_diagnosis", "fdr_diagnosis_BH", "diagnosis_fdr_lt_0_10",
      "diagnosis_raw_p_lt_0_05")
  )
  enrichment_results <- ensure_columns(
    enrichment_results,
    c("region", "module", "gene_set", "gene_set_family", "overlap_count",
      "p_value", "fdr_BH", "enrichment_fdr_lt_0_10", "enrichment_raw_p_lt_0_05")
  )
  ahba_results <- ensure_columns(
    ahba_results,
    c("region", "mdd_module", "ahba_module", "highlighted_ahba_module",
      "overlap_count", "p_value", "fdr_BH", "ahba_highlight_fdr_lt_0_10",
      "ahba_highlight_raw_p_lt_0_05")
  )
  important_modules <- ensure_columns(
    important_modules,
    c("region", "module", "important_module", "interpretation", "top_hub_genes")
  )

  usable <- summary_table %>% filter(usable_solution)
  best_regions <- unique(summary_table$region)
  diagnosis_hits <- trait_results %>% filter(!is_grey_module(module), !is.na(fdr_diagnosis_BH), fdr_diagnosis_BH < 0.10)
  diagnosis_nominal <- trait_results %>% filter(!is_grey_module(module), diagnosis_raw_p_lt_0_05, !diagnosis_fdr_lt_0_10)
  tier_hits <- enrichment_results %>% filter(!is_grey_module(module), !is.na(fdr_BH), fdr_BH < 0.10)
  ahba_hits <- ahba_results %>% filter(!is_grey_module(mdd_module), highlighted_ahba_module, !is.na(fdr_BH), fdr_BH < 0.10)
  run_scope <- if (is_limited_run) {
    paste0(
      "Limited validation run. Active STEP6 controls: ",
      paste(names(run_control)[nzchar(unlist(run_control))], unlist(run_control)[nzchar(unlist(run_control))], sep = "=", collapse = "; ")
    )
  } else {
    "Full default run. No STEP6 limiting environment variables were set; all available normalized region matrices, default gene filters, full parameter grid, and unadjusted plus residualized-if-available variants were attempted."
  }

  lines <- c(
    ifelse(full_report, "# MDD WGCNA Module Report FULL", "# MDD WGCNA Module Report"),
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    paste0("Run scope: ", run_scope),
    "",
    "## Limited Validation vs Full Run",
    "",
    "The earlier Frontal Cortex-only result was a smoke test run with STEP6_MAX_REGIONS=1, STEP6_GENE_FILTER_TARGETS=5000, STEP6_PARAM_GRID_LIMIT=1, and STEP6_VARIANTS=unadjusted_only. That limited run checked that the workflow could execute and should not be interpreted as the final biological result. This full report is the intended biological summary after running without those limiting controls.",
    "",
    "## What was done",
    "",
    "This analysis built WGCNA coexpression modules directly from the normalized MDD/control GeTMM-log2 expression matrices.",
    "",
    "The modules were not built from DEG-only genes, Tier 2-only genes, or ML-prioritized-only genes. Genes were filtered by expression quality and variability/MAD, then existing Tier 2 and ML-prioritized gene sets were tested for enrichment after module construction.",
    "",
    "Signed WGCNA networks were used as the main analysis. The script tried robust bicor correlation first and used Pearson only if bicor failed.",
    "",
    "MEgrey is the WGCNA unassigned module. MEgrey rows are retained in raw output files for auditability, but MEgrey is excluded from the main biological interpretation tables.",
    "",
    "## Sample-size and usability overview",
    "",
    paste0("- Regions processed: ", paste(sort(best_regions), collapse = ", ")),
    paste0("- Parameter-grid runs recorded: ", nrow(summary_table)),
    paste0("- Usable parameter-grid runs: ", sum(summary_table$usable_solution, na.rm = TRUE)),
    paste0("- Regions with at least one usable module solution: ", paste(sort(unique(usable$region)), collapse = ", ")),
    "",
    "A solution was marked usable when it had at least 3 non-grey modules, the largest module was below 70% of genes, grey genes were below 50%, and the region had at least 15 samples.",
    "",
    "## Diagnosis-associated modules",
    "",
    if (nrow(diagnosis_hits) == 0) {
      "No module eigengene reached BH-FDR < 0.10 for Diagnosis in the saved best solutions."
    } else {
      paste0(
        "- ", diagnosis_hits$region, " ", diagnosis_hits$module,
        ": beta=", signif(diagnosis_hits$beta_diagnosis, 3),
        ", p=", signif(diagnosis_hits$p_diagnosis, 3),
        ", FDR=", signif(diagnosis_hits$fdr_diagnosis_BH, 3)
      )
    },
    "",
    "Exploratory raw-p Diagnosis modules:",
    "",
    if (nrow(diagnosis_nominal) == 0) {
      "No non-grey module eigengene had raw p < 0.05 without also reaching the FDR threshold."
    } else {
      paste0(
        "- ", diagnosis_nominal$region, " ", diagnosis_nominal$module,
        ": beta=", signif(diagnosis_nominal$beta_diagnosis, 3),
        ", p=", signif(diagnosis_nominal$p_diagnosis, 3),
        ", FDR=", signif(diagnosis_nominal$fdr_diagnosis_BH, 3)
      )
    },
    "",
    "## Tier 2 and ML-prioritized enrichment",
    "",
    if (nrow(tier_hits) == 0) {
      "No Tier 2 or ML-prioritized gene-set enrichment reached BH-FDR < 0.10 in the saved best solutions."
    } else {
      paste0(
        "- ", tier_hits$region, " ", tier_hits$module, " / ", tier_hits$gene_set,
        ": overlap=", tier_hits$overlap_count,
        ", p=", signif(tier_hits$p_value, 3),
        ", FDR=", signif(tier_hits$fdr_BH, 3)
      )
    },
    "",
    "## Published AHBA module overlap",
    "",
    "The AHBA overlap step tests each MDD-derived module against published AHBA modules and highlights M26, M28, M2, M6, and M7.",
    "",
    if (nrow(ahba_hits) == 0) {
      "No highlighted AHBA module overlap reached BH-FDR < 0.10 in the saved best solutions."
    } else {
      paste0(
        "- ", ahba_hits$region, " ", ahba_hits$mdd_module, " vs ", ahba_hits$ahba_module,
        ": overlap=", ahba_hits$overlap_count,
        ", p=", signif(ahba_hits$p_value, 3),
        ", FDR=", signif(ahba_hits$fdr_BH, 3)
      )
    },
    "",
    "## Important Module Summary",
    "",
    if (nrow(important_modules %>% filter(important_module)) == 0) {
      "No non-grey important modules were highlighted by Diagnosis, Tier2/ML enrichment, or highlighted AHBA overlap criteria."
    } else {
      paste0(
        "- ", (important_modules %>% filter(important_module))$region,
        " ", (important_modules %>% filter(important_module))$module,
        ": ", (important_modules %>% filter(important_module))$interpretation,
        " Top hubs: ", (important_modules %>% filter(important_module))$top_hub_genes
      )
    },
    "",
    "## Limitations",
    "",
    "- This is exploratory and should be interpreted as module discovery, not validation.",
    "- Some regions have modest sample sizes, so module stability may be limited.",
    "- Residualized analyses are sensitivity analyses and may use Age/Sex if those covariates are available.",
    "- No external validation cohort has been used yet.",
    "",
    "## Main output files",
    "",
    "- `mdd_wgcna_parameter_sensitivity.csv`",
    "- `best_solution_gene_to_module_assignments.csv`",
    "- `best_solution_module_size_table.csv`",
    "- `best_solution_module_eigengenes.csv`",
    "- `module_trait_association_results.csv`",
    "- `module_enrichment_for_Tier2_and_ML_genes.csv`",
    "- `mdd_modules_vs_published_AHBA_modules_enrichment.csv`",
    "- `module_hub_genes_top50.csv`",
    "- `important_modules_summary.csv`",
    "- `mdd_wgcna_manifest.json`"
  )

  writeLines(lines, output_path)
  output_path
}

## -------------------------------------------------------------------------
## 3. Load shared inputs
## -------------------------------------------------------------------------

metadata <- prepare_metadata(metadata_file)
tier2_genes <- readr::read_csv(tier2_file, show_col_types = FALSE)
ml_genes <- readr::read_csv(ml_file, show_col_types = FALSE)
ahba_modules <- readr::read_csv(ahba_file, show_col_types = FALSE)

message("Metadata samples: ", nrow(metadata))
message("Tier 2 rows:      ", nrow(tier2_genes))
message("ML rows:          ", nrow(ml_genes))
message("AHBA rows:        ", nrow(ahba_modules))

all_parameter_rows <- list()
best_gene_to_module <- list()
best_size_tables <- list()
best_eigengenes <- list()
best_metadata <- list()
all_trait_results <- list()
all_enrichment_results <- list()
all_ahba_results <- list()
all_hub_genes <- list()
plot_files <- character(0)

record_parameter_row <- function(row) {
  all_parameter_rows[[length(all_parameter_rows) + 1]] <<- row
  write_csv_safely(
    bind_rows(all_parameter_rows),
    file.path(output_folder, "mdd_wgcna_parameter_sensitivity.csv")
  )
  invisible(row)
}

## -------------------------------------------------------------------------
## 4. Process each region
## -------------------------------------------------------------------------

for (expression_file in expression_files) {
  region <- infer_region_from_file(expression_file)
  region_safe <- safe_name(region)
  gene_annotation <- read_region_annotation(region)
  message("\n==============================")
  message("Region: ", region)
  message("==============================")

  dat_expr_all <- read_expression_matrix(expression_file)
  n_genes_input <- ncol(dat_expr_all)

  sample_ids <- rownames(dat_expr_all)
  metadata_region <- metadata %>% filter(.sample_id %in% sample_ids)
  matched_ids <- intersect(sample_ids, metadata_region$.sample_id)
  dat_expr_all <- dat_expr_all[matched_ids, , drop = FALSE]
  metadata_region <- metadata_region[match(matched_ids, metadata_region$.sample_id), , drop = FALSE]
  rownames(metadata_region) <- metadata_region$.sample_id

  n_samples <- nrow(dat_expr_all)
  sample_note <- if (n_samples < 15) {
    "Skipped: fewer than 15 matched samples."
  } else if (n_samples < 20) {
    "Low-sample-size exploratory run: 15-19 matched samples."
  } else {
    "Normal sample-size run: at least 20 matched samples."
  }
  message("Matched samples: ", n_samples, ". ", sample_note)

  if (n_samples < 15) {
    record_parameter_row(make_empty_row(
      region, expression_file, "unadjusted_logGeTMM", "not_run",
      n_samples, n_genes_input, sample_note
    ))
    next
  }

  plot_files <- c(plot_files, save_sample_dendrogram(dat_expr_all, paste(region_safe, "unadjusted", sep = "_")))

  analysis_variants <- list(
    unadjusted_logGeTMM = list(
      dat_expr = dat_expr_all,
      metadata = metadata_region,
      notes = sample_note
    )
  )

  if (variant_mode != "unadjusted_only") {
    residualized <- residualize_expression(dat_expr_all, metadata_region)
    if (!is.null(residualized$dat_expr)) {
      analysis_variants$residualized_technical_covariates <- residualized
    } else {
      warning("Skipping residualized branch for ", region, ": ", residualized$notes)
    }
  }

  region_solution_records <- list()

  for (analysis_variant in names(analysis_variants)) {
    variant <- analysis_variants[[analysis_variant]]
    dat_variant <- variant$dat_expr
    metadata_variant <- variant$metadata
    variant_note <- variant$notes

    filter_versions <- build_gene_filter_versions(dat_variant)
    if (nrow(filter_versions) == 0) {
      record_parameter_row(make_empty_row(
        region, expression_file, analysis_variant, "no_filter_version_available",
        nrow(dat_variant), n_genes_input, "No usable gene filter version was available."
      ))
      next
    }

    for (filter_i in seq_len(nrow(filter_versions))) {
      filter_version <- filter_versions$gene_filter_version[filter_i]
      target_size <- filter_versions$target_size[filter_i]
      message("Variant: ", analysis_variant, " | Filter: ", filter_version)

      dat_filtered <- filter_genes_for_wgcna(dat_variant, target_size)
      metadata_filtered <- metadata_variant[rownames(dat_filtered), , drop = FALSE]
      n_genes_after_filter <- ncol(dat_filtered)

      if (nrow(dat_filtered) < 15 || n_genes_after_filter < 100) {
        record_parameter_row(make_empty_row(
          region, expression_file, analysis_variant, filter_version,
          nrow(dat_filtered), n_genes_input,
          paste("Skipped after filtering; samples:", nrow(dat_filtered), "genes:", n_genes_after_filter)
        ))
        next
      }

      soft <- tryCatch(choose_soft_power(dat_filtered, "bicor"), error = function(e) e)
      if (inherits(soft, "error")) {
        record_parameter_row(make_empty_row(
          region, expression_file, analysis_variant, filter_version,
          nrow(dat_filtered), n_genes_input,
          paste("Soft-threshold selection failed:", soft$message)
        ))
        next
      }

      prefix <- paste(region_safe, safe_name(analysis_variant), filter_version, sep = "_")
      plot_files <- c(plot_files, save_soft_threshold_plot(soft$fit, prefix))
      prepared_network <- tryCatch(
        prepare_wgcna_network_once(dat_filtered, soft$power, soft$corType),
        error = function(e) e
      )
      if (inherits(prepared_network, "error") && soft$corType == "bicor") {
        warning("bicor adjacency/TOM preparation failed; retrying this filter with Pearson. Error: ", prepared_network$message)
        soft_pearson <- tryCatch(choose_soft_power(dat_filtered, "pearson"), error = function(e) e)
        if (!inherits(soft_pearson, "error")) {
          soft <- soft_pearson
          prepared_network <- tryCatch(
            prepare_wgcna_network_once(dat_filtered, soft$power, soft$corType),
            error = function(e) e
          )
        }
      }
      if (inherits(prepared_network, "error")) {
        record_parameter_row(make_empty_row(
          region, expression_file, analysis_variant, filter_version,
          nrow(dat_filtered), n_genes_input,
          paste("Adjacency/TOM preparation failed:", prepared_network$message)
        ) %>%
          mutate(
            n_genes_after_filter = n_genes_after_filter,
            power = soft$power,
            power_selection_status = soft$status,
            corType = soft$corType
          ))
        next
      }

      for (grid_i in seq_len(nrow(parameter_grid))) {
        grid_row <- parameter_grid[grid_i, ]
        message(
          "Running dynamic tree cut grid: power=", soft$power,
          ", deepSplit=", grid_row$deepSplit,
          ", minModuleSize=", grid_row$minModuleSize,
          ", mergeCutHeight=", grid_row$mergeCutHeight
        )

        net <- tryCatch(
          run_grid_from_prepared_network(dat_filtered, prepared_network, grid_row),
          error = function(e) e
        )

        if (inherits(net, "error")) {
          row <- make_empty_row(
            region, expression_file, analysis_variant, filter_version,
            nrow(dat_filtered), n_genes_input,
            paste("dynamic tree cut grid failed:", net$message)
          ) %>%
            mutate(
              n_genes_after_filter = n_genes_after_filter,
              power = soft$power,
              power_selection_status = soft$status,
              deepSplit = grid_row$deepSplit,
              minModuleSize = grid_row$minModuleSize,
              mergeCutHeight = grid_row$mergeCutHeight,
              corType = soft$corType
            )
          record_parameter_row(row)
          next
        }

        module_summary <- summarize_modules(net$colors)
        usable <- with(
          module_summary,
          n_non_grey_modules >= 3 &&
            largest_module_fraction < 0.70 &&
            grey_gene_fraction < 0.50 &&
            nrow(dat_filtered) >= 15
        )

        row <- tibble(
          region = region,
          expression_file = basename(expression_file),
          analysis_variant = analysis_variant,
          gene_filter_version = filter_version,
          n_samples = nrow(dat_filtered),
          n_genes_input = n_genes_input,
          n_genes_after_filter = n_genes_after_filter,
          power = soft$power,
          power_selection_status = soft$status,
          deepSplit = grid_row$deepSplit,
          minModuleSize = grid_row$minModuleSize,
          mergeCutHeight = grid_row$mergeCutHeight,
          corType = soft$corType,
          networkType = "signed",
          n_non_grey_modules = module_summary$n_non_grey_modules,
          largest_module_size = module_summary$largest_module_size,
          largest_module_fraction = module_summary$largest_module_fraction,
          grey_gene_count = module_summary$grey_gene_count,
          grey_gene_fraction = module_summary$grey_gene_fraction,
          median_module_size = module_summary$median_module_size,
          usable_solution = usable,
          notes = paste(sample_note, variant_note, sep = " | ")
        )

        record_parameter_row(row)
        region_solution_records[[length(region_solution_records) + 1]] <- list(
          summary = row,
          dat_expr = dat_filtered,
          metadata = metadata_filtered,
          net = net,
          soft = soft,
          prefix = paste(prefix, "ds", grid_row$deepSplit, "min", grid_row$minModuleSize, "merge", grid_row$mergeCutHeight, sep = "_")
        )
      }
      rm(prepared_network)
      gc()
    }
  }

  if (length(region_solution_records) == 0) {
    warning("No WGCNA solutions were created for region: ", region)
    next
  }

  ranking_table <- bind_rows(lapply(region_solution_records, `[[`, "summary")) %>%
    mutate(
      usable_rank = ifelse(usable_solution, 1, 0),
      filter_rank = case_when(
        gene_filter_version == "all_usable_genes" ~ n_genes_after_filter,
        TRUE ~ n_genes_after_filter
      )
    ) %>%
    arrange(desc(usable_rank), desc(filter_rank), desc(n_non_grey_modules), largest_module_fraction)

  best_index <- which(vapply(region_solution_records, function(x) {
    identical(x$summary$analysis_variant, ranking_table$analysis_variant[1]) &&
      identical(x$summary$gene_filter_version, ranking_table$gene_filter_version[1]) &&
      identical(x$summary$deepSplit, ranking_table$deepSplit[1]) &&
      identical(x$summary$minModuleSize, ranking_table$minModuleSize[1]) &&
      identical(x$summary$mergeCutHeight, ranking_table$mergeCutHeight[1])
  }, logical(1)))[1]

  best <- region_solution_records[[best_index]]
  best_prefix <- paste0(region_safe, "_best_solution")
  message("Best solution selected for ", region, ": ", best$summary$analysis_variant, " / ", best$summary$gene_filter_version)

  module_tables <- make_module_tables(region, best$dat_expr, best$net, gene_annotation)
  module_eigengenes <- WGCNA::moduleEigengenes(best$dat_expr, colors = module_tables$gene_to_module$module_color)$eigengenes
  rownames(module_eigengenes) <- rownames(best$dat_expr)
  module_eigengenes <- WGCNA::orderMEs(module_eigengenes)
  rownames(module_eigengenes) <- rownames(best$dat_expr)

  trait_results <- fit_module_trait_models(region, module_eigengenes, best$metadata, module_tables$size_table)
  enrichment_results <- test_enrichment(region, module_tables$gene_to_module, tier2_genes, ml_genes)
  ahba_results <- test_ahba_overlap(region, module_tables$gene_to_module, ahba_modules)
  hub_genes <- calculate_hub_genes(region, best$dat_expr, best$metadata, module_tables$gene_to_module, module_eigengenes)

  plot_files <- c(
    plot_files,
    save_module_dendrogram(best$net, best_prefix),
    save_module_size_plot(module_tables$size_table, best_prefix),
    save_trait_heatmap(trait_results, best_prefix),
    save_dotplot(
      enrichment_results %>% filter(overlap_count > 0),
      "gene_set", "module", "overlap_count", "fdr_BH",
      paste(region, "Tier2/ML enrichment"),
      file.path(plot_folder, paste0(best_prefix, "_tier2_ml_enrichment_dotplot.png"))
    ),
    save_dotplot(
      ahba_results %>% filter(highlighted_ahba_module, overlap_count > 0),
      "ahba_module", "mdd_module", "overlap_count", "fdr_BH",
      paste(region, "highlighted AHBA overlap"),
      file.path(plot_folder, paste0(best_prefix, "_ahba_highlight_overlap_dotplot.png"))
    )
  )

  best_gene_to_module[[region]] <- module_tables$gene_to_module
  best_size_tables[[region]] <- module_tables$size_table
  best_eigengenes[[region]] <- as_tibble(module_eigengenes, rownames = "SampleID") %>%
    mutate(region = region, .before = 1)
  best_metadata[[region]] <- best$metadata %>%
    as_tibble() %>%
    mutate(region = region, .before = 1)
  all_trait_results[[region]] <- trait_results
  all_enrichment_results[[region]] <- enrichment_results
  all_ahba_results[[region]] <- ahba_results
  all_hub_genes[[region]] <- hub_genes
}

## -------------------------------------------------------------------------
## 5. Save combined outputs
## -------------------------------------------------------------------------

parameter_summary <- bind_rows(all_parameter_rows)
trait_results <- bind_rows(all_trait_results)
enrichment_results <- bind_rows(all_enrichment_results)
ahba_results <- bind_rows(all_ahba_results)
hub_genes <- bind_rows(all_hub_genes)
important_modules <- make_important_module_summary(trait_results, enrichment_results, ahba_results, hub_genes)

output_files <- c(
  write_csv_safely(parameter_summary, file.path(output_folder, "mdd_wgcna_parameter_sensitivity.csv")),
  write_csv_safely(bind_rows(best_gene_to_module), file.path(output_folder, "best_solution_gene_to_module_assignments.csv")),
  write_csv_safely(bind_rows(best_size_tables), file.path(output_folder, "best_solution_module_size_table.csv")),
  write_csv_safely(bind_rows(best_eigengenes), file.path(output_folder, "best_solution_module_eigengenes.csv")),
  write_csv_safely(bind_rows(best_metadata), file.path(output_folder, "best_solution_sample_metadata_used.csv")),
  write_csv_safely(trait_results, file.path(output_folder, "module_trait_association_results.csv")),
  write_csv_safely(enrichment_results, file.path(output_folder, "module_enrichment_for_Tier2_and_ML_genes.csv")),
  write_csv_safely(ahba_results, file.path(output_folder, "mdd_modules_vs_published_AHBA_modules_enrichment.csv")),
  write_csv_safely(hub_genes, file.path(output_folder, "module_hub_genes_top50.csv")),
  write_csv_safely(important_modules, file.path(output_folder, "important_modules_summary.csv"))
)

report_path <- file.path(output_folder, "MDD_WGCNA_Module_Report.md")
write_report(parameter_summary, trait_results, enrichment_results, ahba_results, important_modules, report_path, full_report = FALSE)
output_files <- c(output_files, report_path)

full_report_path <- file.path(output_folder, "MDD_WGCNA_Module_Report_FULL.md")
write_report(parameter_summary, trait_results, enrichment_results, ahba_results, important_modules, full_report_path, full_report = TRUE)
output_files <- c(output_files, full_report_path)

docs_report_path <- file.path(project_root, "docs", "MDD_WGCNA_Module_Report.md")
file.copy(full_report_path, docs_report_path, overwrite = TRUE)
output_files <- c(output_files, docs_report_path)

manifest_path <- file.path(output_folder, "mdd_wgcna_manifest.json")
manifest <- list(
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  project_root = project_root,
  inputs = list(
    metadata_file = metadata_file,
    expression_files = expression_files,
    tier2_file = tier2_file,
    ml_file = ml_file,
    ahba_file = ahba_file
    ,
    gene_annotation_folder = gene_annotation_folder
  ),
  outputs = c(output_files, na.omit(plot_files)),
    settings = list(
    run_control = run_control,
    is_limited_run = is_limited_run,
    networkType = "signed",
    TOMType = "signed",
    powers_to_test = powers_to_test,
    target_signed_r2 = target_signed_r2,
    gene_filter_targets = gene_filter_targets,
    parameter_grid = parameter_grid,
    usability_criteria = list(
      n_non_grey_modules_min = 3,
      largest_module_fraction_max = 0.70,
      grey_gene_fraction_max = 0.50,
      n_samples_min = 15
    ),
    ahba_highlight_modules = highlight_ahba_modules
  ),
  sessionInfo = capture.output(sessionInfo())
)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  warning("Optional package 'jsonlite' is not installed. Manifest will be written as plain text.")
  writeLines(capture.output(str(manifest)), manifest_path)
} else {
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE, null = "null")
}
output_files <- c(output_files, manifest_path)

message("\nDone. Main outputs:")
for (path in output_files) {
  message(" - ", path)
}
