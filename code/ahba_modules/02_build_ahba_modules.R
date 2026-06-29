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

ml_tier2_path <- file.path(clean_dir, "10_ML_prioritized_Tier2_extended", "final_ranked_ML_prioritized_Tier2_genes.csv")
ml_tier2_genes <- character(0)
if (file.exists(ml_tier2_path)) {
  ml_tier2_tbl <- read.csv(ml_tier2_path, check.names = FALSE, stringsAsFactors = FALSE)
  if ("gene_symbol" %in% names(ml_tier2_tbl)) {
    ml_tier2_genes <- unique(toupper(trimws(as.character(ml_tier2_tbl$gene_symbol))))
    ml_tier2_genes <- ml_tier2_genes[!is.na(ml_tier2_genes) & ml_tier2_genes != ""]
  }
}

soft_powers <- c(4, 6, 8, 10)
deep_splits <- c(2, 3, 4)
min_module_sizes <- c(30, 50, 100)
merge_cut_heights <- c(0.15, 0.25)
grid <- expand.grid(
  soft_power = soft_powers,
  deepSplit = deep_splits,
  minModuleSize = min_module_sizes,
  mergeCutHeight = merge_cut_heights,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

summarize_modules <- function(colors) {
  color_table <- table(colors)
  grey_count <- if ("grey" %in% names(color_table)) as.integer(color_table[["grey"]]) else 0L
  non_grey <- color_table[names(color_table) != "grey"]
  non_grey_sizes <- as.integer(non_grey)
  list(
    n_modules_excluding_grey = length(non_grey_sizes),
    grey_genes = grey_count,
    largest_module_size = if (length(non_grey_sizes)) max(non_grey_sizes) else 0L,
    median_module_size = if (length(non_grey_sizes)) as.numeric(stats::median(non_grey_sizes)) else 0,
    module_sizes = color_table
  )
}

best_solution <- NULL
audit_rows <- list()
run_index <- 0L

for (soft_power in soft_powers) {
  message("Computing signed adjacency/TOM for soft power ", soft_power)
  power_start <- Sys.time()
  power_state <- tryCatch(
    {
      adjacency_matrix <- adjacency(datExpr, power = soft_power, type = "signed")
      tom <- TOMsimilarity(adjacency_matrix, TOMType = "signed")
      diss_tom <- 1 - tom
      gene_tree <- hclust(as.dist(diss_tom), method = "average")
      list(diss_tom = diss_tom, gene_tree = gene_tree, error = NA_character_)
    },
    error = function(e) {
      list(diss_tom = NULL, gene_tree = NULL, error = conditionMessage(e))
    }
  )

  for (deep_split in deep_splits) {
    for (min_module_size in min_module_sizes) {
      for (merge_cut_height in merge_cut_heights) {
        run_index <- run_index + 1L
        run_start <- Sys.time()
        completed <- FALSE
        error_message <- power_state$error
        colors <- NULL
        eigengenes <- NULL
        module_summary <- list(
          n_modules_excluding_grey = 0L,
          grey_genes = NA_integer_,
          largest_module_size = NA_integer_,
          median_module_size = NA_real_,
          module_sizes = table(character(0))
        )

        if (is.na(error_message)) {
          result <- tryCatch(
            {
              dynamic_mods <- cutreeDynamic(
                dendro = power_state$gene_tree,
                distM = power_state$diss_tom,
                deepSplit = deep_split,
                pamRespectsDendro = FALSE,
                minClusterSize = min_module_size,
                verbose = 0
                      )
                      initial_colors <- labels2colors(dynamic_mods)
                      names(initial_colors) <- colnames(datExpr)
                      initial_summary <- summarize_modules(initial_colors)
                      if (initial_summary$n_modules_excluding_grey < 2) {
                        list(colors = structure(as.character(initial_colors), names = names(initial_colors)))
                      } else {
                        merged <- mergeCloseModules(
                          datExpr,
                          initial_colors,
                          cutHeight = merge_cut_height,
                          verbose = 0
                        )
                        names(merged$colors) <- colnames(datExpr)
                        list(colors = structure(as.character(merged$colors), names = names(merged$colors)))
                      }
            },
            error = function(e) {
              error_message <<- conditionMessage(e)
              NULL
            }
          )
          if (!is.null(result)) {
            colors <- result$colors
            completed <- TRUE
            module_summary <- summarize_modules(colors)
          }
        }

        runtime_seconds <- as.numeric(difftime(Sys.time(), run_start, units = "secs"))
        n_candidate_overlap <- if (!is.null(colors)) {
          sum(toupper(names(colors)) %in% ml_tier2_genes)
        } else {
          0L
        }
        largest_fraction <- if (!is.na(module_summary$largest_module_size)) {
          module_summary$largest_module_size / ncol(datExpr)
        } else {
          NA_real_
        }
        usable <- isTRUE(completed) &&
          module_summary$n_modules_excluding_grey > 1 &&
          is.finite(largest_fraction) &&
          largest_fraction < 0.90 &&
          n_candidate_overlap > 0

        audit_rows[[length(audit_rows) + 1L]] <- data.frame(
          run_index = run_index,
          soft_power = soft_power,
          deepSplit = deep_split,
          minModuleSize = min_module_size,
          mergeCutHeight = merge_cut_height,
          n_genes_used = ncol(datExpr),
          n_modules_excluding_grey = module_summary$n_modules_excluding_grey,
          grey_unassigned_genes = module_summary$grey_genes,
          largest_module_size = module_summary$largest_module_size,
          largest_module_fraction = largest_fraction,
          median_module_size = module_summary$median_module_size,
          ml_prioritized_tier2_overlap_genes = n_candidate_overlap,
          completed = completed,
          usable_solution = usable,
          runtime_seconds = runtime_seconds,
          error_message = ifelse(is.na(error_message), "", error_message),
          stringsAsFactors = FALSE
        )

        if (usable) {
          score <- module_summary$n_modules_excluding_grey -
            largest_fraction -
            (module_summary$grey_genes / ncol(datExpr))
          if (is.null(best_solution) || score > best_solution$score) {
            best_solution <- list(
              score = score,
              colors = colors,
              eigengenes = eigengenes,
              params = list(
                soft_power = soft_power,
                deepSplit = deep_split,
                minModuleSize = min_module_size,
                mergeCutHeight = merge_cut_height
              ),
              summary = module_summary,
              ml_prioritized_tier2_overlap_genes = n_candidate_overlap
            )
          }
        }
      }
    }
  }

  rm(power_state)
  gc()
  message("Finished soft power ", soft_power, " in ", round(as.numeric(difftime(Sys.time(), power_start, units = "mins")), 2), " minutes")
}

sensitivity <- do.call(rbind, audit_rows)
sensitivity_path <- file.path(out_dir, "ahba_wgcna_parameter_sensitivity.csv")
write.csv(sensitivity, sensitivity_path, row.names = FALSE)

params$method <- "WGCNA lower-level adjacency/TOM + cutreeDynamic + mergeCloseModules"
params$sensitivity_grid <- list(
  soft_powers = soft_powers,
  deep_splits = deep_splits,
  min_module_sizes = min_module_sizes,
  merge_cut_heights = merge_cut_heights
)
params$parameter_sensitivity_table <- sensitivity_path

if (is.null(best_solution)) {
  msg <- paste(
    "No usable multi-module AHBA WGCNA solution was found under the tested settings.",
    "The prior one-module result was not promoted as a final biological module solution."
  )
  writeLines(msg, file.path(out_dir, "AHBA_NO_USABLE_MULTIMODULE_SOLUTION.txt"))
  params$selected_solution <- NULL
  params$selection_status <- "no_usable_multi_module_solution"
  write_params(params, file.path(out_dir, "module_parameters.json"))
  report <- c(
    "# AHBA Module Construction Summary",
    "",
    msg,
    "",
    paste("- Input matrix:", expr_path),
    paste("- Samples:", nrow(datExpr)),
    paste("- Genes used:", ncol(datExpr)),
    paste("- Parameter sensitivity table:", sensitivity_path),
    "- No final module assignment table was overwritten with a forced/fake solution."
  )
  writeLines(report, file.path(out_dir, "AHBA_module_construction_summary.md"))
  stop(msg)
}

modules <- data.frame(
  gene_symbol = names(best_solution$colors),
  ahba_module = as.character(best_solution$colors),
  stringsAsFactors = FALSE
)
write.csv(modules, file.path(out_dir, "ahba_gene_to_module_assignments.csv"), row.names = FALSE)

module_sizes <- as.data.frame(table(modules$ahba_module), stringsAsFactors = FALSE)
names(module_sizes) <- c("ahba_module", "n_genes")
module_sizes <- module_sizes[order(module_sizes$n_genes, decreasing = TRUE), ]
write.csv(module_sizes, file.path(out_dir, "ahba_module_size_table.csv"), row.names = FALSE)

final_mEs <- moduleEigengenes(datExpr, best_solution$colors)$eigengenes
eigengenes <- as.data.frame(final_mEs)
eigengenes$ahba_sample_id <- rownames(eigengenes)
write.csv(eigengenes, file.path(out_dir, "ahba_module_eigengenes.csv"), row.names = FALSE)

params$selected_solution <- best_solution$params
params$selection_status <- "usable_multi_module_solution_selected"
params$ml_prioritized_tier2_overlap_genes <- best_solution$ml_prioritized_tier2_overlap_genes
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
  paste("- Parameter sensitivity table:", sensitivity_path),
  paste("- Selected soft power:", best_solution$params$soft_power),
  paste("- Selected deepSplit:", best_solution$params$deepSplit),
  paste("- Selected minimum module size:", best_solution$params$minModuleSize),
  paste("- Selected merge cut height:", best_solution$params$mergeCutHeight),
  paste("- Non-grey modules:", best_solution$summary$n_modules_excluding_grey),
  paste("- Grey/unassigned genes:", best_solution$summary$grey_genes),
  paste("- Largest module size:", best_solution$summary$largest_module_size),
  paste("- ML-prioritized Tier 2 genes overlapping module background:", best_solution$ml_prioritized_tier2_overlap_genes)
)
writeLines(report, file.path(out_dir, "AHBA_module_construction_summary.md"))
