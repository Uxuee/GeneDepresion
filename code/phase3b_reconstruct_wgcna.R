## Phase 3B: independent WGCNA reconstruction into a validation-only directory.
## This script does not read or overwrite manuscript outputs and does not make
## reconstructed assignments canonical.

options(stringsAsFactors = FALSE)
project_root <- "D:/genes/GeneDepresion"
.libPaths(c(file.path(project_root, "environment", "phase3b", "library"), .libPaths()))

suppressPackageStartupMessages({
  library(WGCNA)
  library(dynamicTreeCut)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(limma)
  library(jsonlite)
})
WGCNA::allowWGCNAThreads()

args <- commandArgs(trailingOnly = TRUE)
requested_region <- if (length(args)) paste(args, collapse = " ") else "all"
out_root <- file.path(project_root, "outputs", "length_gt_100_clean_pipeline", "17_WGCNA_reconstruction_validation")
recon_dir <- file.path(out_root, "reconstruction")
manifest_dir <- file.path(out_root, "manifests")
dir.create(recon_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

input_manifest <- read_csv(file.path(out_root, "inputs_manifest", "wgcna_reconstruction_input_manifest.csv"), show_col_types = FALSE)
historical <- read_csv(file.path(out_root, "inputs_manifest", "historical_selected_parameter_rows.csv"), show_col_types = FALSE)

normalize_region <- function(x) {
  y <- gsub("[^a-z0-9]+", " ", tolower(x))
  y <- trimws(gsub(" +", " ", y))
  dplyr::case_when(
    y %in% c("cg25", "cingulate gyrus 25") ~ "Cingulate gyrus 25",
    y %in% c("dlpfc ba8 9", "dlpfc ba8 9 ") ~ "dlPFC BA8/9",
    y == "ofc ba11" ~ "OFC BA11",
    y == "ventral striatum" ~ "Ventral Striatum",
    y == "insula" ~ "Insula",
    y == "hippocampal formation" ~ "Hippocampal Formation",
    TRUE ~ x
  )
}
input_manifest$region <- vapply(input_manifest$region, normalize_region, character(1))

regions <- c("Ventral Striatum", "Insula", "Hippocampal Formation", "Cingulate gyrus 25", "dlPFC BA8/9", "OFC BA11")
if (requested_region != "all") {
  requested_region <- normalize_region(requested_region)
  if (!requested_region %in% regions) stop("Unknown region: ", requested_region)
  regions <- requested_region
}

safe_name <- function(x) trimws(gsub("(^_|_$)", "", gsub("[^A-Za-z0-9]+", "_", x)))
first_col <- function(x, candidates) { y <- candidates[candidates %in% names(x)]; if (length(y)) y[1] else NA_character_ }
`%||%` <- function(x, y) if (is.null(x) || !length(x) || all(is.na(x))) y else x

read_expression <- function(path) {
  raw <- read_csv(path, show_col_types = FALSE, progress = FALSE)
  gc <- first_col(raw, c("gene_symbol", "Gene", "gene", "symbol", "gene_id"))
  if (is.na(gc)) gc <- names(raw)[1]
  genes <- trimws(as.character(raw[[gc]])); genes[genes == ""] <- NA_character_
  expr <- as.data.frame(lapply(raw[setdiff(names(raw), gc)], as.numeric), check.names = FALSE)
  rownames(expr) <- genes
  expr <- expr[!is.na(rownames(expr)) & !duplicated(rownames(expr)), , drop = FALSE]
  x <- t(as.matrix(expr)); storage.mode(x) <- "numeric"; x
}

prepare_metadata <- function(data) {
  sc <- first_col(data, c("sample_id", "SampleID", "sample", "geo_accession", "GSM", "Sample", ".sample_id"))
  if (is.na(sc)) stop("No sample ID column in metadata")
  raw_diag <- if ("Diagnosis" %in% names(data)) data$Diagnosis else if ("phenotype" %in% names(data)) data$phenotype else data[["phenotype:ch1"]]
  data %>% mutate(
    .sample_id = as.character(.data[[sc]]),
    .Diagnosis_raw = as.character(raw_diag),
    Diagnosis_clean = case_when(
      str_detect(str_to_lower(.Diagnosis_raw), "mdd|major|depress") ~ "MDD",
      str_detect(str_to_lower(.Diagnosis_raw), "ctrl|control") ~ "Control",
      TRUE ~ .Diagnosis_raw),
    Diagnosis_numeric = case_when(Diagnosis_clean == "MDD" ~ 1, Diagnosis_clean == "Control" ~ 0, TRUE ~ NA_real_)
  )
}

metadata_for_region <- function(region, sample_ids) {
  if (region %in% c("Ventral Striatum", "Insula", "Hippocampal Formation")) {
    p <- file.path(project_root, "outputs", "length_gt_100_clean_pipeline", "12_MDD_WGCNA_modules", "best_solution_sample_metadata_used.csv")
  } else {
    stub <- c("Cingulate gyrus 25" = "CG25", "dlPFC BA8/9" = "dlPFC_BA8_9", "OFC BA11" = "OFC_BA11")[[region]]
    p <- file.path(project_root, "outputs", "length_gt_100_clean_pipeline", "02_normalized_expression", "frontal_subregion_wgcna_inputs", paste0("normalized_GeTMM_log2_", stub, "_sample_metadata.csv"))
  }
  m <- prepare_metadata(read_csv(p, show_col_types = FALSE, progress = FALSE))
  ids <- intersect(sample_ids, m$.sample_id)
  m <- m[match(ids, m$.sample_id), , drop = FALSE]
  rownames(m) <- m$.sample_id
  list(data = m, path = p, ids = ids)
}

usable_covariates <- function(metadata, candidates) {
  available <- candidates[candidates %in% names(metadata)]
  available[vapply(available, function(z) {
    v <- metadata[[z]]; v <- v[!is.na(v)]; length(v) >= 3 && length(unique(v)) >= 2
  }, logical(1))]
}

residualize_expression <- function(dat_expr, metadata) {
  covars <- usable_covariates(metadata, c("BrainBank", "SeqBatch", "RIN", "rin", "Age", "age", "Sex", "gender"))
  cc <- complete.cases(metadata[, covars, drop = FALSE])
  design <- model.matrix(as.formula(paste("~", paste(covars, collapse = " + "))), metadata[cc, , drop = FALSE])
  fit <- lm.fit(design, dat_expr[cc, , drop = FALSE])
  x <- fit$residuals; rownames(x) <- rownames(dat_expr)[cc]; colnames(x) <- colnames(dat_expr)
  list(expr = x, metadata = metadata[cc, , drop = FALSE], covariates = covars, design_rank = fit$rank, design_columns = ncol(design))
}

filter_mad <- function(dat_expr, n) {
  usable <- dat_expr[, colSums(is.na(dat_expr)) == 0, drop = FALSE]
  gm <- apply(usable, 2, mad, na.rm = TRUE)
  gv <- apply(usable, 2, var, na.rm = TRUE)
  keep <- is.finite(gm) & is.finite(gv) & gv > 0
  usable <- usable[, keep, drop = FALSE]; gm <- gm[keep]
  gsg <- WGCNA::goodSamplesGenes(usable, verbose = 0)
  usable <- usable[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  gm <- gm[colnames(usable)]
  if (ncol(usable) > n) usable <- usable[, names(sort(gm, decreasing = TRUE))[seq_len(n)], drop = FALSE]
  usable
}

choose_power <- function(dat_expr, cor_type = "bicor") {
  cf <- if (cor_type == "bicor") "bicor" else "cor"
  co <- if (cor_type == "bicor") "maxPOutliers = 0.05, use = 'pairwise.complete.obs'" else "use = 'pairwise.complete.obs'"
  fit <- tryCatch(WGCNA::pickSoftThreshold(dat_expr, powerVector = 1:20, networkType = "signed", corFnc = cf, corOptions = co, verbose = 0), error = identity)
  if (inherits(fit, "error") && cor_type == "bicor") return(choose_power(dat_expr, "pearson"))
  if (inherits(fit, "error")) stop(fit)
  fi <- as_tibble(fit$fitIndices)
  rc <- intersect(c("SFT.R.sq", "Rsquared.SFT"), names(fi))[1]
  fi$.signed_r2 <- fi[[rc]]
  reached <- fi %>% filter(.signed_r2 >= 0.80)
  if (nrow(reached)) list(power = reached$Power[1], status = "first_power_reached_0.8", fit = fi, corType = cor_type)
  else { b <- fi %>% arrange(desc(.signed_r2), Power) %>% slice(1); list(power = b$Power, status = "no_power_reached_0.8", fit = fi, corType = cor_type) }
}

summarize_colors <- function(colors) {
  s <- sort(table(colors), decreasing = TRUE); ng <- s[names(s) != "grey"]; grey <- as.integer(s["grey"] %||% 0)
  tibble(n_non_grey_modules = length(ng), largest_module_size = ifelse(length(ng), max(ng), 0),
         largest_module_fraction = ifelse(length(ng), max(ng) / length(colors), 0),
         grey_gene_count = grey, grey_gene_fraction = grey / length(colors),
         median_module_size = ifelse(length(ng), median(as.integer(ng)), NA_real_))
}

fit_traits <- function(region, MEs, metadata, sizes) {
  covars <- usable_covariates(metadata, c("Age", "age", "Sex", "gender", "BrainBank", "SeqBatch", "RIN", "rin", "PMI", "pmi"))
  rows <- lapply(colnames(MEs), function(me) {
    d <- metadata %>% mutate(.ME = MEs[, me]) %>% filter(!is.na(Diagnosis_numeric), !is.na(.ME))
    ftxt <- paste(".ME ~ Diagnosis_numeric", paste(covars, collapse = " + "), sep = ifelse(length(covars), " + ", ""))
    note <- "full_or_available_covariate_model"
    if (nrow(d) < 5 + length(covars) * 2) { ftxt <- ".ME ~ Diagnosis_numeric"; note <- "reduced_model_due_to_sample_size" }
    fit <- lm(as.formula(ftxt), d); ct <- summary(fit)$coefficients
    tibble(region = region, reconstructed_module = me, n_module_genes = sizes$n_genes[match(me, sizes$module)],
           beta_diagnosis = ct["Diagnosis_numeric", "Estimate"], p_diagnosis = ct["Diagnosis_numeric", "Pr(>|t|)"],
           model_formula = ftxt, n_samples = nrow(d), notes = note)
  })
  bind_rows(rows) %>% group_by(region) %>% mutate(fdr_diagnosis_BH = if_else(reconstructed_module == "MEgrey", NA_real_, p.adjust(if_else(reconstructed_module == "MEgrey", NA_real_, p_diagnosis), "BH"))) %>% ungroup()
}

grid <- expand.grid(deepSplit = c(2,3,4), minModuleSize = c(20,30,50), mergeCutHeight = c(.15,.25), KEEP.OUT.ATTRS = FALSE)

for (region in regions) {
  tag <- safe_name(region); region_dir <- file.path(recon_dir, tag); dir.create(region_dir, recursive = TRUE, showWarnings = FALSE)
  done <- file.path(region_dir, "COMPLETE.txt")
  if (file.exists(done)) { message("Already complete: ", region); next }
  existing_assignments <- file.path(region_dir, "reconstructed_gene_to_module_assignments.csv")
  existing_eigengenes <- file.path(region_dir, "reconstructed_module_eigengenes.csv")
  if (file.exists(existing_assignments) && file.exists(existing_eigengenes)) {
    message("Resuming post-network outputs for ", region)
    hist <- historical %>% filter(region == !!region) %>% slice(1)
    im <- input_manifest %>% filter(region == !!region) %>% slice(1)
    dat <- readRDS(file.path(region_dir, "reconstructed_filtered_expression.rds"))
    metadata <- metadata_for_region(region, rownames(dat))$data
    metadata <- metadata[match(rownames(dat), metadata$.sample_id), , drop = FALSE]
    write_csv(as_tibble(metadata), file.path(region_dir, "reconstruction_sample_metadata.csv"))
    assignments <- read_csv(existing_assignments, show_col_types = FALSE) %>% mutate(gene_id = as.character(gene_id))
    sizes <- read_csv(file.path(region_dir, "reconstructed_module_sizes.csv"), show_col_types = FALSE)
    me_tbl <- read_csv(existing_eigengenes, show_col_types = FALSE)
    MEs <- as.matrix(me_tbl[, setdiff(names(me_tbl), "sample_id"), drop = FALSE]); storage.mode(MEs) <- "numeric"; rownames(MEs) <- me_tbl$sample_id
    traits <- fit_traits(region, MEs, metadata, sizes)
    write_csv(traits, file.path(region_dir, "reconstructed_module_trait_associations.csv"))
    hubs <- lapply(setdiff(unique(assignments$reconstructed_module), "MEgrey"), function(me) {
      genes <- assignments$gene_id[assignments$reconstructed_module == me]
      k <- cor(dat[, genes, drop = FALSE], MEs[, me], use = "pairwise.complete.obs")[,1]
      tibble(region = region, reconstructed_module = me, gene_id = genes, kME = as.numeric(k), abs_kME = abs(as.numeric(k))) %>% arrange(desc(abs_kME)) %>% slice_head(n = 50)
    }) %>% bind_rows()
    write_csv(hubs, file.path(region_dir, "reconstructed_module_hub_genes_top50.csv"))
    best <- read_csv(file.path(region_dir, "selected_solution_parameters.csv"), show_col_types = FALSE)
    manifest <- list(region = region, status = "validation_candidate", generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      git_baseline_requested = "a483b4f6267befbb98ef49f23a27ffe6fa48ed69", git_head_at_start = "df46778d1b62d92c8ba2b0a3b9002e3c3f7bdadc",
      expression_path = im$path, expression_sha256 = im$sha256, construction_used_diagnosis = FALSE,
      analysis_variant = hist$analysis_variant, gene_filter = hist$gene_filter_version, powers_tested = 1:20,
      target_signed_r2 = .80, selected_power = best$power, power_selection_status = best$power_selection_status,
      correlation = best$corType, maxPOutliers = ifelse(best$corType == "bicor", .05, NA), networkType = "signed", TOMType = "signed",
      dynamic_cut_grid = grid, selection_rule = "usable first; larger gene universe; more non-grey modules; smaller largest-module fraction",
      selected_solution = best, packages = capture.output(sessionInfo()))
    write_json(manifest, file.path(region_dir, "reconstruction_manifest.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
    writeLines(paste("Completed", region, format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), done)
    next
  }
  message("\n=== Phase 3B reconstruction: ", region, " ===")
  hist <- historical %>% filter(region == !!region) %>% slice(1)
  im <- input_manifest %>% filter(region == !!region) %>% slice(1)
  raw_expr <- read_expression(im$path)
  mi <- metadata_for_region(region, rownames(raw_expr))
  ids <- mi$ids; dat <- raw_expr[ids, , drop = FALSE]; metadata <- mi$data
  variant_info <- list(covariates = character(), design_rank = NA_integer_, design_columns = NA_integer_)
  if (hist$analysis_variant == "residualized_technical_covariates") {
    rr <- residualize_expression(dat, metadata); dat <- rr$expr; metadata <- rr$metadata; variant_info <- rr
  }
  target <- hist$n_genes_after_filter
  dat <- filter_mad(dat, target)
  metadata <- metadata[match(rownames(dat), metadata$.sample_id), , drop = FALSE]
  saveRDS(dat, file.path(region_dir, "reconstructed_filtered_expression.rds"), compress = "xz")
  write_csv(as_tibble(metadata), file.path(region_dir, "reconstruction_sample_metadata.csv"))
  write_csv(tibble(gene_id = colnames(dat), mad_rank = seq_len(ncol(dat))), file.path(region_dir, "mad_filtered_gene_order.csv"))

  soft <- choose_power(dat, "bicor")
  write_csv(soft$fit, file.path(region_dir, "soft_threshold_diagnostics.csv"))
  if (Sys.getenv("PHASE3B_POWER_ONLY", unset = "0") == "1") {
    message("Power-only probe completed for ", region, ": ", soft$power, " (", soft$status, ")")
    quit(save = "no", status = 0)
  }
  cf <- if (soft$corType == "bicor") "bicor" else "cor"
  co <- if (soft$corType == "bicor") "maxPOutliers = 0.05, use = 'pairwise.complete.obs'" else "use = 'pairwise.complete.obs'"
  adjacency_matrix <- WGCNA::adjacency(dat, power = soft$power, type = "signed", corFnc = cf, corOptions = co)
  tom <- WGCNA::TOMsimilarity(adjacency_matrix, TOMType = "signed", verbose = 1)
  rm(adjacency_matrix); gc()
  diss <- 1 - tom; rm(tom); gc()
  tree <- hclust(as.dist(diss), method = "average")

  solutions <- vector("list", nrow(grid)); summaries <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    gr <- grid[i, ]
    dyn <- dynamicTreeCut::cutreeDynamic(dendro = tree, distM = diss, deepSplit = gr$deepSplit,
                                         pamRespectsDendro = FALSE, minClusterSize = gr$minModuleSize, verbose = 0)
    cols <- WGCNA::labels2colors(dyn)
    merged <- WGCNA::mergeCloseModules(dat, cols, cutHeight = gr$mergeCutHeight, verbose = 0)
    cols <- merged$colors
    sm <- summarize_colors(cols)
    usable <- sm$n_non_grey_modules >= 3 && sm$largest_module_fraction < .70 && sm$grey_gene_fraction < .50 && nrow(dat) >= 15
    summaries[[i]] <- bind_cols(tibble(region = region, analysis_variant = hist$analysis_variant,
      gene_filter_version = hist$gene_filter_version, n_samples = nrow(dat), n_genes_after_filter = ncol(dat),
      power = soft$power, power_selection_status = soft$status, deepSplit = gr$deepSplit,
      minModuleSize = gr$minModuleSize, mergeCutHeight = gr$mergeCutHeight, corType = soft$corType,
      networkType = "signed", TOMType = "signed"), sm) %>% mutate(usable_solution = usable)
    solutions[[i]] <- cols
  }
  grid_table <- bind_rows(summaries) %>% arrange(desc(usable_solution), desc(n_genes_after_filter), desc(n_non_grey_modules), largest_module_fraction)
  best <- grid_table[1, ]
  bi <- which(vapply(summaries, function(x) x$deepSplit == best$deepSplit && x$minModuleSize == best$minModuleSize && x$mergeCutHeight == best$mergeCutHeight, logical(1)))[1]
  colors <- solutions[[bi]]
  write_csv(bind_rows(summaries), file.path(region_dir, "dynamic_cut_parameter_grid.csv"))
  write_csv(best, file.path(region_dir, "selected_solution_parameters.csv"))

  assignments <- tibble(region = region, gene_id = colnames(dat), module_color = colors, reconstructed_module = paste0("ME", colors))
  sizes <- assignments %>% count(region, module = reconstructed_module, module_color, name = "n_genes") %>% arrange(desc(n_genes))
  write_csv(assignments, file.path(region_dir, "reconstructed_gene_to_module_assignments.csv"))
  write_csv(sizes, file.path(region_dir, "reconstructed_module_sizes.csv"))
  MEs <- WGCNA::moduleEigengenes(dat, colors = colors)$eigengenes
  rownames(MEs) <- rownames(dat); MEs <- WGCNA::orderMEs(MEs); rownames(MEs) <- rownames(dat)
  write_csv(as_tibble(MEs, rownames = "sample_id"), file.path(region_dir, "reconstructed_module_eigengenes.csv"))
  traits <- fit_traits(region, MEs, metadata, sizes)
  write_csv(traits, file.path(region_dir, "reconstructed_module_trait_associations.csv"))

  hubs <- lapply(setdiff(unique(assignments$reconstructed_module), "MEgrey"), function(me) {
    genes <- assignments$gene_id[assignments$reconstructed_module == me]
    k <- cor(dat[, genes, drop = FALSE], MEs[, me], use = "pairwise.complete.obs")[,1]
    tibble(region = region, reconstructed_module = me, gene_id = genes, kME = as.numeric(k), abs_kME = abs(as.numeric(k))) %>% arrange(desc(abs_kME)) %>% slice_head(n = 50)
  }) %>% bind_rows()
  write_csv(hubs, file.path(region_dir, "reconstructed_module_hub_genes_top50.csv"))

  manifest <- list(region = region, status = "validation_candidate", generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    git_baseline_requested = "a483b4f6267befbb98ef49f23a27ffe6fa48ed69", git_head_at_start = "df46778d1b62d92c8ba2b0a3b9002e3c3f7bdadc",
    expression_path = im$path, expression_sha256 = im$sha256, metadata_path = mi$path,
    construction_used_diagnosis = FALSE, analysis_variant = hist$analysis_variant,
    residualization_covariates = variant_info$covariates, residualization_design_rank = variant_info$design_rank,
    residualization_design_columns = variant_info$design_columns, gene_filter = hist$gene_filter_version,
    powers_tested = 1:20, target_signed_r2 = .80, selected_power = soft$power, power_selection_status = soft$status,
    correlation = soft$corType, maxPOutliers = ifelse(soft$corType == "bicor", .05, NA), networkType = "signed", TOMType = "signed",
    dynamic_cut_grid = grid, selection_rule = "usable first; larger gene universe; more non-grey modules; smaller largest-module fraction",
    selected_solution = best, packages = capture.output(sessionInfo()))
  write_json(manifest, file.path(region_dir, "reconstruction_manifest.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
  writeLines(paste("Completed", region, format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), done)
  rm(raw_expr, dat, diss, tree, solutions, colors, MEs); gc()
  message("Completed ", region, ": power ", soft$power, ", modules ", best$n_non_grey_modules)
}
