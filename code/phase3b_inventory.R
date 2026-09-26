options(stringsAsFactors = FALSE)

project_root <- "D:/genes/GeneDepresion"
lib <- file.path(project_root, "environment", "phase3b", "library")
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(jsonlite)
})

out_root <- file.path(project_root, "outputs", "length_gt_100_clean_pipeline", "17_WGCNA_reconstruction_validation")
dirs <- c("inputs_manifest", "reconstruction", "module_comparison", "eigengene_validation",
          "ml_gene_enrichment_validation", "functional_annotation", "figures", "manifests")
for (d in dirs) dir.create(file.path(out_root, d), recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) unname(tools::md5sum(path))
sha256_real <- function(path) {
  cmd <- sprintf('certutil -hashfile "%s" SHA256', normalizePath(path, winslash = "\\", mustWork = TRUE))
  x <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
  x <- trimws(x)
  x[grepl("^[0-9A-Fa-f]{64}$", x)][1]
}

safe_read <- function(path) readr::read_csv(path, show_col_types = FALSE, progress = FALSE)

broad_dir <- file.path(project_root, "outputs", "step1_length_sensitivity", "Length_gt_100", "02_normalized_expression", "keep_all_samples")
frontal_dir <- file.path(project_root, "outputs", "length_gt_100_clean_pipeline", "02_normalized_expression", "frontal_subregion_wgcna_inputs")
input_rows <- tibble::tribble(
  ~region, ~path, ~role,
  "Ventral Striatum", file.path(broad_dir, "normalized_GeTMM_log2_ventral_striatum.csv"), "expression",
  "Insula", file.path(broad_dir, "normalized_GeTMM_log2_insula.csv"), "expression",
  "Hippocampal Formation", file.path(broad_dir, "normalized_GeTMM_log2_hippocampal_formation.csv"), "expression",
  "CG25", file.path(frontal_dir, "normalized_GeTMM_log2_CG25.csv"), "expression",
  "dlPFC BA8/9", file.path(frontal_dir, "normalized_GeTMM_log2_dlPFC_BA8_9.csv"), "expression",
  "OFC BA11", file.path(frontal_dir, "normalized_GeTMM_log2_OFC_BA11.csv"), "expression"
)

inspect_expression <- function(region, path, role) {
  stopifnot(file.exists(path))
  header <- names(safe_read(path))
  n <- nrow(safe_read(path))
  tibble::tibble(
    region = region, role = role, path = normalizePath(path, winslash = "/"),
    sha256 = sha256_real(path), rows_genes = n, columns_total = length(header),
    sample_count = length(header) - 1L, gene_id_column = header[1],
    expression_representation = "log2(GeTMM + 1)",
    sample_ids = paste(header[-1], collapse = ";")
  )
}

manifest <- bind_rows(Map(inspect_expression, input_rows$region, input_rows$path, input_rows$role))
write_csv(manifest, file.path(out_root, "inputs_manifest", "wgcna_reconstruction_input_manifest.csv"))
write_json(manifest, file.path(out_root, "inputs_manifest", "wgcna_reconstruction_input_manifest.json"), pretty = TRUE)

sample_manifest <- manifest %>%
  select(region, sample_ids) %>%
  tidyr::separate_longer_delim(sample_ids, delim = ";") %>%
  rename(sample_id = sample_ids)
write_csv(sample_manifest, file.path(out_root, "inputs_manifest", "wgcna_reconstruction_sample_ids.csv"))

metadata_rows <- tibble::tribble(
  ~region, ~path,
  "Ventral Striatum", file.path(project_root, "outputs", "length_gt_100_clean_pipeline", "12_MDD_WGCNA_modules", "best_solution_sample_metadata_used.csv"),
  "Insula", file.path(project_root, "outputs", "length_gt_100_clean_pipeline", "12_MDD_WGCNA_modules", "best_solution_sample_metadata_used.csv"),
  "Hippocampal Formation", file.path(project_root, "outputs", "length_gt_100_clean_pipeline", "12_MDD_WGCNA_modules", "best_solution_sample_metadata_used.csv"),
  "CG25", file.path(frontal_dir, "normalized_GeTMM_log2_CG25_sample_metadata.csv"),
  "dlPFC BA8/9", file.path(frontal_dir, "normalized_GeTMM_log2_dlPFC_BA8_9_sample_metadata.csv"),
  "OFC BA11", file.path(frontal_dir, "normalized_GeTMM_log2_OFC_BA11_sample_metadata.csv")
)
metadata_manifest <- metadata_rows %>% rowwise() %>% mutate(
  path = normalizePath(path, winslash = "/", mustWork = TRUE),
  sha256 = sha256_real(path),
  rows_in_file = nrow(safe_read(path)),
  columns_in_file = ncol(safe_read(path))
) %>% ungroup()
write_csv(metadata_manifest, file.path(out_root, "inputs_manifest", "wgcna_reconstruction_metadata_manifest.csv"))

select_best <- function(path) {
  safe_read(path) %>%
    filter(region != "Frontal Cortex") %>%
    mutate(usable_rank = if_else(usable_solution, 1L, 0L), filter_rank = n_genes_after_filter) %>%
    arrange(region, desc(usable_rank), desc(filter_rank), desc(n_non_grey_modules), largest_module_fraction) %>%
    group_by(region) %>% slice(1) %>% ungroup()
}

broad_params <- select_best(file.path(project_root, "outputs", "length_gt_100_clean_pipeline", "12_MDD_WGCNA_modules", "mdd_wgcna_parameter_sensitivity.csv"))
frontal_param_candidates <- c(
  file.path(project_root, "outputs", "length_gt_100_clean_pipeline", "12_MDD_WGCNA_modules", "frontal_subregion_wgcna", "frontal_subregion_wgcna_parameter_sensitivity.csv"),
  file.path(project_root, "outputs", "length_gt_100_clean_pipeline", "mdd_wgcna_modules", "frontal_subregion_wgcna", "frontal_subregion_wgcna_parameter_sensitivity.csv")
)
frontal_param_path <- frontal_param_candidates[file.exists(frontal_param_candidates)][1]
frontal_params <- select_best(frontal_param_path)
best <- bind_rows(broad_params, frontal_params)
write_csv(best, file.path(out_root, "inputs_manifest", "historical_selected_parameter_rows.csv"))

canonical_rows <- tibble::tribble(
  ~role, ~path,
  "canonical_broad_assignments", file.path(project_root, "outputs", "length_gt_100_clean_pipeline", "12_MDD_WGCNA_modules", "best_solution_gene_to_module_assignments.csv"),
  "canonical_frontal_assignments", file.path(project_root, "outputs", "length_gt_100_clean_pipeline", "mdd_wgcna_modules", "frontal_subregion_wgcna", "frontal_subregion_best_solution_gene_to_module_assignments.csv"),
  "canonical_corrected_ml_ranking", file.path(project_root, "outputs", "length_gt_100_clean_pipeline", "10_ML_prioritized_Tier2_corrected", "ML_prioritized_Tier2_genes_all_ranked.csv")
)
canonical_files <- canonical_rows %>% rowwise() %>% mutate(
  exists = file.exists(path), sha256 = if (exists) sha256_real(path) else NA_character_,
  rows = if (exists) nrow(safe_read(path)) else NA_integer_,
  columns = if (exists) ncol(safe_read(path)) else NA_integer_
) %>% ungroup()
write_csv(canonical_files, file.path(out_root, "inputs_manifest", "canonical_reference_file_manifest.csv"))

writeLines(capture.output(sessionInfo()), file.path(out_root, "manifests", "phase3b_sessionInfo.txt"))

cat("Input manifest rows:", nrow(manifest), "\n")
print(manifest %>% select(region, sha256, rows_genes, sample_count))
cat("\nHistorical selected settings:\n")
print(best %>% select(region, analysis_variant, gene_filter_version, n_samples, n_genes_after_filter,
                      power, power_selection_status, deepSplit, minModuleSize, mergeCutHeight,
                      corType, n_non_grey_modules, grey_gene_fraction, largest_module_fraction))
cat("\nCanonical reference files:\n")
print(canonical_files)
