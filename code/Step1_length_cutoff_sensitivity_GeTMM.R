## Step1_length_cutoff_sensitivity_GeTMM.R
##
## Purpose:
##   Rerun the existing GSE102556 GeTMM preprocessing while changing only the
##   minimum gene-length cutoff.
##
## This script is a sensitivity analysis.
## It does not overwrite the baseline Step 1 output folder.
##
## Length cutoffs tested:
##   Length > 20
##   Length > 50
##   Length > 100
##   Length > 200
##   Length > 500
##
## Important assumptions:
##   1. Raw count rows are genes and raw count columns are samples.
##   2. The first raw count column is GeneID.
##   3. The annotation file has GeneID and Length columns.
##   4. The brain-area mapping CSV is the same mapping used by baseline Step 1.
##   5. The per-group expression filter stays fixed:
##        keep genes whose 80th percentile raw count is > 10.
##   6. No differential expression is run in this script.

rm(list = ls())
options(stringsAsFactors = FALSE)

## -------------------------------------------------------------------------
## 0. Packages
## -------------------------------------------------------------------------

required_packages <- c("edgeR", "GEOquery", "Biobase")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

library(edgeR)
library(GEOquery)
library(Biobase)

## -------------------------------------------------------------------------
## 1. User settings
## -------------------------------------------------------------------------

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(script_path)) {
  script_path <- normalizePath("code/Step1_length_cutoff_sensitivity_GeTMM.R", mustWork = FALSE)
}

repo_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
data_dir <- "D:/genes/genedata"

raw_counts_path <- file.path(data_dir, "GSE102556_raw_counts_GRCh38.p13_NCBI.tsv")
annotation_path <- file.path(data_dir, "Human.GRCh38.p13.annot.tsv")

metadata_candidates <- c(
  file.path(data_dir, "GSE102556_metadata_clean.csv"),
  file.path(data_dir, "metadata_clean_matched_all_samples.csv"),
  file.path(data_dir, "Step1_GeTMM_log2_outputs", "metadata_clean_matched_all_samples.csv")
)

brain_area_mapping_path <- file.path(repo_dir, "config", "GSE102556_brain_area_mapping_template.csv")

output_root <- file.path(repo_dir, "outputs", "step1_length_sensitivity")

length_cutoffs <- c(20, 50, 100, 200, 500)

## Keep the same broad analysis groups used in the current project.
target_analysis_groups <- c(
  "Frontal_Cortex",
  "Insula",
  "Ventral_Striatum",
  "Hippocampal_Formation"
)

## Baseline Step 1 writes both sample-handling versions.
## Step 3 uses keep_all_samples, but this script writes both for auditing.
sample_handling_modes <- c(
  "keep_all_samples",
  "keep_highest_RIN_within_analysis_group"
)

analysis_group_column <- "BrainArea"

min_count <- 10
low_expr_exclusion_prop <- 0.80

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------------
## 2. Small helper functions
## -------------------------------------------------------------------------

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

slugify <- function(x) {
  x <- as.character(x)
  x <- gsub("\\+", "plus", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

first_present_col <- function(df, candidates) {
  for (nm in candidates) {
    if (nm %in% colnames(df)) {
      return(nm)
    }
  }
  NA_character_
}

extract_characteristic <- function(pdata, key) {
  characteristic_columns <- grep("^characteristics", colnames(pdata), value = TRUE, ignore.case = TRUE)
  output <- rep(NA_character_, nrow(pdata))

  for (row_index in seq_len(nrow(pdata))) {
    values <- as.character(unlist(pdata[row_index, characteristic_columns, drop = FALSE]))
    values <- values[!is.na(values)]
    hit <- values[grepl(paste0("^", key, "\\s*:"), values, ignore.case = TRUE)]

    if (length(hit) > 0) {
      output[row_index] <- sub("^[^:]+:\\s*", "", hit[1])
    }
  }

  output
}

extract_first_available <- function(pdata, keys) {
  output <- rep(NA_character_, nrow(pdata))

  for (key in keys) {
    temporary_values <- extract_characteristic(pdata, key)
    use_this_key <- (is.na(output) | trimws(output) == "") &
      !is.na(temporary_values) &
      trimws(temporary_values) != ""
    output[use_this_key] <- temporary_values[use_this_key]
  }

  output
}

read_expression_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("Raw count file not found: ", path)
  }

  message("Reading raw count file: ", path)
  raw_table <- read.delim(path, header = TRUE, sep = "\t", check.names = FALSE)

  gene_id <- as.character(raw_table[[1]])
  count_matrix <- as.matrix(raw_table[, -1, drop = FALSE])
  storage.mode(count_matrix) <- "numeric"
  rownames(count_matrix) <- gene_id

  count_matrix
}

collapse_duplicate_gene_ids <- function(count_matrix, output_folder) {
  duplicated_ids <- unique(rownames(count_matrix)[duplicated(rownames(count_matrix))])

  if (length(duplicated_ids) == 0) {
    duplicate_report <- data.frame(
      GeneID = character(),
      n_rows = integer(),
      stringsAsFactors = FALSE
    )
    write.csv(
      duplicate_report,
      file.path(output_folder, "duplicate_gene_ids_in_raw_counts.csv"),
      row.names = FALSE
    )
    return(count_matrix)
  }

  message("Found duplicated GeneID rows: ", length(duplicated_ids))

  duplicate_report <- data.frame(
    GeneID = duplicated_ids,
    n_rows = as.integer(table(rownames(count_matrix))[duplicated_ids]),
    stringsAsFactors = FALSE
  )
  write.csv(
    duplicate_report,
    file.path(output_folder, "duplicate_gene_ids_in_raw_counts.csv"),
    row.names = FALSE
  )

  rowsum(count_matrix, group = rownames(count_matrix), reorder = FALSE, na.rm = TRUE)
}

load_or_fetch_metadata <- function(sample_ids) {
  existing_metadata_files <- metadata_candidates[file.exists(metadata_candidates)]

  if (length(existing_metadata_files) > 0) {
    message("Reading local metadata file: ", existing_metadata_files[1])
    metadata <- read.csv(existing_metadata_files[1], check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    message("No local metadata file found. Fetching GSE102556 metadata with GEOquery.")
    gse_list <- GEOquery::getGEO("GSE102556", GSEMatrix = TRUE)
    gse <- gse_list[[1]]
    metadata <- Biobase::pData(gse)
  }

  if (!("SampleID" %in% colnames(metadata))) {
    if ("geo_accession" %in% colnames(metadata)) {
      metadata$SampleID <- as.character(metadata$geo_accession)
    } else {
      metadata$SampleID <- rownames(metadata)
    }
  }

  ## If the script receives raw GEO pData, extract the useful clean fields.
  if (!("tissue" %in% colnames(metadata)) && any(grepl("^characteristics", colnames(metadata), ignore.case = TRUE))) {
    metadata$tissue <- extract_first_available(metadata, c("tissue", "brain region", "region"))
  }
  if (!("gender" %in% colnames(metadata)) && any(grepl("^characteristics", colnames(metadata), ignore.case = TRUE))) {
    metadata$gender <- extract_first_available(metadata, c("gender", "sex"))
  }
  if (!("phenotype" %in% colnames(metadata)) && any(grepl("^characteristics", colnames(metadata), ignore.case = TRUE))) {
    metadata$phenotype <- extract_first_available(metadata, c("phenotype", "diagnosis"))
  }
  if (!("age" %in% colnames(metadata)) && any(grepl("^characteristics", colnames(metadata), ignore.case = TRUE))) {
    metadata$age <- extract_first_available(metadata, c("age"))
  }
  if (!("pmi" %in% colnames(metadata)) && any(grepl("^characteristics", colnames(metadata), ignore.case = TRUE))) {
    metadata$pmi <- extract_first_available(metadata, c("pmi"))
  }
  if (!("ph" %in% colnames(metadata)) && any(grepl("^characteristics", colnames(metadata), ignore.case = TRUE))) {
    metadata$ph <- extract_first_available(metadata, c("ph", "pH"))
  }
  if (!("rin" %in% colnames(metadata)) && any(grepl("^characteristics", colnames(metadata), ignore.case = TRUE))) {
    metadata$rin <- extract_first_available(metadata, c("rin", "RIN"))
  }

  metadata <- metadata[metadata$SampleID %in% sample_ids, , drop = FALSE]
  metadata
}

prepare_metadata <- function(metadata) {
  if (!("title" %in% colnames(metadata)) && "Title" %in% colnames(metadata)) {
    metadata$title <- metadata$Title
  }

  region_column <- first_present_col(metadata, c("OriginalRegion", "RegionID", "Region", "tissue", "tissue_ch1"))
  if (is.na(region_column)) {
    stop("Could not find a tissue/region column in metadata.")
  }

  metadata$OriginalRegion <- as.character(metadata[[region_column]])
  metadata$RegionID <- metadata$OriginalRegion

  if (!("Diagnosis" %in% colnames(metadata))) {
    diagnosis_column <- first_present_col(metadata, c("phenotype", "phenotype_ch1", "diagnosis", "Diagnosis"))
    metadata$Diagnosis <- if (!is.na(diagnosis_column)) as.character(metadata[[diagnosis_column]]) else NA_character_
  }

  if (!("Sex" %in% colnames(metadata))) {
    sex_column <- first_present_col(metadata, c("gender", "gender_ch1", "sex", "Sex"))
    metadata$Sex <- if (!is.na(sex_column)) as.character(metadata[[sex_column]]) else NA_character_
  }

  if (!("Age" %in% colnames(metadata))) {
    age_column <- first_present_col(metadata, c("age", "age_ch1", "Age"))
    metadata$Age <- if (!is.na(age_column)) safe_num(metadata[[age_column]]) else NA_real_
  }

  if (!("RIN" %in% colnames(metadata))) {
    rin_column <- first_present_col(metadata, c("RIN", "rin", "rin_ch1"))
    metadata$RIN <- if (!is.na(rin_column)) safe_num(metadata[[rin_column]]) else NA_real_
  }

  if (!("PMI" %in% colnames(metadata))) {
    pmi_column <- first_present_col(metadata, c("pmi", "pmi_ch1", "PMI"))
    metadata$PMI <- if (!is.na(pmi_column)) safe_num(metadata[[pmi_column]]) else NA_real_
  }

  if (!("pH" %in% colnames(metadata))) {
    ph_column <- first_present_col(metadata, c("ph", "ph_ch1", "pH"))
    metadata$pH <- if (!is.na(ph_column)) safe_num(metadata[[ph_column]]) else NA_real_
  }

  ## BrainID is usually the donor/person number before ":" in the GEO title.
  if (!("BrainID" %in% colnames(metadata))) {
    if ("title" %in% colnames(metadata)) {
      metadata$BrainID <- sub(":.*$", "", as.character(metadata$title))
    } else {
      metadata$BrainID <- NA_character_
    }
  }

  metadata$BrainID <- as.character(metadata$BrainID)
  bad_brain_id <- is.na(metadata$BrainID) | trimws(metadata$BrainID) == ""
  metadata$BrainID[bad_brain_id] <- paste0("unknown_", metadata$SampleID[bad_brain_id])

  metadata
}

load_brain_area_mapping <- function(original_regions) {
  if (!file.exists(brain_area_mapping_path)) {
    stop("Brain-area mapping file was not found: ", brain_area_mapping_path)
  }

  message("Reading brain-area mapping file: ", brain_area_mapping_path)
  mapping <- read.csv(brain_area_mapping_path, check.names = FALSE, stringsAsFactors = FALSE)

  required_columns <- c("OriginalRegion", "BrainArea", "MajorBrainPart")
  missing_columns <- setdiff(required_columns, colnames(mapping))
  if (length(missing_columns) > 0) {
    stop("Brain-area mapping file is missing columns: ", paste(missing_columns, collapse = ", "))
  }

  mapping$OriginalRegion <- as.character(mapping$OriginalRegion)
  mapping$BrainArea <- as.character(mapping$BrainArea)
  mapping$MajorBrainPart <- as.character(mapping$MajorBrainPart)

  unmapped_regions <- setdiff(sort(unique(original_regions)), mapping$OriginalRegion)
  if (length(unmapped_regions) > 0) {
    stop("These metadata OriginalRegion values are missing from the mapping file: ",
         paste(unmapped_regions, collapse = " | "))
  }

  mapping
}

apply_brain_area_mapping <- function(metadata, mapping) {
  match_index <- match(metadata$OriginalRegion, mapping$OriginalRegion)

  metadata$BrainArea <- mapping$BrainArea[match_index]
  metadata$MajorBrainPart <- mapping$MajorBrainPart[match_index]

  if (!(analysis_group_column %in% colnames(metadata))) {
    stop("analysis_group_column is not present after mapping: ", analysis_group_column)
  }

  metadata$AnalysisGroup <- as.character(metadata[[analysis_group_column]])

  metadata
}

keep_highest_rin_within_analysis_group <- function(metadata) {
  if (nrow(metadata) == 0) {
    return(metadata)
  }

  metadata$.original_order <- seq_len(nrow(metadata))
  metadata$.rin_for_order <- metadata$RIN
  metadata$.rin_for_order[is.na(metadata$.rin_for_order)] <- -Inf

  metadata <- metadata[
    order(metadata$AnalysisGroup, metadata$BrainID, -metadata$.rin_for_order, metadata$.original_order),
    ,
    drop = FALSE
  ]

  metadata <- metadata[!duplicated(paste(metadata$AnalysisGroup, metadata$BrainID, sep = "||")), , drop = FALSE]

  metadata$.original_order <- NULL
  metadata$.rin_for_order <- NULL

  metadata
}

make_group_composition_report <- function(metadata, output_folder) {
  groups <- sort(unique(metadata$AnalysisGroup))

  rows <- lapply(groups, function(group_name) {
    group_metadata <- metadata[metadata$AnalysisGroup == group_name, , drop = FALSE]
    region_table <- sort(table(group_metadata$OriginalRegion), decreasing = TRUE)

    data.frame(
      AnalysisGroup = group_name,
      GroupingColumn = analysis_group_column,
      n_samples = nrow(group_metadata),
      n_unique_BrainID = length(unique(group_metadata$BrainID)),
      n_original_regions = length(unique(group_metadata$OriginalRegion)),
      original_regions = paste(names(region_table), collapse = " | "),
      original_region_counts = paste(paste(names(region_table), as.integer(region_table), sep = "="), collapse = " | "),
      stringsAsFactors = FALSE
    )
  })

  report <- do.call(rbind, rows)
  write.csv(report, file.path(output_folder, "analysis_group_composition_report.csv"), row.names = FALSE)
  report
}

make_model_decision_summary <- function(metadata, output_folder) {
  donor_counts_all <- table(metadata$BrainID)
  donor_analysis_group_counts <- table(paste(metadata$BrainID, metadata$AnalysisGroup, sep = "__"))
  donor_original_region_counts <- table(paste(metadata$BrainID, metadata$OriginalRegion, sep = "__"))

  summary_table <- data.frame(
    Analysis_scope = c("original-region-wise", "brain-area-wise", "all-groups-combined"),
    Repetition_checked = c(
      "same BrainID repeated within the same OriginalRegion",
      "same BrainID repeated within the same broad AnalysisGroup",
      "same BrainID repeated across all samples/groups"
    ),
    n_repeated_BrainID_units = c(
      sum(donor_original_region_counts > 1),
      sum(donor_analysis_group_counts > 1),
      sum(donor_counts_all > 1)
    ),
    max_entries_same_BrainID = c(
      max(as.integer(donor_original_region_counts)),
      max(as.integer(donor_analysis_group_counts)),
      max(as.integer(donor_counts_all))
    ),
    stringsAsFactors = FALSE
  )

  write.csv(summary_table, file.path(output_folder, "model_decision_repetition_summary.csv"), row.names = FALSE)
  summary_table
}

filter_genes_for_group <- function(counts, gene_lengths_bp, length_cutoff) {
  lengths_for_counts <- gene_lengths_bp[rownames(counts)]

  keep_valid_length <- !is.na(lengths_for_counts) & lengths_for_counts > length_cutoff

  ## This is the existing per-group expression filter:
  ## keep genes whose 80th percentile raw count is greater than 10.
  q80 <- apply(counts, 1, quantile, probs = low_expr_exclusion_prop, na.rm = TRUE)
  keep_expression <- q80 > min_count

  keep <- keep_valid_length & keep_expression

  list(
    filtered_counts = counts[keep, , drop = FALSE],
    genes_after_length_filter = sum(keep_valid_length),
    genes_after_expression_filter = sum(keep)
  )
}

compute_getmm <- function(counts, gene_lengths_bp) {
  common_genes <- intersect(rownames(counts), names(gene_lengths_bp))

  counts <- counts[common_genes, , drop = FALSE]
  lengths_bp <- gene_lengths_bp[common_genes]

  keep_length <- !is.na(lengths_bp) & lengths_bp > 0
  counts <- counts[keep_length, , drop = FALSE]
  lengths_kb <- lengths_bp[keep_length] / 1000

  ## Step 1: gene-length correction.
  rpk <- sweep(counts, 1, lengths_kb, "/")
  rpk[!is.finite(rpk)] <- 0

  ## Step 2: TMM correction on length-adjusted expression.
  dge <- edgeR::DGEList(counts = rpk)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")

  ## Step 3: scale by effective library size.
  effective_library_size <- colSums(rpk) * dge$samples$norm.factors
  getmm <- sweep(rpk, 2, effective_library_size, "/") * 1e6

  list(
    getmm = getmm,
    log2_getmm = log2(getmm + 1),
    norm_factors = data.frame(
      SampleID = colnames(counts),
      raw_count_library_size = colSums(counts),
      RPK_library_size = colSums(rpk),
      TMM_factor = dge$samples$norm.factors,
      effective_library_size = effective_library_size,
      row.names = NULL
    )
  )
}

write_simple_donor_report <- function(metadata, output_file) {
  donor_ids <- sort(unique(metadata$BrainID))

  donor_report <- do.call(rbind, lapply(donor_ids, function(brain_id) {
    donor_metadata <- metadata[metadata$BrainID == brain_id, , drop = FALSE]

    data.frame(
      BrainID = brain_id,
      n_samples = nrow(donor_metadata),
      AnalysisGroups = paste(sort(unique(donor_metadata$AnalysisGroup)), collapse = " | "),
      OriginalRegions = paste(sort(unique(donor_metadata$OriginalRegion)), collapse = " | "),
      SampleIDs = paste(donor_metadata$SampleID, collapse = " | "),
      repeated_donor = nrow(donor_metadata) > 1,
      stringsAsFactors = FALSE
    )
  }))

  write.csv(donor_report, output_file, row.names = FALSE)
}

process_one_group <- function(
  group_name,
  mode,
  raw_counts,
  metadata,
  gene_lengths_bp,
  annotation,
  output_paths,
  length_cutoff
) {
  message("")
  message("Processing group: ", group_name)
  message("Sample handling:  ", mode)

  group_metadata_before_mode <- metadata[metadata$AnalysisGroup == group_name, , drop = FALSE]

  if (mode == "keep_all_samples") {
    group_metadata <- group_metadata_before_mode
  } else if (mode == "keep_highest_RIN_within_analysis_group") {
    group_metadata <- keep_highest_rin_within_analysis_group(group_metadata_before_mode)
  } else {
    stop("Unknown sample handling mode: ", mode)
  }

  common_samples <- intersect(colnames(raw_counts), group_metadata$SampleID)
  group_metadata <- group_metadata[match(common_samples, group_metadata$SampleID), , drop = FALSE]
  group_counts <- raw_counts[, common_samples, drop = FALSE]

  message("Samples used: ", nrow(group_metadata))
  message("Genes before filtering: ", nrow(group_counts))

  filter_result <- filter_genes_for_group(group_counts, gene_lengths_bp, length_cutoff)
  filtered_counts <- filter_result$filtered_counts

  message("Genes after Length > ", length_cutoff, ": ", filter_result$genes_after_length_filter)
  message("Genes after length + expression filter: ", nrow(filtered_counts))

  normalization_result <- compute_getmm(filtered_counts, gene_lengths_bp)
  normalized_expression <- normalization_result$log2_getmm

  group_annotation <- annotation[match(rownames(normalized_expression), annotation$GeneID), , drop = FALSE]

  file_stub <- slugify(group_name)

  typed_dirs <- list(
    metadata = file.path(output_paths$metadata, mode),
    normalized = file.path(output_paths$normalized, mode),
    filtered_counts = file.path(output_paths$filtered_counts, mode),
    normalization_factors = file.path(output_paths$normalization_factors, mode),
    gene_annotation = file.path(output_paths$gene_annotation, mode),
    repetition_reports = file.path(output_paths$repetition_reports, mode, file_stub)
  )

  invisible(lapply(typed_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

  write.csv(
    group_metadata,
    file.path(typed_dirs$metadata, paste0("metadata_", file_stub, ".csv")),
    row.names = FALSE
  )

  write.csv(
    filtered_counts,
    file.path(typed_dirs$filtered_counts, paste0("filtered_raw_counts_", file_stub, ".csv")),
    row.names = TRUE
  )

  write.csv(
    normalized_expression,
    file.path(typed_dirs$normalized, paste0("normalized_GeTMM_log2_", file_stub, ".csv")),
    row.names = TRUE
  )

  write.csv(
    normalization_result$norm_factors,
    file.path(typed_dirs$normalization_factors, paste0("GeTMM_normalization_factors_", file_stub, ".csv")),
    row.names = FALSE
  )

  write.csv(
    group_annotation,
    file.path(typed_dirs$gene_annotation, paste0("gene_annotation_", file_stub, ".csv")),
    row.names = FALSE
  )

  write_simple_donor_report(
    group_metadata,
    file.path(typed_dirs$repetition_reports, paste0("donor_repetition_", file_stub, ".csv"))
  )

  data.frame(
    length_cutoff = length_cutoff,
    sample_handling_mode = mode,
    brain_area = group_name,
    genes_after_length_filter = filter_result$genes_after_length_filter,
    genes_after_group_expression_filter = nrow(filtered_counts),
    samples = nrow(group_metadata),
    MDD_samples = sum(group_metadata$Diagnosis == "MDD", na.rm = TRUE),
    Control_samples = sum(group_metadata$Diagnosis == "CTRL", na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

make_output_paths <- function(cutoff_folder) {
  output_paths <- list(
    config = file.path(cutoff_folder, "00_config"),
    metadata = file.path(cutoff_folder, "01_metadata"),
    normalized = file.path(cutoff_folder, "02_normalized_expression"),
    filtered_counts = file.path(cutoff_folder, "03_filtered_raw_counts"),
    normalization_factors = file.path(cutoff_folder, "04_normalization_factors"),
    gene_annotation = file.path(cutoff_folder, "05_gene_annotation"),
    repetition_reports = file.path(cutoff_folder, "06_repetition_reports"),
    qc = file.path(cutoff_folder, "07_qc"),
    rds = file.path(cutoff_folder, "08_rds")
  )

  invisible(lapply(output_paths, dir.create, recursive = TRUE, showWarnings = FALSE))
  output_paths
}

## -------------------------------------------------------------------------
## 3. Load raw counts and annotation once
## -------------------------------------------------------------------------

raw_counts <- read_expression_matrix(raw_counts_path)

if (!file.exists(annotation_path)) {
  stop("Annotation file not found: ", annotation_path)
}

message("Reading annotation file: ", annotation_path)
annotation <- read.delim(annotation_path, header = TRUE, sep = "\t", check.names = FALSE)

if (!("GeneID" %in% colnames(annotation))) {
  stop("Annotation file must contain GeneID.")
}

if (!("Length" %in% colnames(annotation))) {
  stop("Annotation file must contain Length.")
}

annotation$GeneID <- as.character(annotation$GeneID)
annotation$Length <- safe_num(annotation$Length)

temporary_qc_folder <- file.path(output_root, "temporary_qc_before_cutoffs")
dir.create(temporary_qc_folder, recursive = TRUE, showWarnings = FALSE)
raw_counts <- collapse_duplicate_gene_ids(raw_counts, temporary_qc_folder)

gene_lengths_bp <- setNames(annotation$Length, annotation$GeneID)

common_genes <- intersect(rownames(raw_counts), annotation$GeneID)
raw_counts <- raw_counts[common_genes, , drop = FALSE]
annotation <- annotation[match(common_genes, annotation$GeneID), , drop = FALSE]
gene_lengths_bp <- gene_lengths_bp[common_genes]

message("Raw count matrix after annotation matching:")
print(dim(raw_counts))

## -------------------------------------------------------------------------
## 4. Load and prepare metadata once
## -------------------------------------------------------------------------

metadata <- load_or_fetch_metadata(colnames(raw_counts))
metadata <- prepare_metadata(metadata)

required_metadata_columns <- c("SampleID", "OriginalRegion", "BrainID", "Diagnosis", "Age", "Sex", "RIN")
missing_metadata_columns <- setdiff(required_metadata_columns, colnames(metadata))

if (length(missing_metadata_columns) > 0) {
  stop("Metadata is missing required columns: ", paste(missing_metadata_columns, collapse = ", "))
}

metadata$SampleID <- as.character(metadata$SampleID)

common_samples_all <- intersect(colnames(raw_counts), metadata$SampleID)
raw_counts <- raw_counts[, common_samples_all, drop = FALSE]
metadata <- metadata[match(common_samples_all, metadata$SampleID), , drop = FALSE]

mapping <- load_brain_area_mapping(metadata$OriginalRegion)
metadata <- apply_brain_area_mapping(metadata, mapping)

found_groups <- sort(unique(metadata$AnalysisGroup))
message("Analysis groups found:")
print(found_groups)

missing_groups <- setdiff(target_analysis_groups, found_groups)
if (length(missing_groups) > 0) {
  stop("These requested analysis groups were not found: ", paste(missing_groups, collapse = ", "))
}

metadata <- metadata[metadata$AnalysisGroup %in% target_analysis_groups, , drop = FALSE]
raw_counts <- raw_counts[, metadata$SampleID, drop = FALSE]

message("Samples after keeping target groups: ", nrow(metadata))
message("Sample counts by group:")
print(table(metadata$AnalysisGroup))

## -------------------------------------------------------------------------
## 5. Process every length cutoff
## -------------------------------------------------------------------------

all_summary_rows <- list()

for (length_cutoff in length_cutoffs) {
  message("")
  message("============================================================")
  message("Starting Length > ", length_cutoff)
  message("============================================================")

  cutoff_label <- paste0("Length_gt_", length_cutoff)
  cutoff_folder <- file.path(output_root, cutoff_label)
  output_paths <- make_output_paths(cutoff_folder)

  write.csv(mapping, file.path(output_paths$config, "brain_area_mapping_used.csv"), row.names = FALSE)
  write.csv(metadata, file.path(output_paths$metadata, "metadata_clean_matched_all_samples_with_brain_area.csv"), row.names = FALSE)

  composition_report <- make_group_composition_report(metadata, output_paths$qc)
  model_decision_summary <- make_model_decision_summary(metadata, output_paths$qc)

  cutoff_summary_rows <- list()

  for (mode in sample_handling_modes) {
    mode_summary_rows <- list()

    for (group_name in target_analysis_groups) {
      result_row <- process_one_group(
        group_name = group_name,
        mode = mode,
        raw_counts = raw_counts,
        metadata = metadata,
        gene_lengths_bp = gene_lengths_bp,
        annotation = annotation,
        output_paths = output_paths,
        length_cutoff = length_cutoff
      )

      mode_summary_rows[[length(mode_summary_rows) + 1]] <- result_row
      all_summary_rows[[length(all_summary_rows) + 1]] <- result_row
    }

    mode_summary <- do.call(rbind, mode_summary_rows)
    write.csv(
      mode_summary,
      file.path(output_paths$qc, paste0("Step1_QC_summary_by_analysis_group_", mode, ".csv")),
      row.names = FALSE
    )

    cutoff_summary_rows[[length(cutoff_summary_rows) + 1]] <- mode_summary
  }

  cutoff_summary_all_modes <- do.call(rbind, cutoff_summary_rows)
  write.csv(
    cutoff_summary_all_modes,
    file.path(output_paths$qc, "Step1_QC_summary_by_analysis_group_all_modes.csv"),
    row.names = FALSE
  )

  saveRDS(
    list(
      metadata = metadata,
      brain_area_mapping = mapping,
      analysis_group_column = analysis_group_column,
      analysis_group_composition = composition_report,
      annotation = annotation,
      qc_summary_all_modes = cutoff_summary_all_modes,
      model_decision_summary = model_decision_summary,
      normalization = "GeTMM_log2",
      length_cutoff = length_cutoff,
      min_count = min_count,
      low_expr_exclusion_prop = low_expr_exclusion_prop,
      sample_handling_modes = sample_handling_modes
    ),
    file = file.path(output_paths$rds, "Step1_length_cutoff_metadata_QC_and_reports.rds")
  )
}

## -------------------------------------------------------------------------
## 6. Save the requested summary CSV
## -------------------------------------------------------------------------

summary_all_modes <- do.call(rbind, all_summary_rows)

write.csv(
  summary_all_modes,
  file.path(output_root, "Step1_length_cutoff_sensitivity_summary_all_modes.csv"),
  row.names = FALSE
)

summary_keep_all_samples <- summary_all_modes[summary_all_modes$sample_handling_mode == "keep_all_samples", , drop = FALSE]
summary_keep_all_samples <- summary_keep_all_samples[, c(
  "length_cutoff",
  "brain_area",
  "genes_after_length_filter",
  "genes_after_group_expression_filter",
  "samples",
  "MDD_samples",
  "Control_samples"
)]

write.csv(
  summary_keep_all_samples,
  file.path(output_root, "Step1_length_cutoff_sensitivity_summary.csv"),
  row.names = FALSE
)

message("")
message("Finished Step 1 length-cutoff sensitivity preprocessing.")
message("Outputs were saved under:")
message(normalizePath(output_root, mustWork = FALSE))
message("")
message("Requested summary CSV:")
message(file.path(output_root, "Step1_length_cutoff_sensitivity_summary.csv"))
