## Step 1 for GSE102556
## Goal: make one clean normalized expression output per brain region.
## Normalization used here: GeTMM followed by log2(GeTMM + 1).
##
## This script intentionally does NOT run differential expression, PCA,
## WGCNA, or machine learning. It only prepares filtered raw counts,
## metadata, gene annotation, normalization factors, and QC reports.

rm(list = ls())
options(stringsAsFactors = FALSE)

message("Starting Step1.R for GSE102556...")

## -----------------------------
## 1. User settings
## -----------------------------

data_dir <- "D:/genes/genedata"

raw_counts_file <- file.path(data_dir, "GSE102556_raw_counts_GRCh38.p13_NCBI.tsv")
annotation_file <- file.path(data_dir, "Human.GRCh38.p13.annot.tsv")
metadata_clean_file <- file.path(data_dir, "GSE102556_metadata_clean.csv")

output_dir <- file.path(data_dir, "Step1_GeTMM_log2_outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## Low-expression filter:
## Keep genes where the 80th percentile raw count is > 10 within each region.
low_expression_quantile <- 0.80
low_expression_count_cutoff <- 10

## Short-gene filter:
## Keep genes longer than 20 bp, matching your old Step1.R/Ferrari-inspired rule.
minimum_gene_length_bp <- 20

## -----------------------------
## 2. Package setup
## -----------------------------

if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop("The edgeR package is required. Install it with BiocManager::install('edgeR').")
}

## GEOquery is only needed when GSE102556_metadata_clean.csv is not present.
if (!file.exists(metadata_clean_file) && !requireNamespace("GEOquery", quietly = TRUE)) {
  stop("GSE102556_metadata_clean.csv was not found, and GEOquery is not installed. Install it with BiocManager::install('GEOquery').")
}

if (!file.exists(metadata_clean_file) && !requireNamespace("Biobase", quietly = TRUE)) {
  stop("GSE102556_metadata_clean.csv was not found, and Biobase is needed to read GEOquery metadata. Install it with BiocManager::install('Biobase').")
}

## -----------------------------
## 3. Small helper functions
## -----------------------------

clean_column_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

safe_file_part <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nchar(x) == 0, "unknown", x)
}

first_matching_column <- function(df, patterns) {
  column_names <- names(df)
  for (pattern in patterns) {
    hit <- grep(pattern, column_names, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) {
      return(hit[1])
    }
  }
  NA_character_
}

coalesce_character <- function(x, fallback) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- fallback[is.na(x) | trimws(x) == ""]
  x
}

write_report <- function(x, filename) {
  full_path <- file.path(output_dir, filename)
  write.csv(x, full_path, row.names = FALSE)
  message("Wrote: ", full_path)
}

standardize_metadata <- function(metadata) {
  names(metadata) <- clean_column_name(names(metadata))

  sample_col <- first_matching_column(
    metadata,
    c("^geo_accession$", "^sample$", "^sample_id$", "^gsm$", "accession")
  )
  region_col <- first_matching_column(
    metadata,
    c("^RegionID$", "brain.*region", "^region$", "tissue", "source.*name", "brain")
  )
  brain_id_col <- first_matching_column(
    metadata,
    c("^BrainID$", "brain.*id", "donor", "subject", "individual", "person", "case")
  )
  rin_col <- first_matching_column(
    metadata,
    c("^RIN$", "rna.*integrity", "rin")
  )

  if (is.na(sample_col)) {
    if (!is.null(rownames(metadata)) && all(grepl("^GSM", rownames(metadata)))) {
      metadata$SampleID <- rownames(metadata)
    } else {
      stop("Could not find a sample/GSM accession column in metadata.")
    }
  } else {
    metadata$SampleID <- metadata[[sample_col]]
  }

  if (is.na(region_col)) {
    stop("Could not find a brain region column in metadata. Please add one to GSE102556_metadata_clean.csv.")
  }
  metadata$RegionID <- metadata[[region_col]]

  if (is.na(brain_id_col)) {
    if ("title" %in% names(metadata) && any(grepl(":", metadata$title))) {
      message("No explicit BrainID column found. Using the donor/person prefix before ':' in the GEO title field.")
      metadata$BrainID <- trimws(sub(":.*$", "", metadata$title))
      metadata$BrainID <- coalesce_character(metadata$BrainID, metadata$SampleID)
    } else {
      warning("Could not find BrainID/donor/person column. SampleID will be used as BrainID.")
      metadata$BrainID <- metadata$SampleID
    }
  } else {
    metadata$BrainID <- metadata[[brain_id_col]]
    metadata$BrainID <- coalesce_character(metadata$BrainID, metadata$SampleID)
  }

  if (is.na(rin_col)) {
    warning("Could not find a RIN column. RIN will be NA, so deduplication keeps the first sample within each RegionID + BrainID group.")
    metadata$RIN <- NA_real_
  } else {
    metadata$RIN <- suppressWarnings(as.numeric(metadata[[rin_col]]))
  }

  metadata$SampleID <- as.character(metadata$SampleID)
  metadata$RegionID <- as.character(metadata$RegionID)
  metadata$BrainID <- as.character(metadata$BrainID)

  metadata <- metadata[!is.na(metadata$SampleID) & metadata$SampleID != "", , drop = FALSE]
  rownames(metadata) <- metadata$SampleID
  metadata
}

parse_geo_characteristics <- function(pdata) {
  characteristics_cols <- grep("^characteristics_ch", names(pdata), value = TRUE)
  if (length(characteristics_cols) == 0) {
    return(pdata)
  }

  message("Parsing GEO characteristics_ch columns...")
  parsed <- list()

  for (col in characteristics_cols) {
    values <- as.character(pdata[[col]])
    has_key_value <- grepl(":", values)
    keys <- trimws(sub(":.*$", "", values[has_key_value]))

    for (key in unique(keys)) {
      new_col <- clean_column_name(key)
      key_pattern <- paste0("^", gsub("([\\W])", "\\\\\\1", key), "\\s*:")
      extracted <- rep(NA_character_, length(values))
      matched <- grepl(key_pattern, values, ignore.case = TRUE)
      extracted[matched] <- trimws(sub("^[^:]+:\\s*", "", values[matched]))

      if (is.null(parsed[[new_col]])) {
        parsed[[new_col]] <- extracted
      } else {
        old <- parsed[[new_col]]
        old[is.na(old) | old == ""] <- extracted[is.na(old) | old == ""]
        parsed[[new_col]] <- old
      }
    }
  }

  for (new_col in names(parsed)) {
    if (!new_col %in% names(pdata)) {
      pdata[[new_col]] <- parsed[[new_col]]
    }
  }

  pdata
}

load_metadata <- function() {
  if (file.exists(metadata_clean_file)) {
    message("Reading local metadata: ", metadata_clean_file)
    metadata <- read.csv(metadata_clean_file, check.names = FALSE)
    return(standardize_metadata(metadata))
  }

  message("Local metadata not found. Fetching GSE102556 metadata using GEOquery::getGEO()...")
  geo <- GEOquery::getGEO("GSE102556", GSEMatrix = TRUE)
  if (length(geo) > 1) {
    geo <- geo[[1]]
  } else {
    geo <- geo[[1]]
  }

  metadata <- Biobase::pData(geo)
  metadata <- parse_geo_characteristics(metadata)
  standardize_metadata(metadata)
}

make_donor_repetition_overall <- function(metadata) {
  donor_counts <- as.data.frame(table(metadata$BrainID), stringsAsFactors = FALSE)
  names(donor_counts) <- c("BrainID", "number_of_samples")
  donor_counts <- donor_counts[order(-donor_counts$number_of_samples, donor_counts$BrainID), ]
  donor_counts$repeated <- donor_counts$number_of_samples > 1
  donor_counts
}

make_donor_repetition_by_region <- function(metadata) {
  donor_counts <- as.data.frame(table(metadata$RegionID, metadata$BrainID), stringsAsFactors = FALSE)
  names(donor_counts) <- c("RegionID", "BrainID", "number_of_samples")
  donor_counts <- donor_counts[donor_counts$number_of_samples > 0, , drop = FALSE]
  donor_counts <- donor_counts[order(donor_counts$RegionID, -donor_counts$number_of_samples, donor_counts$BrainID), ]
  donor_counts$repeated_within_region <- donor_counts$number_of_samples > 1
  donor_counts
}

make_gene_entry_repetition_summary <- function(count_matrix, metadata) {
  sample_ids <- colnames(count_matrix)
  metadata <- metadata[sample_ids, , drop = FALSE]

  result <- data.frame(
    GeneID = rownames(count_matrix),
    expression_entries = NA_integer_,
    unique_brain_ids = NA_integer_,
    max_entries_from_same_person = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(count_matrix))) {
    present <- !is.na(count_matrix[i, ])
    brain_ids <- metadata$BrainID[present]
    brain_ids <- brain_ids[!is.na(brain_ids) & brain_ids != ""]
    counts_by_person <- table(brain_ids)

    result$expression_entries[i] <- sum(present)
    result$unique_brain_ids[i] <- length(unique(brain_ids))
    result$max_entries_from_same_person[i] <- if (length(counts_by_person) == 0) 0 else max(counts_by_person)
  }

  result
}

make_model_decision_summary <- function(metadata) {
  overall <- make_donor_repetition_overall(metadata)
  by_region <- make_donor_repetition_by_region(metadata)

  repeated_overall <- any(overall$repeated)
  repeated_within_region <- any(by_region$repeated_within_region)

  recommendation <- if (repeated_within_region) {
    "Repeated BrainID appears within at least one region; consider duplicate handling or a random intercept."
  } else if (repeated_overall) {
    "Repeated BrainID appears across all regions; an all-region model should use LMM/dream with (1 | BrainID). Region-wise linear models are probably OK if each region has no repeated BrainID."
  } else {
    "No repeated BrainID was detected; a linear model is probably OK for a region-wise analysis."
  }

  data.frame(
    repeated_brain_id_overall = repeated_overall,
    repeated_brain_id_within_region = repeated_within_region,
    recommendation = recommendation,
    stringsAsFactors = FALSE
  )
}

keep_highest_rin_within_region <- function(metadata) {
  metadata$dedup_order <- seq_len(nrow(metadata))
  metadata$RIN_for_order <- metadata$RIN
  metadata$RIN_for_order[is.na(metadata$RIN_for_order)] <- -Inf

  metadata <- metadata[order(metadata$RegionID, metadata$BrainID, -metadata$RIN_for_order, metadata$dedup_order), ]
  metadata <- metadata[!duplicated(paste(metadata$RegionID, metadata$BrainID, sep = "||")), ]
  metadata$dedup_order <- NULL
  metadata$RIN_for_order <- NULL
  metadata
}

run_getmm_log2 <- function(counts, gene_annotation) {
  length_kb <- gene_annotation$Length / 1000
  names(length_kb) <- gene_annotation$GeneID
  length_kb <- length_kb[rownames(counts)]

  rpk <- counts / length_kb
  rpk[!is.finite(rpk)] <- NA
  rpk[is.na(rpk)] <- 0

  dge <- edgeR::DGEList(counts = rpk)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")

  effective_library_size <- colSums(rpk) * dge$samples$norm.factors
  getmm <- t(t(rpk) / effective_library_size) * 1e6
  getmm[!is.finite(getmm)] <- NA

  normalized <- log2(getmm + 1)

  factors <- data.frame(
    SampleID = colnames(counts),
    raw_count_library_size = colSums(counts),
    rpk_library_size = colSums(rpk),
    TMM_factor = dge$samples$norm.factors,
    effective_library_size = effective_library_size,
    stringsAsFactors = FALSE
  )

  list(normalized = normalized, factors = factors)
}

process_region <- function(region_name, mode_name, metadata_mode, counts_collapsed, annotation_clean) {
  message("Processing mode = ", mode_name, " | region = ", region_name)

  region_metadata <- metadata_mode[metadata_mode$RegionID == region_name, , drop = FALSE]
  region_sample_ids <- region_metadata$SampleID
  region_counts <- counts_collapsed[, region_sample_ids, drop = FALSE]

  counts_before_filter_n <- nrow(region_counts)

  keep_by_expression <- apply(
    region_counts,
    1,
    function(x) as.numeric(stats::quantile(x, probs = low_expression_quantile, na.rm = TRUE)) > low_expression_count_cutoff
  )
  region_counts <- region_counts[keep_by_expression, , drop = FALSE]

  genes_after_expression_filter_n <- nrow(region_counts)

  region_annotation <- annotation_clean[match(rownames(region_counts), annotation_clean$GeneID), , drop = FALSE]
  rownames(region_annotation) <- region_annotation$GeneID

  normalized <- run_getmm_log2(region_counts, region_annotation)

  region_safe <- safe_file_part(region_name)
  mode_dir <- file.path(output_dir, mode_name)
  dir.create(mode_dir, showWarnings = FALSE, recursive = TRUE)

  write.csv(normalized$normalized, file.path(mode_dir, paste0("normalized_GeTMM_log2_", region_safe, ".csv")), row.names = TRUE)
  write.csv(region_counts, file.path(mode_dir, paste0("filtered_raw_counts_", region_safe, ".csv")), row.names = TRUE)
  write.csv(region_metadata, file.path(mode_dir, paste0("metadata_", region_safe, ".csv")), row.names = FALSE)
  write.csv(region_annotation, file.path(mode_dir, paste0("gene_annotation_", region_safe, ".csv")), row.names = FALSE)
  write.csv(normalized$factors, file.path(mode_dir, paste0("GeTMM_normalization_factors_", region_safe, ".csv")), row.names = FALSE)

  write.csv(make_donor_repetition_overall(region_metadata), file.path(mode_dir, paste0("region_", region_safe, "_", mode_name, "_donor_repetition_overall.csv")), row.names = FALSE)
  write.csv(make_donor_repetition_by_region(region_metadata), file.path(mode_dir, paste0("region_", region_safe, "_", mode_name, "_donor_repetition_by_region.csv")), row.names = FALSE)
  write.csv(make_gene_entry_repetition_summary(region_counts, region_metadata), file.path(mode_dir, paste0("region_", region_safe, "_", mode_name, "_gene_entry_repetition_summary.csv")), row.names = FALSE)

  data.frame(
    mode = mode_name,
    region = region_name,
    samples = nrow(region_metadata),
    unique_brain_ids = length(unique(region_metadata$BrainID)),
    genes_before_region_expression_filter = counts_before_filter_n,
    genes_after_region_expression_filter = genes_after_expression_filter_n,
    genes_removed_by_region_expression_filter = counts_before_filter_n - genes_after_expression_filter_n,
    stringsAsFactors = FALSE
  )
}

## -----------------------------
## 4. Read counts and annotation
## -----------------------------

message("Reading raw counts: ", raw_counts_file)
raw_counts_df <- read.delim(raw_counts_file, check.names = FALSE)

if (!"GeneID" %in% names(raw_counts_df)) {
  stop("Raw counts file must contain a GeneID column.")
}

raw_counts_df$GeneID <- as.character(raw_counts_df$GeneID)
raw_counts_df <- raw_counts_df[!is.na(raw_counts_df$GeneID) & trimws(raw_counts_df$GeneID) != "", , drop = FALSE]
message("Raw count rows after removing missing/blank GeneID: ", nrow(raw_counts_df))

sample_columns <- setdiff(names(raw_counts_df), "GeneID")
raw_counts_df[sample_columns] <- lapply(raw_counts_df[sample_columns], function(x) suppressWarnings(as.numeric(x)))

duplicate_gene_ids <- raw_counts_df$GeneID[duplicated(raw_counts_df$GeneID) | duplicated(raw_counts_df$GeneID, fromLast = TRUE)]
if (length(duplicate_gene_ids) == 0) {
  duplicate_report <- data.frame(
    GeneID = character(),
    number_of_rows_in_raw_counts = integer(),
    stringsAsFactors = FALSE
  )
} else {
  duplicate_report <- data.frame(table(duplicate_gene_ids), stringsAsFactors = FALSE)
  names(duplicate_report) <- c("GeneID", "number_of_rows_in_raw_counts")
  duplicate_report <- duplicate_report[order(-duplicate_report$number_of_rows_in_raw_counts, duplicate_report$GeneID), ]
}
write_report(duplicate_report, "duplicate_gene_ids_in_raw_counts.csv")

message("Collapsing duplicated GeneID rows by summing raw counts...")
counts_collapsed_df <- aggregate(. ~ GeneID, data = raw_counts_df, FUN = sum, na.rm = TRUE)
rownames(counts_collapsed_df) <- counts_collapsed_df$GeneID
counts_collapsed <- as.matrix(counts_collapsed_df[, sample_columns, drop = FALSE])
mode(counts_collapsed) <- "numeric"

message("Reading gene annotation: ", annotation_file)
annotation <- read.delim(annotation_file, check.names = FALSE)
if (!all(c("GeneID", "Length") %in% names(annotation))) {
  stop("Annotation file must contain GeneID and Length columns.")
}

annotation$GeneID <- as.character(annotation$GeneID)
annotation$Length <- suppressWarnings(as.numeric(annotation$Length))
annotation_clean <- annotation[
  !is.na(annotation$GeneID) &
    trimws(annotation$GeneID) != "" &
    !is.na(annotation$Length) &
    annotation$Length > minimum_gene_length_bp,
  ,
  drop = FALSE
]

annotation_clean <- annotation_clean[!duplicated(annotation_clean$GeneID), , drop = FALSE]

common_genes <- intersect(rownames(counts_collapsed), annotation_clean$GeneID)
counts_collapsed <- counts_collapsed[common_genes, , drop = FALSE]
annotation_clean <- annotation_clean[match(common_genes, annotation_clean$GeneID), , drop = FALSE]
rownames(annotation_clean) <- annotation_clean$GeneID

message("Genes matched to annotation and Length > ", minimum_gene_length_bp, " bp: ", length(common_genes))

## -----------------------------
## 5. Load and match metadata
## -----------------------------

metadata <- load_metadata()
message("Metadata rows loaded: ", nrow(metadata))

matched_samples <- intersect(colnames(counts_collapsed), metadata$SampleID)
if (length(matched_samples) == 0) {
  stop("No samples matched between raw counts and metadata.")
}

counts_collapsed <- counts_collapsed[, matched_samples, drop = FALSE]
metadata_matched_all <- metadata[matched_samples, , drop = FALSE]
metadata_matched_all <- metadata_matched_all[!is.na(metadata_matched_all$RegionID) & metadata_matched_all$RegionID != "", , drop = FALSE]
counts_collapsed <- counts_collapsed[, metadata_matched_all$SampleID, drop = FALSE]

write.csv(metadata_matched_all, file.path(output_dir, "metadata_clean_matched_all_samples.csv"), row.names = FALSE)
message("Matched samples with non-missing region: ", nrow(metadata_matched_all))

## -----------------------------
## 6. Reports before deduplication
## -----------------------------

write_report(make_donor_repetition_overall(metadata_matched_all), "all_matched_samples_before_region_dedup_donor_repetition_overall.csv")
write_report(make_donor_repetition_by_region(metadata_matched_all), "all_matched_samples_before_region_dedup_donor_repetition_by_region.csv")
write_report(make_gene_entry_repetition_summary(counts_collapsed, metadata_matched_all), "all_matched_samples_before_region_dedup_gene_entry_repetition_summary.csv")
write_report(make_model_decision_summary(metadata_matched_all), "model_decision_repetition_summary.csv")

## -----------------------------
## 7. Run the two requested sample-handling modes
## -----------------------------

modes <- list(
  keep_all_samples = metadata_matched_all,
  keep_highest_RIN_within_region = keep_highest_rin_within_region(metadata_matched_all)
)

all_qc <- data.frame()

for (mode_name in names(modes)) {
  message("Running sample mode: ", mode_name)
  metadata_mode <- modes[[mode_name]]
  counts_mode <- counts_collapsed[, metadata_mode$SampleID, drop = FALSE]

  mode_qc <- data.frame()
  regions <- sort(unique(metadata_mode$RegionID))

  for (region_name in regions) {
    region_qc <- process_region(
      region_name = region_name,
      mode_name = mode_name,
      metadata_mode = metadata_mode,
      counts_collapsed = counts_mode,
      annotation_clean = annotation_clean
    )
    mode_qc <- rbind(mode_qc, region_qc)
  }

  mode_qc_file <- file.path(output_dir, paste0("Step1_QC_summary_by_region_", mode_name, ".csv"))
  write.csv(mode_qc, mode_qc_file, row.names = FALSE)
  message("Wrote: ", mode_qc_file)

  all_qc <- rbind(all_qc, mode_qc)
}

all_qc_file <- file.path(output_dir, "Step1_QC_summary_by_region_all_modes.csv")
write.csv(all_qc, all_qc_file, row.names = FALSE)
message("Wrote: ", all_qc_file)

message("Step 1 complete. Outputs are in: ", output_dir)
