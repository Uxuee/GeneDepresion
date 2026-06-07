############################################################
## Step 1 — GSE102556 preprocessing + GeTMM by broad brain areas
## Project: depression biomarker discovery from brain RNA-seq
##
## Why this version exists:
##   The raw GEO metadata contains specific tissue labels such as:
##     - Orbitofrontal (OFC; BA11)
##     - Dorsolateral prefrontal cortex (dlPFC; BA8/9)
##     - Cingulate gyrus 25 (Cg25)
##     - Anterior Insula (aINS)
##     - Nucleus Accumbens (Nac)
##     - Subiculum (Sub)
##
##   Elisa-style scripts usually split the transcriptomic data into broader
##   brain areas before normalization, instead of treating every tiny anatomical
##   label as a fully independent analysis unit.
##
## This script therefore creates:
##   - OriginalRegion: the exact GEO tissue label
##   - BrainArea: a broader user-editable grouping, e.g. Frontal_Cortex
##   - MajorBrainPart: an even coarser grouping, e.g. Cerebral_Cortex
##   - AnalysisGroup: whichever grouping column you choose below
##
## Two sample-handling versions are produced:
##   A) keep_all_samples
##   B) keep_highest_RIN_within_analysis_group
##
## Main normalization:
##   filtered raw counts -> GeTMM -> log2(GeTMM + 1)
##
## GeTMM handles:
##   - gene length: counts / gene_length_kb
##   - sequencing depth: length-adjusted library total
##   - RNA composition: TMM normalization factors
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

############################
## 0. User configuration
############################

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(script_path)) {
  script_path <- normalizePath("code/Step1_GSE102556_GeTMM_brain_area_groups.R", mustWork = FALSE)
}
repo_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
data_dir <- "D:/genes/genedata"

output_dir <- file.path(data_dir, "Step1_GSE102556_GeTMM_brain_area_groups")
clean_output_dir_before_run <- TRUE

if (clean_output_dir_before_run && dir.exists(output_dir)) {
  message("Cleaning previous Step 1 output folder: ", output_dir)
  unlink(output_dir, recursive = TRUE, force = TRUE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_paths <- list(
  config = file.path(output_dir, "00_config"),
  metadata = file.path(output_dir, "01_metadata"),
  normalized = file.path(output_dir, "02_normalized_expression"),
  filtered_counts = file.path(output_dir, "03_filtered_raw_counts"),
  normalization_factors = file.path(output_dir, "04_normalization_factors"),
  gene_annotation = file.path(output_dir, "05_gene_annotation"),
  repetition_reports = file.path(output_dir, "06_repetition_reports"),
  qc = file.path(output_dir, "07_qc"),
  rds = file.path(output_dir, "08_rds")
)
invisible(lapply(output_paths, dir.create, recursive = TRUE, showWarnings = FALSE))

raw_counts_path <- file.path(data_dir, "GSE102556_raw_counts_GRCh38.p13_NCBI.tsv")
annot_path      <- file.path(data_dir, "Human.GRCh38.p13.annot.tsv")

## The script first looks for these local metadata files.
## If neither exists, it tries to fetch metadata from GEO using GEOquery.
metadata_candidates <- c(
  file.path(data_dir, "GSE102556_metadata_clean.csv"),
  file.path(data_dir, "metadata_clean_matched_all_samples.csv"),
  file.path(data_dir, "Step1_GeTMM_log2_outputs", "metadata_clean_matched_all_samples.csv")
)

## Editable mapping file. If it does not exist, the script creates a default one.
brain_area_mapping_path <- file.path(repo_dir, "GSE102556_brain_area_mapping_template.csv")

## Choose the grouping used for normalization/output:
##   "BrainArea"      = recommended first version, moderately broad
##   "MajorBrainPart" = coarser version
##
## You can edit GSE102556_brain_area_mapping_template.csv before running
## if you want a different grouping.
analysis_group_column <- "BrainArea"

## Keep NULL to process every available broad brain area.
target_analysis_groups <- NULL

## Ferrari/thesis-inspired filters:
##   - remove genes too short
##   - remove genes with < min_count counts in at least low_expr_exclusion_prop samples
## This keeps genes whose 80th percentile count is > 10.
min_gene_length_bp <- 20
min_count <- 10
low_expr_exclusion_prop <- 0.80

## Two versions. Do not edit unless you only want one.
sample_handling_modes <- c(
  "keep_all_samples",
  "keep_highest_RIN_within_analysis_group"
)

############################
## 1. Package helpers
############################

install_if_missing_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

install_if_missing_bioc <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

install_if_missing_bioc("edgeR")
install_if_missing_bioc("Biobase")
install_if_missing_bioc("GEOquery")

library(edgeR)
library(GEOquery)
library(Biobase)

############################
## 2. Small helper functions
############################

read_expression_matrix <- function(path) {
  if (!file.exists(path)) stop("Raw count file not found: ", path)
  message("Reading raw count file: ", path)

  x <- read.delim(path, header = TRUE, sep = "\t", check.names = FALSE)
  gene_id <- as.character(x[[1]])

  mat <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- gene_id
  mat
}

slugify <- function(x) {
  x <- as.character(x)
  x <- gsub("\\+", "plus", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x <- tolower(x)
  ifelse(nchar(x) == 0, "unknown", x)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

first_present_col <- function(df, candidates) {
  for (nm in candidates) {
    if (nm %in% colnames(df)) return(nm)
  }
  NA_character_
}

coalesce_nonempty <- function(x, fallback) {
  x <- as.character(x)
  fallback <- as.character(fallback)
  bad <- is.na(x) | trimws(x) == ""
  x[bad] <- fallback[bad]
  x
}

extract_characteristic <- function(pdata, key) {
  char_cols <- grep("^characteristics", colnames(pdata), value = TRUE, ignore.case = TRUE)
  out <- rep(NA_character_, nrow(pdata))

  for (i in seq_len(nrow(pdata))) {
    vals <- as.character(unlist(pdata[i, char_cols, drop = FALSE]))
    vals <- vals[!is.na(vals)]
    hit <- vals[grepl(paste0("^", key, "\\s*:"), vals, ignore.case = TRUE)]
    if (length(hit) > 0) {
      out[i] <- sub("^[^:]+:\\s*", "", hit[1])
    }
  }
  out
}

extract_first_available <- function(pdata, keys) {
  out <- rep(NA_character_, nrow(pdata))
  for (key in keys) {
    tmp <- extract_characteristic(pdata, key)
    use <- (is.na(out) | trimws(out) == "") & !is.na(tmp) & trimws(tmp) != ""
    out[use] <- tmp[use]
  }
  out
}

collapse_duplicate_gene_ids <- function(mat, out_dir) {
  duplicated_ids <- unique(rownames(mat)[duplicated(rownames(mat))])
  if (length(duplicated_ids) == 0) {
    duplicate_report <- data.frame(
      GeneID = character(),
      n_rows = integer(),
      stringsAsFactors = FALSE
    )
    write.csv(
      duplicate_report,
      file.path(out_dir, "duplicate_gene_ids_in_raw_counts.csv"),
      row.names = FALSE
    )
    return(mat)
  }

  message("Found duplicated GeneID rows: ", length(duplicated_ids))
  duplicate_report <- data.frame(
    GeneID = duplicated_ids,
    n_rows = as.integer(table(rownames(mat))[duplicated_ids]),
    stringsAsFactors = FALSE
  )
  write.csv(
    duplicate_report,
    file.path(out_dir, "duplicate_gene_ids_in_raw_counts.csv"),
    row.names = FALSE
  )

  ## Collapse duplicated GeneID rows by summing raw counts.
  rowsum(mat, group = rownames(mat), reorder = FALSE, na.rm = TRUE)
}

load_or_fetch_metadata <- function(sample_ids) {
  existing <- metadata_candidates[file.exists(metadata_candidates)]
  if (length(existing) > 0) {
    message("Reading local metadata file: ", existing[1])
    meta <- read.csv(existing[1], check.names = FALSE)
  } else {
    message("No local clean metadata found. Fetching GSE102556 metadata with GEOquery...")
    gse_list <- GEOquery::getGEO("GSE102556", GSEMatrix = TRUE)
    gse <- gse_list[[1]]
    pdata <- Biobase::pData(gse)
    meta <- pdata
  }

  if (!("SampleID" %in% colnames(meta))) {
    if ("geo_accession" %in% colnames(meta)) {
      meta$SampleID <- as.character(meta$geo_accession)
    } else {
      meta$SampleID <- rownames(meta)
    }
  }

  ## If we got raw GEO pData, extract useful clean columns.
  if (!("tissue" %in% colnames(meta)) && any(grepl("^characteristics", colnames(meta), ignore.case = TRUE))) {
    meta$tissue <- extract_first_available(meta, c("tissue", "brain region", "region"))
  }
  if (!("gender" %in% colnames(meta)) && any(grepl("^characteristics", colnames(meta), ignore.case = TRUE))) {
    meta$gender <- extract_first_available(meta, c("gender", "sex"))
  }
  if (!("phenotype" %in% colnames(meta)) && any(grepl("^characteristics", colnames(meta), ignore.case = TRUE))) {
    meta$phenotype <- extract_first_available(meta, c("phenotype", "diagnosis"))
  }
  if (!("age" %in% colnames(meta)) && any(grepl("^characteristics", colnames(meta), ignore.case = TRUE))) {
    meta$age <- extract_first_available(meta, c("age"))
  }
  if (!("pmi" %in% colnames(meta)) && any(grepl("^characteristics", colnames(meta), ignore.case = TRUE))) {
    meta$pmi <- extract_first_available(meta, c("pmi"))
  }
  if (!("ph" %in% colnames(meta)) && any(grepl("^characteristics", colnames(meta), ignore.case = TRUE))) {
    meta$ph <- extract_first_available(meta, c("ph", "pH"))
  }
  if (!("rin" %in% colnames(meta)) && any(grepl("^characteristics", colnames(meta), ignore.case = TRUE))) {
    meta$rin <- extract_first_available(meta, c("rin", "RIN"))
  }
  if (!("alcool" %in% colnames(meta)) && any(grepl("^characteristics", colnames(meta), ignore.case = TRUE))) {
    meta$alcool <- extract_first_available(meta, c("alcool", "alcohol"))
  }
  if (!("medication" %in% colnames(meta)) && any(grepl("^characteristics", colnames(meta), ignore.case = TRUE))) {
    meta$medication <- extract_first_available(meta, c("medication"))
  }
  if (!("drugs" %in% colnames(meta)) && any(grepl("^characteristics", colnames(meta), ignore.case = TRUE))) {
    meta$drugs <- extract_first_available(meta, c("drugs"))
  }
  if (!("smoking" %in% colnames(meta)) && any(grepl("^characteristics", colnames(meta), ignore.case = TRUE))) {
    meta$smoking <- extract_first_available(meta, c("smoking"))
  }

  ## Keep only samples in the count matrix.
  meta <- meta[meta$SampleID %in% sample_ids, , drop = FALSE]
  meta
}

prepare_metadata <- function(meta) {
  if (!("title" %in% colnames(meta)) && "Title" %in% colnames(meta)) {
    meta$title <- meta$Title
  }

  ## Original fine anatomical region.
  region_col <- first_present_col(meta, c("OriginalRegion", "RegionID", "Region", "tissue", "tissue_ch1"))
  if (is.na(region_col)) {
    stop("Could not find a tissue/region column in metadata.")
  }
  meta$OriginalRegion <- as.character(meta[[region_col]])
  meta$RegionID <- meta$OriginalRegion

  ## Clean diagnosis/sex/covariate aliases.
  if (!("Diagnosis" %in% colnames(meta))) {
    pheno_col <- first_present_col(meta, c("phenotype", "phenotype_ch1", "diagnosis", "Diagnosis"))
    meta$Diagnosis <- if (!is.na(pheno_col)) as.character(meta[[pheno_col]]) else NA_character_
  }
  if (!("Sex" %in% colnames(meta))) {
    sex_col <- first_present_col(meta, c("gender", "gender_ch1", "sex", "Sex"))
    meta$Sex <- if (!is.na(sex_col)) as.character(meta[[sex_col]]) else NA_character_
  }
  if (!("Age" %in% colnames(meta))) {
    age_col <- first_present_col(meta, c("age", "age_ch1", "Age"))
    meta$Age <- if (!is.na(age_col)) safe_num(meta[[age_col]]) else NA_real_
  }
  if (!("PMI" %in% colnames(meta))) {
    pmi_col <- first_present_col(meta, c("pmi", "pmi_ch1", "PMI"))
    meta$PMI <- if (!is.na(pmi_col)) safe_num(meta[[pmi_col]]) else NA_real_
  }
  if (!("pH" %in% colnames(meta))) {
    ph_col <- first_present_col(meta, c("ph", "ph_ch1", "pH"))
    meta$pH <- if (!is.na(ph_col)) safe_num(meta[[ph_col]]) else NA_real_
  }
  if (!("RIN" %in% colnames(meta))) {
    rin_col <- first_present_col(meta, c("RIN", "rin", "rin_ch1"))
    meta$RIN <- if (!is.na(rin_col)) safe_num(meta[[rin_col]]) else NA_real_
  }

  ## BrainID/donor ID: usually the number before ":" in the GEO title.
  if (!("BrainID" %in% colnames(meta))) {
    if ("title" %in% colnames(meta)) {
      meta$BrainID <- sub(":.*$", "", as.character(meta$title))
    } else {
      meta$BrainID <- NA_character_
    }
  }
  meta$BrainID <- as.character(meta$BrainID)
  meta$BrainID[is.na(meta$BrainID) | trimws(meta$BrainID) == ""] <- paste0("unknown_", meta$SampleID[is.na(meta$BrainID) | trimws(meta$BrainID) == ""])

  meta
}

default_brain_area_mapping <- function(regions) {
  brain_area <- rep("Other", length(regions))
  major <- rep("Other", length(regions))
  comment <- rep("Default fallback mapping. Please check manually.", length(regions))

  is_frontal <- grepl("Orbitofrontal|Dorsolateral|prefrontal|Cingulate|Cg25|BA11|BA8|BA9|BA25", regions, ignore.case = TRUE)
  brain_area[is_frontal] <- "Frontal_Cortex"
  major[is_frontal] <- "Cerebral_Cortex"
  comment[is_frontal] <- "Default: OFC/dlPFC/Cg25 grouped as frontal cortical samples."

  is_insula <- grepl("Insula|aINS", regions, ignore.case = TRUE)
  brain_area[is_insula] <- "Insula"
  major[is_insula] <- "Cerebral_Cortex"
  comment[is_insula] <- "Default: anterior insula kept as its own broad cortical area."

  is_nacc <- grepl("Nucleus Accumbens|Accumbens|Nac", regions, ignore.case = TRUE)
  brain_area[is_nacc] <- "Ventral_Striatum"
  major[is_nacc] <- "Subcortical_Striatum"
  comment[is_nacc] <- "Default: nucleus accumbens grouped as ventral striatum/subcortical."

  is_sub <- grepl("Subiculum|\\bSub\\b", regions, ignore.case = TRUE)
  brain_area[is_sub] <- "Hippocampal_Formation"
  major[is_sub] <- "Temporal_Limbic_Hippocampal"
  comment[is_sub] <- "Default: subiculum grouped as hippocampal formation / temporal-limbic."

  data.frame(
    OriginalRegion = regions,
    BrainArea = brain_area,
    MajorBrainPart = major,
    MappingComment = comment,
    stringsAsFactors = FALSE
  )
}

load_or_create_brain_area_mapping <- function(original_regions) {
  if (file.exists(brain_area_mapping_path)) {
    message("Reading brain-area mapping file: ", brain_area_mapping_path)
    mapping <- read.csv(brain_area_mapping_path, check.names = FALSE)
  } else {
    message("Creating default brain-area mapping file: ", brain_area_mapping_path)
    mapping <- default_brain_area_mapping(sort(unique(original_regions)))
    write.csv(mapping, brain_area_mapping_path, row.names = FALSE)
  }

  required <- c("OriginalRegion", "BrainArea", "MajorBrainPart")
  missing <- setdiff(required, colnames(mapping))
  if (length(missing) > 0) {
    stop("Brain-area mapping file is missing columns: ", paste(missing, collapse = ", "))
  }

  mapping$OriginalRegion <- as.character(mapping$OriginalRegion)
  mapping$BrainArea <- as.character(mapping$BrainArea)
  mapping$MajorBrainPart <- as.character(mapping$MajorBrainPart)
  mapping
}

apply_brain_area_mapping <- function(meta, mapping) {
  idx <- match(meta$OriginalRegion, mapping$OriginalRegion)
  meta$BrainArea <- mapping$BrainArea[idx]
  meta$MajorBrainPart <- mapping$MajorBrainPart[idx]

  meta$BrainArea[is.na(meta$BrainArea) | trimws(meta$BrainArea) == ""] <- "Unmapped"
  meta$MajorBrainPart[is.na(meta$MajorBrainPart) | trimws(meta$MajorBrainPart) == ""] <- "Unmapped"

  if (!(analysis_group_column %in% colnames(meta))) {
    stop("analysis_group_column='", analysis_group_column, "' is not present in metadata after mapping.")
  }

  meta$AnalysisGroup <- as.character(meta[[analysis_group_column]])
  meta$AnalysisGroup[is.na(meta$AnalysisGroup) | trimws(meta$AnalysisGroup) == ""] <- "Unmapped"
  meta
}

make_group_composition_report <- function(meta, out_dir) {
  groups <- sort(unique(meta$AnalysisGroup))
  rows <- lapply(groups, function(g) {
    m <- meta[meta$AnalysisGroup == g, , drop = FALSE]
    tab <- sort(table(m$OriginalRegion), decreasing = TRUE)
    data.frame(
      AnalysisGroup = g,
      GroupingColumn = analysis_group_column,
      n_samples = nrow(m),
      n_unique_BrainID = length(unique(m$BrainID)),
      n_original_regions = length(unique(m$OriginalRegion)),
      original_regions = paste(names(tab), collapse = " | "),
      original_region_counts = paste(paste(names(tab), as.integer(tab), sep = "="), collapse = " | "),
      stringsAsFactors = FALSE
    )
  })
  report <- do.call(rbind, rows)
  write.csv(report, file.path(out_dir, "analysis_group_composition_report.csv"), row.names = FALSE)
  report
}

keep_highest_rin_within_analysis_group <- function(meta) {
  if (nrow(meta) == 0) return(meta)

  meta$.original_order <- seq_len(nrow(meta))
  rin_for_order <- meta$RIN
  rin_for_order[is.na(rin_for_order)] <- -Inf
  meta$.rin_for_order <- rin_for_order

  ## One sample per BrainID inside each broad analysis group.
  meta <- meta[order(meta$AnalysisGroup, meta$BrainID, -meta$.rin_for_order, meta$.original_order), , drop = FALSE]
  meta <- meta[!duplicated(paste(meta$AnalysisGroup, meta$BrainID, sep = "||")), , drop = FALSE]

  meta$.original_order <- NULL
  meta$.rin_for_order <- NULL
  meta
}

filter_genes_for_group <- function(counts, gene_lengths_bp) {
  gene_lengths_group <- gene_lengths_bp[rownames(counts)]

  keep_valid_length <- !is.na(gene_lengths_group) & gene_lengths_group > min_gene_length_bp

  ## Exclude genes with <10 counts in >=80% of samples.
  ## Equivalent: keep if 80th percentile count is >10.
  q80 <- apply(counts, 1, quantile, probs = low_expr_exclusion_prop, na.rm = TRUE)
  keep_expression <- q80 > min_count

  keep <- keep_valid_length & keep_expression
  counts[keep, , drop = FALSE]
}

compute_getmm <- function(counts, gene_lengths_bp) {
  common <- intersect(rownames(counts), names(gene_lengths_bp))
  counts <- counts[common, , drop = FALSE]
  lengths_bp <- gene_lengths_bp[common]
  keep <- !is.na(lengths_bp) & lengths_bp > 0

  counts <- counts[keep, , drop = FALSE]
  lengths_kb <- lengths_bp[keep] / 1000

  ## 1) Gene-length correction.
  rpk <- sweep(counts, 1, lengths_kb, "/")
  rpk[!is.finite(rpk)] <- 0

  ## 2) RNA composition correction.
  dge <- edgeR::DGEList(counts = rpk)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")

  ## 3) Sequencing-depth / effective library scaling.
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

write_repetition_reports <- function(metadata_subset, genes, out_dir, prefix) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  meta <- metadata_subset
  if (nrow(meta) == 0) return(NULL)

  ## Overall donor repetition table for this subset.
  brain_ids <- unique(meta$BrainID)
  donor_overall <- do.call(rbind, lapply(brain_ids, function(id) {
    m <- meta[meta$BrainID == id, , drop = FALSE]
    data.frame(
      BrainID = id,
      n_samples_total = nrow(m),
      n_analysis_groups = length(unique(m$AnalysisGroup)),
      analysis_groups = paste(unique(m$AnalysisGroup), collapse = " | "),
      n_original_regions = length(unique(m$OriginalRegion)),
      original_regions = paste(unique(m$OriginalRegion), collapse = " | "),
      sample_ids = paste(m$SampleID, collapse = " | "),
      RIN_values = paste(m$RIN, collapse = " | "),
      repeated_person = nrow(m) > 1,
      stringsAsFactors = FALSE
    )
  }))
  write.csv(
    donor_overall,
    file.path(out_dir, paste0(prefix, "_donor_repetition_overall.csv")),
    row.names = FALSE
  )

  ## Donor-by-analysis-group repetition table.
  combos <- unique(meta[, c("BrainID", "AnalysisGroup"), drop = FALSE])
  donor_group <- do.call(rbind, lapply(seq_len(nrow(combos)), function(i) {
    id <- combos$BrainID[i]
    grp <- combos$AnalysisGroup[i]
    m <- meta[meta$BrainID == id & meta$AnalysisGroup == grp, , drop = FALSE]
    rin_for_order <- m$RIN
    rin_for_order[is.na(rin_for_order)] <- -Inf
    best_idx <- which.max(rin_for_order)

    data.frame(
      BrainID = id,
      AnalysisGroup = grp,
      n_samples_same_person_same_analysis_group = nrow(m),
      n_original_regions_inside_group = length(unique(m$OriginalRegion)),
      original_regions_inside_group = paste(unique(m$OriginalRegion), collapse = " | "),
      sample_ids = paste(m$SampleID, collapse = " | "),
      RIN_values = paste(m$RIN, collapse = " | "),
      highest_RIN_sample = if (length(best_idx) > 0) m$SampleID[best_idx] else NA_character_,
      highest_RIN = if (length(best_idx) > 0) m$RIN[best_idx] else NA_real_,
      n_samples_that_would_be_dropped_by_highest_RIN_rule = max(0, nrow(m) - 1),
      repeated_within_analysis_group = nrow(m) > 1,
      stringsAsFactors = FALSE
    )
  }))
  write.csv(
    donor_group,
    file.path(out_dir, paste0(prefix, "_donor_repetition_by_analysis_group.csv")),
    row.names = FALSE
  )

  ## Donor-by-original-region repetition table.
  combos2 <- unique(meta[, c("BrainID", "OriginalRegion"), drop = FALSE])
  donor_original_region <- do.call(rbind, lapply(seq_len(nrow(combos2)), function(i) {
    id <- combos2$BrainID[i]
    reg <- combos2$OriginalRegion[i]
    m <- meta[meta$BrainID == id & meta$OriginalRegion == reg, , drop = FALSE]
    data.frame(
      BrainID = id,
      OriginalRegion = reg,
      n_samples_same_person_same_original_region = nrow(m),
      sample_ids = paste(m$SampleID, collapse = " | "),
      RIN_values = paste(m$RIN, collapse = " | "),
      repeated_within_original_region = nrow(m) > 1,
      stringsAsFactors = FALSE
    )
  }))
  write.csv(
    donor_original_region,
    file.path(out_dir, paste0(prefix, "_donor_repetition_by_original_region.csv")),
    row.names = FALSE
  )

  ## Gene-level repetition summary. This answers:
  ## "For one gene, how many expression entries can come from the same person?"
  donor_counts <- table(meta$BrainID)
  n_repeated_donors <- sum(donor_counts > 1)
  max_entries_same_person <- max(as.integer(donor_counts))

  gene_entry <- data.frame(
    GeneID = genes,
    n_expression_entries = nrow(meta),
    n_unique_people_or_BrainID = length(donor_counts),
    n_people_with_multiple_entries = n_repeated_donors,
    max_entries_from_same_person_for_this_gene = max_entries_same_person,
    repeated_measure_present_for_this_gene = max_entries_same_person > 1,
    grouping_column = analysis_group_column,
    analysis_groups_present = paste(unique(meta$AnalysisGroup), collapse = " | "),
    original_regions_present = paste(unique(meta$OriginalRegion), collapse = " | "),
    stringsAsFactors = FALSE
  )
  write.csv(
    gene_entry,
    file.path(out_dir, paste0(prefix, "_gene_entry_repetition_summary.csv")),
    row.names = FALSE
  )

  ## Compact one-row version for quick inspection.
  gene_entry_compact <- unique(gene_entry[, setdiff(colnames(gene_entry), "GeneID"), drop = FALSE])
  write.csv(
    gene_entry_compact,
    file.path(out_dir, paste0(prefix, "_gene_entry_repetition_summary_compact.csv")),
    row.names = FALSE
  )

  invisible(list(
    donor_overall = donor_overall,
    donor_group = donor_group,
    donor_original_region = donor_original_region,
    gene_entry = gene_entry
  ))
}

process_analysis_group <- function(group_name, mode, raw_counts, metadata, gene_lengths_bp, annotation) {
  message("\nProcessing analysis group: ", group_name, " | mode: ", mode)

  group_meta_before <- metadata[metadata$AnalysisGroup == group_name, , drop = FALSE]
  n_samples_before <- nrow(group_meta_before)

  if (mode == "keep_highest_RIN_within_analysis_group") {
    group_meta <- keep_highest_rin_within_analysis_group(group_meta_before)
  } else if (mode == "keep_all_samples") {
    group_meta <- group_meta_before
  } else {
    stop("Unknown sample-handling mode: ", mode)
  }

  common_samples <- intersect(colnames(raw_counts), group_meta$SampleID)
  group_meta <- group_meta[match(common_samples, group_meta$SampleID), , drop = FALSE]
  counts_group <- raw_counts[, common_samples, drop = FALSE]

  counts_filtered <- filter_genes_for_group(counts_group, gene_lengths_bp)
  norm <- compute_getmm(counts_filtered, gene_lengths_bp)
  expr_norm <- norm$log2_getmm

  annot_group <- annotation[match(rownames(expr_norm), annotation$GeneID), , drop = FALSE]

  tag <- slugify(group_name)
  typed_dirs <- list(
    normalized = file.path(output_paths$normalized, mode),
    filtered_counts = file.path(output_paths$filtered_counts, mode),
    metadata = file.path(output_paths$metadata, mode),
    gene_annotation = file.path(output_paths$gene_annotation, mode),
    normalization_factors = file.path(output_paths$normalization_factors, mode),
    repetition_reports = file.path(output_paths$repetition_reports, mode, tag)
  )
  invisible(lapply(typed_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

  write.csv(expr_norm,
            file.path(typed_dirs$normalized, paste0("normalized_GeTMM_log2_", tag, ".csv")),
            row.names = TRUE)
  write.csv(counts_filtered,
            file.path(typed_dirs$filtered_counts, paste0("filtered_raw_counts_", tag, ".csv")),
            row.names = TRUE)
  write.csv(group_meta,
            file.path(typed_dirs$metadata, paste0("metadata_", tag, ".csv")),
            row.names = FALSE)
  write.csv(annot_group,
            file.path(typed_dirs$gene_annotation, paste0("gene_annotation_", tag, ".csv")),
            row.names = FALSE)
  write.csv(norm$norm_factors,
            file.path(typed_dirs$normalization_factors, paste0("GeTMM_normalization_factors_", tag, ".csv")),
            row.names = FALSE)

  write_repetition_reports(
    metadata_subset = group_meta,
    genes = rownames(expr_norm),
    out_dir = typed_dirs$repetition_reports,
    prefix = paste0("analysis_group_", tag, "_", mode)
  )

  donor_counts_group <- table(group_meta$BrainID)
  original_region_counts <- table(group_meta$OriginalRegion)

  data.frame(
    Mode = mode,
    GroupingColumn = analysis_group_column,
    AnalysisGroup = group_name,
    n_samples_before_mode_filter = n_samples_before,
    n_samples_after_mode_filter = nrow(group_meta),
    n_samples_dropped_by_highest_RIN_rule = n_samples_before - nrow(group_meta),
    n_unique_BrainID = length(unique(group_meta$BrainID)),
    n_BrainID_repeated_within_analysis_group = sum(donor_counts_group > 1),
    max_samples_per_BrainID_within_analysis_group = if (length(donor_counts_group) > 0) max(as.integer(donor_counts_group)) else NA_integer_,
    n_original_regions_in_group = length(unique(group_meta$OriginalRegion)),
    original_regions_in_group = paste(names(original_region_counts), collapse = " | "),
    original_region_counts = paste(paste(names(original_region_counts), as.integer(original_region_counts), sep = "="), collapse = " | "),
    n_genes_before_filter = nrow(counts_group),
    n_genes_after_length_and_expression_filter = nrow(counts_filtered),
    n_CTRL = if ("Diagnosis" %in% colnames(group_meta)) sum(group_meta$Diagnosis == "CTRL", na.rm = TRUE) else NA_integer_,
    n_MDD  = if ("Diagnosis" %in% colnames(group_meta)) sum(group_meta$Diagnosis == "MDD",  na.rm = TRUE) else NA_integer_,
    n_missing_RIN = sum(is.na(group_meta$RIN)),
    n_missing_Age = sum(is.na(group_meta$Age)),
    n_missing_Sex = sum(is.na(group_meta$Sex) | trimws(group_meta$Sex) == ""),
    n_missing_Diagnosis = sum(is.na(group_meta$Diagnosis) | trimws(group_meta$Diagnosis) == ""),
    min_gene_length_bp = min_gene_length_bp,
    min_count = min_count,
    low_expr_exclusion_prop = low_expr_exclusion_prop,
    normalization = "GeTMM_log2",
    stringsAsFactors = FALSE
  )
}

make_model_decision_summary <- function(metadata, out_dir) {
  donor_counts_all <- table(metadata$BrainID)
  donor_analysis_group_counts <- table(paste(metadata$BrainID, metadata$AnalysisGroup, sep = "__"))
  donor_original_region_counts <- table(paste(metadata$BrainID, metadata$OriginalRegion, sep = "__"))

  summary <- data.frame(
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
    Suggested_first_model = c(
      ifelse(sum(donor_original_region_counts > 1) == 0,
             "linear model is probably OK inside each original fine region",
             "consider duplicate handling or random intercept inside original fine region"),
      ifelse(sum(donor_analysis_group_counts > 1) == 0,
             "linear model is probably OK inside each broad brain area",
             "if using keep_all_samples, consider LMM/dream with (1 | BrainID), or use highest-RIN version"),
      ifelse(sum(donor_counts_all > 1) == 0,
             "linear model is probably OK",
             "LMM/dream with (1 | BrainID) is advisable if combining all groups")
    ),
    stringsAsFactors = FALSE
  )

  write.csv(summary, file.path(out_dir, "model_decision_repetition_summary.csv"), row.names = FALSE)
  summary
}

############################
## 3. Load raw counts + annotation
############################

raw_counts <- read_expression_matrix(raw_counts_path)
raw_counts <- collapse_duplicate_gene_ids(raw_counts, output_paths$qc)

if (!file.exists(annot_path)) stop("Annotation file not found: ", annot_path)
annotation <- read.delim(annot_path, header = TRUE, sep = "\t", check.names = FALSE)

if (!("GeneID" %in% colnames(annotation))) stop("Annotation file must contain GeneID.")
if (!("Length" %in% colnames(annotation))) stop("Annotation file must contain Length.")

annotation$GeneID <- as.character(annotation$GeneID)
annotation$Length <- safe_num(annotation$Length)

gene_lengths_bp <- setNames(annotation$Length, annotation$GeneID)

common_genes <- intersect(rownames(raw_counts), annotation$GeneID)
raw_counts <- raw_counts[common_genes, , drop = FALSE]
annotation <- annotation[match(common_genes, annotation$GeneID), , drop = FALSE]
gene_lengths_bp <- gene_lengths_bp[common_genes]

message("Raw count matrix after annotation matching: ",
        nrow(raw_counts), " genes x ", ncol(raw_counts), " samples")

############################
## 4. Load/harmonize metadata + broad brain-area mapping
############################

metadata <- load_or_fetch_metadata(colnames(raw_counts))
metadata <- prepare_metadata(metadata)

required_meta <- c("SampleID", "OriginalRegion", "BrainID")
missing_required <- setdiff(required_meta, colnames(metadata))
if (length(missing_required) > 0) {
  stop("Metadata is missing required columns: ", paste(missing_required, collapse = ", "))
}

metadata$SampleID <- as.character(metadata$SampleID)

common_samples_all <- intersect(colnames(raw_counts), metadata$SampleID)
raw_counts <- raw_counts[, common_samples_all, drop = FALSE]
metadata <- metadata[match(common_samples_all, metadata$SampleID), , drop = FALSE]

mapping <- load_or_create_brain_area_mapping(metadata$OriginalRegion)
metadata <- apply_brain_area_mapping(metadata, mapping)

write.csv(mapping, file.path(output_paths$config, "brain_area_mapping_used.csv"), row.names = FALSE)
write.csv(metadata, file.path(output_paths$metadata, "metadata_clean_matched_all_samples_with_brain_area.csv"), row.names = FALSE)

message("After metadata matching: ", nrow(metadata), " samples")
message("Original fine regions detected:")
print(sort(table(metadata$OriginalRegion)))
message("Broad analysis groups detected using column '", analysis_group_column, "':")
print(sort(table(metadata$AnalysisGroup)))

composition_report <- make_group_composition_report(metadata, output_paths$qc)

## Reports before any group-level deduplication.
write_repetition_reports(
  metadata_subset = metadata,
  genes = rownames(raw_counts),
  out_dir = file.path(output_paths$repetition_reports, "all_matched_samples_before_group_dedup"),
  prefix = "all_matched_samples_before_group_dedup"
)
model_summary <- make_model_decision_summary(metadata, output_paths$qc)

############################
## 5. Process both versions by broad analysis group
############################

all_groups <- sort(unique(metadata$AnalysisGroup))
if (!is.null(target_analysis_groups)) {
  all_groups <- intersect(all_groups, target_analysis_groups)
}
if (length(all_groups) == 0) stop("No analysis groups selected/found.")

all_qc <- list()

for (mode in sample_handling_modes) {
  qc_mode <- do.call(rbind, lapply(
    all_groups,
    process_analysis_group,
    mode = mode,
    raw_counts = raw_counts,
    metadata = metadata,
    gene_lengths_bp = gene_lengths_bp,
    annotation = annotation
  ))

  write.csv(qc_mode,
            file.path(output_paths$qc, paste0("Step1_QC_summary_by_analysis_group_", mode, ".csv")),
            row.names = FALSE)
  all_qc[[mode]] <- qc_mode
}

qc_summary_all_modes <- do.call(rbind, all_qc)
write.csv(qc_summary_all_modes,
          file.path(output_paths$qc, "Step1_QC_summary_by_analysis_group_all_modes.csv"),
          row.names = FALSE)

############################
## 6. Save compact RDS for Step 2
############################

saveRDS(
  list(
    metadata = metadata,
    brain_area_mapping = mapping,
    analysis_group_column = analysis_group_column,
    analysis_group_composition = composition_report,
    annotation = annotation,
    qc_summary_all_modes = qc_summary_all_modes,
    model_decision_summary = model_summary,
    normalization = "GeTMM_log2",
    min_gene_length_bp = min_gene_length_bp,
    min_count = min_count,
    low_expr_exclusion_prop = low_expr_exclusion_prop,
    sample_handling_modes = sample_handling_modes
  ),
  file = file.path(output_paths$rds, "Step1_metadata_QC_and_repetition_reports_by_brain_area.rds")
)

message("\nDone. Files saved in: ", normalizePath(output_dir))
message("Main files to inspect:")
message("  - 00_config/brain_area_mapping_used.csv")
message("  - 01_metadata/metadata_clean_matched_all_samples_with_brain_area.csv")
message("  - 02_normalized_expression/<mode>/normalized_GeTMM_log2_<group>.csv")
message("  - 06_repetition_reports/<mode>/<group>/")
message("  - 07_qc/analysis_group_composition_report.csv")
message("  - 07_qc/Step1_QC_summary_by_analysis_group_all_modes.csv")
message("  - 07_qc/model_decision_repetition_summary.csv")
message("\nIf the grouping is not what you want, edit GSE102556_brain_area_mapping_template.csv and rerun.")
