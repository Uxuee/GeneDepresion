suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) >= 1) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else normalizePath(".", winslash = "/")
root <- file.path(repo, "outputs", "length_gt_100_clean_pipeline", "21_Statistical_method_benchmark")
input_dir <- file.path(root, "matched_inputs")
result_dir <- file.path(root, "region_results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

regions <- data.frame(
  region = c("Ventral Striatum", "Insula", "Hippocampal Formation", "CG25", "dlPFC BA8/9", "OFC BA11"),
  slug = c("ventral_striatum", "insula", "hippocampal_formation", "cg25", "dlpfc_ba8_9", "ofc_ba11"),
  stringsAsFactors = FALSE
)

rank_primary <- function(p, gene) {
  order_index <- order(p, gene, na.last = TRUE, method = "radix")
  answer <- rep(NA_integer_, length(p))
  answer[order_index] <- seq_along(order_index)
  answer
}

limma_all <- list()
deseq_all <- list()

for (i in seq_len(nrow(regions))) {
  region <- regions$region[[i]]
  slug <- regions$slug[[i]]
  message("Fitting ", region)
  counts <- read.csv(gzfile(file.path(input_dir, paste0(slug, "_raw_counts.csv.gz"))), row.names = 1, check.names = FALSE)
  counts <- as.matrix(counts)
  storage.mode(counts) <- "integer"
  metadata <- read.csv(file.path(input_dir, paste0(slug, "_metadata.csv")), stringsAsFactors = FALSE, check.names = FALSE)
  rownames(metadata) <- metadata$SampleID
  stopifnot(identical(colnames(counts), metadata$SampleID))
  stopifnot(all(counts >= 0L), all(counts == floor(counts)))
  metadata$Diagnosis <- factor(metadata$Diagnosis, levels = c("CTRL", "MDD"))
  metadata$Sex <- factor(metadata$Sex)
  metadata$Age <- as.numeric(metadata$Age)
  metadata$RIN <- as.numeric(metadata$RIN)
  stopifnot(!anyNA(metadata[, c("Diagnosis", "Age", "Sex", "RIN")]))

  design_limma <- model.matrix(~ Diagnosis + Age + Sex + RIN, data = metadata)
  stopifnot("DiagnosisMDD" %in% colnames(design_limma))
  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  voomed <- limma::voom(dge, design = design_limma, plot = FALSE)
  fit <- limma::eBayes(limma::lmFit(voomed, design_limma))
  tab <- limma::topTable(fit, coef = "DiagnosisMDD", number = Inf, sort.by = "none", adjust.method = "BH")
  limma_result <- data.frame(
    region = region,
    gene_id = rownames(tab),
    gene_symbol = NA_character_,
    logFC = tab$logFC,
    t_statistic = tab$t,
    raw_p = tab$P.Value,
    BH = tab$adj.P.Val,
    average_expression = tab$AveExpr,
    stringsAsFactors = FALSE
  )
  limma_result$primary_rank <- rank_primary(limma_result$raw_p, limma_result$gene_id)
  limma_result$signed_evidence <- sign(limma_result$logFC) * -log10(pmax(limma_result$raw_p, .Machine$double.xmin))
  write.csv(limma_result, gzfile(file.path(result_dir, paste0(slug, "_limma_voom.csv.gz"))), row.names = FALSE)
  limma_all[[i]] <- limma_result

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ Age + Sex + RIN + Diagnosis
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  deseq_tab <- as.data.frame(DESeq2::results(dds, contrast = c("Diagnosis", "MDD", "CTRL"), independentFiltering = TRUE))
  deseq_result <- data.frame(
    region = region,
    gene_id = rownames(deseq_tab),
    gene_symbol = NA_character_,
    baseMean = deseq_tab$baseMean,
    log2FoldChange = deseq_tab$log2FoldChange,
    lfcSE = deseq_tab$lfcSE,
    statistic = deseq_tab$stat,
    raw_p = deseq_tab$pvalue,
    native_padj = deseq_tab$padj,
    stringsAsFactors = FALSE
  )
  deseq_result$BH_same_universe <- p.adjust(deseq_result$raw_p, method = "BH")
  deseq_result$primary_rank <- rank_primary(deseq_result$raw_p, deseq_result$gene_id)
  deseq_result$signed_evidence <- sign(deseq_result$log2FoldChange) * -log10(pmax(deseq_result$raw_p, .Machine$double.xmin))
  write.csv(deseq_result, gzfile(file.path(result_dir, paste0(slug, "_deseq2.csv.gz"))), row.names = FALSE)
  deseq_all[[i]] <- deseq_result
}

write.csv(do.call(rbind, limma_all), file.path(root, "limma_voom_results_all_regions.csv"), row.names = FALSE)
write.csv(do.call(rbind, deseq_all), file.path(root, "deseq2_results_all_regions.csv"), row.names = FALSE)

versions <- data.frame(
  component = c("R", "limma", "edgeR", "DESeq2", "limma_normalization", "DESeq2_size_factors", "DESeq2_dispersion_fit", "contrast"),
  value = c(
    R.version.string,
    as.character(packageVersion("limma")),
    as.character(packageVersion("edgeR")),
    as.character(packageVersion("DESeq2")),
    "TMM via edgeR::calcNormFactors; standard voom",
    "DESeq2 median-of-ratios estimateSizeFactors",
    "standard DESeq() defaults",
    "Diagnosis MDD - CTRL"
  ),
  stringsAsFactors = FALSE
)
write.csv(versions, file.path(root, "software_versions.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(root, "R_sessionInfo.txt"))
message("Completed limma-voom and DESeq2 for six regions")
