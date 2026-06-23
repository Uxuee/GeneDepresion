## run_Tier2_enrichment_length_gt_100.R
##
## Purpose:
##   Run over-representation enrichment analysis for Length > 100 Tier 2
##   exploratory MDD DE candidate genes, using the tested genes from the same
##   region/subregion as the enrichment background.
##
## Important:
##   This script does not rerun LMMs. It uses the already packaged Tier 2
##   candidate gene sets and background universes.

rm(list = ls())
options(stringsAsFactors = FALSE)

local_d_library <- "D:/Rlibs"
if (dir.exists(local_d_library)) {
  .libPaths(unique(c(local_d_library, .libPaths())))
}

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(script_path)) {
  script_path <- normalizePath("code/run_Tier2_enrichment_length_gt_100.R", mustWork = FALSE)
}
repo_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

root <- file.path(repo_dir, "outputs", "length_gt_100_elisa_style")
input_dir <- file.path(root, "DE_candidate_gene_sets")
out_dir <- file.path(root, "enrichment_Tier2")
tables_dir <- file.path(out_dir, "tables")
plots_dir <- file.path(out_dir, "plots")
unmapped_dir <- file.path(out_dir, "unmapped_genes")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(unmapped_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  file.path(input_dir, "Tier2_working_exploratory"),
  file.path(input_dir, "background_genes_four_group_by_region.csv"),
  file.path(input_dir, "background_genes_frontal_subregion_by_subregion.csv")
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Missing required inputs:\n", paste(missing_inputs, collapse = "\n"))
}

required_packages <- c("clusterProfiler", "org.Hs.eg.db")
missing_required_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_required_packages) > 0) {
  stop(
    "Missing required enrichment packages: ",
    paste(missing_required_packages, collapse = ", "),
    "\nInstall them with BiocManager before running this script."
  )
}

has_reactome <- requireNamespace("ReactomePA", quietly = TRUE)
has_msigdb <- requireNamespace("msigdbr", quietly = TRUE)
has_enrichplot <- requireNamespace("enrichplot", quietly = TRUE)
has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)

library(clusterProfiler)
library(org.Hs.eg.db)

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

read_csv <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE)
}

parse_tier2_file <- function(path) {
  name <- basename(path)
  name <- sub("\\.csv$", "", name)
  parts <- strsplit(name, "__", fixed = TRUE)[[1]]
  if (length(parts) < 2) {
    stop("Could not parse Tier 2 file name: ", path)
  }
  list(
    label = parts[1],
    direction = parts[2]
  )
}

standardize_direction <- function(x) {
  x <- sub("^Tier2_", "", x)
  x
}

map_entrez <- function(df, id_col = "gene_id", symbol_col = "gene_symbol") {
  ids <- unique(as.character(df[[id_col]]))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  suppressMessages({
    mapped <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = ids,
      keytype = "ENTREZID",
      columns = c("ENTREZID", "SYMBOL")
    )
  })
  mapped <- mapped[!is.na(mapped$ENTREZID), , drop = FALSE]
  mapped <- mapped[!duplicated(mapped$ENTREZID), , drop = FALSE]
  unmapped <- df[!(as.character(df[[id_col]]) %in% mapped$ENTREZID), , drop = FALSE]
  list(mapped = mapped, unmapped = unmapped)
}

result_to_table <- function(enrich_result) {
  if (is.null(enrich_result)) {
    return(data.frame())
  }
  result <- tryCatch(as.data.frame(enrich_result), error = function(e) data.frame())
  if (nrow(result) == 0) {
    return(result)
  }
  result
}

save_enrichment_tables <- function(result, prefix) {
  full <- result_to_table(result)
  full_path <- file.path(tables_dir, paste0(prefix, "__full.csv"))
  sig_path <- file.path(tables_dir, paste0(prefix, "__padj_lt_005.csv"))
  relaxed_path <- file.path(tables_dir, paste0(prefix, "__padj_lt_010.csv"))

  write_csv(full, full_path)
  if (nrow(full) == 0 || !("p.adjust" %in% colnames(full))) {
    write_csv(full[0, , drop = FALSE], sig_path)
    write_csv(full[0, , drop = FALSE], relaxed_path)
  } else {
    write_csv(full[full$p.adjust < 0.05, , drop = FALSE], sig_path)
    write_csv(full[full$p.adjust < 0.10, , drop = FALSE], relaxed_path)
  }

  full
}

make_plots <- function(result, prefix, title_text) {
  if (!has_enrichplot || !has_ggplot2 || is.null(result)) {
    return(invisible(FALSE))
  }
  table <- result_to_table(result)
  if (nrow(table) == 0 || !("p.adjust" %in% colnames(table)) || !any(table$p.adjust < 0.10, na.rm = TRUE)) {
    return(invisible(FALSE))
  }

  dot_path <- file.path(plots_dir, paste0(prefix, "__dotplot.png"))
  bar_path <- file.path(plots_dir, paste0(prefix, "__barplot.png"))

  show_n <- min(15, sum(table$p.adjust < 0.10, na.rm = TRUE))

  tryCatch({
    png(dot_path, width = 1200, height = 900)
    print(enrichplot::dotplot(result, showCategory = show_n, title = title_text))
    dev.off()
  }, error = function(e) {
    if (names(dev.cur()) != "null device") dev.off()
  })

  if ("barplot" %in% getNamespaceExports("enrichplot")) {
    tryCatch({
      png(bar_path, width = 1200, height = 900)
      print(enrichplot::barplot(result, showCategory = show_n, title = title_text))
      dev.off()
    }, error = function(e) {
      if (names(dev.cur()) != "null device") dev.off()
    })
  }

  if (nrow(table) >= 5 && show_n >= 5) {
    sim_result <- tryCatch(enrichplot::pairwise_termsim(result), error = function(e) NULL)
    if (!is.null(sim_result)) {
      emap_path <- file.path(plots_dir, paste0(prefix, "__emapplot.png"))
      tryCatch({
        png(emap_path, width = 1300, height = 1000)
        print(enrichplot::emapplot(sim_result, showCategory = min(20, show_n)))
        dev.off()
      }, error = function(e) {
        if (names(dev.cur()) != "null device") dev.off()
      })
    }
  }

  invisible(TRUE)
}

run_enrich_go <- function(genes, universe, ont) {
  tryCatch(
    clusterProfiler::enrichGO(
      gene = genes,
      universe = universe,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    ),
    error = function(e) NULL
  )
}

run_enrich_kegg <- function(genes, universe) {
  tryCatch(
    clusterProfiler::enrichKEGG(
      gene = genes,
      universe = universe,
      organism = "hsa",
      keyType = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1
    ),
    error = function(e) NULL
  )
}

run_enrich_reactome <- function(genes, universe) {
  if (!has_reactome) {
    return(NULL)
  }
  tryCatch(
    ReactomePA::enrichPathway(
      gene = genes,
      universe = universe,
      organism = "human",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    ),
    error = function(e) NULL
  )
}

get_msigdb_term2gene <- function(collection, subcollection = NULL) {
  if (!has_msigdb) {
    return(NULL)
  }
  args <- list(species = "Homo sapiens")
  formals_names <- names(formals(msigdbr::msigdbr))
  if ("collection" %in% formals_names) {
    args$collection <- collection
    if (!is.null(subcollection) && "subcollection" %in% formals_names) {
      args$subcollection <- subcollection
    }
  } else {
    args$category <- collection
    if (!is.null(subcollection)) {
      args$subcategory <- subcollection
    }
  }
  df <- tryCatch(do.call(msigdbr::msigdbr, args), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }
  gene_col <- intersect(c("ncbi_gene", "entrez_gene", "entrez_gene_id"), colnames(df))[1]
  term_col <- intersect(c("gs_name", "gs_exact_source"), colnames(df))[1]
  if (is.na(gene_col) || is.na(term_col)) {
    return(NULL)
  }
  term2gene <- unique(data.frame(
    term = df[[term_col]],
    gene = as.character(df[[gene_col]]),
    stringsAsFactors = FALSE
  ))
  term2gene <- term2gene[!is.na(term2gene$gene) & nzchar(term2gene$gene), , drop = FALSE]
  term2gene
}

run_enrich_msigdb <- function(genes, universe, term2gene) {
  if (is.null(term2gene) || nrow(term2gene) == 0) {
    return(NULL)
  }
  tryCatch(
    clusterProfiler::enricher(
      gene = genes,
      universe = universe,
      TERM2GENE = term2gene,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1
    ),
    error = function(e) NULL
  )
}

top_term <- function(result_table, threshold) {
  if (nrow(result_table) == 0 || !("p.adjust" %in% colnames(result_table))) {
    return(NA_character_)
  }
  hit <- result_table[result_table$p.adjust < threshold, , drop = FALSE]
  if (nrow(hit) == 0) {
    return(NA_character_)
  }
  hit <- hit[order(hit$p.adjust, hit$pvalue, na.last = TRUE), , drop = FALSE]
  hit$Description[1]
}

make_summary_row <- function(meta, database, input_n, mapped_input_n, bg_n, mapped_bg_n, result_table, status) {
  data.frame(
    analysis_level = meta$analysis_level,
    region_or_subregion = meta$region_or_subregion,
    direction = meta$direction,
    database = database,
    n_input_genes = input_n,
    n_mapped_input_genes = mapped_input_n,
    n_background_genes = bg_n,
    n_mapped_background_genes = mapped_bg_n,
    n_terms_padj005 = if (nrow(result_table) && "p.adjust" %in% colnames(result_table)) sum(result_table$p.adjust < 0.05, na.rm = TRUE) else 0,
    n_terms_padj010 = if (nrow(result_table) && "p.adjust" %in% colnames(result_table)) sum(result_table$p.adjust < 0.10, na.rm = TRUE) else 0,
    top_term_padj005 = top_term(result_table, 0.05),
    top_term_padj010 = top_term(result_table, 0.10),
    status = status,
    stringsAsFactors = FALSE
  )
}

load_backgrounds <- function() {
  four <- read_csv(file.path(input_dir, "background_genes_four_group_by_region.csv"))
  frontal <- read_csv(file.path(input_dir, "background_genes_frontal_subregion_by_subregion.csv"))
  rbind(four, frontal)
}

build_gene_sets <- function() {
  tier2_dir <- file.path(input_dir, "Tier2_working_exploratory")
  files <- list.files(tier2_dir, pattern = "\\.csv$", full.names = TRUE)
  files <- files[grepl("Tier2_(all|up_in_MDD|down_in_MDD)\\.csv$", files)]

  sets <- lapply(files, function(path) {
    parsed <- parse_tier2_file(path)
    data <- read_csv(path)
    if (nrow(data) == 0) {
      analysis_level <- if (grepl("broad_Frontal|Insula|Ventral|Hippocampal", parsed$label)) "four_group" else "frontal_subregion"
      region <- gsub("_", " ", parsed$label)
    } else {
      analysis_level <- data$analysis_level[1]
      region <- data$region_or_subregion[1]
    }
    list(
      path = path,
      label = parsed$label,
      direction = standardize_direction(parsed$direction),
      analysis_level = analysis_level,
      region_or_subregion = region,
      data = data
    )
  })

  wanted_regions <- c(
    "Ventral Striatum", "Insula", "Hippocampal Formation", "broad Frontal Cortex",
    "Cingulate gyrus 25", "dlPFC BA8/9", "OFC BA11"
  )
  wanted_directions <- c("all", "up_in_MDD", "down_in_MDD")
  sets <- sets[vapply(sets, function(x) {
    x$region_or_subregion %in% wanted_regions && x$direction %in% wanted_directions
  }, logical(1))]
  sets
}

msigdb_collections <- list(
  MSigDB_Hallmark = list(collection = "H", subcollection = NULL),
  MSigDB_C2_curated = list(collection = "C2", subcollection = NULL),
  MSigDB_C5_GO = list(collection = "C5", subcollection = NULL),
  MSigDB_C8_cell_type = list(collection = "C8", subcollection = NULL)
)

message("Loading backgrounds and Tier 2 gene sets...")
backgrounds <- load_backgrounds()
gene_sets <- build_gene_sets()

message("Preloading MSigDB collections when available...")
msigdb_term2gene <- list()
if (has_msigdb) {
  for (db in names(msigdb_collections)) {
    spec <- msigdb_collections[[db]]
    message("  ", db)
    msigdb_term2gene[[db]] <- get_msigdb_term2gene(spec$collection, spec$subcollection)
  }
}

summary_rows <- list()
top_rows <- list()

for (set in gene_sets) {
  message("Running enrichment: ", set$region_or_subregion, " | ", set$direction)
  input_df <- set$data
  bg_df <- backgrounds[
    backgrounds$analysis_level == set$analysis_level &
      backgrounds$region_or_subregion == set$region_or_subregion,
    ,
    drop = FALSE
  ]

  meta <- list(
    analysis_level = set$analysis_level,
    region_or_subregion = set$region_or_subregion,
    direction = set$direction
  )

  input_map <- map_entrez(input_df)
  bg_map <- map_entrez(bg_df)
  input_genes <- unique(input_map$mapped$ENTREZID)
  bg_genes <- unique(bg_map$mapped$ENTREZID)
  input_genes <- intersect(input_genes, bg_genes)

  prefix_base <- paste(safe_name(set$analysis_level), safe_name(set$region_or_subregion), safe_name(set$direction), sep = "__")
  write_csv(input_map$unmapped, file.path(unmapped_dir, paste0(prefix_base, "__input_unmapped.csv")))
  write_csv(bg_map$unmapped, file.path(unmapped_dir, paste0(prefix_base, "__background_unmapped.csv")))

  databases <- list()
  databases[["GO_BP"]] <- function() run_enrich_go(input_genes, bg_genes, "BP")
  databases[["GO_MF"]] <- function() run_enrich_go(input_genes, bg_genes, "MF")
  databases[["GO_CC"]] <- function() run_enrich_go(input_genes, bg_genes, "CC")
  databases[["KEGG"]] <- function() run_enrich_kegg(input_genes, bg_genes)
  databases[["Reactome"]] <- function() run_enrich_reactome(input_genes, bg_genes)
  for (db in names(msigdb_term2gene)) {
    local({
      local_db <- db
      databases[[local_db]] <<- function() {
        run_enrich_msigdb(input_genes, bg_genes, msigdb_term2gene[[local_db]])
      }
    })
  }

  for (database in names(databases)) {
    result_table <- data.frame()
    status <- "ok"
    result <- NULL

    if (set$region_or_subregion == "broad Frontal Cortex") {
      status <- "record_only_broad_frontal_cortex"
    }
    if (length(input_genes) < 5) {
      status <- ifelse(
        status == "ok",
        "skipped_fewer_than_5_mapped_input_genes",
        paste(status, "skipped_fewer_than_5_mapped_input_genes", sep = "; ")
      )
    } else if (length(bg_genes) < 5) {
      status <- "skipped_fewer_than_5_mapped_background_genes"
    } else {
      result <- databases[[database]]()
      if (is.null(result)) {
        status <- ifelse(status == "ok", "not_available_or_failed", paste(status, "not_available_or_failed", sep = "; "))
      }
      result_table <- save_enrichment_tables(result, paste(prefix_base, safe_name(database), sep = "__"))
      if (!is.null(result)) {
        make_plots(result, paste(prefix_base, safe_name(database), sep = "__"), paste(set$region_or_subregion, set$direction, database))
      }
    }

    if (length(input_genes) < 5 || length(bg_genes) < 5) {
      save_enrichment_tables(NULL, paste(prefix_base, safe_name(database), sep = "__"))
    }

    summary_rows[[length(summary_rows) + 1]] <- make_summary_row(
      meta = meta,
      database = database,
      input_n = nrow(input_df),
      mapped_input_n = length(input_genes),
      bg_n = nrow(bg_df),
      mapped_bg_n = length(bg_genes),
      result_table = result_table,
      status = status
    )

    if (nrow(result_table) > 0) {
      top <- result_table[order(result_table$p.adjust, result_table$pvalue, na.last = TRUE), , drop = FALSE]
      top <- head(top, 20)
      wanted_top_columns <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
      for (column in wanted_top_columns) {
        if (!(column %in% colnames(top))) {
          top[[column]] <- NA
        }
      }
      top <- top[, wanted_top_columns, drop = FALSE]
      top$analysis_level <- set$analysis_level
      top$region_or_subregion <- set$region_or_subregion
      top$direction <- set$direction
      top$database <- database
      top_rows[[length(top_rows) + 1]] <- top
    }
  }
}

summary_df <- do.call(rbind, summary_rows)
write_csv(summary_df, file.path(out_dir, "enrichment_summary_by_gene_set.csv"))

if (length(top_rows) > 0) {
  top_df <- do.call(rbind, top_rows)
} else {
  top_df <- data.frame()
}
write_csv(top_df, file.path(out_dir, "top_enrichment_terms_all_sets.csv"))

write_report <- function(summary_df, top_df) {
  input_counts <- unique(summary_df[, c("analysis_level", "region_or_subregion", "direction", "n_input_genes", "n_mapped_input_genes", "n_background_genes", "n_mapped_background_genes")])
  top_sig <- summary_df[summary_df$n_terms_padj010 > 0, , drop = FALSE]
  top_sig <- top_sig[order(-top_sig$n_terms_padj010, top_sig$region_or_subregion, top_sig$direction), , drop = FALSE]
  top_sig <- head(top_sig, 30)

  region_sections <- c("Ventral Striatum", "Insula", "Hippocampal Formation", "Cingulate gyrus 25", "dlPFC BA8/9", "OFC BA11")
  section_text <- character()
  for (region in region_sections) {
    region_summary <- summary_df[summary_df$region_or_subregion == region & summary_df$n_terms_padj010 > 0, , drop = FALSE]
    if (nrow(region_summary) == 0) {
      section_text <- c(section_text, paste0("### ", region), "", "No relaxed adjusted p < 0.10 enrichment terms were detected for this region/subregion.", "")
    } else {
      region_summary <- region_summary[order(-region_summary$n_terms_padj010, region_summary$direction, region_summary$database), , drop = FALSE]
      section_text <- c(
        section_text,
        paste0("### ", region),
        "",
        paste(capture.output(print(head(region_summary[, c("direction", "database", "n_terms_padj005", "n_terms_padj010", "top_term_padj010", "status")], 10), row.names = FALSE)), collapse = "\n"),
        ""
      )
    }
  }

  up_down <- summary_df[summary_df$n_terms_padj010 > 0 & summary_df$direction %in% c("up_in_MDD", "down_in_MDD"), , drop = FALSE]
  if (nrow(up_down) == 0) {
    up_down_text <- "No direction-specific gene set produced relaxed adjusted p < 0.10 terms."
  } else {
    up_down_text <- paste(capture.output(print(head(up_down[order(up_down$region_or_subregion, up_down$direction, -up_down$n_terms_padj010), c("region_or_subregion", "direction", "database", "n_terms_padj010", "top_term_padj010")], 30), row.names = FALSE)), collapse = "\n")
  }

  report <- c(
    "# Tier 2 Enrichment Report",
    "",
    "## Scope",
    "",
    "This enrichment analysis used the Length > 100 Tier 2 exploratory DE candidate gene sets.",
    "",
    "Tier 2 definition:",
    "",
    "```text",
    "p_value_diagnosis < 0.01 and abs_standardized_effect_size_diagnosis > 0.8",
    "```",
    "",
    "Tier 2 genes are exploratory DE candidate genes, not FDR-significant DE genes.",
    "",
    "The enrichment background for each region/subregion was the tested-gene universe from that same analysis, not all genes in the genome.",
    "",
    "## Input And Mapping Summary",
    "",
    paste(capture.output(print(input_counts, row.names = FALSE)), collapse = "\n"),
    "",
    "Gene sets with fewer than 5 mapped input genes were skipped. Broad Frontal Cortex was retained for documentation only and is not useful for enrichment because it has only 1 Tier 2 gene.",
    "",
    "## Main Enrichment Summary",
    "",
    if (nrow(top_sig) > 0) paste(capture.output(print(top_sig[, c("analysis_level", "region_or_subregion", "direction", "database", "n_terms_padj005", "n_terms_padj010", "top_term_padj010", "status")], row.names = FALSE)), collapse = "\n") else "No relaxed adjusted p < 0.10 terms were detected.",
    "",
    "## Region/Subregion Interpretation",
    "",
    section_text,
    "## Up-In-MDD Versus Down-In-MDD Patterns",
    "",
    up_down_text,
    "",
    "## Cautious Interpretation",
    "",
    "These enrichment results are exploratory and depend on the Tier 2 candidate definition. They should be used to prioritize biological hypotheses and downstream analyses, not as final proof of pathway-level disease mechanisms.",
    "",
    "The p-value correction sensitivity audit showed that larger single-gene counts mainly appear under effect-size-gated or direction-specific correction scopes. Those remain sensitivity results and do not replace Tier 2.",
    "",
    "## Recommended Next Step",
    "",
    "After this ORA enrichment, run ranked GSEA using the Tier 4 signed standardized log-p signatures.",
    "",
    "## Outputs",
    "",
    "- `tables/`: full, adjusted p < 0.05, and adjusted p < 0.10 enrichment tables.",
    "- `plots/`: dotplots, barplots, and enrichment maps when enough terms were available.",
    "- `unmapped_genes/`: unmapped input and background genes for each gene set.",
    "- `enrichment_summary_by_gene_set.csv`: summary of mapping and term counts.",
    "- `top_enrichment_terms_all_sets.csv`: top 20 terms from each successful enrichment analysis."
  )

  writeLines(report, file.path(out_dir, "Tier2_Enrichment_Report.md"))
}

write_readme <- function() {
  readme <- c(
    "# Tier 2 Enrichment Analysis",
    "",
    "This folder contains over-representation enrichment analysis for the Length > 100 Tier 2 exploratory DE candidate gene sets.",
    "",
    "Analyses include GO Biological Process, GO Molecular Function, GO Cellular Component, KEGG, Reactome when available, and MSigDB collections when available.",
    "",
    "Tables are in `tables/`.",
    "Plots are in `plots/`.",
    "Unmapped input and background genes are in `unmapped_genes/`.",
    "",
    "Backgrounds were defined as all genes tested in the same region/subregion analysis, using the packaged background universe files from `DE_candidate_gene_sets/`.",
    "",
    "To rerun:",
    "",
    "```powershell",
    "& \"C:\\Program Files\\R\\R-4.4.2\\bin\\Rscript.exe\" code\\run_Tier2_enrichment_length_gt_100.R",
    "```",
    "",
    "Tier 2 genes are exploratory DE candidate genes and should not be described as FDR-significant DE genes."
  )
  writeLines(readme, file.path(out_dir, "README.md"))
}

write_report(summary_df, top_df)
write_readme()

message("Saved Tier 2 enrichment outputs to: ", out_dir)
message("Summary rows: ", nrow(summary_df))
message("Top-term rows: ", nrow(top_df))
