#!/usr/bin/env Rscript
# Step7_Annotate_MDD_WGCNA_Modules.R
# Biological annotation of prioritized exploratory MDD WGCNA modules.
# Uses saved Step 12 outputs only. Does not rerun WGCNA.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(jsonlite)
})

project_root <- "D:/genes/GeneDepresion"
pipeline_root <- file.path(project_root, "outputs", "length_gt_100_clean_pipeline")
wgcna_root <- file.path(pipeline_root, "12_MDD_WGCNA_modules")
ahba_root <- file.path(pipeline_root, "11_AHBA_modules_enrichment")
out_dir <- file.path(wgcna_root, "biological_annotation")
cards_dir <- file.path(out_dir, "module_cards")
plots_dir <- file.path(out_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cards_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

focus_modules <- tibble::tribble(
  ~region, ~module,
  "Ventral Striatum", "MElavenderblush3",
  "Ventral Striatum", "MElightcyan1",
  "Ventral Striatum", "MEdarkred",
  "Ventral Striatum", "MElightgreen",
  "Insula", "MEorange",
  "Hippocampal Formation", "MEgreen",
  "Frontal Cortex", "MEdarkgreen"
)
secondary_modules <- tibble::tribble(
  ~region, ~module,
  "Frontal Cortex", "MEblack",
  "Frontal Cortex", "MEsalmon",
  "Insula", "MEdarkorange"
)
modules_to_annotate <- bind_rows(focus_modules, secondary_modules) %>% distinct()

input_files <- list(
  top_priority = file.path(wgcna_root, "mdd_wgcna_top_priority_modules_for_report.csv"),
  triage = file.path(wgcna_root, "mdd_wgcna_module_triage_table.csv"),
  assignments = file.path(wgcna_root, "best_solution_gene_to_module_assignments.csv"),
  hubs = file.path(wgcna_root, "module_hub_genes_top50.csv"),
  trait = file.path(wgcna_root, "module_trait_association_results.csv"),
  enrichment = file.path(wgcna_root, "module_enrichment_for_Tier2_and_ML_genes.csv"),
  ahba_overlap = file.path(wgcna_root, "mdd_modules_vs_published_AHBA_modules_enrichment.csv"),
  tier2 = file.path(pipeline_root, "07_Tier2_candidate_gene_sets", "Tier2_candidate_genes_all_regions.csv"),
  ml = file.path(pipeline_root, "10_ML_prioritized_Tier2_corrected", "ML_prioritized_Tier2_genes_all_ranked.csv"),
  ahba_annotation_summary = file.path(ahba_root, "published_modules", "annotations", "published_ahba_module_annotation_summary.csv"),
  ahba_corrected_annotations = file.path(ahba_root, "published_modules", "annotations", "published_ahba_corrected_suggestive_hits_with_annotations.csv")
)

read_csv_quiet <- function(path) readr::read_csv(path, show_col_types = FALSE)
na0 <- function(x) ifelse(is.na(x), "", x)
fmt_num <- function(x) {
  trimws(ifelse(is.na(x), "NA", ifelse(abs(x) < 0.001, formatC(x, format = "e", digits = 2), formatC(x, format = "fg", digits = 3))))
}

trim_markdown_file <- function(path) {
  if (!file.exists(path)) return(invisible(FALSE))
  lines <- readLines(path, warn = FALSE)
  lines <- sub("[[:space:]]+$", "", lines)
  while (length(lines) > 0 && lines[[length(lines)]] == "") {
    lines <- lines[-length(lines)]
  }
  writeLines(lines, path, useBytes = TRUE)
  invisible(TRUE)
}
trim_markdown_outputs <- function(root) {
  md_files <- list.files(root, pattern = "\\.md$", recursive = TRUE, full.names = TRUE)
  invisible(lapply(md_files, trim_markdown_file))
}

clean_term <- function(x) {
  x %>% str_replace("^(GOBP|GOCC|GOMF|HALLMARK|REACTOME|KEGG|WP)_", "") %>% str_replace_all("_", " ") %>% str_to_sentence()
}
safe_name <- function(x) str_replace_all(x, "[^A-Za-z0-9]+", "_")

packages <- c("AnnotationDbi", "org.Hs.eg.db", "clusterProfiler", "ReactomePA", "msigdbr", "enrichplot")
package_available <- setNames(vapply(packages, requireNamespace, logical(1), quietly = TRUE), packages)
skipped_databases <- tibble::tibble(
  database = character(), status = character(), note = character()
)
if (!package_available[["org.Hs.eg.db"]]) {
  skipped_databases <- bind_rows(skipped_databases, tibble(database="org.Hs.eg.db", status="skipped", note="Package unavailable; gene descriptions and formal Entrez/Ensembl mapping from org.Hs.eg.db were not run. MSigDB-derived mappings were used when available."))
}
if (!package_available[["clusterProfiler"]]) {
  skipped_databases <- bind_rows(skipped_databases, tibble(database="clusterProfiler", status="skipped", note="Package unavailable; GO/KEGG enrichment used custom hypergeometric tests over MSigDB gene sets when available."))
}
if (!package_available[["ReactomePA"]]) {
  skipped_databases <- bind_rows(skipped_databases, tibble(database="ReactomePA", status="skipped", note="Package unavailable; Reactome-like pathway enrichment used MSigDB C2 Reactome gene sets when available."))
}

triage <- read_csv_quiet(input_files$triage)
assignments <- read_csv_quiet(input_files$assignments) %>% mutate(gene_id = as.character(gene_id), gene_symbol = as.character(gene_symbol))
hubs <- read_csv_quiet(input_files$hubs) %>% mutate(gene_id = as.character(gene_id), gene_symbol = as.character(gene_symbol))
trait <- read_csv_quiet(input_files$trait)
module_enrich <- read_csv_quiet(input_files$enrichment)
ahba_overlap <- read_csv_quiet(input_files$ahba_overlap)
tier2 <- read_csv_quiet(input_files$tier2) %>% mutate(gene_id = as.character(gene_id), gene_symbol = as.character(gene_symbol))
ml <- read_csv_quiet(input_files$ml) %>% mutate(gene_id = as.character(gene_id), gene_symbol = as.character(gene_symbol))
ahba_annotation <- read_csv_quiet(input_files$ahba_annotation_summary)
ahba_corrected <- read_csv_quiet(input_files$ahba_corrected_annotations)

module_gene_sets <- assignments %>%
  semi_join(modules_to_annotate, by = c("region", "module")) %>%
  filter(module != "MEgrey") %>%
  distinct(region, module, gene_id, gene_symbol)

region_backgrounds <- assignments %>%
  filter(module != "MEgrey", !is.na(gene_symbol), gene_symbol != "") %>%
  distinct(region, gene_symbol)

all_module_symbols <- unique(module_gene_sets$gene_symbol)

# MSigDB mappings and gene sets.
msig_available <- package_available[["msigdbr"]]
msig_mapping <- tibble(gene_symbol = character(), entrez_gene = character(), ensembl_gene = character())
msig_sets <- list()
if (msig_available) {
  suppressPackageStartupMessages(library(msigdbr))
  safe_msig <- function(...) {
    tryCatch(msigdbr::msigdbr(species = "Homo sapiens", ...), error = function(e) {
      skipped_databases <<- bind_rows(skipped_databases, tibble(database="msigdbr", status="partial_skip", note=paste("msigdbr call failed:", conditionMessage(e))))
      tibble()
    })
  }
  msig_h <- safe_msig(collection = "H")
  msig_c2_cp <- safe_msig(collection = "C2", subcollection = "CP")
  msig_c5 <- safe_msig(collection = "C5")
  msig_all <- bind_rows(msig_h, msig_c2_cp, msig_c5)
  if (nrow(msig_all) > 0) {
    entrez_col <- if ("entrez_gene" %in% names(msig_all)) "entrez_gene" else if ("ncbi_gene" %in% names(msig_all)) "ncbi_gene" else NA_character_
    ensembl_col <- if ("ensembl_gene" %in% names(msig_all)) "ensembl_gene" else if ("db_ensembl_gene" %in% names(msig_all)) "db_ensembl_gene" else NA_character_
    msig_mapping <- msig_all %>%
      mutate(
        entrez_gene = if (!is.na(entrez_col)) as.character(.data[[entrez_col]]) else NA_character_,
        ensembl_gene = if (!is.na(ensembl_col)) as.character(.data[[ensembl_col]]) else NA_character_
      ) %>%
      transmute(gene_symbol = as.character(gene_symbol), entrez_gene, ensembl_gene) %>%
      filter(!is.na(gene_symbol), gene_symbol != "") %>%
      distinct(gene_symbol, .keep_all = TRUE)
  }
  msig_sets$GO_BP <- msig_c5 %>% filter(gs_subcollection == "GO:BP")
  msig_sets$GO_CC <- msig_c5 %>% filter(gs_subcollection == "GO:CC")
  msig_sets$GO_MF <- msig_c5 %>% filter(gs_subcollection == "GO:MF")
  msig_sets$HALLMARK <- msig_h
  msig_sets$C2_CP <- msig_c2_cp
  msig_sets$REACTOME <- msig_c2_cp %>% filter(str_detect(gs_name, "^REACTOME_"))
  msig_sets$KEGG <- msig_c2_cp %>% filter(str_detect(gs_name, "^KEGG_"))
} else {
  skipped_databases <- bind_rows(skipped_databases, tibble(database="msigdbr", status="skipped", note="Package unavailable; MSigDB-derived enrichments were not run."))
}

mapping_report <- module_gene_sets %>%
  left_join(msig_mapping, by = "gene_symbol") %>%
  mutate(
    mapping_source = if_else(!is.na(entrez_gene) | !is.na(ensembl_gene), "msigdbr", "unmapped"),
    mapping_note = if_else(mapping_source == "unmapped", "No org.Hs.eg.db available and gene was not found in loaded MSigDB mapping.", "Mapped from MSigDB gene metadata.")
  )
write_csv(mapping_report, file.path(out_dir, "module_gene_id_mapping_report.csv"))

unmapped <- mapping_report %>% filter(mapping_source == "unmapped") %>% arrange(region, module, gene_symbol)
write_csv(unmapped, file.path(out_dir, "unmapped_module_genes.csv"))

hyper_enrich <- function(query_symbols, background_symbols, gene_sets_df, source_name, min_size = 5, max_size = 500) {
  query_symbols <- unique(na.omit(query_symbols))
  background_symbols <- unique(na.omit(background_symbols))
  if (length(query_symbols) == 0 || length(background_symbols) == 0 || nrow(gene_sets_df) == 0) return(tibble())
  query_symbols <- intersect(query_symbols, background_symbols)
  sets <- gene_sets_df %>%
    transmute(gs_name, gs_description = ifelse("gs_description" %in% names(gene_sets_df), gs_description, NA_character_), gene_symbol = as.character(gene_symbol)) %>%
    filter(gene_symbol %in% background_symbols) %>%
    distinct(gs_name, gene_symbol, .keep_all = TRUE) %>%
    group_by(gs_name, gs_description) %>%
    summarise(set_genes = list(unique(gene_symbol)), set_size = n_distinct(gene_symbol), .groups = "drop") %>%
    filter(set_size >= min_size, set_size <= max_size)
  if (nrow(sets) == 0) return(tibble())
  N <- length(background_symbols)
  n <- length(query_symbols)
  out <- sets %>%
    rowwise() %>%
    mutate(
      overlap_genes = list(intersect(query_symbols, set_genes)),
      overlap_count = length(overlap_genes),
      p_value = phyper(overlap_count - 1, set_size, N - set_size, n, lower.tail = FALSE)
    ) %>%
    ungroup() %>%
    filter(overlap_count > 0) %>%
    mutate(
      fdr_BH = p.adjust(p_value, method = "BH"),
      source = source_name,
      term_label = clean_term(gs_name),
      overlap_genes = vapply(overlap_genes, paste, character(1), collapse = ";")
    ) %>%
    arrange(fdr_BH, p_value, desc(overlap_count))
  out
}

module_queries <- module_gene_sets %>%
  group_by(region, module) %>%
  summarise(query_symbols = list(unique(gene_symbol)), module_size = n_distinct(gene_symbol), .groups = "drop")

run_module_enrichment <- function(collection_name, output_file) {
  gs <- msig_sets[[collection_name]]
  if (is.null(gs) || nrow(gs) == 0) {
    skipped_databases <<- bind_rows(skipped_databases, tibble(database=collection_name, status="skipped", note="No gene sets available from installed resources."))
    empty <- tibble(region=character(), module=character(), source=character(), gs_name=character(), term_label=character(), overlap_count=integer(), set_size=integer(), module_size=integer(), p_value=double(), fdr_BH=double(), overlap_genes=character())
    write_csv(empty, file.path(out_dir, output_file))
    return(empty)
  }
  rows <- list()
  for (i in seq_len(nrow(module_queries))) {
    reg <- module_queries$region[[i]]
    mod <- module_queries$module[[i]]
    bg <- region_backgrounds %>% filter(region == reg) %>% pull(gene_symbol)
    enr <- hyper_enrich(module_queries$query_symbols[[i]], bg, gs, collection_name)
    if (nrow(enr) > 0) {
      enr <- enr %>% mutate(region = reg, module = mod, module_size = module_queries$module_size[[i]], .before = 1)
      rows[[length(rows) + 1]] <- enr
    }
  }
  res <- bind_rows(rows)
  if (nrow(res) == 0) {
    res <- tibble(region=character(), module=character(), source=character(), gs_name=character(), term_label=character(), overlap_count=integer(), set_size=integer(), module_size=integer(), p_value=double(), fdr_BH=double(), overlap_genes=character())
  }
  write_csv(res, file.path(out_dir, output_file))
  res
}

go_bp <- run_module_enrichment("GO_BP", "go_BP_enrichment_by_module.csv")
go_cc <- run_module_enrichment("GO_CC", "go_CC_enrichment_by_module.csv")
go_mf <- run_module_enrichment("GO_MF", "go_MF_enrichment_by_module.csv")
hallmark <- run_module_enrichment("HALLMARK", "msigdb_hallmark_enrichment_by_module.csv")
c2_cp <- run_module_enrichment("C2_CP", "msigdb_c2_cp_enrichment_by_module.csv")
reactome <- run_module_enrichment("REACTOME", "reactome_enrichment_by_module.csv")
kegg <- run_module_enrichment("KEGG", "kegg_enrichment_by_module.csv")

# Cell-type marker sanity-check list. Used only because no local marker files were present.
cell_markers <- tibble::tribble(
  ~cell_type, ~gene_symbol,
  "neuron_synaptic", "SNAP25", "neuron_synaptic", "SYT1", "neuron_synaptic", "SYN1", "neuron_synaptic", "RBFOX3", "neuron_synaptic", "MAP2", "neuron_synaptic", "GRIN1", "neuron_synaptic", "CAMK2A",
  "astrocyte", "GFAP", "astrocyte", "AQP4", "astrocyte", "ALDH1L1", "astrocyte", "SLC1A3", "astrocyte", "SLC1A2", "astrocyte", "GJA1",
  "oligodendrocyte_myelin", "MBP", "oligodendrocyte_myelin", "MOG", "oligodendrocyte_myelin", "PLP1", "oligodendrocyte_myelin", "MAG", "oligodendrocyte_myelin", "MOBP", "oligodendrocyte_myelin", "OLIG2",
  "microglia_immune", "AIF1", "microglia_immune", "CX3CR1", "microglia_immune", "CSF1R", "microglia_immune", "TYROBP", "microglia_immune", "C1QA", "microglia_immune", "TREM2", "microglia_immune", "PTPRC",
  "endothelial_vascular", "CLDN5", "endothelial_vascular", "PECAM1", "endothelial_vascular", "VWF", "endothelial_vascular", "KDR", "endothelial_vascular", "FLT1", "endothelial_vascular", "RGS5", "endothelial_vascular", "PDGFRB",
  "ependymal_ciliary", "FOXJ1", "ependymal_ciliary", "DNAH5", "ependymal_ciliary", "PIFO", "ependymal_ciliary", "TPPP3", "ependymal_ciliary", "TUBB4B", "ependymal_ciliary", "HYDIN",
  "inhibitory_neuron", "GAD1", "inhibitory_neuron", "GAD2", "inhibitory_neuron", "SLC6A1", "inhibitory_neuron", "SST", "inhibitory_neuron", "PVALB", "inhibitory_neuron", "VIP", "inhibitory_neuron", "RELN",
  "excitatory_neuron", "SLC17A7", "excitatory_neuron", "SLC17A6", "excitatory_neuron", "CAMK2A", "excitatory_neuron", "SATB2", "excitatory_neuron", "TBR1", "excitatory_neuron", "GRIA2"
) %>% distinct()

cell_rows <- list()
for (i in seq_len(nrow(module_queries))) {
  reg <- module_queries$region[[i]]
  mod <- module_queries$module[[i]]
  bg <- region_backgrounds %>% filter(region == reg) %>% pull(gene_symbol)
  res <- hyper_enrich(module_queries$query_symbols[[i]], bg, cell_markers %>% transmute(gs_name = cell_type, gs_description = cell_type, gene_symbol), "built_in_minimal_cell_marker_sanity_list", min_size = 2, max_size = 50)
  if (nrow(res) > 0) cell_rows[[length(cell_rows)+1]] <- res %>% mutate(region=reg, module=mod, module_size=module_queries$module_size[[i]], .before=1)
}
cell_enrich <- bind_rows(cell_rows)
if (nrow(cell_enrich) == 0) cell_enrich <- tibble(region=character(), module=character(), source=character(), gs_name=character(), term_label=character(), overlap_count=integer(), set_size=integer(), module_size=integer(), p_value=double(), fdr_BH=double(), overlap_genes=character())
write_csv(cell_enrich, file.path(out_dir, "cell_type_marker_enrichment_by_module.csv"))
cell_notes <- c(
  "# Cell-Type Annotation Notes",
  "",
  "No local marker files were found under data/reference/cell_type_markers/.",
  "Cell-type marker enrichment therefore used a small built-in sanity-check marker list for broad classes: neuron/synaptic, astrocyte, oligodendrocyte/myelin, microglia/immune, endothelial/vascular, ependymal/ciliary, inhibitory neuron, and excitatory neuron.",
  "Because this fallback marker list is intentionally small, these results should be interpreted only as weak supportive clues, not formal cell-type deconvolution or validated cell-type enrichment.",
  ""
)
writeLines(cell_notes, file.path(out_dir, "cell_type_annotation_notes.md"))

# Module overlaps and hub annotations.
tier2_ids <- unique(tier2$gene_id)
tier2_symbols <- unique(tier2$gene_symbol)
ml_high <- ml %>% filter(str_detect(tolower(ML_priority_class), "high-confidence"))
ml_high_mod <- ml %>% filter(str_detect(tolower(ML_priority_class), "high-confidence|moderate-confidence"))

same_region_tier2 <- function(reg) {
  tier2 %>% filter(region_or_subregion == reg | str_replace_all(region_or_subregion, " ", "_") == str_replace_all(reg, " ", "_"))
}

hub_annotations <- hubs %>%
  semi_join(modules_to_annotate, by = c("region", "module")) %>%
  filter(module != "MEgrey") %>%
  mutate(
    is_Tier2 = gene_id %in% tier2_ids | gene_symbol %in% tier2_symbols,
    is_high_confidence_ML = gene_id %in% ml_high$gene_id | gene_symbol %in% ml_high$gene_symbol,
    is_high_moderate_ML = gene_id %in% ml_high_mod$gene_id | gene_symbol %in% ml_high_mod$gene_symbol,
    known_gene_description = NA_character_,
    annotation_category = case_when(
      is_high_confidence_ML ~ "high-confidence ML-prioritized Tier2 hub",
      is_high_moderate_ML ~ "high/moderate ML-prioritized Tier2 hub",
      is_Tier2 ~ "Tier2 hub",
      TRUE ~ "module hub"
    )
  ) %>%
  transmute(region, module, hub_gene = gene_symbol, gene_id, kME = kME, abs_kME, gene_significance_diagnosis, abs_gene_significance_diagnosis, whether_Tier2 = is_Tier2, whether_high_confidence_ML = is_high_confidence_ML, whether_high_moderate_ML = is_high_moderate_ML, known_gene_name_description = known_gene_description, short_annotation_category = annotation_category)
write_csv(hub_annotations, file.path(out_dir, "prioritized_module_hub_gene_annotations.csv"))

# AHBA context.
ahba_best <- ahba_overlap %>%
  filter(mdd_module != "MEgrey") %>%
  semi_join(modules_to_annotate, by = c("region", "mdd_module" = "module")) %>%
  arrange(region, mdd_module, fdr_BH, p_value, desc(overlap_count)) %>%
  group_by(region, module = mdd_module) %>%
  slice(1) %>%
  ungroup()
ahba_highlight <- ahba_overlap %>%
  filter(ahba_module %in% c("M26", "M28", "M2", "M6", "M7"), mdd_module != "MEgrey") %>%
  semi_join(modules_to_annotate, by = c("region", "mdd_module" = "module")) %>%
  mutate(hit = paste0(ahba_module, " overlap=", overlap_count, " FDR=", fmt_num(fdr_BH))) %>%
  group_by(region, module = mdd_module) %>%
  summarise(highlighted_AHBA_overlaps = paste(hit, collapse = "; "), .groups = "drop")

ahba_context <- modules_to_annotate %>%
  left_join(ahba_best %>% select(region, module, best_AHBA_module = ahba_module, best_AHBA_overlap = overlap_count, best_AHBA_p = p_value, best_AHBA_FDR = fdr_BH), by = c("region", "module")) %>%
  left_join(ahba_highlight, by = c("region", "module")) %>%
  left_join(ahba_annotation %>% select(AHBA_module_ID, top_annotation_categories, top_annotation_terms, best_annotation_q_value, short_biological_interpretation), by = c("best_AHBA_module" = "AHBA_module_ID"))
write_csv(ahba_context, file.path(out_dir, "mdd_modules_with_AHBA_annotation_context.csv"))

# Helper summaries.
top_terms <- function(df, reg, mod, n = 3) {
  x <- df %>% filter(region == reg, module == mod) %>% arrange(fdr_BH, p_value) %>% slice_head(n = n) %>% pull(term_label)
  paste(x, collapse = "; ")
}
top_path_terms <- function(reg, mod, n = 3) {
  bind_rows(reactome, kegg, hallmark, c2_cp) %>%
    filter(region == reg, module == mod) %>%
    arrange(fdr_BH, p_value) %>%
    slice_head(n = n) %>%
    pull(term_label) %>% paste(collapse = "; ")
}
cell_theme <- function(reg, mod) {
  x <- cell_enrich %>% filter(region == reg, module == mod) %>% arrange(fdr_BH, p_value) %>% slice_head(n=1)
  if (nrow(x) == 0) "no clear marker signal" else paste0(x$gs_name, " (FDR=", fmt_num(x$fdr_BH), ")")
}
biol_label <- function(bp, cc, path, cell, ahba_terms, hubs_txt) {
  combined <- tolower(paste(bp, cc, path, cell, ahba_terms, hubs_txt, collapse = " "))
  case_when(
    str_detect(combined, "synap|neuron|axon|dendrit|glutamate|gaba|ion channel|calcium") ~ "neuronal/synaptic",
    str_detect(combined, "cilium|ciliary|cilia|flagell") ~ "ciliary/ependymal-like",
    str_detect(combined, "immune|microglia|cytokine|interferon|complement|inflamm") ~ "immune/glial",
    str_detect(combined, "myelin|oligodendro|ensheath|axon ensheath") ~ "myelin/oligodendrocyte",
    str_detect(combined, "vascular|endothelial|angiogenesis|blood vessel") ~ "vascular/endothelial",
    TRUE ~ "mixed/unclear"
  )
}

summary_table <- modules_to_annotate %>%
  left_join(triage %>% select(region, module, module_size, diagnosis_beta, diagnosis_p, diagnosis_FDR, best_Tier2_or_ML_enrichment_set, best_Tier2_or_ML_FDR, top_10_hub_genes, priority_score), by = c("region", "module")) %>%
  left_join(ahba_context %>% select(region, module, best_AHBA_module, best_AHBA_FDR, highlighted_AHBA_context = short_biological_interpretation, highlighted_AHBA_overlaps), by = c("region", "module")) %>%
  rowwise() %>%
  mutate(
    top_GO_BP_terms = top_terms(go_bp, region, module),
    top_GO_CC_terms = top_terms(go_cc, region, module),
    top_pathway_terms = top_path_terms(region, module),
    cell_type_theme = cell_theme(region, module),
    biological_theme_label = biol_label(top_GO_BP_terms, top_GO_CC_terms, top_pathway_terms, cell_type_theme, highlighted_AHBA_context, top_10_hub_genes),
    confidence_label = case_when(
      (!is.na(best_Tier2_or_ML_FDR) & best_Tier2_or_ML_FDR < 0.05) | (!is.na(best_AHBA_FDR) & best_AHBA_FDR < 0.05) ~ "high-priority exploratory",
      priority_score >= 3 ~ "candidate/supportive",
      TRUE ~ "technical/unclear"
    ),
    caution = "Exploratory WGCNA module; no module eigengene survived Diagnosis BH-FDR < 0.10 in the full Step6 run."
  ) %>%
  ungroup() %>%
  transmute(region, module, module_size, diagnosis_beta, diagnosis_p, diagnosis_fdr = diagnosis_FDR, best_Tier2_ML_set = best_Tier2_or_ML_enrichment_set, best_Tier2_ML_FDR = best_Tier2_or_ML_FDR, best_AHBA_module, best_AHBA_FDR, highlighted_AHBA_context = paste(na0(highlighted_AHBA_overlaps), na0(highlighted_AHBA_context), sep = " | "), top_hub_genes = top_10_hub_genes, top_GO_BP_terms, top_GO_CC_terms, top_pathway_terms, cell_type_theme, biological_theme_label, confidence_label, caution)
write_csv(summary_table, file.path(out_dir, "prioritized_modules_biological_summary.csv"))

# Module cards.
for (i in seq_len(nrow(summary_table))) {
  r <- summary_table[i,]
  hub_top <- hub_annotations %>% filter(region == r$region, module == r$module) %>% arrange(desc(as.numeric(abs_kME))) %>% slice_head(n=10) %>% pull(hub_gene) %>% paste(collapse = "; ")
  tier_summary <- module_enrich %>% filter(region == r$region, module == r$module) %>% arrange(fdr_BH, p_value) %>% slice_head(n=5) %>% mutate(txt = paste0(gene_set, " overlap=", overlap_count, " FDR=", fmt_num(fdr_BH))) %>% pull(txt) %>% paste(collapse = "; ")
  ahba_sum <- ahba_context %>% filter(region == r$region, module == r$module) %>% transmute(txt = paste0("best ", best_AHBA_module, " FDR=", fmt_num(best_AHBA_FDR), "; ", na0(highlighted_AHBA_overlaps), "; ", na0(short_biological_interpretation))) %>% pull(txt)
  interpretation <- paste0(r$region, " ", r$module, " is a prioritized exploratory module. Its annotation is consistent with ", r$biological_theme_label, " biology based on enrichment, AHBA context, and hub-gene patterns. However, its module eigengene did not survive BH-FDR correction for Diagnosis, so it should be treated as a follow-up candidate rather than a validated MDD-associated module.")
  lines <- c(
    paste0("# ", r$region, " ", r$module), "",
    paste0("- Module size: ", r$module_size),
    paste0("- Diagnosis beta/p/FDR: ", fmt_num(r$diagnosis_beta), " / ", fmt_num(r$diagnosis_p), " / ", fmt_num(r$diagnosis_fdr)),
    paste0("- Tier2/ML enrichment summary: ", ifelse(tier_summary == "", "none detected", tier_summary)),
    paste0("- AHBA overlap summary: ", ifelse(length(ahba_sum)==0 || is.na(ahba_sum), "none detected", ahba_sum)),
    paste0("- Top hub genes: ", hub_top),
    paste0("- Top GO BP terms: ", r$top_GO_BP_terms),
    paste0("- Top GO CC terms: ", r$top_GO_CC_terms),
    paste0("- Top pathway terms: ", r$top_pathway_terms),
    paste0("- Cell-type marker summary: ", r$cell_type_theme), "",
    "## Careful Interpretation", "", interpretation, "",
    "## Caution", "", "Exploratory module annotation only. No module eigengene survived Diagnosis BH-FDR < 0.10 in the full Step6 run."
  )
  writeLines(lines, file.path(cards_dir, paste0(safe_name(paste(r$region, r$module)), ".md")))
}

# Plots.
save_plot <- function(p, name, width = 10, height = 7) {
  ggsave(file.path(plots_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 220)
  ggsave(file.path(plots_dir, paste0(name, ".pdf")), p, width = width, height = height)
}

plot_go <- bind_rows(
  go_bp %>% mutate(ontology="BP"),
  go_cc %>% mutate(ontology="CC"),
  go_mf %>% mutate(ontology="MF")
) %>%
  group_by(region, module, ontology) %>%
  arrange(fdr_BH, p_value) %>%
  slice_head(n=4) %>%
  ungroup() %>%
  mutate(module_label = paste(region, module, sep="\n"), neglog10_fdr = -log10(pmax(fdr_BH, 1e-300)), term_label = str_trunc(term_label, 55))
if (nrow(plot_go) > 0) {
  p <- ggplot(plot_go, aes(x = neglog10_fdr, y = reorder(term_label, neglog10_fdr), size = overlap_count, color = ontology)) +
    geom_point(alpha = 0.85) + facet_wrap(~module_label, scales="free_y") +
    labs(title="Top GO/MSigDB C5 terms for prioritized MDD WGCNA modules", x="-log10(BH-FDR)", y="GO term", size="Overlap") + theme_bw(base_size=10)
  save_plot(p, "top_GO_terms_dotplot_by_module", width=14, height=10)
}

category_df <- summary_table %>%
  mutate(module_label = paste(region, module, sep="\n")) %>%
  select(module_label, biological_theme_label, confidence_label)
if (nrow(category_df) > 0) {
  p <- ggplot(category_df, aes(x = biological_theme_label, y = module_label, fill = confidence_label)) + geom_tile(color="white") +
    labs(title="Biological theme labels across prioritized modules", x="Theme", y="Module") + theme_bw(base_size=11) + theme(axis.text.x = element_text(angle=35, hjust=1))
  save_plot(p, "biological_category_heatmap_prioritized_modules", width=10, height=6)
}

if (nrow(cell_enrich) > 0) {
  cell_plot <- cell_enrich %>% group_by(region, module) %>% arrange(fdr_BH, p_value) %>% slice_head(n=5) %>% ungroup() %>% mutate(module_label = paste(region, module, sep="\n"), neglog10_fdr = -log10(pmax(fdr_BH, 1e-300)))
  p <- ggplot(cell_plot, aes(x=neglog10_fdr, y=gs_name, size=overlap_count, color=fdr_BH < 0.10)) + geom_point(alpha=0.85) + facet_wrap(~module_label) + theme_bw(base_size=10) + labs(title="Cell-type marker sanity-check enrichment", x="-log10(BH-FDR)", y="Marker class", color="FDR < 0.10")
  save_plot(p, "cell_type_enrichment_dotplot", width=13, height=8)
}

ahba_plot <- ahba_context %>% mutate(module_label = paste(region, module, sep="\n"), neglog10_fdr = -log10(pmax(best_AHBA_FDR, 1e-300)))
if (nrow(ahba_plot) > 0) {
  p <- ggplot(ahba_plot, aes(x=best_AHBA_module, y=module_label, size=best_AHBA_overlap, color=neglog10_fdr)) + geom_point(alpha=0.85) + scale_color_viridis_c() + theme_bw(base_size=11) + labs(title="Best AHBA overlap context for prioritized MDD WGCNA modules", x="Best AHBA module", y="MDD WGCNA module", color="-log10(FDR)", size="Overlap")
  save_plot(p, "AHBA_overlap_annotation_summary_plot", width=9, height=6)
}

hub_plot <- hub_annotations %>% group_by(region, module) %>% arrange(desc(as.numeric(abs_kME))) %>% slice_head(n=12) %>% ungroup() %>% mutate(module_label = paste(region, module, sep="\n"), hub_gene = str_trunc(hub_gene, 18))
if (nrow(hub_plot) > 0) {
  p <- ggplot(hub_plot, aes(x=as.numeric(abs_kME), y=reorder(hub_gene, as.numeric(abs_kME)), fill=short_annotation_category)) + geom_col() + facet_wrap(~module_label, scales="free_y") + theme_bw(base_size=9) + labs(title="Top hub genes and candidate status", x="abs(kME)", y="Hub gene", fill="Category")
  save_plot(p, "hub_gene_annotation_plot_top_modules", width=14, height=10)
}

# Reports.
report_lines <- c(
  "# MDD WGCNA Biological Annotation Report", "",
  "This report annotates prioritized non-grey WGCNA modules from the full Step6 MDD/control run. No module eigengene survived BH-FDR < 0.10 for Diagnosis, so all biological interpretations are exploratory and intended for follow-up prioritization.", "",
  "## Methods and Databases", "",
  paste0("- Modules annotated: ", nrow(summary_table), " prioritized/secondary modules."),
  "- Background: region-specific WGCNA gene universe from the saved best-solution module assignments.",
  "- GO-like enrichment: custom hypergeometric tests over MSigDB C5 GO:BP, GO:CC, and GO:MF gene sets because clusterProfiler was not installed.",
  "- Pathways: MSigDB Hallmark, C2/CP, and available Reactome/KEGG-style MSigDB sets.",
  "- Cell type: built-in minimal marker sanity-check list because no local marker resource was found.",
  "- AHBA context: integrated from Step 11 published AHBA module annotation tables.", "",
  "## Top Biological Themes by Module", "",
  paste0("- ", summary_table$region, " ", summary_table$module, ": ", summary_table$biological_theme_label, "; confidence: ", summary_table$confidence_label, "; top GO BP: ", ifelse(summary_table$top_GO_BP_terms == "", "none", summary_table$top_GO_BP_terms)), "",
  "## Strong Tier2/ML or AHBA Evidence", "",
  paste0("- ", summary_table$region, " ", summary_table$module, ": best Tier2/ML FDR=", fmt_num(summary_table$best_Tier2_ML_FDR), "; best AHBA=", summary_table$best_AHBA_module, " FDR=", fmt_num(summary_table$best_AHBA_FDR)), "",
  "## Limitations", "",
  "- These are prioritized exploratory modules, not validated disease modules.",
  "- No module eigengene survived Diagnosis BH-FDR < 0.10.",
  "- org.Hs.eg.db, clusterProfiler, ReactomePA, and enrichplot were unavailable in this environment; the script used MSigDB and custom enrichment fallbacks where possible.",
  "- The cell-type marker analysis is a minimal sanity check and should not be overinterpreted.",
  "- Enrichment terms suggest biological context but do not establish causality or diagnostic utility.", ""
)
writeLines(report_lines, file.path(out_dir, "MDD_WGCNA_Biological_Annotation_Report.md"))

write_csv(skipped_databases %>% distinct(), file.path(out_dir, "annotation_database_skip_notes.csv"))

manifest <- list(
  generated_at = as.character(Sys.time()),
  project_root = project_root,
  inputs = input_files,
  outputs = as.list(normalizePath(list.files(out_dir, recursive = TRUE, full.names = TRUE), winslash = "/", mustWork = FALSE)),
  packages_available = as.list(package_available),
  databases_unavailable_or_skipped = skipped_databases,
  interpretation_caution = "No module eigengene survived Diagnosis BH-FDR < 0.10; annotations are exploratory.",
  sessionInfo = paste(capture.output(sessionInfo()), collapse = "\n")
)
write_json(manifest, file.path(out_dir, "mdd_wgcna_biological_annotation_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
trim_markdown_outputs(out_dir)

message("Biological annotation complete: ", out_dir)
message("Modules annotated: ", nrow(summary_table))
message("GO BP rows: ", nrow(go_bp), "; GO CC rows: ", nrow(go_cc), "; GO MF rows: ", nrow(go_mf))
message("Reactome rows: ", nrow(reactome), "; KEGG rows: ", nrow(kegg), "; Hallmark rows: ", nrow(hallmark), "; C2/CP rows: ", nrow(c2_cp))
