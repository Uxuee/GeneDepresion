options(stringsAsFactors = FALSE)

project_root <- "D:/genes/GeneDepresion"
.libPaths(c(file.path(project_root, "environment", "phase3b", "library"), .libPaths()))
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(GO.db)
})

source_path <- file.path(
  project_root, "outputs", "length_gt_100_clean_pipeline",
  "17_WGCNA_reconstruction_validation", "functional_annotation",
  "functional_enrichment_GO_Reactome_Hallmark.csv"
)
output_dir <- file.path(
  project_root, "outputs", "length_gt_100_clean_pipeline",
  "18_WGCNA_enrichment_interpretation_audit", "functional_annotation"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

terms <- read_csv(source_path, show_col_types = FALSE) %>%
  filter(database == "GO", BH_q < 0.10) %>%
  distinct(term_id, ontology)

parent_map <- list(BP = GOBPPARENTS, MF = GOMFPARENTS, CC = GOCCPARENTS)
rows <- list()
for (i in seq_len(nrow(terms))) {
  child <- terms$term_id[i]
  ontology <- terms$ontology[i]
  parents <- mget(child, parent_map[[ontology]], ifnotfound = NA)[[1]]
  if (!length(parents) || all(is.na(parents))) next
  relationship <- names(parents)
  if (is.null(relationship)) relationship <- rep(NA_character_, length(parents))
  rows[[length(rows) + 1]] <- tibble(
    child_term_id = child,
    parent_term_id = unname(parents),
    relationship = relationship,
    ontology = ontology
  )
}

out <- if (length(rows)) bind_rows(rows) else tibble(
  child_term_id = character(), parent_term_id = character(),
  relationship = character(), ontology = character()
)
write_csv(out, file.path(output_dir, "phase3c_go_immediate_parent_relations.csv"))
cat("GO terms audited:", nrow(terms), "\n")
cat("Immediate parent relations exported:", nrow(out), "\n")
