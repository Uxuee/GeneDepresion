options(repos = c(CRAN = "https://cloud.r-project.org"))

project_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
library_path <- file.path(project_root, "environment", "phase3b", "library")
dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library_path, .libPaths()))

cran_packages <- c(
  "renv", "WGCNA", "dynamicTreeCut", "readr", "dplyr", "tidyr",
  "stringr", "jsonlite", "mclust", "clue", "msigdbr"
)
missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  install.packages(missing_cran, lib = library_path, dependencies = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = library_path)
}

bioc_packages <- c(
  "limma", "AnnotationDbi", "org.Hs.eg.db", "clusterProfiler",
  "ReactomePA", "enrichplot", "qvalue"
)
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, lib = library_path, ask = FALSE, update = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(project_root, "environment", "phase3b", "sessionInfo_after_install.txt"))
