options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 1200)
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
library_path <- file.path(project_root, "environment", "phase3b", "library")
.libPaths(c(library_path, .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = library_path)
}

ordered <- c(
  "BiocGenerics", "impute", "preprocessCore", "minet", "KEGGREST",
  "GO.db", "org.Hs.eg.db", "DOSE", "clusterProfiler", "reactome.db",
  "ReactomePA", "enrichplot"
)
for (package in ordered) {
  if (!requireNamespace(package, quietly = TRUE)) {
    message("Installing ", package)
    BiocManager::install(package, lib = library_path, ask = FALSE, update = FALSE)
  }
  if (!requireNamespace(package, quietly = TRUE)) {
    warning("Package remains unavailable after retry: ", package)
  }
}

writeLines(capture.output(sessionInfo()), file.path(project_root, "environment", "phase3b", "sessionInfo_after_install.txt"))
