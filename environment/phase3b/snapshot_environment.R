options(stringsAsFactors = FALSE)
project_root <- "D:/genes/GeneDepresion"
phase3b_root <- file.path(project_root, "environment", "phase3b")
phase3b_lib <- file.path(phase3b_root, "library")
.libPaths(c(phase3b_lib, .libPaths()))

required <- c(
  "WGCNA", "dynamicTreeCut", "readr", "dplyr", "tidyr", "stringr", "limma",
  "matrixStats", "data.table", "jsonlite", "mclust", "clue", "clusterProfiler",
  "org.Hs.eg.db", "AnnotationDbi", "ReactomePA", "msigdbr", "enrichplot", "qvalue", "renv"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing Phase-3B packages: ", paste(missing, collapse = ", "))

versions <- data.frame(
  package = required,
  version = vapply(required, function(x) as.character(packageVersion(x)), character(1)),
  stringsAsFactors = FALSE
)
write.csv(versions, file.path(phase3b_root, "package_versions.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(phase3b_root, "sessionInfo.txt"))

renv::snapshot(
  project = phase3b_root,
  library = phase3b_lib,
  lockfile = file.path(phase3b_root, "renv.lock"),
  type = "explicit",
  packages = required,
  prompt = FALSE,
  force = TRUE
)

cat("Phase-3B environment captured;", nrow(versions), "required packages recorded.\n")
