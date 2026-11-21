packages <- c("edgeR", "limma", "DESeq2")
checks <- lapply(packages, function(package) tryCatch({
  loadNamespace(package)
  list(available = TRUE, error = "")
}, error = function(error) list(available = FALSE, error = conditionMessage(error))))
status <- vapply(checks, `[[`, logical(1), "available")
print(data.frame(package = packages, available = status, error = vapply(checks, `[[`, character(1), "error")))
if (!all(status)) quit(status = 2)
print(vapply(packages, function(package) as.character(packageVersion(package)), character(1)))
