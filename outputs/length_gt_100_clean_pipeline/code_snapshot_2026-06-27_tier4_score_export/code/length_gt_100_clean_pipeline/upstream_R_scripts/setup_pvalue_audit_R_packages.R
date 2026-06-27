## Setup helper for the Length > 100 p-value correction audit.
##
## This script installs the optional Bioconductor packages used by the fixed
## audit when they are missing:
##   - qvalue: Storey q-value correction
##   - IHW: Independent Hypothesis Weighting

options(stringsAsFactors = FALSE)

packages_needed <- c("qvalue", "IHW")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("Installing BiocManager...")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

missing_packages <- packages_needed[
  !vapply(packages_needed, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  message("Installing missing Bioconductor packages: ", paste(missing_packages, collapse = ", "))
  BiocManager::install(missing_packages, ask = FALSE, update = FALSE)
} else {
  message("All requested packages were already installed.")
}

availability <- data.frame(
  package = packages_needed,
  available_after_setup = vapply(packages_needed, requireNamespace, logical(1), quietly = TRUE),
  stringsAsFactors = FALSE
)

print(availability)

if (!all(availability$available_after_setup)) {
  warning(
    "Some optional p-value audit packages are still unavailable after setup: ",
    paste(availability$package[!availability$available_after_setup], collapse = ", ")
  )
}
