# workflow/scripts/install_orgdb.R
#
# Purpose:
# This script installs the organism-specific OrgDb package required for
# enrichment analysis. It determines the package name, checks if it's
# already installed, and installs it from Bioconductor or a local source
# if needed. This script is intended to be run by a dedicated Snakemake
# rule before any analysis rules that depend on the package.
#
# Inputs: None
# Outputs: None (creates a flag file via the Snakemake rule's 'touch' command)
#
# Parameters:
#   - A single 'enrichment' object containing all enrichment-related params.

# --- 1. Setup and Logging ---
log <- base::file(snakemake@log[[1]], open = "wt")
base::sink(log)
base::sink(log, type = "message")
base::date()

base::library(package = "devtools", character.only = TRUE)
devtools::session_info()

# --- 2. Determine Package Name and Install if Necessary ---
enrichment_params <- snakemake@params[["enrichment"]]
install_method <- enrichment_params$install_method
install_source <- enrichment_params$install_source

# Determine the actual package name
org_db_pkg <- if (install_method == "local") {
  # When built locally, the package name is derived inside AnnotationForge.
  # We find the directory name it created to get the actual package name.
  pkg_dirs <- base::list.dirs(path = install_source, full.names = FALSE, recursive = FALSE)
  pkg_name <- pkg_dirs[grepl("^org\\..+\\.db$", pkg_dirs)]
  if (length(pkg_name) == 0) base::stop("Could not find locally built OrgDb package directory.")
  pkg_name[1]
} else {
  enrichment_params$org_db_pkg
}
base::message("Target OrgDb package: ", org_db_pkg)

# Install if necessary
if (!base::require(org_db_pkg, character.only = TRUE)) {
  base::message(base::paste("Package", org_db_pkg, "not found, installing..."))
  base::message(base::paste("Method:", install_method, "| Source:", install_source))

  if (install_method == "local") {
    devtools::install(file.path(install_source, org_db_pkg), upgrade = "never")
  } else {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(install_source, update = FALSE, ask = FALSE)
  }

  # Verify that the package was installed successfully
  if (!base::require(org_db_pkg, character.only = TRUE)) {
    base::stop(paste("Failed to install package:", org_db_pkg))
  }
  base::message("Package successfully installed.")
} else {
  base::message("Package already installed.")
}

base::message("Installation script finished successfully.")
base::date()
