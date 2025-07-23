# workflow/scripts/enrichment_utils.R
#
# Purpose:
# This script provides a helper function to safely apply clusterProfiler::setReadable,
# ensuring that the function is only called when the provided OrgDb package
# contains SYMBOL mappings. This avoids errors with custom-built or
# non-standard annotation packages.
#
# Functions:
#   - safely_set_readable: Checks for SYMBOL support before mapping gene IDs.

#' Safely apply setReadable to an enrichment result
#'
#' This function checks if the provided OrgDb object contains a 'SYMBOL' column.
#' If it does, it runs clusterProfiler::setReadable to map gene IDs to symbols.
#' If not, it returns the original enrichment result object with a warning message.
#'
#' @param enrich_result A `gseaResult` or `enrichResult` object from clusterProfiler.
#' @param orgdb_obj The OrgDb annotation package object (e.g., org.Hs.eg.db).
#' @param has_symbol_support A logical value indicating if the OrgDb contains SYMBOL mappings.
#'
#' @return An enrichment result object, with gene IDs mapped to symbols if possible.
safely_set_readable <- function(enrich_result, orgdb_obj, has_symbol_support) {
  if (has_symbol_support) {
    clusterProfiler::setReadable(
      enrich_result,
      OrgDb = orgdb_obj,
      keyType = "ENTREZID"
    )
  } else {
    base::message(
      "The OrgDb package does not contain SYMBOL mappings. ",
      "Returning results with original IDs."
    )
    enrich_result
  }
}

#' Set up logging and load common libraries
#'
#' This function initializes the log file sink and loads all the R packages
#' that are required by both the ORA and GSEA scripts. It also logs the
#' current session information for reproducibility.
#'
#' @param log_file The path to the log file from the Snakemake object.
#'
setup_logging_and_libs <- function(log_file) {
  # Setup logging
  log <- base::file(log_file, open = "wt")
  base::sink(log)
  base::sink(log, type = "message")
  base::date()
  base::message("--- Logging setup complete ---")

  # Load libraries and log session info
  base::message("--- Loading libraries and logging session info ---")
  base::library(package = "devtools", character.only = TRUE)
  devtools::session_info()

  base::library(package = "clusterProfiler", character.only = TRUE)
  base::library(package = "AnnotationDbi", character.only = TRUE)
  base::library(package = "magrittr", character.only = TRUE)
  base::library(package = "dplyr", character.only = TRUE)
  base::library(package = "readr", character.only = TRUE)
  base::message("--- Library loading complete ---")
}

#' Determine the name of and load the OrgDb package
#'
#' This function inspects the Snakemake parameters to determine if the OrgDb
#' package should be sourced from a local build or from Bioconda, determines
#' its name, and loads it into the session.
#'
#' @param enrichment_params The `enrichment` section of the Snakemake params.
#'
#' @return The name of the loaded OrgDb package as a string.
get_and_load_orgdb <- function(enrichment_params) {
  base::message("--- Determining and loading OrgDb package ---")
  install_method <- enrichment_params$install_method
  install_source <- enrichment_params$install_source

  org_db_pkg <- if (install_method == "local") {
    pkg_dirs <- base::list.dirs(
      path = install_source, full.names = FALSE, recursive = FALSE
    )
    pkg_name <- pkg_dirs[grepl("^org\\..+\\.db$", pkg_dirs)]
    if (length(pkg_name) == 0) {
      base::stop("Could not find locally built OrgDb package directory.")
    }
    pkg_name[1]
  } else {
    enrichment_params$org_db_pkg
  }

  base::message("Target OrgDb package: ", org_db_pkg)
  base::library(org_db_pkg, character.only = TRUE)
  base::message("--- OrgDb package loaded ---")

  return(org_db_pkg)
}

#' Load DGE results and map ENSEMBL IDs to Entrez IDs
#'
#' This function reads the TSV file of differential gene expression results,
#' and uses the provided OrgDb package to map the ENSEMBL gene IDs from the
#' first column to Entrez IDs, which are required for enrichment analysis.
#'
#' @param dge_path The path to the DGE results TSV file.
#' @param orgdb_pkg_name The name of the OrgDb package to use for ID mapping.
#'
#' @return A tibble of the DGE results with an `entrez_id` column added.
load_and_map_dge_results <- function(dge_path, orgdb_pkg_name) {
  base::message("--- Loading DGE results and mapping IDs ---")
  base::message("Loading DGE results from: ", dge_path)
  res_tb <- readr::read_tsv(dge_path, show_col_types = FALSE) %>%
    dplyr::rename(feature_id = 1)

  base::message("Mapping feature IDs to Entrez IDs...")
  orgdb_obj <- base::get(orgdb_pkg_name)
  res_tb$entrez_id <- AnnotationDbi::mapIds(
    orgdb_obj,
    keys = res_tb$feature_id,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  base::message("--- DGE loading and mapping complete ---")
  return(res_tb)
}

#' Log a script completion message
#'
#' This function prints a standardized success message to the log, including
#' the name of the script that finished and the current date and time.
#'
#' @param script_name The name of the script that has completed.
log_script_completion <- function(script_name) {
  base::message(script_name, " finished successfully.")
  base::date()
}
