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
  base::library(package = "msigdbr", character.only = TRUE)
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

#' Validate MSigDB species support
#'
#' This function checks if the specified species is supported by MSigDB
#' using msigdbr::msigdbr_show_species().
#'
#' @param species The species name in format "Genus_species"
#'
#' @return A logical indicating if the species is supported
validate_msigdb_species <- function(species) {
  base::message("Validating MSigDB species support for: ", species)

  # Convert underscore to space for MSigDB format
  species_formatted <- base::gsub("_", " ", species)

  # Get supported species from MSigDB
  supported_species <- msigdbr::msigdbr_show_species()

  is_supported <- species_formatted %in% supported_species$species_name

  if (is_supported) {
    base::message("Species '", species_formatted, "' is supported by MSigDB.")
  } else {
    base::message(
      "Species '", species_formatted, "' is not supported by MSigDB. ",
      "Available species: ", base::paste(supported_species$species_name[1:5], collapse = ", "), "..."
    )
  }

  return(is_supported)
}

#' Load MSigDB gene sets for specified collections
#'
#' This function loads gene sets from MSigDB collections using msigdbr.
#'
#' @param species The species name in format "Genus_species"
#' @param collections Vector of MSigDB collection codes (H, C1, C2, etc.)
#'
#' @return A data.frame with columns: gs_name, gs_description, entrez_gene
load_msigdb_genesets <- function(species, collections) {
  base::message("Loading MSigDB gene sets for collections: ", base::paste(collections, collapse = ", "))

  # Convert underscore to space for MSigDB format
  species_formatted <- base::gsub("_", " ", species)

  # Load gene sets for each collection
  all_genesets <- base::list()

  for (collection in collections) {
    base::message("Loading collection: ", collection)
    genesets <- msigdbr::msigdbr(species = species_formatted, category = collection)
    all_genesets[[collection]] <- genesets
  }

  # Combine all collections
  combined_genesets <- base::do.call(rbind, all_genesets)

  base::message(
    "Loaded ", base::nrow(combined_genesets), " gene sets from ",
    base::length(collections), " collections."
  )

  return(combined_genesets)
}

#' Load custom GMT files
#'
#' This function loads custom gene set files in GMT format using clusterProfiler.
#'
#' @param gmt_files Vector of paths to GMT files
#'
#' @return A list of gene set data.frames, one per file
load_custom_gmt_files <- function(gmt_files) {
  base::message("Loading custom GMT files: ", base::paste(gmt_files, collapse = ", "))

  gmt_results <- base::list()

  for (gmt_file in gmt_files) {
    if (base::file.exists(gmt_file)) {
      base::message("Loading GMT file: ", gmt_file)
      gmt_data <- clusterProfiler::read.gmt(gmt_file)
      gmt_results[[gmt_file]] <- gmt_data
      base::message("Loaded ", base::length(gmt_data), " gene sets from ", gmt_file)
    } else {
      base::message("Warning: GMT file not found: ", gmt_file)
    }
  }

  return(gmt_results)
}

#' Prepare TERM2GENE and TERM2NAME for universal enrichment analysis
#'
#' This function prepares the TERM2GENE and TERM2NAME data.frames required
#' for clusterProfiler's enricher() and GSEA() functions.
#'
#' @param genesets_data Data.frame from load_msigdb_genesets or load_custom_gmt_files
#' @param source_type Either "msigdb" or "gmt" to indicate the data source
#'
#' @return A list with TERM2GENE and TERM2NAME data.frames
prepare_enrichment_data <- function(genesets_data, source_type = "msigdb") {
  if (source_type == "msigdb") {
    # For MSigDB data
    term2gene <- genesets_data[, c("gs_name", "entrez_gene")]
    term2name <- genesets_data[, c("gs_name", "gs_description")]
    # Remove duplicates
    term2gene <- base::unique(term2gene)
    term2name <- base::unique(term2name)
  } else if (source_type == "gmt") {
    # For GMT data - convert list of gene sets to data.frame
    term2gene_list <- base::list()
    term2name_list <- base::list()

    for (geneset_name in base::names(genesets_data)) {
      genes <- genesets_data[[geneset_name]]
      term2gene_list[[geneset_name]] <- base::data.frame(
        term = geneset_name,
        gene = genes
      )
      term2name_list[[geneset_name]] <- base::data.frame(
        term = geneset_name,
        name = geneset_name # Use name as description for GMT files
      )
    }

    term2gene <- base::do.call(rbind, term2gene_list)
    term2name <- base::do.call(rbind, term2name_list)
  }

  return(base::list(
    TERM2GENE = term2gene,
    TERM2NAME = term2name
  ))
}

#' Process KEGGLINK and convert to human-readable gene symbols
#'
#' This function processes KEGGLINK URLs from SPIA results and converts
#' ENTREZ gene IDs to human-readable gene symbols. It handles cases where
#' the OrgDb doesn't have SYMBOL mappings by falling back to ENTREZ IDs.
#'
#' @param link A KEGGLINK URL string from SPIA results.
#' @param orgdb_obj The OrgDb annotation package object (e.g., org.Hs.eg.db).
#'
#' @return A string containing gene symbols (or ENTREZ IDs if symbols unavailable)
#'         separated by '/'.
processKEGGLINK <- function(link, orgdb_obj) {
  if (base::is.na(link) || link == "") {
    return("")
  }

  # Extract the part of the URL after the last "?" character
  queryString <- stringr::str_remove(link, ".*\\?")

  # Split the string by '+' into parts
  parts <- base::strsplit(queryString, "\\+")[[1]]

  # Assuming the first elem is always the pathway ID, we skip it for gene IDs
  deGenes <- parts[-1] # Extract gene IDs

  if (base::length(deGenes) == 0) {
    return("")
  }

  # Check if the OrgDb has SYMBOL mappings
  has_symbol_support <- "SYMBOL" %in% AnnotationDbi::columns(orgdb_obj)

  if (has_symbol_support) {
    # Convert ENTREZ IDs to gene symbols
    deGenesSym <- AnnotationDbi::select(
      orgdb_obj,
      keys = deGenes,
      columns = "SYMBOL",
      keytype = "ENTREZID"
    )

    # Filter out NA symbols and combine
    valid_symbols <- deGenesSym$SYMBOL[!base::is.na(deGenesSym$SYMBOL)]
    if (base::length(valid_symbols) == 0) {
      return("")
    }

    # Combine gene symbols into a single string separated by '/'
    geneString <- base::paste(valid_symbols, collapse = "/")

    return(geneString)
  } else {
    # If no SYMBOL support, return the original ENTREZ IDs
    base::message(
      "The OrgDb package does not contain SYMBOL mappings. ",
      "Using ENTREZ IDs for gene names."
    )
    return(base::paste(deGenes, collapse = "/"))
  }
}
