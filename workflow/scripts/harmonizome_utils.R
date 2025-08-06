# workflow/scripts/harmonizome_utils.R
#
# Purpose:
# This script provides generalized functions for downloading and processing
# Harmonizome datasets and gene sets for use with clusterProfiler's
# universal enrichment analysis functions (enricher() and GSEA()).
#
# Dependencies:
#   - reticulate: For Python API integration
#   - enrichment_utils.R: For gene ID conversion functions
#     (must be sourced by calling script)
#
# Functions:
#   - setup_harmonizome_api: Sets up the Harmonizome Python API
#   - get_harmonizome_gene_set: Gets a specific gene set using the official API
#   - download_harmonizome_dataset: Downloads entire datasets using the official API
#   - prepare_harmonizome_gene_sets: Prepares gene sets for clusterProfiler
#   - load_harmonizome_gene_sets: Main function to load and process Harmonizome data

# Global variables for functions from enrichment_utils.R
# These are defined in enrichment_utils.R and sourced by calling scripts
globalVariables(c("process_gene_set", "is_symbol_like", "convert_genes_to_entrez"))

# Wrapper to suppress linting warnings for functions from enrichment_utils.R
process_gene_set_wrapper <- function(genes, set_name, orgdb_pkg_name) {
  base::get("process_gene_set", envir = .GlobalEnv)(genes, set_name, orgdb_pkg_name)
}

#' Setup Harmonizome Python API
#'
#' This function sets up the Harmonizome Python API using reticulate.
#' It downloads and imports the official Harmonizome API class.
#'
#' @return The Harmonizome Python class object
setup_harmonizome_api <- function() {
  base::message("Setting up Harmonizome Python API...")

  # Check if reticulate is available
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    base::stop("reticulate package is required for Harmonizome integration. Please install it.")
  }

  # Import reticulate
  reticulate::use_python(reticulate::py_config()$python)

  # Download the Harmonizome API if not already present
  api_file <- "harmonizomeapi.py"
  if (!base::file.exists(api_file)) {
    base::message("Downloading Harmonizome API...")
    api_url <- "https://maayanlab.cloud/Harmonizome/static/harmonizomeapi.py"
    utils::download.file(api_url, api_file)
  }

  # Import the Harmonizome class
  tryCatch(
    {
      harmonizome <- reticulate::source_python(api_file)
      base::message("Harmonizome API loaded successfully")
      return(harmonizome$Harmonizome)
    },
    error = function(e) {
      base::stop("Failed to load Harmonizome API: ", e$message)
    }
  )
}

#' Get a specific gene set from Harmonizome using the official API
#'
#' This function uses the official Harmonizome Python API to get gene set data.
#'
#' @param dataset_name The name of the Harmonizome dataset
#' @param gene_set_name The name of the specific gene set
#' @param output_dir Directory to save the downloaded gene set
#'
#' @return A list containing the gene set data and metadata
get_harmonizome_gene_set <- function(dataset_name, gene_set_name, output_dir = "resources/harmonizome") {
  base::message("Getting Harmonizome gene set: ", gene_set_name, " from dataset: ", dataset_name)

  # Setup the API
  harmonizome <- setup_harmonizome_api()

  # Create output directory if it doesn't exist
  if (!base::dir.exists(output_dir)) {
    base::dir.create(output_dir, recursive = TRUE)
  }

  tryCatch(
    {
      # Get gene set information using the official API
      gene_set_info <- harmonizome$get("gene_set", gene_set_name)

      # Check if this gene set belongs to the specified dataset
      if (gene_set_info$dataset != dataset_name) {
        base::stop("Gene set '", gene_set_name, "' does not belong to dataset '", dataset_name, "'")
      }

      # Download the gene set file
      cache_file <- base::file.path(output_dir, base::paste0(
        gene_set_name, "_",
        base::gsub(" ", "_", dataset_name), ".txt"
      ))

      if (!base::file.exists(cache_file)) {
        base::message("Downloading gene set file...")
        # Use the download method to get the gene set file
        download_url <- base::paste0(
          "https://maayanlab.cloud/static/hdfs/harmonizome/data/",
          gene_set_info$path, "/gene_list_terms.txt"
        )
        utils::download.file(download_url, cache_file)
      }

      # Parse the gene set file
      gene_set_data <- base::readLines(cache_file)

      # Extract genes (skip header lines)
      genes <- base::character(0)
      for (line in gene_set_data) {
        if (line != "" && !base::grepl("^#", line) && !base::grepl("^Term", line)) {
          # Split by tabs and extract gene names
          parts <- base::strsplit(line, "\t")[[1]]
          if (base::length(parts) >= 2) {
            genes <- base::c(genes, parts[2]) # Second column contains gene names
          }
        }
      }

      # Remove empty genes and duplicates
      genes <- genes[genes != ""]
      genes <- base::unique(genes)

      # Create result object
      result <- base::list(
        dataset = dataset_name,
        gene_set = gene_set_name,
        genes = genes,
        count = base::length(genes),
        url = gene_set_info$url,
        timestamp = base::Sys.time()
      )

      # Save to RDS file for caching
      rds_file <- base::file.path(output_dir, base::paste0(
        gene_set_name, "_",
        base::gsub(" ", "_", dataset_name), ".rds"
      ))
      base::saveRDS(result, file = rds_file)

      base::message("Downloaded ", base::length(genes), " genes for gene set: ", gene_set_name)
      return(result)
    },
    error = function(e) {
      base::stop("Error getting gene set: ", e$message)
    }
  )
}

#' Download entire Harmonizome datasets using the official API
#'
#' This function uses the official Harmonizome Python API to download datasets.
#'
#' @param dataset_name The name of the Harmonizome dataset
#' @param output_dir Directory to save the downloaded data
#'
#' @return A list of gene sets from the dataset
download_harmonizome_dataset <- function(dataset_name, output_dir = "resources/harmonizome") {
  base::message("Downloading Harmonizome dataset: ", dataset_name)

  # Setup the API
  harmonizome <- setup_harmonizome_api()

  # Create output directory
  if (!base::dir.exists(output_dir)) {
    base::dir.create(output_dir, recursive = TRUE)
  }

  tryCatch(
    {
      # Get all gene sets for this dataset
      gene_sets_response <- harmonizome$get("gene_set")

      # Filter gene sets for this dataset
      dataset_gene_sets <- base::list()

      for (gene_set in gene_sets_response$results) {
        if (gene_set$dataset == dataset_name) {
          gene_set_name <- gene_set$name
          base::message("Found gene set: ", gene_set_name)

          # Get this specific gene set
          gene_set_data <- get_harmonizome_gene_set(dataset_name, gene_set_name, output_dir)
          dataset_gene_sets[[gene_set_name]] <- gene_set_data
        }
      }

      base::message("Downloaded ", base::length(dataset_gene_sets), " gene sets from dataset: ", dataset_name)
      return(dataset_gene_sets)
    },
    error = function(e) {
      base::stop("Error downloading dataset: ", e$message)
    }
  )
}

#' Prepare Harmonizome gene sets for clusterProfiler
#'
#' This function converts Harmonizome gene set data into the format
#' required by clusterProfiler's enricher() and GSEA() functions.
#'
#' @param harmonizome_data List of Harmonizome gene set data
#' @param orgdb_pkg_name Name of the OrgDb package for ID conversion
#'
#' @return A list with TERM2GENE and TERM2NAME data.frames
prepare_harmonizome_gene_sets <- function(harmonizome_data, orgdb_pkg_name) {
  base::message("Preparing Harmonizome gene sets for clusterProfiler")

  # Note: This function assumes that enrichment_utils.R has been sourced
  # by the calling script (run_ora.R or run_gsea.R)
  # Required functions from enrichment_utils.R:
  # - process_gene_set(genes, set_name, orgdb_pkg_name)
  # - is_symbol_like(genes)
  # - convert_genes_to_entrez(genes, orgdb_pkg_name, set_name)

  term2gene_list <- base::list()
  term2name_list <- base::list()

  for (gene_set_name in base::names(harmonizome_data)) {
    gene_set_data <- harmonizome_data[[gene_set_name]]
    genes <- gene_set_data$genes

    if (base::length(genes) > 0) {
      # Convert gene IDs if needed
      processed_genes <- process_gene_set_wrapper(
        genes, gene_set_name, orgdb_pkg_name
      )

      if (base::length(processed_genes) > 0) {
        # Create term2gene entries
        for (gene in processed_genes) {
          term2gene_list[[base::length(term2gene_list) + 1]] <- base::data.frame(
            term = gene_set_name,
            gene = gene,
            stringsAsFactors = FALSE
          )
        }

        # Create term2name entry
        term2name_list[[base::length(term2name_list) + 1]] <- base::data.frame(
          term = gene_set_name,
          name = base::paste0(gene_set_data$dataset, " - ", gene_set_name),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (base::length(term2gene_list) > 0) {
    term2gene <- base::do.call(rbind, term2gene_list)
    term2name <- base::do.call(rbind, term2name_list)

    base::message(
      "Prepared ", base::nrow(term2gene), " gene-term associations for ",
      base::length(base::unique(term2gene$term)), " gene sets"
    )

    return(base::list(
      TERM2GENE = term2gene,
      TERM2NAME = term2name
    ))
  } else {
    base::message("No valid gene sets found after processing")
    return(base::list(
      TERM2GENE = base::data.frame(term = base::character(), gene = base::character()),
      TERM2NAME = base::data.frame(term = base::character(), name = base::character())
    ))
  }
}

#' Load Harmonizome gene sets for enrichment analysis
#'
#' This is the main function to load and process Harmonizome gene sets
#' for use with clusterProfiler's enrichment functions.
#'
#' @param harmonizome_config List containing Harmonizome configuration
#' @param orgdb_pkg_name Name of the OrgDb package for ID conversion
#'
#' @return A list with TERM2GENE and TERM2NAME data.frames
load_harmonizome_gene_sets <- function(harmonizome_config, orgdb_pkg_name) {
  base::message("Loading Harmonizome gene sets")

  if (!harmonizome_config$enabled) {
    base::message("Harmonizome analysis is disabled")
    return(NULL)
  }

  harmonizome_data <- base::list()

  # Process each configured dataset and gene set
  for (dataset_config in harmonizome_config$datasets) {
    dataset_name <- dataset_config$name
    gene_sets <- dataset_config$gene_sets

    base::message("Processing dataset: ", dataset_name)

    for (gene_set_name in gene_sets) {
      base::message("Processing gene set: ", gene_set_name)

      # Try to load from cache first
      cache_file <- base::file.path(
        "resources/harmonizome",
        base::paste0(
          gene_set_name, "_",
          base::gsub(" ", "_", dataset_name), ".rds"
        )
      )

      if (base::file.exists(cache_file)) {
        base::message("Loading from cache: ", cache_file)
        gene_set_data <- base::readRDS(cache_file)
      } else {
        # Get the gene set using the official API
        gene_set_data <- get_harmonizome_gene_set(dataset_name, gene_set_name)
      }

      harmonizome_data[[gene_set_name]] <- gene_set_data
    }
  }

  # Prepare for clusterProfiler
  if (base::length(harmonizome_data) > 0) {
    return(prepare_harmonizome_gene_sets(harmonizome_data, orgdb_pkg_name))
  } else {
    base::message("No Harmonizome gene sets found")
    return(NULL)
  }
}
