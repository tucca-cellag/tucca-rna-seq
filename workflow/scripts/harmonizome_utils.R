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

#' Create mock Harmonizome data for testing
#'
#' This function creates mock Harmonizome gene sets using MSigDB data
#' for testing purposes, particularly useful for CI/CD pipelines.
#'
#' @param species The species name (e.g., "Bos taurus", "Homo sapiens")
#' @param collections MSigDB collections to use (default: "H" for Hallmarks)
#' @param max_genes_per_set Maximum genes per gene set (default: 50)
#'
#' @return A list of mock Harmonizome gene set data
create_mock_harmonizome_data <- function(species = "Saccharomyces cerevisiae",
                                         collections = "H",
                                         max_genes_per_set = 50) {
  base::message("Creating mock Harmonizome data for species: ", species)

  # Check if msigdbr is available
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    base::stop("msigdbr package is required for mock Harmonizome data")
  }

  # Get MSigDB gene sets for the species
  msig_data <- msigdbr::msigdbr(species = species, category = collections)

  if (base::nrow(msig_data) == 0) {
    base::stop("No MSigDB gene sets found for species: ", species)
  }

  # Create mock Harmonizome format
  mock_data <- base::list()

  # Get unique gene sets
  unique_sets <- base::unique(msig_data$gs_name)

  for (set_name in unique_sets[1:min(5, base::length(unique_sets))]) { # Limit to 5 sets for testing
    set_genes <- msig_data$gene_symbol[msig_data$gs_name == set_name]

    # Limit genes per set for testing
    if (base::length(set_genes) > max_genes_per_set) {
      set_genes <- set_genes[1:max_genes_per_set]
    }

    # Create mock Harmonizome format
    mock_data[[set_name]] <- base::list(
      dataset = "Mock MSigDB Dataset",
      gene_set = set_name,
      genes = set_genes,
      count = base::length(set_genes),
      url = "https://www.gsea-msigdb.org/",
      timestamp = base::Sys.time()
    )
  }

  base::message("Created ", base::length(mock_data), " mock gene sets")
  return(mock_data)
}

#' Create custom yeast gene sets for testing
#'
#' This function creates custom gene sets using common yeast genes
#' that are likely to be found in test data.
#'
#' @param max_genes_per_set Maximum genes per gene set (default: 20)
#'
#' @return A list of custom yeast gene set data
create_custom_yeast_gene_sets <- function(max_genes_per_set = 20) {
  base::message("Creating custom yeast gene sets for testing")

  # Common yeast genes that are likely to be in test data
  yeast_gene_sets <- base::list(
    "HALLMARK_GLYCOLYSIS" = base::c(
      "TEF1", "TEF2", "PGK1", "ENO1", "ENO2", "PYK1", "PYK2",
      "TDH1", "TDH2", "TDH3", "GPM1", "GPM2", "GPM3", "CDC19"
    ),
    "HALLMARK_AMINO_ACID_METABOLISM" = base::c(
      "HIS1", "HIS2", "HIS3", "HIS4", "HIS5", "HIS6", "HIS7",
      "LEU1", "LEU2", "LEU4", "LEU9", "ARG1", "ARG3", "ARG4"
    ),
    "HALLMARK_RIBOSOME" = base::c(
      "RPL1A", "RPL1B", "RPL2A", "RPL2B", "RPL3", "RPL4A", "RPL4B",
      "RPS1A", "RPS1B", "RPS2", "RPS3", "RPS4A", "RPS4B", "RPS5"
    ),
    "HALLMARK_STRESS_RESPONSE" = base::c(
      "HSP12", "HSP26", "HSP30", "HSP42", "HSP78", "HSP82", "HSP104",
      "CTT1", "SOD1", "SOD2", "GPX1", "GPX2", "TRX1", "TRX2"
    ),
    "HALLMARK_CELL_CYCLE" = base::c(
      "CLN1", "CLN2", "CLN3", "CLB1", "CLB2", "CLB3", "CLB4",
      "CDC28", "CDC20", "CDC14", "SWI4", "SWI6", "MBP1", "SWE1"
    )
  )

  # Create mock Harmonizome format
  mock_data <- base::list()

  for (set_name in base::names(yeast_gene_sets)) {
    genes <- yeast_gene_sets[[set_name]]

    # Limit genes per set for testing
    if (base::length(genes) > max_genes_per_set) {
      genes <- genes[1:max_genes_per_set]
    }

    # Create mock Harmonizome format
    mock_data[[set_name]] <- base::list(
      dataset = "Custom Yeast Gene Sets",
      gene_set = set_name,
      genes = genes,
      count = base::length(genes),
      url = "https://www.yeastgenome.org/",
      timestamp = base::Sys.time()
    )
  }

  base::message("Created ", base::length(mock_data), " custom yeast gene sets")
  return(mock_data)
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

  # Check if we should use mock data for testing
  use_mock_data <- harmonizome_config$use_mock_data %||% FALSE

  if (use_mock_data) {
    base::message("Using mock Harmonizome data for testing")
    species <- harmonizome_config$mock_species %||% "Saccharomyces cerevisiae"

    # Use custom yeast gene sets for better testing
    if (species == "Saccharomyces cerevisiae") {
      harmonizome_data <- create_custom_yeast_gene_sets()
    } else {
      harmonizome_data <- create_mock_harmonizome_data(species = species)
    }
  } else {
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
  }

  # Prepare for clusterProfiler
  if (base::length(harmonizome_data) > 0) {
    return(prepare_harmonizome_gene_sets(harmonizome_data, orgdb_pkg_name))
  } else {
    base::message("No Harmonizome gene sets found")
    return(NULL)
  }
}
