# workflow/scripts/make_spia_data.R
#
# Purpose:
# This script generates SPIA data files for a specific organism using
# SPIA::makeSPIAdata(). First, it downloads KEGG pathway XML files using
# the KEGG REST API, then processes them into the required data format
# for SPIA analysis.
#
# Inputs:
#   - None (downloads KEGG data automatically)
#
# Outputs:
#   - spia_data_dir: Directory containing the generated SPIA data files
#
# Parameters:
#   - kegg_organism: The 3-letter KEGG organism code (e.g., "sce" for yeast).

# --- 1. Setup ---
# Source the helper functions
print(snakemake)
source(base::file.path(snakemake@scriptdir, "enrichment_utils.R"))

# Unpack snakemake object
log_file <- snakemake@log[[1]]
output_dir <- snakemake@output$spia_data_dir
enrichment_params <- snakemake@params[["enrichment"]]

# --- 2. Logging and Library Loading ---
setup_logging_and_libs(log_file)

# Load SPIA library
base::library(package = "SPIA", character.only = TRUE)
base::message("--- SPIA library loaded ---")

# --- 3. Generate SPIA Data ---
kegg_organism <- enrichment_params$kegg_organism
base::message("Generating SPIA data for organism: ", kegg_organism)

# Create output directory
base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set working directory to output directory for makeSPIAdata
original_wd <- base::getwd()
base::setwd(output_dir)

tryCatch(
  {
    # Step 1: Download pathway list for the organism
    base::message("Downloading pathway list for organism: ", kegg_organism)
    pathway_list_url <- base::paste0("https://rest.kegg.jp/list/pathway/", kegg_organism)

    pathway_list <- base::tryCatch(
      {
        base::readLines(pathway_list_url)
      },
      error = function(e) {
        base::message("Failed to download pathway list: ", e$message)
        base::stop("Cannot proceed without pathway list")
      }
    )

    if (base::length(pathway_list) == 0) {
      base::stop("No pathways found for organism: ", kegg_organism)
    }

    base::message("Found ", base::length(pathway_list), " pathways for ", kegg_organism)

    # Step 2: Download each pathway XML file
    base::message("Downloading pathway XML files...")
    successful_downloads <- 0
    for (pathway_line in pathway_list) {
      # Parse pathway ID from the line (format: "pathway_id\tpathway_name")
      pathway_parts <- base::strsplit(pathway_line, "\t")[[1]]
      pathway_id <- pathway_parts[1]

      # Download the XML file
      xml_url <- base::paste0("https://rest.kegg.jp/get/", pathway_id, "/kgml")
      xml_file <- base::paste0(pathway_id, ".xml")

      base::message("Downloading ", pathway_id, "...")
      download_result <- base::tryCatch(
        {
          utils::download.file(xml_url, xml_file, quiet = TRUE)
          base::message("Successfully downloaded ", pathway_id)
          successful_downloads <- successful_downloads + 1
        },
        error = function(e) {
          base::message("Failed to download ", pathway_id, ": ", e$message)
        }
      )

      # Add a small delay to be respectful to the KEGG server
      base::Sys.sleep(0.4)
    }

    base::message("Successfully downloaded ", successful_downloads, " out of ", base::length(pathway_list), " pathways")

    if (successful_downloads == 0) {
      base::stop("No pathway files were downloaded successfully")
    }

    # Step 3: Run makeSPIAdata on the downloaded files
    base::message("Running makeSPIAdata for organism: ", kegg_organism)
    SPIA::makeSPIAdata(
      organism = kegg_organism,
      out.path = "."
    )
    base::message("SPIA data generation completed successfully.")
  },
  error = function(e) {
    base::message("SPIA data generation failed: ", e$message)
    base::message("Creating empty data directory.")
  }
)

# Restore original working directory
base::setwd(original_wd)

log_script_completion("SPIA data generation script")
