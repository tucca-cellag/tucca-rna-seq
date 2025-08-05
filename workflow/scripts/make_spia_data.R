# workflow/scripts/make_spia_data.R
#
# Purpose:
# This script generates SPIA data files for a specific organism using
# SPIA::makeSPIAdata(). This function processes KEGG XML files into
# the required data format for SPIA analysis.
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

# Check if the organism is supported by SPIA
supported_organisms <- SPIA::getSPIAOrganisms()

if (!kegg_organism %in% supported_organisms) {
  base::message(
    "Organism '", kegg_organism, "' not found in SPIA supported organisms. ",
    "Available organisms: ", base::paste(supported_organisms, collapse = ", ")
  )
  base::message("Creating empty data directory.")
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
} else {
  base::message("Running makeSPIAdata for organism: ", kegg_organism)

  # Create output directory
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Set working directory to output directory for makeSPIAdata
  original_wd <- base::getwd()
  base::setwd(output_dir)

  tryCatch(
    {
      # Run makeSPIAdata
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
}

log_script_completion("SPIA data generation script")
