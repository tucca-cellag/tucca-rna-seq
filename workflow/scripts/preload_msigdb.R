# workflow/scripts/preload_msigdb.R
#
# Purpose:
# This script pre-loads MSigDB gene sets to prevent race conditions
# during concurrent enrichment analysis. It downloads and caches the
# MSigDB data before the main enrichment analysis runs.
#
# Inputs:
#   - None (uses parameters from Snakemake)
#
# Outputs:
#   - A flag file indicating successful pre-loading
#
# Parameters:
#   - species: The species name in format "Genus_species"
#   - collections: Vector of MSigDB collection codes (H, C1, C2, etc.)

# --- 1. Setup, Logging, and Library Loading ---
# Setup logging
log <- base::file(snakemake@log[[1]], open = "wt")
base::sink(log)
base::sink(log, type = "message")
base::date()

# Log session info for reproducibility
base::library(package = "devtools", character.only = TRUE)
devtools::session_info()

# Load required libraries
base::library(package = "msigdbr", character.only = TRUE)
base::library(package = "dplyr", character.only = TRUE)

# --- 2. Extract Parameters ---
species <- snakemake@params$species
collections <- snakemake@params$collections

base::message("Pre-loading MSigDB data for species: ", species)
base::message("Collections to load: ", base::paste(collections, collapse = ", "))

# --- 3. Pre-load MSigDB Data ---
# Convert underscore to space for MSigDB format
species_formatted <- base::gsub("_", " ", species)

# Load gene sets for each collection to trigger download and caching
for (collection in collections) {
  base::message("Pre-loading collection: ", collection)

  # This will trigger the download and caching of MSigDB data
  genesets <- msigdbr::msigdbr(species = species_formatted, category = collection)

  base::message("Loaded ", base::nrow(genesets), " gene sets from collection: ", collection)
}

# --- 4. Create Output Directory and Flag ---
output_dir <- "resources/msigdb"
if (!base::dir.exists(output_dir)) {
  base::dir.create(output_dir, recursive = TRUE)
}

# Create flag file to indicate successful pre-loading
flag_file <- snakemake@output[[1]]
base::file.create(flag_file)

base::message("MSigDB data pre-loaded successfully.")
base::message("Flag file created: ", flag_file)

# --- 5. Cleanup ---
base::message("Script finished successfully.")
base::date()
base::sink()
base::sink(type = "message")
