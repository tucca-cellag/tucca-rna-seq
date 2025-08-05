# workflow/scripts/run_spia.R
#
# Purpose:
# This script performs Signaling Pathway Impact Analysis (SPIA) for a single
# DESeq2 contrast using the SPIA package. It generates pathway impact analysis
# results for KEGG pathways, identifying pathways that are significantly
# perturbed based on both gene expression changes and pathway topology.
#
# Inputs:
#   - dge_tsv: Path to the differential gene expression results from DESeq2.
#   - spia_data: Directory containing pre-generated SPIA data files.
#
# Outputs:
#   - spia_rds: An RDS file containing SPIA analysis results.
#
# Parameters:
#   - org_db_pkg: The name of the organism-specific annotation package.
#   - species: The scientific name of the species (e.g., "Saccharomyces_cerevisiae").
#   - kegg_organism: The 3-letter KEGG organism code (e.g., "sce" for yeast).
#   - extra: Additional parameters to pass to SPIA::spia() function.

# --- 1. Setup ---
# Source the helper functions
print(snakemake)
source(base::file.path(snakemake@scriptdir, "enrichment_utils.R"))

# Unpack snakemake object
log_file <- snakemake@log[[1]]
dge_path <- snakemake@input$dge_tsv
spia_data_path <- snakemake@input$spia_data
output_path <- snakemake@output$spia_rds
enrichment_params <- snakemake@params[["enrichment"]]

# --- 2. Logging and Library Loading ---
setup_logging_and_libs(log_file)

# Load SPIA library
base::library(package = "SPIA", character.only = TRUE)
base::message("--- SPIA library loaded ---")

# --- 3. Load Organism-Specific Database ---
org_db_pkg <- get_and_load_orgdb(enrichment_params)

# --- 4. Load Data and Prepare Gene Lists ---
res_tb <- load_and_map_dge_results(dge_path, org_db_pkg)

res_tb_filtered <- res_tb %>%
  dplyr::filter(!is.na(entrez_id) & !is.na(log2FoldChange)) %>%
  dplyr::distinct(entrez_id, .keep_all = TRUE)

# Prepare gene list for SPIA (Entrez IDs with fold changes)
genelist_fc <- res_tb_filtered$log2FoldChange
names(genelist_fc) <- res_tb_filtered$entrez_id

base::message(base::paste(
  "Created gene list for SPIA with", base::length(genelist_fc),
  "genes."
))

if (base::length(genelist_fc) == 0) {
  base::stop(
    "The gene list for SPIA is empty. Check gene ID mapping and DGE results."
  )
}

# --- 5. Perform SPIA Analysis ---
base::set.seed(123)
base::message("Running SPIA analysis...")

# Check if SPIA data directory exists and contains required files
if (!base::dir.exists(spia_data_path)) {
  base::message("SPIA data directory does not exist: ", spia_data_path)
  base::message("Creating empty results.")
  spia_results <- base::list()
} else {
  kegg_organism <- enrichment_params$kegg_organism
  base::message("Running SPIA for organism: ", kegg_organism)
  base::message("Using SPIA data from: ", spia_data_path)

  # Run SPIA analysis with pre-generated data
  spia_defaults <- base::paste0(
    "de = genelist_fc, all = names(genelist_fc), organism = '",
    kegg_organism, "', data.dir = '", spia_data_path, "'"
  )
  spia_final_args <- base::paste(
    spia_defaults, enrichment_params$spia$extra,
    sep = ", "
  )
  spia_cmd <- base::paste0("SPIA::spia(", spia_final_args, ")")
  base::message("Command: ", spia_cmd)

  tryCatch(
    {
      spia_results <- base::eval(base::parse(text = spia_cmd))
      base::message("SPIA analysis completed successfully.")
    },
    error = function(e) {
      base::message("SPIA analysis failed: ", e$message)
      base::message("Creating empty results.")
      spia_results <- base::list()
    }
  )
}

# --- 6. Save Results ---
base::message("Saving SPIA results to: ", output_path)
base::saveRDS(spia_results, file = output_path)

log_script_completion("SPIA script")
