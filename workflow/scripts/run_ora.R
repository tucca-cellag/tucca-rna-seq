# workflow/scripts/run_ora.R
#
# Purpose:
# This script performs Over-Representation Analysis (ORA) for a single
# DESeq2 contrast using clusterProfiler. It generates enrichResult objects for
# Gene Ontology (GO) and KEGG pathways.
#
# Inputs:
#   - dge_tsv: Path to the differential gene expression results from DESeq2.
#
# Outputs:
#   - ora_rds: An RDS file containing a named list of enrichResult objects.
#
# Parameters:
#   - org_db_pkg: The name of the organism-specific annotation package.
#   - species: The scientific name of the species (e.g., "Saccharomyces_cerevisiae").
#   - kegg_organism: The 3-letter KEGG organism code (e.g., "sce" for yeast).
#   - alpha_pathway: The significance threshold for pathway enrichment.
#   - padj_cutoff: The adjusted p-value cutoff to define significant genes.

# --- 1. Setup ---
# Source the helper functions
print(snakemake)
source(base::file.path(snakemake@scriptdir, "enrichment_utils.R"))

# Unpack snakemake object
log_file <- snakemake@log[[1]]
dge_path <- snakemake@input$dge_tsv
output_path <- snakemake@output$ora_rds
enrichment_params <- snakemake@params[["enrichment"]]

# --- 2. Logging and Library Loading ---
setup_logging_and_libs(log_file)

# --- 3. Load Organism-Specific Database ---
org_db_pkg <- get_and_load_orgdb(enrichment_params)

# --- 4. Load Data and Prepare Gene Lists ---
res_tb <- load_and_map_dge_results(dge_path, org_db_pkg)

res_tb_filtered <- res_tb %>%
  dplyr::filter(!is.na(entrez_id)) %>%
  dplyr::distinct(entrez_id, .keep_all = TRUE)

significant_genes <- res_tb_filtered %>%
  dplyr::filter(padj < enrichment_params$padj_cutoff) %>%
  dplyr::pull(entrez_id)
base::message(base::paste(
  "Found", base::length(significant_genes), "significant genes for ORA."
))

universe_genes <- res_tb_filtered$entrez_id
base::message(base::paste(
  "Universe contains", base::length(universe_genes), "genes."
))

if (base::length(significant_genes) == 0) {
  base::message(
    "No significant genes found. ORA will not be performed.",
    "Creating empty output."
  )
  base::saveRDS(base::list(), file = output_path)
  base::quit(save = "no", status = 0)
}

# --- 5. Perform ORA ---
base::set.seed(123)
ora_results <- base::list()

# Check once if the OrgDb has SYMBOL mappings to avoid repeated lookups
orgdb_obj <- base::get(org_db_pkg)
has_symbol_support <- "SYMBOL" %in% AnnotationDbi::columns(orgdb_obj)


# ORA for Gene Ontology (GO)
base::message("Running ORA for GO (BP)...")
enrichgo_defaults <- "gene = significant_genes, universe = universe_genes, OrgDb = get(org_db_pkg), ont = 'BP', keyType = 'ENTREZID'"
enrichgo_final_args <- base::paste(
  enrichgo_defaults, enrichment_params$clusterprofiler$ora$enrichGO$extra,
  sep = ", "
)
enrichgo_cmd <- base::paste0(
  "clusterProfiler::enrichGO(", enrichgo_final_args, ")"
)
base::message("Command: ", enrichgo_cmd)
go_res <- base::eval(base::parse(text = enrichgo_cmd))
ora_results$GO <- safely_set_readable(
  go_res, orgdb_obj,
  has_symbol_support
)

# ORA for KEGG Pathways
base::message("Running ORA for KEGG...")
enrichkegg_defaults <- base::paste0(
  "gene = significant_genes, universe = universe_genes, organism = '",
  enrichment_params$kegg_organism, "'"
)
enrichkegg_final_args <- base::paste(
  enrichkegg_defaults, enrichment_params$clusterprofiler$ora$enrichKEGG$extra,
  sep = ", "
)
enrichkegg_cmd <- base::paste0(
  "clusterProfiler::enrichKEGG(", enrichkegg_final_args, ")"
)
base::message("Command: ", enrichkegg_cmd)
kegg_res <- base::eval(base::parse(text = enrichkegg_cmd))
ora_results$KEGG <- safely_set_readable(
  kegg_res, orgdb_obj,
  has_symbol_support
)

# ORA for KEGG Modules (MKEGG)
if (enrichment_params$clusterprofiler$kegg_module$enabled) {
  base::message("Running ORA for KEGG Modules (MKEGG)...")
  enrichmkegg_defaults <- base::paste0(
    "gene = significant_genes, universe = universe_genes, organism = '",
    enrichment_params$kegg_organism, "'"
  )
  enrichmkegg_final_args <- base::paste(
    enrichmkegg_defaults, enrichment_params$clusterprofiler$kegg_module$enrichMKEGG$extra,
    sep = ", "
  )
  enrichmkegg_cmd <- base::paste0(
    "clusterProfiler::enrichMKEGG(", enrichmkegg_final_args, ")"
  )
  base::message("Command: ", enrichmkegg_cmd)
  mkegg_res <- base::eval(base::parse(text = enrichmkegg_cmd))
  ora_results$MKEGG <- safely_set_readable(
    mkegg_res, orgdb_obj,
    has_symbol_support
  )
}

# ORA for WikiPathways
if (enrichment_params$clusterprofiler$wikipathways$enabled) {
  # Replaces underscores with spaces in species name for matching
  wp_species <- base::gsub("_", " ", enrichment_params$species)
  supported_wp_species <- clusterProfiler::get_wp_organisms()

  if (!wp_species %in% supported_wp_species) {
    base::message(
      "Skipping WikiPathways ORA: species '", wp_species,
      "' not found in clusterProfiler::get_wp_organisms()."
    )
  } else {
    base::message("Running ORA for WikiPathways...")
    enrichwp_defaults <- base::paste0(
      "gene = significant_genes, universe = universe_genes, organism = '",
      wp_species, "'"
    )
    enrichwp_final_args <- base::paste(
      enrichwp_defaults,
      enrichment_params$clusterprofiler$wikipathways$enrichWP$extra,
      sep = ", "
    )
    enrichwp_cmd <- base::paste0(
      "clusterProfiler::enrichWP(", enrichwp_final_args, ")"
    )
    base::message("Command: ", enrichwp_cmd)
    wp_res <- base::eval(base::parse(text = enrichwp_cmd))
    ora_results$WP <- safely_set_readable(
      wp_res, orgdb_obj,
      has_symbol_support
    )
  }
}

# --- 6. Save Results ---
base::message("Saving ORA results to: ", output_path)
base::saveRDS(ora_results, file = output_path)

log_script_completion("ORA script")
