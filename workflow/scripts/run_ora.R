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
    enrichmkegg_defaults,
    enrichment_params$clusterprofiler$kegg_module$enrichMKEGG$extra,
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

# ORA for MSigDB
if (enrichment_params$msigdb$enabled) {
  base::message("Running ORA for MSigDB...")

  # Validate species support
  species_supported <- validate_msigdb_species(enrichment_params$species)

  if (!species_supported) {
    base::message(
      "Skipping MSigDB ORA: species '", enrichment_params$species,
      "' not supported by MSigDB."
    )
  } else {
    # Load MSigDB gene sets
    msigdb_genesets <- load_msigdb_genesets(
      enrichment_params$species,
      enrichment_params$msigdb$collections
    )

    # Prepare enrichment data
    msigdb_data <- prepare_enrichment_data(msigdb_genesets, "msigdb")

    # Run ORA for each collection
    for (collection in enrichment_params$msigdb$collections) {
      base::message("Running ORA for MSigDB collection: ", collection)

      # Filter for this collection
      collection_genesets <- msigdb_genesets[msigdb_genesets$gs_cat == collection, ]
      collection_data <- prepare_enrichment_data(collection_genesets, "msigdb")

      # Run enricher
      enricher_defaults <- "gene = significant_genes, universe = universe_genes"
      enricher_final_args <- base::paste(
        enricher_defaults,
        "TERM2GENE = collection_data$TERM2GENE",
        "TERM2NAME = collection_data$TERM2NAME",
        enrichment_params$msigdb$ora$extra,
        sep = ", "
      )
      enricher_cmd <- base::paste0(
        "clusterProfiler::enricher(", enricher_final_args, ")"
      )
      base::message("Command: ", enricher_cmd)

      collection_res <- base::eval(base::parse(text = enricher_cmd))
      ora_results[[paste0("MSigDB_", collection)]] <- safely_set_readable(
        collection_res, orgdb_obj,
        has_symbol_support
      )
    }

    # Load and process custom GMT files if specified
    if (base::length(enrichment_params$msigdb$custom_gmt_files) > 0) {
      base::message("Processing custom GMT files for ORA...")

      for (gmt_file in enrichment_params$msigdb$custom_gmt_files) {
        if (base::file.exists(gmt_file)) {
          base::message("Running ORA for custom GMT file: ", gmt_file)

          # Load GMT file
          gmt_data <- clusterProfiler::read.gmt(gmt_file)

          # Convert to TERM2GENE format
          term2gene_list <- base::list()
          term2name_list <- base::list()

          for (geneset_name in base::names(gmt_data)) {
            genes <- gmt_data[[geneset_name]]
            term2gene_list[[geneset_name]] <- base::data.frame(
              term = geneset_name,
              gene = genes
            )
            term2name_list[[geneset_name]] <- base::data.frame(
              term = geneset_name,
              name = geneset_name
            )
          }

          term2gene <- base::do.call(rbind, term2gene_list)
          term2name <- base::do.call(rbind, term2name_list)

          # Run enricher
          enricher_defaults <- "gene = significant_genes, universe = universe_genes"
          enricher_final_args <- base::paste(
            enricher_defaults,
            "TERM2GENE = term2gene",
            "TERM2NAME = term2name",
            enrichment_params$msigdb$ora$extra,
            sep = ", "
          )
          enricher_cmd <- base::paste0(
            "clusterProfiler::enricher(", enricher_final_args, ")"
          )
          base::message("Command: ", enricher_cmd)

          gmt_res <- base::eval(base::parse(text = enricher_cmd))

          # Create a safe name for the result
          gmt_name <- base::gsub("^.*/", "", base::gsub("\\.gmt$", "", gmt_file))
          ora_results[[paste0("CustomGMT_", gmt_name)]] <- safely_set_readable(
            gmt_res, orgdb_obj,
            has_symbol_support
          )
        } else {
          base::message("Warning: Custom GMT file not found: ", gmt_file)
        }
      }
    }
  }
}

# --- 6. Save Results ---
base::message("Saving ORA results to: ", output_path)
base::saveRDS(ora_results, file = output_path)

log_script_completion("ORA script")
