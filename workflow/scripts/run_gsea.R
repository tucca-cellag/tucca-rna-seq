# workflow/scripts/run_gsea.R
#
# Purpose:
# This script performs Gene Set Enrichment Analysis (GSEA) for a single
# DESeq2 contrast using clusterProfiler. It generates gseaResult objects for
# Gene Ontology (GO) and KEGG pathways.
#
# Inputs:
#   - dge_tsv: Path to the differential gene expression results from DESeq2.
#
# Outputs:
#   - gsea_rds: An RDS file containing a named list of gseaResult objects.
#
# Parameters:
#   - org_db_pkg: The name of the organism-specific annotation package.
#   - species: The scientific name of the species (e.g., "Saccharomyces_cerevisiae").
#   - kegg_organism: The 3-letter KEGG organism code (e.g., "sce" for yeast).
#   - alpha_pathway: The significance threshold for pathway enrichment.

# --- 1. Setup ---
# Source the helper functions
print(snakemake)
source(base::file.path(snakemake@scriptdir, "enrichment_utils.R"))

# Unpack snakemake object
log_file <- snakemake@log[[1]]
dge_path <- snakemake@input$dge_tsv
output_path <- snakemake@output$gsea_rds
enrichment_params <- snakemake@params[["enrichment"]]

# --- 2. Logging and Library Loading ---
setup_logging_and_libs(log_file)

# --- 3. Load Organism-Specific Database ---
org_db_pkg <- get_and_load_orgdb(enrichment_params)

# --- 4. Load Data and Prepare Gene List ---
res_tb <- load_and_map_dge_results(dge_path, org_db_pkg)

res_tb_filtered <- res_tb %>%
  dplyr::filter(!is.na(entrez_id) & !is.na(log2FoldChange)) %>%
  dplyr::distinct(entrez_id, .keep_all = TRUE)

genelist_fc_sort <- res_tb_filtered$log2FoldChange
names(genelist_fc_sort) <- res_tb_filtered$entrez_id
genelist_fc_sort <- base::sort(genelist_fc_sort, decreasing = TRUE)
base::message(base::paste(
  "Created ranked gene list for GSEA with", base::length(genelist_fc_sort),
  "genes."
))

if (base::length(genelist_fc_sort) == 0) {
  base::stop(
    "The gene list for GSEA is empty. Check gene ID mapping and DGE results."
  )
}

# --- 5. Perform GSEA ---
base::set.seed(123)
gsea_results <- base::list()

# Check once if the OrgDb has SYMBOL mappings to avoid repeated lookups
orgdb_obj <- base::get(org_db_pkg)
has_symbol_support <- "SYMBOL" %in% AnnotationDbi::columns(orgdb_obj)

# GSEA for Gene Ontology (GO)
base::message("Running GSEA for GO (BP)...")
gsego_defaults <- "geneList = genelist_fc_sort, OrgDb = get(org_db_pkg), keyType = 'ENTREZID'"
gsego_final_args <- base::paste(
  gsego_defaults, enrichment_params$clusterprofiler$gsea$gseGO$extra,
  sep = ", "
)
gsego_cmd <- base::paste0("clusterProfiler::gseGO(", gsego_final_args, ")")
base::message("Command: ", gsego_cmd)
go_res <- base::eval(base::parse(text = gsego_cmd))
gsea_results$GO <- safely_set_readable(
  go_res, orgdb_obj,
  has_symbol_support
)

# GSEA for KEGG Pathways
base::message("Running GSEA for KEGG...")
gsekegg_defaults <- base::paste0(
  "geneList = genelist_fc_sort, organism = '", enrichment_params$kegg_organism,
  "', keyType = 'ncbi-geneid'"
)
gsekegg_final_args <- base::paste(
  gsekegg_defaults, enrichment_params$clusterprofiler$gsea$gseKEGG$extra,
  sep = ", "
)
gsekegg_cmd <- base::paste0(
  "clusterProfiler::gseKEGG(", gsekegg_final_args, ")"
)
base::message("Command: ", gsekegg_cmd)
kegg_res <- base::eval(base::parse(text = gsekegg_cmd))
gsea_results$KEGG <- safely_set_readable(
  kegg_res, orgdb_obj,
  has_symbol_support
)

# GSEA for KEGG Modules (MKEGG)
if (enrichment_params$clusterprofiler$kegg_module$enabled) {
  base::message("Running GSEA for KEGG Modules (MKEGG)...")
  gsemkegg_defaults <- base::paste0(
    "geneList = genelist_fc_sort, organism = '",
    enrichment_params$kegg_organism, "', keyType = 'ncbi-geneid'"
  )
  gsemkegg_final_args <- base::paste(
    gsemkegg_defaults,
    enrichment_params$clusterprofiler$kegg_module$gseMKEGG$extra,
    sep = ", "
  )
  gsemkegg_cmd <- base::paste0(
    "clusterProfiler::gseMKEGG(", gsemkegg_final_args, ")"
  )
  base::message("Command: ", gsemkegg_cmd)
  mkegg_res <- base::eval(base::parse(text = gsemkegg_cmd))
  gsea_results$MKEGG <- safely_set_readable(
    mkegg_res, orgdb_obj,
    has_symbol_support
  )
}

# GSEA for WikiPathways
if (enrichment_params$clusterprofiler$wikipathways$enabled) {
  # Replaces underscores with spaces in species name for matching
  wp_species <- base::gsub("_", " ", enrichment_params$species)
  supported_wp_species <- clusterProfiler::get_wp_organisms()

  if (!wp_species %in% supported_wp_species) {
    base::message(
      "Skipping WikiPathways GSEA: species '", wp_species,
      "' not found in clusterProfiler::get_wp_organisms()."
    )
  } else {
    base::message("Running GSEA for WikiPathways...")
    gsewp_defaults <- base::paste0(
      "geneList = genelist_fc_sort, organism = '", wp_species, "'"
    )
    gsewp_final_args <- base::paste(
      gsewp_defaults,
      enrichment_params$clusterprofiler$wikipathways$gseWP$extra,
      sep = ", "
    )
    gsewp_cmd <- base::paste0(
      "clusterProfiler::gseWP(", gsewp_final_args, ")"
    )
    base::message("Command: ", gsewp_cmd)
    wp_res <- base::eval(base::parse(text = gsewp_cmd))
    gsea_results$WP <- safely_set_readable(
      wp_res, orgdb_obj,
      has_symbol_support
    )
  }
}

# GSEA for MSigDB
if (enrichment_params$msigdb$enabled) {
  base::message("Running GSEA for MSigDB...")

  # Validate species support
  species_supported <- validate_msigdb_species(enrichment_params$species)

  if (!species_supported) {
    base::message(
      "Skipping MSigDB GSEA: species '", enrichment_params$species,
      "' not supported by MSigDB."
    )
  } else {
    # Load MSigDB gene sets
    msigdb_genesets <- load_msigdb_genesets(
      enrichment_params$species,
      enrichment_params$msigdb$collections
    )

    # Run GSEA for each collection
    for (collection in enrichment_params$msigdb$collections) {
      base::message("Running GSEA for MSigDB collection: ", collection)

      # Filter for this collection
      collection_genesets <- msigdb_genesets[msigdb_genesets$gs_cat == collection, ]
      collection_data <- prepare_enrichment_data(collection_genesets, "msigdb")

      # Run GSEA
      gsea_defaults <- "geneList = genelist_fc_sort"
      gsea_final_args <- base::paste(
        gsea_defaults,
        "TERM2GENE = collection_data$TERM2GENE",
        "TERM2NAME = collection_data$TERM2NAME",
        enrichment_params$msigdb$gsea$extra,
        sep = ", "
      )
      gsea_cmd <- base::paste0(
        "clusterProfiler::GSEA(", gsea_final_args, ")"
      )
      base::message("Command: ", gsea_cmd)

      collection_res <- base::eval(base::parse(text = gsea_cmd))
      gsea_results[[paste0("MSigDB_", collection)]] <- safely_set_readable(
        collection_res, orgdb_obj,
        has_symbol_support
      )
    }

    # Load and process custom GMT files if specified
    if (base::length(enrichment_params$msigdb$custom_gmt_files) > 0) {
      base::message("Processing custom GMT files for GSEA...")

      for (gmt_file in enrichment_params$msigdb$custom_gmt_files) {
        if (base::file.exists(gmt_file)) {
          base::message("Running GSEA for custom GMT file: ", gmt_file)

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

          # Run GSEA
          gsea_defaults <- "geneList = genelist_fc_sort"
          gsea_final_args <- base::paste(
            gsea_defaults,
            "TERM2GENE = term2gene",
            "TERM2NAME = term2name",
            enrichment_params$msigdb$gsea$extra,
            sep = ", "
          )
          gsea_cmd <- base::paste0(
            "clusterProfiler::GSEA(", gsea_final_args, ")"
          )
          base::message("Command: ", gsea_cmd)

          gmt_res <- base::eval(base::parse(text = gsea_cmd))

          # Create a safe name for the result
          gmt_name <- base::gsub("^.*/", "", base::gsub("\\.gmt$", "", gmt_file))
          gsea_results[[paste0("CustomGMT_", gmt_name)]] <- safely_set_readable(
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
base::message("Saving GSEA results to: ", output_path)
base::saveRDS(gsea_results, file = output_path)

log_script_completion("GSEA script")
