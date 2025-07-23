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

# --- 1. Setup, Logging, and Library Loading ---
log <- base::file(snakemake@log[[1]], open = "wt")
base::sink(log)
base::sink(log, type = "message")
base::date()

base::library(package = "devtools", character.only = TRUE)
devtools::session_info()

base::library(package = "clusterProfiler", character.only = TRUE)
base::library(package = "AnnotationDbi", character.only = TRUE)
base::library(package = "magrittr", character.only = TRUE) # For the %>% pipe
base::library(package = "dplyr", character.only = TRUE)
base::library(package = "readr", character.only = TRUE)

# --- 2. Load Organism-Specific Database ---
enrichment_params <- snakemake@params[["enrichment"]]
install_method <- enrichment_params$install_method
install_source <- enrichment_params$install_source

# Determine the actual package name
org_db_pkg <- if (install_method == "local") {
  # When built locally, the package name is derived inside AnnotationForge.
  # We find the directory name it created to get the actual package name.
  pkg_dirs <- base::list.dirs(
    path = install_source, full.names = FALSE, recursive = FALSE
  )
  pkg_name <- pkg_dirs[grepl("^org\\..+\\.db$", pkg_dirs)]
  if (length(pkg_name) == 0) {
    base::stop("Could not find locally built OrgDb package directory.")
  }
  pkg_name[1]
} else {
  enrichment_params[["org_db_pkg"]]
}
base::message("Target OrgDb package: ", org_db_pkg)

# Load the package (it should be pre-installed by the 'install_orgdb' rule)
base::library(org_db_pkg, character.only = TRUE)

# --- 3. Load Data and Prepare Gene List ---
base::message("Loading DGE results from: ", snakemake@input$dge_tsv)
res_tb <- readr::read_tsv(snakemake@input$dge_tsv, show_col_types = FALSE) %>%
  dplyr::rename(feature_id = 1)

base::message("Mapping feature IDs to Entrez IDs...")
res_tb$entrez_id <- AnnotationDbi::mapIds(
  base::get(org_db_pkg),
  keys = res_tb$feature_id,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

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

# --- 4. Perform GSEA ---
base::set.seed(123)
gsea_results <- base::list()

# GSEA for Gene Ontology (GO)
base::message("Running GSEA for GO (BP)...")
gsego_defaults <- "geneList = genelist_fc_sort, OrgDb = get(org_db_pkg), keyType = 'ENTREZID'"
gsego_final_args <- base::paste(
  gsego_defaults, enrichment_params$clusterprofiler$gsea$gseGO$extra,
  sep = ", "
)
gsego_cmd <- base::paste0("clusterProfiler::gseGO(", gsego_final_args, ")")
base::message("Command: ", gsego_cmd)
gsea_results$GO <- base::eval(base::parse(text = gsego_cmd))


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
gsea_results$KEGG <- base::eval(base::parse(text = gsekegg_cmd))

# GSEA for KEGG Modules (MKEGG)
if (enrichment_params$clusterprofiler$kegg_module$enabled) {
  base::message("Running GSEA for KEGG Modules (MKEGG)...")
  gsemkegg_defaults <- base::paste0(
    "geneList = genelist_fc_sort, organism = '",
    enrichment_params$kegg_organism, "', keyType = 'ncbi-geneid'"
  )
  gsemkegg_final_args <- base::paste(
    gsemkegg_defaults, enrichment_params$clusterprofiler$kegg_module$gseMKEGG$extra,
    sep = ", "
  )
  gsemkegg_cmd <- base::paste0(
    "clusterProfiler::gseMKEGG(", gsemkegg_final_args, ")"
  )
  base::message("Command: ", gsemkegg_cmd)
  gsea_results$MKEGG <- base::eval(base::parse(text = gsemkegg_cmd))
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
      "geneList = genelist_fc_sort, organism = '", wp_species,
      "', keyType = 'ncbi-geneid'"
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
    gsea_results$WP <- base::eval(base::parse(text = gsewp_cmd))
  }
}

# --- 5. Save Results ---
base::message("Saving GSEA results to: ", snakemake@output$gsea_rds)
base::saveRDS(gsea_results, file = snakemake@output$gsea_rds)

base::message("GSEA script finished successfully.")
base::date()
