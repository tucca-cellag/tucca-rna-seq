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
  pkg_dirs <- base::list.dirs(path = install_source, full.names = FALSE, recursive = FALSE)
  pkg_name <- pkg_dirs[grepl("^org\\..+\\.db$", pkg_dirs)]
  if (length(pkg_name) == 0) base::stop("Could not find locally built OrgDb package directory.")
  pkg_name[1]
} else {
  enrichment_params$org_db_pkg
}
base::message("Target OrgDb package: ", org_db_pkg)

# Load the package (it should be pre-installed by the 'install_orgdb' rule)
base::library(org_db_pkg, character.only = TRUE)

# --- 3. Load Data and Prepare Gene Lists ---
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
  base::message("No significant genes found. ORA will not be performed. Creating empty output.")
  base::saveRDS(base::list(), file = snakemake@output$ora_rds)
  base::quit(save = "no", status = 0)
}

# --- 4. Perform ORA ---
base::set.seed(123)
ora_results <- base::list()

# ORA for Gene Ontology (GO)
base::message("Running ORA for GO (BP)...")
enrichgo_defaults <- "gene = significant_genes, universe = universe_genes, OrgDb = get(org_db_pkg), ont = 'BP', keyType = 'ENTREZID'"
enrichgo_final_args <- base::paste(
  enrichgo_defaults, enrichment_params$enrichgo_extra,
  sep = ", "
)
enrichgo_cmd <- base::paste0(
  "clusterProfiler::enrichGO(", enrichgo_final_args, ")"
)
base::message("Command: ", enrichgo_cmd)
ora_results$GO <- base::eval(base::parse(text = enrichgo_cmd))

# ORA for KEGG Pathways
base::message("Running ORA for KEGG...")
enrichkegg_defaults <- base::paste0(
  "gene = significant_genes, universe = universe_genes, organism = '",
  enrichment_params$kegg_organism, "'"
)
enrichkegg_final_args <- base::paste(
  enrichkegg_defaults, enrichment_params$enrichkegg_extra,
  sep = ", "
)
enrichkegg_cmd <- base::paste0(
  "clusterProfiler::enrichKEGG(", enrichkegg_final_args, ")"
)
base::message("Command: ", enrichkegg_cmd)
ora_results$KEGG <- base::eval(base::parse(text = enrichkegg_cmd))

# --- 5. Save Results ---
base::message("Saving ORA results to: ", snakemake@output$ora_rds)
base::saveRDS(ora_results, file = snakemake@output$ora_rds)

base::message("ORA script finished successfully.")
base::date()
