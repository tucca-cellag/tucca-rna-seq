# -----------------------------------------------------------------------------
# Snakemake-script for Pathway Enrichment Analysis
#
# This script performs GSEA and SPIA for a single contrast.
#
# Inputs:
#   - dge_tsv: Path to the differential gene expression results table (TSV).
# Outputs:
#   - RDS objects for each enrichment analysis (gseaGO, gseaKegg, etc.)
#   - CSV summaries of significant results.
#   - RDS object for SPIA results.
#   - CSV of pathways containing specific target genes.
# Params:
#   - alpha_pathway, species, org_db_pkg, targets
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 1. Setup, Logging, and Library Loading
# -----------------------------------------------------------------------------
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")
date()

library(devtools)
session_info()

library(clusterProfiler)
library(SPIA)
library(ReactomePA)
library(msigdbr)
library(tidyverse)
# Load the organism-specific database package dynamically
library(snakemake@params$org_db_pkg, character.only = TRUE)

# -----------------------------------------------------------------------------
# 2. Load Data and Prepare Gene Lists
# -----------------------------------------------------------------------------
message("Loading DGE results from: ", snakemake@input$dge_tsv)
res_tb <- read_tsv(snakemake@input$dge_tsv)

# --- Create Gene Lists ---
# Ensure 'entrez_id' and 'log2FoldChange' columns exist and handle NAs
if (!"entrez_id" %in% colnames(res_tb) || !"log2FoldChange" %in% colnames(res_tb)) {
  stop("Input DGE table must contain 'entrez_id' and 'log2FoldChange' columns.")
}
res_tb_filtered <- res_tb %>%
  filter(!is.na(entrez_id) & !is.na(log2FoldChange)) %>%
  distinct(entrez_id, .keep_all = TRUE) # Remove duplicate Entrez IDs

# Ranked gene list for GSEA
genelist_fc_sort <- res_tb_filtered$log2FoldChange
names(genelist_fc_sort) <- res_tb_filtered$entrez_id
genelist_fc_sort <- sort(genelist_fc_sort, decreasing = TRUE)
message(paste("Created ranked gene list for GSEA with", length(genelist_fc_sort), "genes."))

# Gene list for SPIA (significant genes)
res_sig <- res_tb_filtered %>% filter(padj < snakemake@params$alpha_pathway)
genelist_sig_fc <- res_sig$log2FoldChange
names(genelist_sig_fc) <- res_sig$entrez_id
message(paste("Created significant gene list for SPIA with", length(genelist_sig_fc), "genes."))

# Universe for SPIA
universe_by_entrez <- names(genelist_fc_sort)
message(paste("Created universe for SPIA with", length(universe_by_entrez), "genes."))

# -----------------------------------------------------------------------------
# 3. Perform GSEA (GO, KEGG, Reactome, WikiPathways, MSigDB)
# -----------------------------------------------------------------------------

# --- GSEA against GO (Gene Ontology) ---
message("Running GSEA against GO...")
set.seed(12345)
gseaGO <- gseGO(
  geneList = genelist_fc_sort,
  OrgDb = snakemake@params$org_db_pkg,
  ont = "BP",
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = snakemake@params$alpha_pathway,
  verbose = FALSE
)
saveRDS(gseaGO, file = snakemake@output$gsea_go_rds)
if (nrow(gseaGO) > 0) {
  write_csv(as.data.frame(gseaGO@result), file = snakemake@output$gsea_go_csv)
}

# --- GSEA against KEGG ---
message("Running GSEA against KEGG...")
set.seed(12345)
gseaKegg <- gseKEGG(
  geneList = genelist_fc_sort,
  organism = "bta",
  keyType = "ncbi-geneid",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = snakemake@params$alpha_pathway,
  verbose = FALSE
)
saveRDS(gseaKegg, file = snakemake@output$gsea_kegg_rds)
if (nrow(gseaKegg) > 0) {
  write_csv(as.data.frame(gseaKegg@result), file = snakemake@output$gsea_kegg_csv)
}

# --- GSEA against Reactome ---
message("Running GSEA against Reactome...")
set.seed(12345)
gseaReact <- gsePathway(
  genelist_fc_sort,
  organism = "bovine",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = snakemake@params$alpha_pathway,
  verbose = FALSE
)
saveRDS(gseaReact, file = snakemake@output$gsea_reactome_rds)
if (nrow(gseaReact) > 0) {
  write_csv(as.data.frame(gseaReact@result), file = snakemake@output$gsea_reactome_csv)
}

# --- GSEA against WikiPathways ---
message("Running GSEA against WikiPathways...")
set.seed(12345)
gseaWP <- gseWP(
  genelist_fc_sort,
  organism = snakemake@params$species,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = snakemake@params$alpha_pathway,
  verbose = FALSE
)
saveRDS(gseaWP, file = snakemake@output$gsea_wp_rds)
if (nrow(gseaWP) > 0) {
  write_csv(as.data.frame(gseaWP@result), file = snakemake@output$gsea_wp_csv)
}


# --- GSEA against MSigDB Hallmark ---
message("Running GSEA against MSigDB Hallmark gene sets...")
H_t2g <- msigdbr(species = snakemake@params$species, category = "H") %>%
  dplyr::select(gs_name, entrez_gene)
set.seed(12345)
gseaMSigH <- GSEA(
  genelist_fc_sort,
  TERM2GENE = H_t2g,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = snakemake@params$alpha_pathway,
  verbose = FALSE
)
saveRDS(gseaMSigH, file = snakemake@output$gsea_msigdb_h_rds)
if (nrow(gseaMSigH) > 0) {
  write_csv(as.data.frame(gseaMSigH@result), file = snakemake@output$gsea_msigdb_h_csv)
}


# -----------------------------------------------------------------------------
# 4. Perform Signaling Pathway Impact Analysis (SPIA)
# -----------------------------------------------------------------------------
message("Running SPIA...")
# Generate SPIA data if it doesn't exist. Assumes a writable 'workflow' directory.
spia_data_dir <- "workflow/resources/SPIAdata"
if (!dir.exists(spia_data_dir)) dir.create(spia_data_dir, recursive = TRUE)
if (!file.exists(file.path(spia_data_dir, "btaSPIA.RData"))) {
  makeSPIAdata(kgml.path = file.path(spia_data_dir, "kgml"), organism = "bta", out.path = spia_data_dir)
}

spia_result <- spia(
  de = genelist_sig_fc,
  all = universe_by_entrez,
  organism = "bta",
  data.dir = spia_data_dir,
  verbose = FALSE
)
saveRDS(spia_result, file = snakemake@output$spia_rds)
if (!is.null(spia_result) && nrow(spia_result) > 0) {
  spia_result_sig <- spia_result %>% filter(pGFdr < snakemake@params$alpha_pathway)
  write_csv(spia_result_sig, file = snakemake@output$spia_csv)
}


# -----------------------------------------------------------------------------
# 5. Find Pathways with Target Genes
# -----------------------------------------------------------------------------
message("Searching for target genes in significant pathways...")
targets_vect <- snakemake@params$targets
search_pattern <- paste(targets_vect, collapse = "|")

# Function to safely bind rows even if some results are NULL or empty
bind_enrichment_results <- function(...) {
  results_list <- list(...)
  # Filter out NULL or 0-row dataframes
  valid_results <- Filter(function(x) !is.null(x) && nrow(x) > 0, results_list)
  if (length(valid_results) == 0) {
    return(tibble())
  }
  # Select and rename columns for consistency, then bind
  processed_list <- lapply(valid_results, function(df) {
    df %>%
      as_tibble() %>%
      select(ID, Description, p.adjust, core_enrichment)
  })
  bind_rows(processed_list)
}

all_gsea_results <- bind_enrichment_results(
  gseaGO@result, gseaKegg@result, gseaReact@result, gseaWP@result, gseaMSigH@result
)

if (nrow(all_gsea_results) > 0) {
  target_pathways <- all_gsea_results %>%
    filter(str_detect(core_enrichment, regex(search_pattern, ignore_case = TRUE)) |
      str_detect(Description, regex(search_pattern, ignore_case = TRUE))) %>%
    mutate(targets_found = map_chr(str_extract_all(core_enrichment, regex(search_pattern, ignore_case = TRUE)), ~ paste(unique(.x), collapse = ", ")))

  write_csv(target_pathways, file = snakemake@output$target_pathways_csv)
}

message("Enrichment script finished successfully.")
date()
