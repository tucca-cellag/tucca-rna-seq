rule run_enrichment:
    input:
        dge_tsv="resources/deseq2/{analysis}/{contrast}/dge.tsv",
    output:
        # RDS objects for clusterProfiler visualization
        gsea_go_rds="resources/enrichment/{analysis}/{contrast}/gsea_go.rds",
        gsea_kegg_rds="resources/enrichment/{analysis}/{contrast}/gsea_kegg.rds",
        gsea_reactome_rds="resources/enrichment/{analysis}/{contrast}/gsea_reactome.rds",
        gsea_wp_rds="resources/enrichment/{analysis}/{contrast}/gsea_wp.rds",
        gsea_msigdb_h_rds="resources/enrichment/{analysis}/{contrast}/gsea_msigdb_h.rds",
        spia_rds="resources/enrichment/{analysis}/{contrast}/spia.rds",
        # CSV summaries for easy viewing
        gsea_go_csv="resources/enrichment/{analysis}/{contrast}/gsea_go_summary.csv",
        gsea_kegg_csv="resources/enrichment/{analysis}/{contrast}/gsea_kegg_summary.csv",
        gsea_reactome_csv="resources/enrichment/{analysis}/{contrast}/gsea_reactome_summary.csv",
        gsea_wp_csv="resources/enrichment/{analysis}/{contrast}/gsea_wp_summary.csv",
        gsea_msigdb_h_csv="resources/enrichment/{analysis}/{contrast}/gsea_msigdb_h_summary.csv",
        spia_csv="resources/enrichment/{analysis}/{contrast}/spia_summary.csv",
        target_pathways_csv="resources/enrichment/{analysis}/{contrast}/target_pathways.csv",
    params:
        alpha_pathway=config["enrichment"]["alpha_pathway"],
        species=config["enrichment"]["species"],
        org_db_pkg=config["enrichment"]["org_db_pkg"],
        targets=config["enrichment"]["targets"],
    log:
        "logs/enrichment/{analysis}/{contrast}.log",
    conda:
        ENRICHMENT_ENV_FILE
    script:
        "workflow/scripts/run_enrichment.R"


rule all_enrichment:
    input:
        expand(
            [
                "resources/enrichment/{contrast_path}/gsea_go.rds",
                "resources/enrichment/{contrast_path}/gsea_kegg.rds",
                "resources/enrichment/{contrast_path}/spia.rds",
                "resources/enrichment/{contrast_path}/target_pathways.csv",
            ],
            contrast_path=ALL_CONTRASTS,
        ),
    output:
        touch("resources/enrichment/enrichment_analyses_complete.done"),
