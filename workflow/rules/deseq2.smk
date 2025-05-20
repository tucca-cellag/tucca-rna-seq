# Rule 1: Create DESeqDataSet for each analysis configuration
# DDS is parameterized by the analysis name from config.yaml
rule DESeqDataSet_from_ranged_se_per_analysis:
    input:
        se="resources/tximeta/tximeta_gse.RDS",
    output:
        dds="resources/deseq2/{analysis_name}/dds.RDS",
    log:
        "logs/deseq2/{analysis_name}/DESeqDataSet_se.log",
    threads: get_dds_threads
    params:
        formula=get_dds_formula,
        min_counts=get_dds_min_counts,
        extra=get_dds_extra,
    wrapper:
        "v6.2.0/bio/deseq2/deseqdataset"


# Rule 2: Run DESeq2 Wald test for each specified contrast within each analysis
# DDS comes from the DESeqDataSet_from_ranged_se_per_analysis rule, specific to the {analysis_name}
rule deseq2_wald_per_analysis:
    input:
        dds="resources/deseq2/{analysis_name}/dds.RDS",
    output:
        # dds after DESeq()
        wald_rds="resources/deseq2/{analysis_name}/{contrast_name}/wald.RDS",
        # res for a given contrast [generated from dds after DESeq()]
        wald_tsv="resources/deseq2/{analysis_name}/{contrast_name}/dge.tsv",
        # counts(dds) as a TSV
        normalized_counts_table="resources/deseq2/{analysis_name}/{contrast_name}/counts.tsv",
        # counts(dds) as an RDS
        normalized_counts_rds="resources/deseq2/{analysis_name}/{contrast_name}/counts.RDS",
    log:
        "logs/deseq2/{analysis_name}/{contrast_name}/wald.log",
    threads: get_wald_threads
    params:
        # get user extra parameters from config.yaml
        deseq_extra=get_wald_deseq_extra,
        shrink_extra=get_wald_shrink_extra,
        results_extra=get_wald_results_extra,
        # get contrast to be evaluated
        contrast=get_wald_contrast_elements,
    wrapper:
        "v6.2.0/bio/deseq2/wald"


rule get_results_from_all_deseq_analyses:
    input:
        # Request all DDS files (one per analysis_name)
        expand(
            "resources/deseq2/{analysis_name}/dds.RDS",
            analysis_name=DESEQ_ANALYSES_NAMES,
        ),
        # Request all Wald test outputs (for each contrast in each analysis)
        # Use zip to correctly pair analysis_name and contrast_name
        expand(
            [
                "resources/deseq2/{analysis_name}/{contrast_name}/wald.RDS",
                "resources/deseq2/{analysis_name}/{contrast_name}/dge.tsv",
                "resources/deseq2/{analysis_name}/{contrast_name}/counts.tsv",
                "resources/deseq2/{analysis_name}/{contrast_name}/counts.RDS",
            ],
            zip,
            analysis_name=[job["analysis_name"] for job in CONTRAST_JOBS],
            contrast_name=[job["contrast_name"] for job in CONTRAST_JOBS],
        ),
    output:
        touch("resources/deseq2/deseq2_analyses_complete.done"),
    log:
        "logs/deseq2/get_results_from_all_deseq_analyses.log",
    shell:
        "touch {output}"
