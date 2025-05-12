rule DESeqDataSet_from_ranged_se:
    input:
        se="resources/tximeta/tximeta_gse.RDS",  # Output from your tximeta.smk rule
    output:
        dds="resources/deseq2/dds.RDS",
    threads: config["diffexp"]["deseq2"]["deseqdataset"]["threads"]
    log:
        "logs/deseq2/DESeqDataSet_se.log",
    params:
        # Required: design formula from config
        formula=config["diffexp"]["deseq2"]["deseqdataset"]["formula"],
        # Optional params for deseqdataset wrapper
        factor=config["diffexp"]["deseq2"]["deseqdataset"]["factor"],
        reference_level=config["diffexp"]["deseq2"]["deseqdataset"]["reference_level"],
        tested_level=config["diffexp"]["deseq2"]["deseqdataset"]["tested_level"],
        min_counts=config["diffexp"]["deseq2"]["deseqdataset"]["min_counts"],
        extra=config["diffexp"]["deseq2"]["deseqdataset"]["extra"],
    wrapper:
        "v5.10.0/bio/deseq2/deseqdataset"


# Rule 2: Run DESeq2 Wald test for each specified contrast using the snakemake-wrapper
rule deseq2_wald_per_contrast:
    input:
        dds="resources/deseq2/dds.RDS",  # Output from deseqdataset_from_se rule
    output:
        wald_rds="resources/deseq2/wald.RDS",
        wald_tsv="resources/deseq2/dge.tsv",
        deseq2_result_dir=directory("results/deseq2")
        normalized_counts_table="resources/deseq2/counts.tsv",
        normalized_counts_rds="resources/deseq2/counts.RDS",
    params:
        deseq_extra=config["diffexp"]["deseq2"]["wald"]["deseq_extra"],
        shrink_extra=config["diffexp"]["deseq2"]["wald"]["shrink_extra"],
        results_extra=config["diffexp"]["deseq2"]["wald"]["results_extra"],
        contrast=config["diffexp"]["deseq2"]["wald"]["contrast"],
    threads: config["diffexp"]["deseq2"]["deseqdataset"]["threads"]
    log:
        "logs/deseq2/DESeq2_wald.log",
    wrapper:
        # Using the commit from your viewed files for wald wrapper
        # Recommended: Change to a release tag like "vX.Y.Z/bio/deseq2/wald"
        "v5.10.0/bio/deseq2/wald"
