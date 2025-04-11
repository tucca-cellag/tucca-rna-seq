# TODO: requires OrgDb for the species
rule tximport:
    input:
        quant="results/salmon/{sample}_{unit}/quant.sf",
        lib="results/salmon/{sample}_{unit}/lib_format_counts.json",
        aux_info="results/salmon/{sample}_{unit}/aux_info",
        cmd_info="results/salmon/{sample}_{unit}/cmd_info.json",
        libparams="results/salmon/{sample}_{unit}/libParams",
        logs="results/salmon/{sample}_{unit}/logs",
        tx_to_gene="resources/tx2gene.tsv",
    output:
        txi="results/tximport/SummarizedExperimentObject.RDS",
    params:
        extra="type='salmon', countsFromAbundance='lengthScaledTPM', importer=read.delim, ignoreTxVersion=TRUE",
    log:
        "logs/tximport/tximport.log",
    wrapper:
        "v5.10.0/bio/tximport"
