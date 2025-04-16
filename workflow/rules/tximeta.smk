rule tximeta:
    input:
        salmonDir="results/salmon/*/*_quant.sf",
        linkedTxome="results/salmon/transcriptome_index.json",
    output:
        se="se.RDS",
    params:
        extra=config["params"]["tximeta"]["extra"],
    conda:
        "../envs/tximeta.yaml"
    log:
        "logs/tximeta/tximeta.log",
    script:
        "../scripts/tximeta.R"
