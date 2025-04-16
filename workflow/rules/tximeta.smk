rule tximeta:
    input:
        quants_paths=lambda wildcards: expand(
            "results/salmon/{sample}_{unit}/{sample}_{unit}_quant.sf"
        ),
        indexDir="",
        source="",
        organism="",
        release="",
        genome="",
        fasta="",
        gtf="",
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
