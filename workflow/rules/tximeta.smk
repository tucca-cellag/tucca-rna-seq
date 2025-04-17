rule tximeta:
    input:
        files=expand(
            "results/salmon/{sample_unit}/quant.sf",
            zip,
            sample_unit=units.sample_unit.values.tolist(),
        ),
        linkedTxome="results/salmon/transcriptome_index.json",
    output:
        se="resources/tximeta/se.RDS",
    params:
        extra=config["params"]["tximeta"]["extra"],
    conda:
        "../envs/tximeta.yaml"
    log:
        "logs/tximeta/tximeta.log",
    script:
        "../scripts/tximeta.R"
