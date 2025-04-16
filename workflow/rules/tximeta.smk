rule tximeta:
    input:
        salmonDir=expand(
            os.path.join(
                "results/salmon", "{sample}_{unit}", "{sample}_{unit}_quant.sf"
            ),
            sample=samples.names.values.tolist(),
            unit=units.names.values.tolist(),
        ),
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
