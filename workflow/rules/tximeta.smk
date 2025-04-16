rule tximeta:
    input:
        expand(
            os.path.join(
                "results/salmon", "{sample}_{unit}", "{sample}_{unit}_quant.sf"
            ),
            sample=samples.sample_name.values.tolist(),
            unit=units.unit_name.values.tolist(),
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
