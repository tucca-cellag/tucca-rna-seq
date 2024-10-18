rule salmon_quant:
    input:
        transcriptome="results/salmon/transcriptome_index",
        reads=get_paired_reads,
    output:
        multiext(
            "results/salmon/{sample}_{unit}.salmon",
            "cmd_info.json",
            "lib_format_counts.json",
            "quant.sf",
        ),
        directory("results/salmon/{sample}_{unit}.salmon/aux_info"),
        directory("results/salmon/{sample}_{unit}.salmon/libParams"),
        directory("results/salmon/{sample}_{unit}.salmon/logs"),
    params:
        lib_type=config["params"]["salmon_quant"]["lib_type"],
        mapping_strategy=config["params"]["salmon_quant"]["mapping_strategy"],
        bias_correction=config["params"]["salmon_quant"]["bias_correction"],
        extra=config["params"]["salmon_quant"]["extra"],
    threads: 12
    conda:
        "../envs/salmon.yaml"
    log:
        "logs/salmon/salmon_quant_{sample}_{unit}.log",
    shell:
        """
        (salmon quant -i {input.transcriptome} \
        -p {threads} \
        -l {params.lib_type} \
        -1 {input.reads[0]} \
        -2 {input.reads[1]} \
        -o results/salmon/{sample}_{unit}.salmon \
        {input.mapping_strategy} \
        {input.bias_correction} \
        {input.extra}) &> {log}
        """
