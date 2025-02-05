# workflow/rules/salmon_quant.smk


rule salmon_quant:
    input:
        transcriptome=multiext(
            "results/salmon/transcriptome_index/",
            "complete_ref_lens.bin",
            "ctable.bin",
            "ctg_offsets.bin",
            "duplicate_clusters.tsv",
            "info.json",
            "mphf.bin",
            "pos.bin",
            "pre_indexing.log",
            "rank.bin",
            "refAccumLengths.bin",
            "ref_indexing.log",
            "reflengths.bin",
            "refseq.bin",
            "seq.bin",
            "versionInfo.json",
        ),
        reads=get_paired_reads,
    output:
        multiext(
            "results/salmon/{sample}_{unit}.salmon/",
            "cmd_info.json",
            "lib_format_counts.json",
            "quant.sf",
        ),
        directory("results/salmon/{sample}_{unit}.salmon/aux_info"),
        directory("results/salmon/{sample}_{unit}.salmon/libParams"),
        directory("results/salmon/{sample}_{unit}.salmon/logs"),
    params:
        transcriptome_dir="results/salmon/transcriptome_index",
        lib_type=config["params"]["salmon_quant"]["lib_type"],
        mapping_strategy=config["params"]["salmon_quant"]["mapping_strategy"],
        bias_correction=config["params"]["salmon_quant"]["bias_correction"],
        extra=config["params"]["salmon_quant"]["extra"],
    threads: 12
    conda:
        "../envs/salmon.yaml"
    log:
        "logs/salmon/salmon_quant_{sample}_{unit}.log",
    message:
        """
        Running Salmon Quant for:
            sample = {wildcards.sample},
            unit = {wildcards.unit}
        Running Salmon Quant with the inputs:
            {input.reads[0]}
            {input.reads[1]}
        """
    shell:
        """
        (salmon quant -i {params.transcriptome_dir} \
        -p {threads} \
        -l {params.lib_type} \
        -1 {input.reads[0]} \
        -2 {input.reads[1]} \
        -o results/salmon/{wildcards.sample}_{wildcards.unit}.salmon \
        {params.mapping_strategy} \
        {params.bias_correction} \
        {params.extra}) &> {log}
        """
