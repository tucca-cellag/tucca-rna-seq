# workflow/rules/qualimap.smk


rule create_qualimap_flags:
    output:
        branch(
            is_qualimap_enabled(),
            then="resources/qualimap/qualimap_enabled.flag",
            otherwise="resources/qualimap/qualimap_disabled.flag",
        ),
    log:
        "logs/qualimap/create_flags.log",
    shell:
        """
        touch {output}
        """


rule qualimap_rnaseq:
    input:
        "resources/qualimap/qualimap_enabled.flag",
        bam="resources/star/{sample_unit}/Aligned.sortedByCoord.out.bam",
        genome_gtf=branch(
            config["ref_assembly"]["source"] == "RefSeq",
            then="resources/datasets/ncbi_dataset/data/{genome_asc}/genomic.gtf".format(
                genome_asc=config["ref_assembly"]["accession"],
            ),
            otherwise="resources/ensembl/{species}.{genome_name}.gtf".format(
                species=config["ref_assembly"]["species"],
                genome_name=config["ref_assembly"]["name"],
            ),
        ),
    output:
        multiext(
            "results/qualimap/{sample_unit}.qualimap/",
            "qualimapReport.html",
            "rnaseq_qc_results.txt",
        ),
        directory("results/qualimap/{sample_unit}.qualimap/css"),
        directory("results/qualimap/{sample_unit}.qualimap/images_qualimapReport"),
        directory("results/qualimap/{sample_unit}.qualimap/raw_data_qualimapReport"),
    params:
        counting_alg=lambda wildcards: config["params"]["qualimap_rnaseq"][
            "counting_alg"
        ],
        sequencing_protocol=lambda wildcards: config["params"]["qualimap_rnaseq"][
            "sequencing_protocol"
        ],
        extra=lambda wildcards: config["params"]["qualimap_rnaseq"]["extra"],
    conda:
        "../envs/qualimap.yaml"
    log:
        "logs/qualimap/qualimap_rnaseq_{sample_unit}.log",
    message:
        "Running Qualimap RNA-Seq on Bam for {wildcards.sample_unit}"
    shell:
        """
        unset DISPLAY

        (qualimap rnaseq \
        -outdir results/qualimap/{wildcards.sample_unit}.qualimap \
        -a {params.counting_alg} \
        -bam {input.bam} \
        -gtf {input.genome_gtf} \
        --sequencing-protocol {params.sequencing_protocol} \
        {params.extra}) &> {log}
        """
