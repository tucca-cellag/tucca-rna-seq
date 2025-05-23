# workflow/rules/qualimap.smk


rule qualimap_rnaseq:
    input:
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
        counting_alg=config["params"]["qualimap_rnaseq"]["counting_alg"],
        sequencing_protocol=config["params"]["qualimap_rnaseq"]["sequencing_protocol"],
        extra=config["params"]["qualimap_rnaseq"]["extra"],
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
