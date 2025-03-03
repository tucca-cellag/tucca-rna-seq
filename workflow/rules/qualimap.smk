# workflow/rules/qualimap.smk


rule qualimap_rnaseq:
    input:
        bam="results/star/{sample}_{unit}_Aligned.sortedByCoord.out.bam",
        genome_gtf="results/datasets/ncbi_dataset/data/{genome}/genomic.gtf".format(
            genome=config["ref"]["ncbi_genome_accession"]
        ),
    output:
        multiext(
            "results/qualimap/{sample}_{unit}.qualimap/",
            "qualimapReport.html",
            "rnaseq_qc_results.txt",
        ),
        directory("results/qualimap/{sample}_{unit}.qualimap/css"),
        directory("results/qualimap/{sample}_{unit}.qualimap/images_qualimapReport"),
        directory("results/qualimap/{sample}_{unit}.qualimap/raw_data_qualimapReport"),
    params:
        counting_alg=config["params"]["qualimap_rnaseq"]["counting_alg"],
        sequencing_protocol=config["params"]["qualimap_rnaseq"]["sequencing_protocol"],
        extra=config["params"]["qualimap_rnaseq"]["extra"],
    container:
        config["containers"]["qualimap"]
    log:
        "logs/qualimap/qualimap_rnaseq_{sample}_{unit}.log",
    message:
        "Running Qualimap RNA-Seq on Bam for {wildcards.sample} {wildcards.unit}"
    shell:
        """
        unset DISPLAY

        (qualimap rnaseq \
        -outdir results/qualimap/{wildcards.sample}_{wildcards.unit}.qualimap \
        -a {params.counting_alg} \
        -bam {input.bam} \
        -gtf {input.genome_gtf} \
        --sequencing-protocol {params.sequencing_protocol} \
        {params.extra}) &> {log}
        """
