rule create_gtf:
    input:
        genome_fna=("results/datasets/ncbi_dataset/data/{genome}/genomic.gff").format(
            genome=config["ncbi_genome_accession"]
        ),
    output:
        "results/datasets/ncbi_dataset/data/{genome}/genomic.gtf".format(
            genome=config["ncbi_genome_accession"]
        ),
    conda:
        "../envs/gffread.yaml"
    log:
        "logs/gffread/gffread.log",
    shell:
        """
        gffread {input} -T -o {output}
        """
