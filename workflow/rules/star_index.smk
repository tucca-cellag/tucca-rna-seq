import glob


rule star_index:
    input:
        genome_fna=glob.glob(
            (
                "results/datasets/ncbi_dataset/data/{genome}/{genome}_"
                + "*"
                + "_genomic.fna"
            ).format(genome=config["ref"]["ncbi_genome_accession"])
        ),
        genome_gtf="results/datasets/ncbi_dataset/data/{genome}/genomic.gtf".format(
            genome=config["ref"]["ncbi_genome_accession"]
        ),
    output:
        directory(
            "results/star/{genome}_index/".format(
                genome=config["ref"]["ncbi_genome_accession"]
            )
        ),
    threads: 12
    conda:
        "../envs/star.yaml"
    log:
        "logs/star/star_index.log",
    shell:
        """
        (STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {output} \
        --genomeFastaFiles {input.genome_fna} \
        --sjdbGTFfile {input.genome_gtf} \
        --sjdbOverhang {config["params"]["star_index"]["sjdbOverhang"]}) &> {log}
        """
