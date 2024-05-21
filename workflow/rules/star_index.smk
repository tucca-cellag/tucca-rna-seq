import glob


rule star_index:
    input:
        genome_fna=glob.glob(
            (
                "results/datasets/ncbi_dataset/data/{genome}/{genome}_"
                + "*"
                + "_genomic.fna"
            ).format(genome=config["ncbi_genome_accession"])
        ),
        genome_gtf="results/datasets/ncbi_dataset/data/{genome}/genomic.gtf".format(
            genome=config["ncbi_genome_accession"]
        ),
    output:
        directory("results/star/{genome}_index").format(
            genome=config["ncbi_genome_accession"]
        ),
    threads: 12
    conda:
        "../envs/star.yaml"
    log:
        "logs/star_index_{genome}.log".format(genome=config["ncbi_genome_accession"]),
    params:
        sjdbOverhang=149,
    resources:
        mem_mb=64000,
    shell:
        """
        (STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {output} \
        --genomeFastaFiles {input.genome_fna} \
        --sjdbGTFfile {input.genome_gtf} \
        --sjdbOverhang {sjdbOverhang}) &> {log}
        """
