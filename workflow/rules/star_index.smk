# workflow/rules/star_index.smk

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
        multiext(
            "results/star/{genome}_index/".format(
                genome=config["ref"]["ncbi_genome_accession"]
            ),
            "chrLength.txt",
            "chrName.txt",
            "chrNameLength.txt",
            "chrStart.txt",
            "exonGeTrInfo.tab",
            "exonInfo.tab",
            "geneInfo.tab",
            "Genome",
            "genomeParameters.txt",
            "Log.out",
            "SA",
            "SAindex",
            "sjdbInfo.txt",
            "sjdbList.fromGTF.out.tab",
            "sjdbList.out.tab",
            "transcriptInfo.tab",
        ),
    params:
        genome=config["ref"]["ncbi_genome_accession"],
        sjdb_overhang=config["params"]["star_index"]["sjdbOverhang"],
        extra=config["params"]["star_index"]["extra"],
    threads: 12
    conda:
        "../envs/star.yaml"
    log:
        "logs/star/star_index.log",
    shell:
        """
        (STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir results/star/{params.genome}_index \
        --genomeFastaFiles {input.genome_fna} \
        --sjdbGTFfile {input.genome_gtf} \
        --sjdbOverhang {params.sjdb_overhang} \
        {params.extra}) &> {log}
        """
