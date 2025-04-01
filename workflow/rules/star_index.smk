# workflow/rules/star_index.smk

import glob


rule star_index:
    input:
        genome_fna="results/datasets/ncbi_dataset/data/{genome_asc}/{genome_asc}_{genome_name}_genomic.fna".format(
            genome_asc=config["genome"]["assembly_accession"],
            genome_name=config["genome"]["assembly_name"],
        ),
        genome_gtf="results/datasets/ncbi_dataset/data/{genome_asc}/genomic.gtf".format(
            genome_asc=config["genome"]["assembly_accession"],
        ),
    output:
        multiext(
            "results/star/{genome_asc}_index/".format(
                genome_asc=config["genome"]["assembly_accession"],
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
        genome_asc=config["genome"]["assembly_accession"],
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
        --genomeDir results/star/{params.genome_asc}_index \
        --genomeFastaFiles {input.genome_fna} \
        --sjdbGTFfile {input.genome_gtf} \
        --sjdbOverhang {params.sjdb_overhang} \
        {params.extra}) &> {log}
        """
