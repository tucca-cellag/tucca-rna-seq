# workflow/rules/datasets_download_genome.smk

import glob


rule datasets_download_genome:
    output:
        "ncbi_dataset_{genome}.zip".format(
            genome=config["ref"]["ncbi_genome_accession"]
        ),
    params:
        genome_accession=config["ref"]["ncbi_genome_accession"],
        api_key=config["api_keys"]["ncbi"],
    container:
        "docker://staphb/ncbi-datasets:16.41.0"
    log:
        "logs/datasets/datasets_download_genome.log",
    message:
        "Downloading genome for NCBI genome accession: {genome}".format(
            genome=config["ref"]["ncbi_genome_accession"]
        )
    shell:
        """
        (datasets download genome accession {params.genome_accession} \
            --include gtf,gff3,rna,cds,protein,genome,seq-report \
            --api-key {params.api_key} \
            --filename {output}) &> {log}
        """


rule unzip_genome:
    input:
        "ncbi_dataset_{genome}.zip".format(
            genome=config["ref"]["ncbi_genome_accession"]
        ),
    output:
        multiext(
            "results/datasets/",
            "README.md",
        ),
        multiext(
            "results/datasets/ncbi_dataset/data/",
            "assembly_data_report.jsonl",
            "dataset_catalog.json",
        ),
        multiext(
            "results/datasets/ncbi_dataset/data/{genome}/".format(
                genome=config["ref"]["ncbi_genome_accession"]
            ),
            "genomic.gtf",
            "genomic.gff",
            "rna.fna",
            "cds_from_genomic.fna",
            "protein.faa",
            "sequence_report.jsonl",
        ),
        glob.glob(
            (
                "results/datasets/ncbi_dataset/data/{genome}/{genome}_"
                + "*"
                + "_genomic.fna"
            ).format(genome=config["ref"]["ncbi_genome_accession"])
        ),
    container:
        "docker://quay.io/biocontainers/p7zip:16.02"
    log:
        "logs/datasets/unzip_genome.log",
    message:
        "Unzipping genome for NCBI genome accession: {genome}".format(
            genome=config["ref"]["ncbi_genome_accession"]
        )
    shell:
        """
        (mkdir -p results/datasets && 7z x {input} -oresults/datasets) &> {log}
        """
