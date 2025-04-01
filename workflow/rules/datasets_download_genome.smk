# workflow/rules/datasets_download_genome.smk

import glob


rule datasets_download_genome:
    output:
        "ncbi_dataset_{genome_asc}.zip".format(
            genome_asc=config["genome"]["assembly_accession"]
        ),
    params:
        genome_asc=config["genome"]["assembly_accession"],
        api_key=config["api_keys"]["ncbi"],
    conda:
        "../envs/ncbi_datasets.yaml"
    log:
        "logs/datasets/datasets_download_genome.log",
    message:
        "Downloading genome for NCBI genome accession: {genome_asc}".format(
            genome_asc=config["genome"]["assembly_accession"]
        )
    shell:
        """
        (datasets download genome accession {params.genome_asc} \
            --include gtf,gff3,rna,cds,protein,genome,seq-report \
            --api-key {params.api_key} \
            --filename {output}) &> {log}
        """


rule unzip_genome:
    input:
        "ncbi_dataset_{genome_asc}.zip".format(
            genome_asc=config["genome"]["assembly_accession"]
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
            "results/datasets/ncbi_dataset/data/{genome_asc}/".format(
                genome_asc=config["genome"]["assembly_accession"]
            ),
            "genomic.gtf",
            "genomic.gff",
            "rna.fna",
            "cds_from_genomic.fna",
            "protein.faa",
            "sequence_report.jsonl",
            "{genome_asc}_{genome_name}_genomic.fna".format(
                genome_asc=config["genome"]["assembly_accession"],
                genome_name=config["genome"]["assembly_name"],
            ),
        ),
    conda:
        "../envs/p7zip.yaml"
    log:
        "logs/datasets/unzip_genome.log",
    message:
        "Unzipping genome for NCBI genome accession: {genome_asc}".format(
            genome_asc=config["genome"]["assembly_accession"]
        )
    shell:
        """
        (mkdir -p results/datasets && \
        7z x {input} -oresults/datasets && \
        echo "Directory structure in results/datasets:" && \
        ls -laR results/datasets) &> {log}
        """
