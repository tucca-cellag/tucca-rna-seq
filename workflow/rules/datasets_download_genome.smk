import glob


rule datasets_download_genome:
    output:
        "ncbi_dataset_{genome}.zip".format(
            genome=config["ref"]["ncbi_genome_accession"]
        ),
    params:
        genome_accession=config["ref"]["ncbi_genome_accession"],
        api_key=config["api_keys"]["ncbi"],
    conda:
        "../envs/ncbi_datasets_cli.yaml"
    log:
        "logs/datasets/datasets_download_genome.log",
    shell:
        """
        (datasets download genome accession {params.genome_accession} \
            --include gtf,gff3,rna,cds,protein,genome,seq-report \
            --api-key {params.api_key} \
            --filename {output} > {output}) &> {log}
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
    log:
        "logs/datasets/unzip_genome.log",
    shell:
        """
        (mkdir -p results/datasets
        unzip -o {input} -d results/datasets) &> {log}
        """
