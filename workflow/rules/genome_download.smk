# workflow/rules/genome_download.smk

if config["genome"]["assembly_db"] == "ensembl":

    rule get_transcriptome:
        output:
            "resources/ensembl/transcriptome.cdna.fasta",
        params:
            species=config["resources"]["ref"]["species"],
            datatype="cdna",
            build=config["resources"]["ref"]["build"],
            release=config["resources"]["ref"]["release"],
        log:
            "logs/ensembl/get-transcriptome/cdna.log",
        wrapper:
            "v5.10.0/bio/reference/ensembl-sequence"

    rule get_annotation:
        output:
            "resources/ensembl/genome.gtf",
        params:
            species=config["resources"]["ref"]["species"],
            release=config["resources"]["ref"]["release"],
            build=config["resources"]["ref"]["build"],
            fmt="gtf",
        log:
            "logs/ensembl/get-annotation.log",
        wrapper:
            "v5.10.0/bio/reference/ensembl-annotation"


if config["genome"]["assembly_db"] == "ncbi_refseq":

    # https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/
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
                --include gtf,rna,genome,seq-report \
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
                "resources/datasets/",
                "README.md",
            ),
            multiext(
                "resources/datasets/ncbi_dataset/data/",
                "assembly_data_report.jsonl",
                "dataset_catalog.json",
            ),
            multiext(
                "resources/datasets/ncbi_dataset/data/{genome_asc}/".format(
                    genome_asc=config["genome"]["assembly_accession"]
                ),
                "genomic.gtf",
                "rna.fna",
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
            (mkdir -p resources/datasets && \
            7z x {input} -oresources/datasets && \
            echo "Directory structure in resources/datasets:" && \
            ls -laR resources/datasets) &> {log}
            """
