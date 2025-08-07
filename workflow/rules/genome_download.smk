# workflow/rules/genome_download.smk

if config["ref_assembly"]["source"] in ("Ensembl", "GENCODE"):

    rule get_genome:
        output:
            "resources/ensembl/{species}.{genome_name}.dna.fa".format(
                species=config["ref_assembly"]["species"],
                genome_name=config["ref_assembly"]["name"],
            ),
        params:
            species=config["ref_assembly"]["species"],
            datatype="dna",
            build=config["ref_assembly"]["name"],
            release=config["ref_assembly"]["release"],
            chromosome=lambda wildcards: get_chromosome_param(),
        log:
            "logs/ensembl/get_genome.log",
        wrapper:
            "v5.10.0/bio/reference/ensembl-sequence"

    rule get_transcriptome:
        output:
            "resources/ensembl/{species}.{genome_name}.cdna.fa".format(
                species=config["ref_assembly"]["species"],
                genome_name=config["ref_assembly"]["name"],
            ),
        params:
            species=config["ref_assembly"]["species"],
            datatype="cdna",
            build=config["ref_assembly"]["name"],
            release=config["ref_assembly"]["release"],
        log:
            "logs/ensembl/get_transcriptome_cdna.log",
        wrapper:
            "v5.10.0/bio/reference/ensembl-sequence"

    rule get_annotation:
        output:
            "resources/ensembl/{species}.{genome_name}.gtf".format(
                species=config["ref_assembly"]["species"],
                genome_name=config["ref_assembly"]["name"],
            ),
        params:
            species=config["ref_assembly"]["species"],
            release=config["ref_assembly"]["release"],
            build=config["ref_assembly"]["name"],
            fmt="gtf",
        log:
            "logs/ensembl/get_annotation.log",
        wrapper:
            "v5.10.0/bio/reference/ensembl-annotation"


if config["ref_assembly"]["source"] in ("RefSeq"):

    # TODO: Add more warnings for datasets_download_genome and make sure inputs
    # are properly sanitized

    # https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/
    rule datasets_download_genome:
        output:
            "ncbi_dataset_{genome_asc}.zip".format(
                genome_asc=config["ref_assembly"]["accession"]
            ),
        params:
            genome_asc=config["ref_assembly"]["accession"],
            api_key=config["api_keys"]["ncbi"],
            chromosome=lambda wildcards: get_chromosome_param(),
        conda:
            "../envs/ncbi_datasets.yaml"
        log:
            "logs/datasets/datasets_download_genome.log",
        message:
            "Downloading genome for NCBI genome accession: {genome_asc}".format(
                genome_asc=config["ref_assembly"]["accession"]
            )
        shell:
            """
            # RefSeq downloads always get full genome (single chromosome not supported for alignment)
            echo "Downloading full genome for RefSeq"
            (datasets download genome accession {params.genome_asc} \
                --include gtf,rna,genome,seq-report \
                --api-key {params.api_key} \
                --filename {output}) &> {log}
            """

    rule unzip_genome:
        input:
            "ncbi_dataset_{genome_asc}.zip".format(
                genome_asc=config["ref_assembly"]["accession"]
            ),
        output:
            multiext("resources/datasets/", "README.md", "md5sum.txt"),
            multiext(
                "resources/datasets/ncbi_dataset/data/",
                "assembly_data_report.jsonl",
                "dataset_catalog.json",
            ),
            multiext(
                "resources/datasets/ncbi_dataset/data/{genome_asc}/".format(
                    genome_asc=config["ref_assembly"]["accession"]
                ),
                "genomic.gtf",
                "rna.fna",
                "sequence_report.jsonl",
                "{genome_asc}_{genome_name}_genomic.fna".format(
                    genome_asc=config["ref_assembly"]["accession"],
                    genome_name=config["ref_assembly"]["name"],
                ),
            ),
        conda:
            "../envs/p7zip.yaml"
        log:
            "logs/datasets/unzip_genome.log",
        message:
            "Unzipping genome for NCBI genome accession: {genome_asc}".format(
                genome_asc=config["ref_assembly"]["accession"]
            )
        shell:
            """
            (mkdir -p resources/datasets && \
            7z x {input} -oresources/datasets -y && \
            echo "Directory structure in resources/datasets:" && \
            ls -laR resources/datasets) &> {log}
            """
