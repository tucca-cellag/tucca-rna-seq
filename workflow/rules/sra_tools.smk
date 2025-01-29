# sra_tools.smk


rule download_sra:
    output:
        # the wildcard name must be accession, pointing to an SRA number
        "data/pe/{accession}_1.fastq",
        "data/pe/{accession}_2.fastq",
    log:
        "logs/pe/{accession}.log",
    params:
        extra="--skip-technical",
    threads: 6  # defaults to 6
    wrapper:
        "v5.2.1/bio/sra-tools/fasterq-dump"
