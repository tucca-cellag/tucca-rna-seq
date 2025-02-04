# sra_tools.smk


rule download_sra:
    output:
        # the wildcard name must be accession, pointing to an SRA number
        "data/pe/{accession}_1.fastq.gz",
        "data/pe/{accession}_2.fastq.gz",
    log:
        "logs/pe/download_{accession}.gz.log",
    params:
        extra="--split-files --skip-technical",
    threads: 6  # defaults to 6
    wrapper:
        "v5.2.1/bio/sra-tools/fasterq-dump"
