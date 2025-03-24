# workflow/rules/sra_tools.smk


rule configure_sra_tools:
    output:
        "results/sra_tools/sra_config_completed.txt",
    params:
        vdb_config_ra_path=config["params"]["sra_tools"]["vdb_config_ra_path"],
    log:
        "logs/sra_tools/configure_sra_tools.log",
    container:
        config["containers"]["sra_tools"]
    shell:
        """
        (vdb-config --set {params.vdb_config_ra_path} --verbose
        touch {output}) &> {log}
        """


rule prefetch_sra:
    input:
        "results/sra_tools/sra_config_completed.txt",
    output:
        multiext("data/sra_cache/{accession}/", "{accession}.sra"),
    log:
        "logs/sra_tools/prefetch/prefetch_{accession}.log",
    threads: 6
    container:
        config["containers"]["sra_tools"]
    shell:
        """
        (prefetch {wildcards.accession} -O ./data/sra_cache --verbose) &> {log}
        """


rule download_sra_pe_reads:
    input:
        "results/sra_tools/sra_config_completed.txt",
        "data/sra_cache/{accession}/{accession}.sra",
    output:
        "data/sra_reads/{accession}_1.fastq",
        "data/sra_reads/{accession}_2.fastq",
    log:
        "logs/sra_tools/fasterq_dump/fasterq_dump_{accession}.log",
    threads: 6
    container:
        config["containers"]["sra_tools"]
    shell:
        """
        (fasterq-dump ./data/sra_cache/{wildcards.accession} \
        -O ./data/sra_reads -e {threads} --split-files --verbose) &> {log}
        """
