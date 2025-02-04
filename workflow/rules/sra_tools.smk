# workflow/rules/sra_tools.smk


rule configure_sra_tools:
    output:
        "sra_config_completed.txt",
    log:
        "logs/sra_tools/configure_sra_tools.log",
    conda:
        "../envs/sra_tools.yaml"
    shell:
        """
        (vdb-config --set "/repository/user/main/remote_access=true"
        vdb-config --set "/repository/user/main/user_repository=$(pwd)/sra_cache"
        touch {output}) &> {log}
        """


rule prefetch:
    input:
        "sra_config_completed.txt",
    output:
        directory("sra_cache/{accession}"),
    log:
        "logs/sra_tools/prefetch/prefetch_{accession}.log",
    threads: 6
    conda:
        "../envs/sra_tools.yaml"
    shell:
        """
        (prefetch {wildcards.accession} --output-directory {output} --verbose \
        ) &> {log}
        """


rule download_sra_pe_reads:
    input:
        "sra_config_completed.txt",
        "sra_cache/{accession}.sra",
    output:
        "data/sra_reads/{accession}_1.fastq",
        "data/sra_reads/{accession}_2.fastq",
    log:
        "logs/sra_tools/fasterq_dump/download_{accession}.log",
    threads: 6
    conda:
        "../envs/sra_tools.yaml"
    shell:
        """
        (fasterq-dump {wildcards.accession} -e {threads} --split-files \
        --verbose) &> {log}
        """
