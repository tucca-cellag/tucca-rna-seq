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
        (vdb-config --set "/repository/user/main/remote_access=true" \
            --set "/repository/user/main/user_repository=$(pwd)/sra_cache" \
            --save && \
        touch {output}) &> {log}
        """


rule download_sra_pe_reads:
    input:
        "sra_config_completed.txt",
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
        (fasterq-dump {wildcards.accession} --split-files) &> {log}
        """
