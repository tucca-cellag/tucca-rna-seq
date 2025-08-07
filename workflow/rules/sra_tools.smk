# workflow/rules/sra_tools.smk


# Rule order to resolve ambiguity between download and subsample rules
# When subsampling is enabled, subsample_sra_pe_reads should take precedence
ruleorder: subsample_sra_pe_reads > download_sra_pe_reads


rule configure_sra_tools:
    output:
        touch("resources/sra_tools/sra_config_completed.done"),
    params:
        vdb_config_ra_path=config["params"]["sra_tools"]["vdb_config_ra_path"],
    log:
        "logs/sra_tools/configure_sra_tools.log",
    conda:
        "../envs/sra_tools.yaml"
    shell:
        """
        # Check SRA-tools version for compatibility
        (fastq-dump --version) &> {log}
        (vdb-config --set {params.vdb_config_ra_path} --verbose) &> {log}
        """


rule prefetch_sra:
    input:
        "resources/sra_tools/sra_config_completed.done",
    output:
        multiext("data/sra_cache/{accession}/", "{accession}.sra"),
    log:
        "logs/sra_tools/prefetch/prefetch_{accession}.log",
    threads: 6
    conda:
        "../envs/sra_tools.yaml"
    shell:
        """
        (prefetch {wildcards.accession} -O ./data/sra_cache --verbose) &> {log}
        """


rule download_sra_pe_reads:
    input:
        "resources/sra_tools/sra_config_completed.done",
        "data/sra_cache/{accession}/{accession}.sra",
    output:
        "data/sra_reads/{accession}_1.fastq",
        "data/sra_reads/{accession}_2.fastq",
    log:
        "logs/sra_tools/fasterq_dump/fasterq_dump_{accession}.log",
    threads: 6
    conda:
        "../envs/sra_tools.yaml"
    shell:
        """
        (fasterq-dump ./data/sra_cache/{wildcards.accession} \
        -O ./data/sra_reads -e {threads} --split-files --verbose) &> {log}
        """


rule subsample_sra_pe_reads:
    input:
        "resources/sra_tools/sra_config_completed.done",
        "data/sra_cache/{accession}/{accession}.sra",
    output:
        "data/sra_reads/{accession}_1.fastq",
        "data/sra_reads/{accession}_2.fastq",
    params:
        num_reads=lambda wildcards: get_sra_subsample_params()["num_reads"],
        skip_reads=lambda wildcards: get_sra_subsample_params()["skip_reads"],
    log:
        "logs/sra_tools/subsample/subsample_{accession}.log",
    threads: 2
    conda:
        "../envs/sra_tools.yaml"
    shell:
        """
        (fastq-dump ./data/sra_cache/{wildcards.accession} \
        --minSpotId {params.skip_reads} \
        --readids \
        --read-filter pass \
        --dumpbase \
        --clip \
        --split-files \
        --outdir ./data/sra_reads \
        -N {params.num_reads} \
        -v \
        --force) &> {log}
        """


# Aggregates all SRA-based paired-end FASTQ files
#   - Used for targeting download of SRA reads for testing purposes, because
#     targeting download_sra_pe_reads directly is not possible, due to
#     wildcards in input
#   - Target via: snakemake resources/sra_tools/sra_pe_aggregate.done
# TODO: Refactor CI testing to use '--omit-from' tag rather than this rule
rule aggregate_sra_pe_reads:
    # Using a lambda func so the list of SRA accessions is computed at runtime
    input:
        lambda wildcards: expand(
            "data/sra_reads/{accession}_1.fastq",
            accession=[r.sra for _, r in units.iterrows() if is_sra_read(r)],
        )
        + expand(
            "data/sra_reads/{accession}_2.fastq",
            accession=[r.sra for _, r in units.iterrows() if is_sra_read(r)],
        ),
    output:
        touch("resources/sra_tools/sra_pe_aggregate.done"),
