# workflow/rules/sra_tools.smk


# Rule order to resolve ambiguity between download and subsample rules
# When subsampling is enabled, subsample_sra_pe_reads should take precedence
# and prefetch_sra should be skipped to avoid storage issues


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


rule create_sra_flags:
    input:
        "resources/sra_tools/sra_config_completed.done",
    output:
        branch(
            is_sra_subsampling_enabled(),
            then="resources/sra_tools/subsampling_enabled.flag",
            otherwise="resources/sra_tools/prefetch_enabled.flag",
        ),
    log:
        "logs/sra_tools/create_flags.log",
    shell:
        """
        touch {output}
        """


rule prefetch_sra:
    input:
        "resources/sra_tools/prefetch_enabled.flag",
    output:
        multiext("data/sra_cache/{accession}/", "{accession}.sra"),
    log:
        "logs/sra_tools/prefetch/prefetch_{accession}.log",
    threads: 6
    conda:
        "../envs/sra_tools.yaml"
    shell:
        """
        # Only prefetch when subsampling is disabled to avoid storage issues
        (prefetch {wildcards.accession} -O ./data/sra_cache --verbose) &> {log}
        """


rule download_sra_pe_reads:
    input:
        "resources/sra_tools/prefetch_enabled.flag",
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
        "resources/sra_tools/subsampling_enabled.flag",
    output:
        "data/sra_reads/{accession}_1.fastq",
        "data/sra_reads/{accession}_2.fastq",
    params:
        min_spot_id=lambda wildcards: get_sra_subsample_params()["min_spot_id"],
        max_spot_id=lambda wildcards: get_sra_subsample_params()["max_spot_id"],
    log:
        "logs/sra_tools/subsample/subsample_{accession}.log",
    threads: 2
    conda:
        "../envs/sra_tools.yaml"
    shell:
        """
        (fastq-dump {wildcards.accession} \
        --minSpotId {params.min_spot_id} \
        --maxSpotId {params.max_spot_id} \
        --split-3 \
        --outdir ./data/sra_reads \
        -v) &> {log}
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
