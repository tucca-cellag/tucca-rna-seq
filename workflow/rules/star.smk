# workflow/rules/star.smk


rule star_pe_multi:
    input:
        fq1=lambda wildcards: [get_paired_reads(wildcards)[0]],
        fq2=lambda wildcards: [get_paired_reads(wildcards)[1]],
        # path to STAR reference genome index
        idx="resources/star/{genome_asc}_index/".format(
            genome_asc=config["ref_assembly"]["accession"]
        ),
    output:
        aln="resources/star/{sample_unit}_Aligned.sortedByCoord.out.bam",
        log="resources/star/{sample_unit}_Log.out",
        log_progress="resources/star/{sample_unit}_Log.progress.out",
        log_final="resources/star/{sample_unit}_Log.final.out",
        sj="resources/star/{sample_unit}_SJ.out.tab",
    log:
        "logs/star/star_{sample_unit}.log",
    params:
        extra=config["params"]["star"]["extra"],
    threads: 12
    wrapper:
        "v6.2.0/bio/star/align"
