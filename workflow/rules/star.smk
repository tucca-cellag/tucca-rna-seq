# workflow/rules/star.smk


rule star_pe_multi:
    input:
        "resources/qualimap/qualimap_enabled.flag",
        fq1=lambda wildcards: [get_paired_reads(wildcards)[0]],
        fq2=lambda wildcards: [get_paired_reads(wildcards)[1]],
        # path to STAR reference genome index
        idx="resources/star/{genome_asc}_index/".format(
            genome_asc=config["ref_assembly"]["accession"]
        ),
    output:
        aln="resources/star/{sample_unit}/Aligned.sortedByCoord.out.bam",
        log="logs/star/{sample_unit}/Log.out",
        log_progress="logs/star/{sample_unit}/Log.progress.out",
        log_final="logs/star/{sample_unit}/Log.final.out",
        sj="resources/star/{sample_unit}/SJ.out.tab",
    log:
        "logs/star/{sample_unit}/star_pe_multi.log",
    params:
        extra=config["params"]["star"]["extra"],
    threads: 12
    wrapper:
        "v6.2.0/bio/star/align"
