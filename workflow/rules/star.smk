rule star_pe_multi:
    input:
        get_paired_reads,
        star_index="results/star/{genome}_index".format(
            genome=config["ref"]["ncbi_genome_accession"]
        ),
    output:
        # see STAR manual for additional output files
        aln="results/star/pe/{sample}/pe_aligned.sam",
        log="logs/pe/{sample}/Log.out",
        sj="results/star/pe/{sample}/SJ.out.tab",
    log:
        "logs/pe/{sample}.log",
    params:
        "",
    threads: 12
    conda:
        "../envs/star.yaml"
    shell:
        """
        STAR --runThreadN {threads} \
        --genomeDir {input.star_index} \
        --readFilesIn {wildcards.fq1} {wildcards.fq2} \
        --outFileNamePrefix results/star/{sample}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --readFilesCommand zcat \
        --outFilterMultimapNmax 1 \
        --outFilterScoreMinOverLread 0 \
        --outFilterMatchNminOverLread 0 \
        --alignIntronMin 1 \
        --alignIntronMax 2500
        """
