rule star_pe_multi:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        # paired end reads needs to be ordered so each item in the two lists match
        get_paired_reads(),
        # path to STAR reference genome index
        idx="results/star/{genome}_index".format(
            genome=config["ref"]["ncbi_genome_accession"]
        ),
    output:
        # see STAR manual for additional output files
        aln="results/star/pe/{sample}/pe_aligned.sam",
        log="logs/pe/{sample}/Log.out",
        sj="results/star/pe/{sample}/SJ.out.tab",
        unmapped=[
            "results/star/pe/{sample}/unmapped.1.fastq.gz",
            "results/star/pe/{sample}/unmapped.2.fastq.gz",
        ],
    log:
        "logs/pe/{sample}.log",
    params:
        # optional parameters
        extra="--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --readFilesCommand zcat --outFilterMultimapNmax 1 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignIntronMin 1 --alignIntronMax 2500",
    threads: 8
    wrapper:
        "v3.13.3/bio/star/align"


""" rule star_se:
    input:
        fq1=get_fq_files(),
        # path to STAR reference genome index
        idx="results/star/{genome}_index".format(
            genome=config["ref"]["ncbi_genome_accession"]
        ),
    output:
        # see STAR manual for additional output files
        aln="star/se/{sample}/se_aligned.bam",
        log="logs/se/{sample}/Log.out",
        log_final="logs/se/{sample}/Log.final.out",
        unmapped="star/se/{sample}/unmapped.fastq",
    log:
        "logs/se/{sample}.log",
    params:
        # optional parameters
        extra="--outSAMtype BAM Unsorted",
    threads: 8
    wrapper:
        "v3.13.3/bio/star/align" """
