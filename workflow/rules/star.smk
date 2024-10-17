rule star:
    input:
        reads=lambda wildcards: get_paired_reads(wildcards),
        star_index="results/star/{genome}_index".format(
            genome=config["ref"]["ncbi_genome_accession"]
        ),
    output:
        multiext(
            "results/star/{sample}_{unit}_",
            "Aligned.sortedByCoord.out.bam",
            "Log.final.out",
            "Log.out",
            "Log.progress.out",
            "SJ.out.tab",
        ),
        temp(directory("results/star/{sample}_{unit}__STARtmp")),
    log:
        "logs/star/star_{sample}_{unit}.log",
    params:
        outSAMtype=config["params"]["star"]["outSAMtype"],
        outSAMunmapped=config["params"]["star"]["outSAMunmapped"],
        outSAMattributes=config["params"]["star"]["outSAMattributes"],
        readFilesCommand=config["params"]["star"]["readFilesCommand"],
        outFilterMultimapNmax=config["params"]["star"]["outFilterMultimapNmax"],
        outFilterScoreMinOverLread=config["params"]["star"][
            "outFilterScoreMinOverLread"
        ],
        outFilterMatchNminOverLread=config["params"]["star"][
            "outFilterMatchNminOverLread"
        ],
        alignIntronMin=config["params"]["star"]["alignIntronMin"],
        alignIntronMax=config["params"]["star"]["alignIntronMax"],
        extra=config["params"]["star"]["extra"],
    threads: 12
    conda:
        "../envs/star.yaml"
    shell:
        """
        (set -x # activate debugging
        echo "Running STAR alignment for sample={wildcards.sample}, unit={wildcards.unit}\n"
        echo "inputs: {input.reads[0]['fq1']} {input.reads[0]['fq2']}"
        STAR --runThreadN {threads} \
        --genomeDir {input.star_index} \
        --readFilesIn {input.reads[0]['fq1']} {input.reads[0]['fq2']} \
        --outFileNamePrefix results/star/{sample}_{unit}_ \
        --outSAMtype {params.outSAMtype} \
        --outSAMunmapped {params.outSAMunmapped} \
        --outSAMattributes {params.outSAMattributes} \
        --readFilesCommand {params.readFilesCommand} \
        --outFilterMultimapNmax {params.outFilterMultimapNmax} \
        --outFilterScoreMinOverLread {params.outFilterScoreMinOverLread} \
        --outFilterMatchNminOverLread {params.outFilterMatchNminOverLread} \
        --alignIntronMin {params.alignIntronMin} \
        --alignIntronMax {params.alignIntronMax} \
        {params.extra} \
        set -x # deactivate debugging) &> {log}
        """
