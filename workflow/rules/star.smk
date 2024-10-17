rule star:
    input:
        reads=get_paired_reads,
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
    message:
        "Mapping {wildcards.sample} {wildcards.unit} reads to genome"
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
        echo "Running STAR alignment for sample={wildcards.sample}, unit={wildcards.unit}" \
        echo "Running STAR with the inputs: {input.reads[0]} {input.reads[1]}"
        (STAR --runThreadN {threads} \
        --genomeDir {input.star_index} \
        --readFilesIn {input.reads[0]} {input.reads[1]} \
        --outFileNamePrefix results/star/{wildcards.sample}_{wildcards.unit}_ \
        --outSAMtype {params.outSAMtype} \
        --outSAMunmapped {params.outSAMunmapped} \
        --outSAMattributes {params.outSAMattributes} \
        --readFilesCommand {params.readFilesCommand} \
        --outFilterMultimapNmax {params.outFilterMultimapNmax} \
        --outFilterScoreMinOverLread {params.outFilterScoreMinOverLread} \
        --outFilterMatchNminOverLread {params.outFilterMatchNminOverLread} \
        --alignIntronMin {params.alignIntronMin} \
        --alignIntronMax {params.alignIntronMax} \
        {params.extra}) &> {log}
        """
