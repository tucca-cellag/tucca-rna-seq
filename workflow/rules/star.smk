# workflow/rules/star.smk


rule star:
    input:
        reads=get_paired_reads,
        star_index=multiext(
            "results/star/{genome_asc}_index/".format(
                genome_asc=config["ref_assembly"]["accession"]
            ),
            "chrLength.txt",
            "chrName.txt",
            "chrNameLength.txt",
            "chrStart.txt",
            "exonGeTrInfo.tab",
            "exonInfo.tab",
            "geneInfo.tab",
            "Genome",
            "genomeParameters.txt",
            "Log.out",
            "SA",
            "SAindex",
            "sjdbInfo.txt",
            "sjdbList.fromGTF.out.tab",
            "sjdbList.out.tab",
            "transcriptInfo.tab",
        ),
    output:
        multiext(
            "results/star/{sample_unit}_",
            "Aligned.sortedByCoord.out.bam",
            "Log.final.out",
            "Log.out",
            "Log.progress.out",
            "SJ.out.tab",
        ),
        directory("results/star/{sample_unit}__STARtmp"),
    log:
        "logs/star/star_{sample_unit}.log",
    params:
        star_index_dir=lambda wildcards, input: os.path.dirname(input.star_index[0]),
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
    message:
        """
        Running STAR alignment for:
            sample_unit = {wildcards.sample_unit},
        Running STAR with the inputs:
            {input.reads[0]}
            {input.reads[1]}
        """
    shell:
        """
        (echo "Running STAR alignment for sample_unit={wildcards.sample_unit}"
        echo "Running STAR with the inputs: {input.reads[0]} {input.reads[1]}"

        # Dynamically set readFilesCommand: use the provided command if gz,
        # else use 'cat'
        readFilesCmd=""
        if [[ "{input.reads[0]}" == *.gz ]]; then
            readFilesCmd="{params.readFilesCommand}"
        else
            readFilesCmd="cat"
        fi

        STAR --runThreadN {threads} \
        --genomeDir {params.star_index_dir} \
        --readFilesIn {input.reads[0]} {input.reads[1]} \
        --readFilesCommand $readFilesCmd \
        --outFileNamePrefix results/star/{wildcards.sample_unit}_ \
        --outSAMtype {params.outSAMtype} \
        --outSAMunmapped {params.outSAMunmapped} \
        --outSAMattributes {params.outSAMattributes} \
        --outFilterMultimapNmax {params.outFilterMultimapNmax} \
        --outFilterScoreMinOverLread {params.outFilterScoreMinOverLread} \
        --outFilterMatchNminOverLread {params.outFilterMatchNminOverLread} \
        --alignIntronMin {params.alignIntronMin} \
        --alignIntronMax {params.alignIntronMax} \
        {params.extra}) &> {log}
        """
