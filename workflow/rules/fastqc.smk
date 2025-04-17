# workflow/rules/fastqc.smk


rule fastqc:
    input:
        get_fq_files,
    output:
        htmls="results/fastqc/{sample_unit}_{read}.html",
        zips="results/fastqc/{sample_unit}_{read}_fastqc.zip",
    params:
        extra=config["params"]["fastqc"]["extra"],
        memory=config["params"]["fastqc"]["memory"],
    threads: 12
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/fastqc/{sample_unit}_{read}.log",
    message:
        """
        Generating FastQC report for:
            sample_unit = {wildcards.sample_unit}
            read = {wildcards.read}
        """
    shell:
        """
        (# Perform fastqc on each read
        fastqc --threads {threads} --memory {params.memory} \
        {params.extra} --outdir results/fastqc/ {input}) &> {log}
        """
