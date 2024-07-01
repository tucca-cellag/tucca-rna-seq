rule fastqc:
    input:
        get_fq_files,
    output:
        htmls=["results/fastqc/{sample}_{unit}_{read}.html"],
        zips=["results/fastqc/{sample}_{unit}_{read}_fastqc.zip"],
    params:
        extra="--quiet",
    threads: 1
    resources:
        mem_mb=1024,
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/fastqc/{sample}_{unit}_{read}.log",
    run:
        fq = input[0]
        shell(
            f"""
            fastqc --threads {threads} --memory {resources.mem_mb} \
            {params.extra} --outdir results/fastqc/ {fq} &> {log}
            """
        )
