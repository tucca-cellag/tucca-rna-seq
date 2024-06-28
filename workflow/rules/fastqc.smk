rule fastqc:
    input:
        unpack(get_fq),
    output:
        html="results/fastqc/{sample}-{unit}.html",
        zip="results/fastqc/{sample}-{unit}_fastqc.zip",
    params:
        extra="--quiet",
    threads: 1
    resources:
        mem_mb=1024,
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/fastqc/{sample}.log",
    shell:
        """
        (fastqc --threads {threads} --memory {resources.mem_mb} {params.extra} \
        --outdir results/fastqc/ {input}) &> {log}
        """
