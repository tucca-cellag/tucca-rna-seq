rule fastqc:
    input:
        unpack(get_fq),
    output:
        "results/fastqc/{sample}{unit}_R1_001.html",
        "results/fastqc/{sample}{unit}_R1_001_fastqc.zip",
        "results/fastqc/{sample}{unit}_R2_001.html",
        "results/fastqc/{sample}{unit}_R2_001_fastqc.zip",
    params:
        extra="--quiet",
    threads: 1
    resources:
        mem_mb=1024,
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/fastqc/{sample}{unit}.log",
    shell:
        """
        (fastqc --threads {threads} --memory {resources.mem_mb} {params.extra} \
        --outdir results/fastqc/ {input}) &> {log}
        """
