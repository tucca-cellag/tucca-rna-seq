rule fastqc:
    input:
        unpack(get_fq),
    output:
        expand(
            [
                "results/fastqc/{sample}-{unit}_R{rep}_001.html",
                "results/fastqc/{sample}-{unit}_R{rep}_001_fastqc.zip",
            ],
            rep=[1, 2],
        ),
    params:
        extra="--quiet",
    threads: 1
    resources:
        mem_mb=1024,
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/fastqc/{sample}-{unit}.log",
    shell:
        """
        (fastqc --threads {threads} --memory {resources.mem_mb} {params.extra} \
        --outdir results/fastqc/ {input}) &> {log}
        """
