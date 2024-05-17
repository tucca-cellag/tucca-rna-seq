rule fastqc:
    input:
        expand("data/{sample}.fastq", sample=samples),
    output:
        html="results/{sample}.html",
    shell:
        "fastqc {input} -o results/"
