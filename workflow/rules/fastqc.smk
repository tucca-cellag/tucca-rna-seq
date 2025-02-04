# workflow/rules/fastqc.smk


rule fastqc:
    input:
        get_fq_files,
    output:
        htmls="results/fastqc/{sample}_{unit}_{read}.html",
        zips="results/fastqc/{sample}_{unit}_{read}_fastqc.zip",
    message:
        "Running FASTQC on {wildcards.sample} {wildcards.unit} {wildcards.read}"
    params:
        extra=config["params"]["fastqc"]["extra"],
    threads: 1
    resources:
        mem_mb=1024,
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/fastqc/{sample}_{unit}_{read}.log",
    shell:
        """
        # Perform fastqc on each read
        fastqc --threads {threads} --memory {resources.mem_mb} \
        {params.extra} --outdir results/fastqc/ {input} &> {log}
        
        # Extract the base name to handle both .fq.gz and .fastq.gz extensions
        base_name=$(basename "{input}")
        if [[ "$base_name" == *.fq.gz ]]; then
            base_name=$(basename "{input}" .fq.gz)
        elif [[ "$base_name" == *.fastq.gz ]]; then
            base_name=$(basename "{input}" .fastq.gz)
        fi
        
        html_output_name=${{base_name}}_fastqc.html
        zip_output_name=${{base_name}}_fastqc.zip
        
        mv results/fastqc/$html_output_name {output.htmls}
        mv results/fastqc/$zip_output_name {output.zips}
        """
