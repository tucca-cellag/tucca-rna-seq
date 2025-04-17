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
        {params.extra} --outdir results/fastqc/ {input}
        
        # Determine the base name by removing known fastq extensions
        base_name=$(basename "{input}")
        if [[ "$base_name" == *.fq.gz ]]; then
            base_name=$(basename "$base_name" .fq.gz)
        elif [[ "$base_name" == *.fastq.gz ]]; then
            base_name=$(basename "$base_name" .fastq.gz)
        elif [[ "$base_name" == *.fq ]]; then
            base_name=$(basename "$base_name" .fq)
        elif [[ "$base_name" == *.fastq ]]; then
            base_name=$(basename "$base_name" .fastq)
        fi
        
        html_output_name=${{base_name}}_fastqc.html
        zip_output_name=${{base_name}}_fastqc.zip
        
        mv results/fastqc/$html_output_name {output.htmls}
        mv results/fastqc/$zip_output_name {output.zips}) &> {log}
        """
