# workflow/rules/multiqc.smk


rule multiqc:
    input:
        lambda wc: get_final_output()[:-2],
        config="config/multiqc_config.yaml",
    output:
        report(
            "results/multiqc/multiqc.html",
            category="Quality control",
            caption="../report/multiqc.rst",
        ),
    params:
        extra=config["params"]["multiqc"]["extra"],
    log:
        "logs/multiqc/multiqc.log",
    wrapper:
        "v6.2.0/bio/multiqc"
