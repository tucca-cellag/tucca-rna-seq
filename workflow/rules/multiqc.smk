# workflow/rules/multiqc.smk


rule multiqc:
    input:
        lambda wc: get_final_output()[:-1],
    output:
        "results/multiqc/{report_name}.html".format(
            report_name=config["params"]["multiqc"]["report_name"]
        ),
        directory(
            "results/multiqc/{report_name}_data".format(
                report_name=config["params"]["multiqc"]["report_name"]
            )
        ),
    params:
        report_name=config["params"]["multiqc"]["report_name"],
        overwrite_existing=config["params"]["multiqc"]["overwrite_existing"],
        multiqc_config_path=config["params"]["multiqc"]["multiqc_config_path"],
        extra=config["params"]["multiqc"]["extra"],
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/multiqc/multiqc.log",
    message:
        "Compiling MultiQC report summary"
    shell:
        """
        unset DISPLAY

        (multiqc --verbose -n {params.report_name} \
        -o results/multiqc/ \
        {params.overwrite_existing} \
        -c {params.multiqc_config_path} \
        {params.extra} \
        {input}) &> {log}
        """
