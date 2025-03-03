# workflow/rules/multiqc.smk


rule multiqc:
    input:
        lambda wc: get_final_output()[:-1],
    output:
        "results/multiqc/{report_name}.html".format(
            report_name=config["params"]["multiqc"]["report_name"]
        ),
        temp(
            directory(
                "results/multiqc/{report_name}_data".format(
                    report_name=config["params"]["multiqc"]["report_name"]
                )
            )
        ),
    params:
        report_name=config["params"]["multiqc"]["report_name"],
        overwrite_existing=config["params"]["multiqc"]["overwrite_existing"],
        extra=config["params"]["multiqc"]["extra"],
    container:
        config["containers"]["multiqc"]
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
        -c config/multiqc_config.yaml \
        {params.extra} \
        {input}) &> {log}
        """
