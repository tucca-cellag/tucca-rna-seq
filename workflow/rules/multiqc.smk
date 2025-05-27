# workflow/rules/multiqc.smk


# rule multiqc:
#     input:
#
#     output:
#         "results/multiqc/{report_name}.html".format(
#             report_name=config["params"]["multiqc"]["report_name"]
#         ),
#         directory(
#             "results/multiqc/{report_name}_data".format(
#                 report_name=config["params"]["multiqc"]["report_name"]
#             )
#         ),
#     params:
#         report_name=config["params"]["multiqc"]["report_name"],
#         overwrite_existing=config["params"]["multiqc"]["overwrite_existing"],
#         multiqc_config_path=config["params"]["multiqc"]["multiqc_config_path"],
#         extra=config["params"]["multiqc"]["extra"],
#     conda:
#         "../envs/multiqc.yaml"
#     log:
#         "logs/multiqc/multiqc.log",
#     message:
#         "Compiling MultiQC report summary"
#     shell:
#         """
#         unset DISPLAY
#
#         (multiqc --verbose -n {params.report_name} \
#         -o results/multiqc/ \
#         {params.overwrite_existing} \
#         -c {params.multiqc_config_path} \
#         {params.extra} \
#         {input}) &> {log}
#         """


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
