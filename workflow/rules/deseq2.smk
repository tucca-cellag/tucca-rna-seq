# Rule 1: Create DESeqDataSet for each analysis configuration
# DDS is parameterized by the analysis name from config.yaml
rule DESeqDataSet_from_ranged_se_per_analysis:
    input:
        se="resources/tximeta/tximeta_gse.RDS",
    output:
        dds="resources/deseq2/{analysis_name}/dds.RDS",
    log:
        "logs/deseq2/{analysis_name}/DESeqDataSet_se.log",
    threads:
        lambda wildcards: get_analysis_config_by_name(wildcards.analysis_name)["deseqdataset"]["threads"]
    params:
        formula=lambda wildcards: get_analysis_config_by_name(wildcards.analysis_name)["deseqdataset"]["formula"],
        min_counts=lambda wildcards: get_analysis_config_by_name(wildcards.analysis_name)["deseqdataset"]["min_counts"],
        extra=lambda wildcards: get_analysis_config_by_name(wildcards.analysis_name)["deseqdataset"]["extra"],
    wrapper:
        "v6.2.0/bio/deseq2/deseqdataset"


# Rule 2: Run DESeq2 Wald test for each specified contrast within each analysis
# DDS comes from the DESeqDataSet_from_ranged_se_per_analysis rule, specific to the {analysis_name}
rule deseq2_wald_per_analysis:
    input:
        dds="resources/deseq2/{analysis_name}/dds.RDS",
    output:
        # dds after DESeq()
        wald_rds="resources/deseq2/{analysis_name}/{contrast_name}.wald.RDS",
        # res for a given contrast [generated from dds after DESeq()]
        wald_tsv="resources/deseq2/{analysis_name}/{contrast_name}.dge.tsv",
        # counts(dds) as a TSV
        normalized_counts_table="resources/deseq2/{analysis_name}/{contrast_name}.counts.tsv",
        # counts(dds) as an RDS
        normalized_counts_rds="resources/deseq2/{analysis_name}/{contrast_name}.counts.RDS",
    log:
        "logs/deseq2/{analysis_name}/wald_{contrast_name}.log",
    threads:
        lambda wildcards: get_analysis_config_by_index(
            get_contrast_job_details(
                wildcards.analysis_name, wildcards.contrast_name
            )["config_index"]
        )["wald"]["threads"]
    params:
        # get user extra parameters from config.yaml
        deseq_extra=lambda wildcards: get_analysis_config_by_index(
            get_contrast_job_details(
                wildcards.analysis_name, wildcards.contrast_name
            )["config_index"]
        )["wald"]["deseq_extra"],
        shrink_extra=lambda wildcards: get_analysis_config_by_index(
            get_contrast_job_details(
                wildcards.analysis_name, wildcards.contrast_name
            )["config_index"]
        )["wald"]["shrink_extra"],
        results_extra=lambda wildcards: get_analysis_config_by_index(
            get_contrast_job_details(
                wildcards.analysis_name, wildcards.contrast_name
            )["config_index"]
        )["wald"]["results_extra"],
        # get contrast to be evaluated
        contrast=lambda wildcards: get_contrast_job_details(
            wildcards.analysis_name, wildcards.contrast_name
        )["elements"],
    wrapper:
        "v6.2.0/bio/deseq2/wald"


rule get_results_from_all_deseq_analyses:
    input:
        # Request all DDS files (one per analysis_name)
        expand("resources/deseq2/{analysis_name}/dds.RDS", analysis_name=DESEQ_ANALYSES_NAMES),
        # Request all Wald test outputs (for each contrast in each analysis)
        expand(
            [
                "resources/deseq2/{analysis_name}/{contrast_name}.wald.RDS",
                "resources/deseq2/{analysis_name}/{contrast_name}.dge.tsv",
                "resources/deseq2/{analysis_name}/{contrast_name}.counts.tsv",
                "resources/deseq2/{analysis_name}/{contrast_name}.counts.RDS",
            ],
            zip, # Use zip to correctly pair analysis_name and contrast_name
            analysis_name=[job["analysis_name"] for job in CONTRAST_JOBS],
            contrast_name=[job["contrast_name"] for job in CONTRAST_JOBS]
        )
    output:
        # A dummy output to mark completion
        touch("workflow/deseq2_analyses_complete.done")

# **Explanation of Changes:**
# 
# *   **`config/config.yaml`**:
#     *   The `diffexp.deseq2.analyses` list allows you to define multiple, independent DESeq2 runs.
#     *   Each entry in `analyses` has a `name` (for file organization), a `formula`, and an optional list of `contrasts`.
#     *   Each `contrast` has a `name` (for its specific output files) and `elements` (the 3-part list: factor, numerator, denominator).
#     *   Global defaults for `min_counts`, `threads`, and various `extra` parameters can be set, and then optionally overridden within each specific analysis block.
#     *   `deseqdataset_factor` and `deseqdataset_reference_level` are added to analysis blocks to allow for factor releveling during DDS creation if needed for a specific formula.
# 
# *   **`workflow/rules/deseq2.smk`**:
#     *   Helper Python code is added at the top to parse the `analyses` and `contrasts` from the config. This makes it easier to use with Snakemake's `expand` and parameter functions.
#     *   `DESeqDataSet_from_ranged_se_multi`:
#         *   Generates a `dds.RDS` file for each `analysis_name` (e.g., `resources/deseq2/main_treatment_effects/dds.RDS`).
#         *   Pulls `formula`, `min_counts`, `extra`, `factor`, and `reference_level` from the corresponding analysis block in the config.
#     *   `deseq2_wald_multi`:
#         *   Takes the `dds.RDS` specific to its `{analysis_name}`.
#         *   Generates output files (Wald RDS, DGE table, counts tables) that include both `{analysis_name}` and `{contrast_name}` in their paths (e.g., `resources/deseq2/main_treatment_effects/treatmentA_vs_control.dge.tsv`).
#         *   Pulls the `contrast` elements and `extra` parameters from the config, allowing for overrides at the analysis level.
#     *   `all_deseq_analyses`: An optional rule that you can use as a target (`snakemake all_deseq_analyses`) to run all configured DESeq2 steps.
# 
# **To Use:**
# 
# 1.  Replace the content of your `config/config.yaml` under `diffexp.deseq2` with the new structure, populating it with your desired formulas and contrasts.
# 2.  Replace your existing `workflow/rules/deseq2.smk` with the new content.
# 3.  Run Snakemake, potentially targeting the `all_deseq_analyses` rule or specific output files you need.
# 
# This setup provides a flexible way to manage multiple DESeq2 analyses within your existing workflow structure. Remember to adjust the `name` fields in the config to be descriptive for your experiments, as these will be used in file and directory names.
# 
# **Key Changes and Assumptions:**
# 
# *   **Explicit Definitions:** The `config.yaml` must now explicitly define all specified parameters (`threads`, `min_counts`, `deseqdataset_extra`, `deseqdataset_factor`, `deseqdataset_reference_level`, `wald_deseq_extra`, `wald_shrink_extra`, `wald_results_extra`, and `contrasts` list) within each block under `diffexp.deseq2.analyses`.
# *   **Empty Strings/Lists for Non-Applicable:** If a parameter like `deseqdataset_factor` isn't used for an analysis, or if an analysis has no contrasts, the config should still include the key with an appropriate empty value (e.g., `deseqdataset_factor: ""` or `contrasts: []`) to satisfy schema validation.
# *   **Direct Key Access:** The Snakemake rules now use direct dictionary access (e.g., `config_block["threads"]`). This relies on your Snakemake YAML schema validation to ensure these keys are present before the workflow logic runs, preventing `KeyError` exceptions.
# *   **`threads` Parameter:** The current setup assumes a single `threads` parameter per analysis block is used for both the `DESeqDataSet_from_ranged_se_multi` rule and the `deseq2_wald_multi` rule. If you need separate thread counts for these rules within the same analysis, you'd need to add distinct parameters to your `config.yaml` (e.g., `deseqdataset_threads` and `wald_threads`) and adjust the rules to access the correct one.
# 
# This approach aligns with your requirement for a fully explicit configuration, making it easier to validate with a schema and ensuring all settings for an analysis are clearly recorded in the `config.yaml`.