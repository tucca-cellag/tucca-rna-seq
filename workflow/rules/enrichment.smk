rule prepare_enrichment_env:
    input:
        template="workflow/envs/enrichment_template.yaml",
        # Add config to inputs to ensure rule re-runs if config changes
        config="config/config.yaml",
    output:
        env_file="workflow/envs/enrichment.yaml",
    params:
        org_db_pkg=get_orgdb_pkg_name(config),
    log:
        "logs/enrichment/prepare_enrichment_env.log",
    script:
        "workflow/scripts/prepare_enrichment_env.py"


rule run_enrichment:
    input:
        dge_tsv="resources/deseq2/{analysis}/{contrast}/dge.tsv",
        env_yaml="workflow/envs/enrichment.yaml",
    output:
        directory("resources/enrichment/{analysis}/{contrast}"),
    params:
        alpha_pathway=config["enrichment"]["alpha_pathway"],
        org_db_pkg=get_orgdb_pkg_name(config),
        targets=config["enrichment"]["targets"],
    wildcard_constraints:
        analysis="[^/]+",
        contrast="[^/]+",
    log:
        "logs/enrichment/{analysis}/{contrast}/enrichment.log",
    conda:
        "../envs/enrichment.yaml"
    script:
        "workflow/scripts/run_enrichment.R"


rule all_enrichment:
    input:
        get_all_enrichment_dirs(config),
    output:
        touch("resources/enrichment/enrichment_analyses_complete.done"),
