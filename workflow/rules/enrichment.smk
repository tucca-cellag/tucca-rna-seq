# This rule generates the dynamic post-deployment script
rule prepare_enrichment_env:
    output:
        script="workflow/envs/enrichment.post-deploy.sh",
    params:
        org_db_pkg=get_orgdb_pkg_name(config),
    log:
        "logs/enrichment/prepare_enrichment_env.log",
    script:
        "../scripts/prepare_enrichment_env.py"


rule run_enrichment:
    input:
        post_deploy_script="workflow/envs/enrichment.post-deploy.sh",
        dge_tsv="resources/deseq2/{analysis}/{contrast}/dge.tsv",
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
        "../scripts/run_enrichment.R"


rule all_enrichment:
    input:
        get_all_enrichment_dirs(config),
    output:
        touch("resources/enrichment/enrichment_analyses_complete.done"),
