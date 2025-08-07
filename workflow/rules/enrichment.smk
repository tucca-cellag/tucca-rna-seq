# workflow/rules/enrichment.smk


# Rule to get the taxonomy ID using the assembly accession
rule get_tax_id:
    output:
        tax_id_file="resources/enrichment/tax_id.txt",
    params:
        accession=config["ref_assembly"]["accession"],
        api_key=config["api_keys"]["ncbi"],
    log:
        "logs/enrichment/get_tax_id.log",
    conda:
        "../envs/ncbi_datasets.yaml"
    shell:
        """
        (datasets summary genome accession {params.accession} --api-key {params.api_key} | \
        jq -r '.assemblies[0].assembly.org.tax_id' > {output.tax_id_file}) 2> {log}
        """


# Rule to build the OrgDb package locally if it's not on Bioconda
rule build_local_orgdb:
    input:
        tax_id_file="resources/enrichment/tax_id.txt",
        gtf=branch(
            config["ref_assembly"]["source"] == "RefSeq",
            then="resources/datasets/ncbi_dataset/data/{genome_asc}/genomic.gtf".format(
                genome_asc=config["ref_assembly"]["accession"],
            ),
            otherwise="resources/ensembl/{species}.{genome_name}.gtf".format(
                species=config["ref_assembly"]["species"],
                genome_name=config["ref_assembly"]["name"],
            ),
        ),
    output:
        out_dir=directory("resources/enrichment/local_orgdb_build"),
    params:
        species=config["ref_assembly"]["species"].replace("_", " "),
        genus=config["ref_assembly"]["species"].split("_")[0],
        version=config["enrichment"]["annotationforge"]["version"],
        author=config["enrichment"]["annotationforge"]["author"],
        extra=config["enrichment"]["annotationforge"]["extra"],
    log:
        "logs/enrichment/build_local_orgdb.log",
    conda:
        "../envs/annotationforge.yaml"
    script:
        "../scripts/build_local_orgdb.R"


# Rule to install the appropriate OrgDb package
rule install_orgdb:
    output:
        flag=touch(get_orgdb_install_flag(config)),
        params_rds="resources/enrichment/enrichment_params.RDS",
    params:
        enrichment=get_enrichment_params,
    log:
        "logs/enrichment/install_orgdb.log",
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/install_orgdb.R"


# SPIA rules are conditionally executed based on config setting
if config["enrichment"]["spia"]["enabled"]:

    rule make_spia_data:
        output:
            spia_data_dir=directory("resources/enrichment/spia_data"),
        params:
            enrichment=get_enrichment_params,
        log:
            "logs/enrichment/make_spia_data.log",
        conda:
            "../envs/r_env.yaml"
        script:
            "../scripts/make_spia_data.R"


rule clusterprofiler_gsea:
    input:
        unpack(get_enrichment_deps),
        msigdb_flag=lambda wildcards: (
            "resources/msigdb/msigdb_loaded.flag"
            if config["enrichment"]["msigdb"]["enabled"]
            else "resources/msigdb/msigdb_disabled.flag"
        ),
    output:
        gsea_rds="resources/enrichment/{analysis}/{contrast}/gsea_results.RDS",
    params:
        enrichment=get_enrichment_params,
    log:
        "logs/enrichment/{analysis}/{contrast}/gsea.log",
    conda:
        "../envs/r_env.yaml"
    wildcard_constraints:
        analysis="[^/]+",
        contrast="[^/]+",
    script:
        "../scripts/run_gsea.R"


rule clusterprofiler_ora:
    input:
        unpack(get_enrichment_deps),
        msigdb_flag=lambda wildcards: (
            "resources/msigdb/msigdb_loaded.flag"
            if config["enrichment"]["msigdb"]["enabled"]
            else "resources/msigdb/msigdb_disabled.flag"
        ),
    output:
        ora_rds="resources/enrichment/{analysis}/{contrast}/ora_results.RDS",
    params:
        enrichment=get_enrichment_params,
    log:
        "logs/enrichment/{analysis}/{contrast}/ora.log",
    conda:
        "../envs/r_env.yaml"
    wildcard_constraints:
        analysis="[^/]+",
        contrast="[^/]+",
    script:
        "../scripts/run_ora.R"


# SPIA analysis rule is conditionally executed based on config setting
if config["enrichment"]["spia"]["enabled"]:

    rule spia:
        input:
            unpack(get_enrichment_deps),
            spia_data="resources/enrichment/spia_data",
        output:
            spia_rds="resources/enrichment/{analysis}/{contrast}/spia_results.RDS",
            spia_readable_rds="resources/enrichment/{analysis}/{contrast}/spia_results_readable.RDS",
        params:
            enrichment=get_enrichment_params,
        log:
            "logs/enrichment/{analysis}/{contrast}/spia.log",
        conda:
            "../envs/r_env.yaml"
        wildcard_constraints:
            analysis="[^/]+",
            contrast="[^/]+",
        script:
            "../scripts/run_spia.R"


# Create flag indicating if MSigDB is enabled
rule create_msigdb_flag:
    output:
        branch(
            config["enrichment"]["msigdb"]["enabled"],
            then="resources/msigdb/msigdb_enabled.flag",
            otherwise="resources/msigdb/msigdb_disabled.flag",
        ),
    log:
        "logs/enrichment/create_msigdb_flag.log",
    shell:
        """
        touch {output}
        """


# Pre-load MSigDB data to prevent race conditions during concurrent enrichment analysis
rule preload_msigdb:
    input:
        "resources/msigdb/msigdb_enabled.flag",
    output:
        "resources/msigdb/msigdb_loaded.flag",
    params:
        species=lambda wildcards: config["ref_assembly"]["species"],
        collections=lambda wildcards: config["enrichment"]["msigdb"]["collections"],
    log:
        "logs/enrichment/preload_msigdb.log",
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/preload_msigdb.R"


rule all_enrichment_analyses:
    input:
        get_enrichment_outputs(),
    output:
        touch("resources/enrichment/enrichment_analyses_complete.done"),
