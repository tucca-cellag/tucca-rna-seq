# tucca-rna-seq/workflow/rules/renv.smk
#
# Purpose:
# This Snakefile contains a rule to integrate the R package manager `renv`
# with the Conda-based environment of the main Snakemake pipeline. This hybrid
# approach uses Conda to ensure pipeline reproducibility while using a
# manually-curated `renv.lock` file to create a portable, user-friendly
# interactive R environment for downstream analysis.
#
#    `renv_restore`: This rule uses the manually-managed `renv.lock` file.
#    It executes `renv::restore()`, which installs all the specified packages
#    into a project-local library at `renv/library/`.
#
# By including this rule in the main `Snakefile`'s `all` rule, the setup of
# the interactive R environment is fully automated. After a successful pipeline
# run, a user can immediately open RStudio or an R session, and `renv` will
# automatically provide access to all the necessary packages without any manual
# setup.


rule renv_restore:
    input:
        "resources/enrichment/enrichment_analyses_complete.done",
    output:
        touch(".renv_restore.done"),
    log:
        "logs/renv/restore.log",
    conda:
        "../envs/renv_restore_env.yaml"
    script:
        "../scripts/renv_restore.R"
