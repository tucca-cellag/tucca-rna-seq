# Snakemake workflow: `tucca-rna-seq`

## `THIS WORKFLOW IS STILL UNDER CONSTRUCTION!`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.27.1-3EB049)](https://snakemake.github.io)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity_≥3.8.4-1d355c)](https://sylabs.io/docs/)
[![run with apptainer](https://img.shields.io/badge/run%20with-apptainer-1d355c)](https://apptainer.org/)
[![GitHub Actions](https://img.shields.io/badge/GitHub_Actions-2088FF?logo=github-actions&logoColor=white)](https://github.com/tucca-cellag/tucca-rna-seq/actions)
[![GitHub actions status](https://github.com/tucca-cellag/tucca-rna-seq/workflows/Tests/badge.svg?branch=main)](https://github.com/tucca-cellag/tucca-rna-seq/actions?query=branch%3Amain+workflow%3ATests)
[![GitHub license](https://img.shields.io/github/license/tucca-cellag/tucca-rna-seq?color=orange)](https://github.com/tucca-cellag/tucca-rna-seq/blob/main/LICENSE)
[![Docusaurus](https://img.shields.io/badge/tucca--cellag.github.io-3EB049?logo=docusaurus&label=Docusaurus)](https://tucca-cellag.github.io/tucca-rna-seq/introduction)
[![Issues](https://img.shields.io/github/issues/tucca-cellag/tucca-rna-seq?style=flat&label=issues&color=3EB049)](https://github.com/tucca-cellag/tucca-rna-seq/issues)
[![Open Issues](https://img.shields.io/badge/GitHub-Open%20Issue-blue?logo=github)](https://github.com/tucca-cellag/tucca-rna-seq/issues/new)
[![Learn about TUCCA](http://img.shields.io/badge/TUCCA-3172AE.svg?label=learn%20about)](https://cellularagriculture.tufts.edu/)
[![Watch on YouTube](http://img.shields.io/badge/TUCCA-FF0000?label=youtube&logo=youtube)](https://www.youtube.com/channel/UC29F8uqsu_K7aRxOgjfG_HQ)
[![Follow on Twitter](http://img.shields.io/badge/tuftscellag-1DA1F2?label=twitter&logo=x)](https://twitter.com/tuftscellag)
[![Connect on LindedIn](https://custom-icon-badges.demolab.com/badge/TUCCA-0077B5?label=LinkedIn&logo=linkedin-white&logoColor=fff)](https://www.linkedin.com/company/tufts-cell-ag/)

This workflow is the
[Tufts University Center for Cellular Agriculture's (TUCCA)][1] RNA-Seq
Snakemake Workflow for Cellular Agriculture Projects.

The usage of this workflow is described in our documentation at
[tucca-cellag.github.io][2].

This workflow is a standardized usage Snakemake workflow that follows the
[best practices][2.5] laid out in the Snakemake documentation (as of v9.3).
This workflow can also be found in the [Snakemake Workflow Catalog][3].

If you use this workflow in a paper, don't forget to give credits to the authors
by citing the URL of this (original) repository and its DOI (above if
available).

## Introduction

<h1>
  <picture>
    <img alt="tucca-rna-seq workflow map" src="images/tucca-rna-seq-workflow.png" width="900">
  </picture>
</h1>

## About the Workflow

The **`tucca-rna-seq`** workflow is designed to provide a seamless and efficient
pipeline for RNA-Seq data analysis, tailored specifically for cellular
agriculture applications. Here's what makes our workflow stand out ⭐:

- **Cell Ag-Specfic Analysis Modes** 🥩🍔:
  - 🚧 More Info Coming Soon! 🚧

- **Automated with Snakemake** 🐍: Utilizes `Snakemake`, a Python-based workflow
management system, to create readable and maintainable pipelines that simplify
complex bioinformatics tasks.
  
- **Comprehensive Data Processing** 📂:
  - **Quality Control** 📋: Implements `FastQC` and `Qualimap` for quality
    assessment.
  - **Salmon for Quantification** 🐟: Employs `Salmon` for fast, accurate
    transcript quantification, while taking into account experimental attributes
    and biases commonly observed in RNA-Seq data.
  - **Meta-Analysis** 📊: Aggregates these results using `MultiQC` to provide a
    unified overview of your data quality and processing metrics.
  
- **Differential Gene Expression Analysis** 🧬:
  - **Robust Statistical Tools** ✖️➗: Leverages `DESeq2` for differential
    expression analysis, ensuring reliable and statistically sound results.
  - **Parallel Processing** ⚙️⏱️: Employs `BiocParallel` and other R
    parallelization packages to efficiently parallelize differential gene
    expression analyses across multiple `DESeq2` contrasts simultaneously,
    significantly reducing computation time.
  - **Pathway and Enrichment Analysis** 🧩: Integrates `ClusterProfiler`, `GO`,
    `KEGG`, `msigdbr`, and `SPIA` to facilitate comprehensive pathway and
    functional enrichment analyses. These analyses include:
    - **Over-Representation Analysis (ORA):** Identifies pathways or gene sets
      that are over-represented in your differentially expressed genes compared
      to a background set.
    - **Functional Class Scoring Analyses (e.g., Gene Set Enrichment Analysis
      [GSEA]):** Assesses whether predefined sets of genes show statistically
      significant differences between two biological states.
    - **Topology-Based Analyses (e.g., SPIA):** Incorporates pathway topology
      information to evaluate the impact of gene expression changes on specific
      biological pathways.
  - **Visualization** 📸: Utilizes `ggplot2`, `EnhancedVolcano`, `pheatmap`,
    `ClusterProfiler`, and many other visualization tools to create insightful
    and publication-ready figures.
  
- **High Reproducibility** 🔄:
  - **Environment Management** 🔧🌐: Employs `singularity` and `renv` to manage and
    replicate computational environments, ensuring consistency across different
    systems and projects.

- **Scalability and Flexibility** 📈🎛️: Designed to handle datasets of varying
  sizes and complexities, making it suitable for both small-scale studies and
  large, high-throughput projects.

By integrating these powerful tools into a cohesive workflow,
**`tucca-rna-seq`** provides a reliable and efficient platform for your RNA-Seq
data analysis needs, allowing you to focus on deriving meaningful biological
insights 🧠 without getting bogged down by technical complexities ⚙️.

## How do I get help?

We're here to help! Please check out our detailed documentation at
[tucca-cellag.github.io][2] and check out [previous issues opened][4] by other
users. If you have additional believe you've found a bug, or there is a feature
that we're missing, lease feel free to [open an issue][5]. You can always reach
out to us through our social channels.

## What is Cellular Agriculture? 🧬🌱

**Cellular Agriculture** is a cutting-edge field that harnesses biotechnology
to produce agricultural products directly from cells. Unlike traditional
farming, which relies on raising and harvesting whole organisms, cellular
agriculture focuses on cultivating animal cells in controlled environments to
create sustainable alternatives for meat, dairy, and other animal-derived
products.

### **Why Cellular Agriculture Matters**

- **Sustainability 🌍:** Reduces the environmental impact associated with
  conventional agriculture, including lower greenhouse gas emissions, reduced
  land and water usage, and minimized waste production.
- **Ethical Considerations 🐮❤️:** Offers humane alternatives by eliminating the
  need for animal slaughter, addressing animal welfare concerns.
- **Food Security 🍽️:** Enhances the ability to produce food in areas with
  limited agricultural resources, contributing to global food security.
- **Innovation and Research 🔬:** Drives advancements in biotechnology,
  genetics, and bioinformatics, fostering interdisciplinary collaboration and
  novel scientific discoveries.

[1]: https://cellularagriculture.tufts.edu/
[2]: https://tucca-cellag.github.io/tucca-rna-seq/introduction
[2.5]: https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html
[3]: https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/tucca-cellag%20tucca-rna-seq.html
[4]: https://github.com/tucca-cellag/tucca-rna-seq/issues
[5]: https://github.com/tucca-cellag/tucca-rna-seq/issues/new