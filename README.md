# Snakemake workflow: `tucca-rna-seq`

# `THIS REPO IS STILL UNDER CONSTRUCTION AND DOES NOT REPRESENT A COMPLETED PIPELINE`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.27.1-brightgreen.svg)](https://snakemake.github.io)
[![Singularity](https://img.shields.io/badge/singularity-≥3.8.4-brightgreen.svg)](https://snakemake.github.io)
[![GitHub license](https://img.shields.io/github/license/tucca-cellag/tucca-rna-seq?color=orange)](https://github.com/tucca-cellag/tucca-rna-seq/blob/main/LICENSE)
[![tucca-rna-seq docs](https://img.shields.io/badge/documentation-tucca--rna--seq_docs-blue)](https://tucca-cellag.github.io/tucca-rna-seq/introduction)
[![GitHub actions status](https://github.com/tucca-cellag/tucca-rna-seq/workflows/Tests/badge.svg?branch=main)](https://github.com/tucca-cellag/tucca-rna-seq/actions?query=branch%3Amain+workflow%3ATests)

This workflow is the
[Tufts University Center for Cellular Agriculture's (TUCCA)](https://cellularagriculture.tufts.edu/) RNA-Seq Snakemake
Workflow for Cellular Agriculture Projects.

The usage of this workflow is described in our documentation at [TUCCA's Bioinformatics Docs for tucca-rna-seq](https://tucca-cellag.github.io/tucca-rna-seq/introduction).

This workflow is a standardized usage Snakemake workflow and can be found in
the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/tucca-cellag%20tucca-rna-seq.html).

If you use this workflow in a paper, don't forget to give credits to the authors
by citing the URL of this (original) repository and its DOI (above if
available).

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

## Connect With Us

We're here to help! If you have questions, feedback, or need assistance, feel
free to reach out through our social channels:

- [Tufts University Center for Cellular Agriculture](https://cellularagriculture.tufts.edu/)
- [TUCCA on LinkedIn](https://www.linkedin.com/company/tufts-cell-ag/)
- [TUCCA on X](https://twitter.com/tuftscellag)
- [TUCCA on YouTube](https://www.youtube.com/channel/UC29F8uqsu_K7aRxOgjfG_HQ)

Alternatively, visit our [GitHub Repository](https://github.com/tucca-cellag)
to explore more of TUCCA's computational projects.