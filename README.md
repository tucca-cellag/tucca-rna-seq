<div align="left">
  <img width="25%" align="left" src="images/tucca-rna-seq-logo-white.png" alt="tucca-rna-seq logo">
</div>

**`tucca-cellag/tucca-rna-seq`** is a modular RNA-Seq workflow developed in the
[Kaplan Lab at TUCCA][tucca]. It provides a complete, end-to-end pipeline that
processes raw FASTQ files through quality control, quantification, differential
expression, and pathway enrichment analysis.

Initially tailored for cellular agriculture research on muscle and fat cell
transcriptomes, the workflow has been designed with modularity and flexibility
at its core. This allows users to easily customize the pipeline by adding,
removing, or modifying steps to suit a broad range of RNA-Seq applications
and experimental designs.

This workflow was developed using the [Snakemake][snakemake] workflow management
system. `tucca-rna-seq` is a standardized usage Snakemake workflow that follows
the [best practices][best-prac] laid out in the Snakemake documentation (as of
Snakemake v9.3) and can be found in the [Snakemake Workflow Catalog][swc].

---

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.27.1-3EB049)](https://snakemake.github.io)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity_≥3.8.4-1d355c)](https://sylabs.io/docs/)
[![run with apptainer](https://img.shields.io/badge/run%20with-apptainer-1d355c)](https://apptainer.org/)
[![GitHub Actions](https://img.shields.io/badge/GitHub_Actions-2088FF?logo=github-actions&logoColor=white)](https://github.com/tucca-cellag/tucca-rna-seq/actions)
[![GitHub actions status](https://github.com/tucca-cellag/tucca-rna-seq/workflows/Tests/badge.svg)](https://github.com/tucca-cellag/tucca-rna-seq/actions/workflows/main.yml?query=branch%3Amain%20workflow%3ATests)
[![GitHub license](https://img.shields.io/github/license/tucca-cellag/tucca-rna-seq?color=orange)](https://github.com/tucca-cellag/tucca-rna-seq/blob/main/LICENSE)
[![Docusaurus](https://img.shields.io/badge/tucca--cellag.github.io-3EB049?logo=docusaurus&label=Docusaurus)](https://tucca-cellag.github.io/tucca-rna-seq/introduction)
[![Issues](https://img.shields.io/github/issues/tucca-cellag/tucca-rna-seq?style=flat&label=issues&color=3EB049)](https://github.com/tucca-cellag/tucca-rna-seq/issues)
[![Open Issues](https://img.shields.io/badge/GitHub-Open%20Issue-blue?logo=github)](https://github.com/tucca-cellag/tucca-rna-seq/issues/new)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15605826.svg)](https://doi.org/10.5281/zenodo.15605826)
<!-- ![GitHub Downloads (all assets, latest release)](https://img.shields.io/github/downloads/tucca-cellag/tucca-rna-seq/total) -->
[![Learn about TUCCA](http://img.shields.io/badge/TUCCA-3172AE.svg?label=learn%20about)](https://cellularagriculture.tufts.edu/)
[![Watch on YouTube](http://img.shields.io/badge/TUCCA-FF0000?label=youtube&logo=youtube)](https://www.youtube.com/channel/UC29F8uqsu_K7aRxOgjfG_HQ)
[![Follow on Twitter](http://img.shields.io/badge/tuftscellag-1DA1F2?label=twitter&logo=x)](https://twitter.com/tuftscellag)
[![Connect on LinkedIn](https://custom-icon-badges.demolab.com/badge/TUCCA-0077B5?label=LinkedIn&logo=linkedin-white&logoColor=fff)](https://www.linkedin.com/company/tufts-cell-ag/)

> **Looking Ahead: AI-Powered Analysis Integration**
>
> Beyond v1.0.0, we're developing an innovative AI-powered analysis platform
> that will integrate directly with this workflow. This platform will feature a
> private Large Language Model (LLM) chat application connected to your RNA-Seq
> results, enabling natural language queries about your data. The system will
> automatically vectorize workflow outputs and provide intelligent,
> context-aware analysis assistance.

## Documentation

The usage of this workflow is described in our documentation at
[tucca-cellag.github.io][our-docs]. If you've found a bug or there is a feature that
we're missing (in the workflow or in our documentation) please
[open an issue][open-issue] to let us know.

### Getting Started

Before running the workflow, we recommend reviewing our
[Data Collection Template][data-collection]
to ensure you have all necessary information organized. This template helps you:

- Organize your raw sequencing data
- Document sample metadata and experimental design
- Specify analysis parameters and preferences
- Plan for quality control and visualization needs

For detailed workflow documentation, visit [tucca-cellag.github.io][our-docs].

## Workflow Overview

<div align="center">
  <img alt="tucca-rna-seq workflow map" src="images/tucca-rna-seq-workflow.png" width="700">
  <p>Created in <a href="https://BioRender.com">https://BioRender.com</a></p>
</div>

## Rulegraph

<div align="center">
  <img alt="tucca-rna-seq workflow map" src="images/rulegraph.png" width="700">
  <p>Created via `snakemake --rulegraph`</p>
</div>

## How do I get help?

For questions or suggestions regarding the workflow, first, check out our
detailed documentation at [tucca-cellag.github.io][our-docs]. If you can't find the
answer to your question in our documentation you can try checking if someone
has [previously opened an issue][see-issues] answering your question. If you still have
a question please [open an issue][open-issue] so we can help! For any other inquiries,
please contact us via [email][contact].

## Citing the Workflow

If you use this workflow in a paper, don't forget to give credits to the
authors by citing the URL of this (original) repository and its DOI (above if
available). You should also cite the individual tools used in this workflow.
See [CITATIONS.md](CITATIONS.md) for details.

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

## Contributing

We welcome your involvement in the development of this workflow via submission
of bug reports, proposing new features, engaging in discussions, or providing
fixes and other code modifications. If you're interested in contributing,
please consult the [contributing guidelines][contrib]. For all interactions within
the `tucca-cellag` community, we ask that you observe our [code of conduct][conduct].

&copy; 2025 [Tufts University Center for Cellular Agriculture][tucca]

[tucca]: https://cellularagriculture.tufts.edu/
[snakemake]: https://snakemake.readthedocs.io/en/stable/
[our-docs]: https://tucca-cellag.github.io/tucca-rna-seq/introduction
[best-prac]: https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html
[swc]: https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/tucca-cellag/tucca-rna-seq.html
[see-issues]: https://github.com/tucca-cellag/tucca-rna-seq/issues
[open-issue]: https://github.com/tucca-cellag/tucca-rna-seq/issues/new/choose
[data-collection]: https://tucca-cellag.github.io/tucca-rna-seq/data-collection/data-collection-why
[contact]: <mailto:benjamin.bromberg@tufts.edu>
[contrib]: .github/CONTRIBUTING.md
[conduct]: CODE_OF_conduct.md
[v0.9.0]: https://github.com/tucca-cellag/tucca-rna-seq/releases/tag/v0.9.0
[pcaExplorer]: https://bioconductor.org/packages/release/bioc/html/pcaExplorer.html
[ideal]: https://bioconductor.org/packages/release/bioc/html/ideal.html
[GeneTonic]: https://bioconductor.org/packages/release/bioc/html/GeneTonic.html
[clusterProfiler]: https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
