<div align="left">
  <img width="25%" align="left" src="images/tucca-rna-seq-logo-white.png" alt="tucca-rna-seq logo">
</div>

**`tucca-cellag/tucca-rna-seq`** is the
[Tufts University Center for Cellular Agriculture's (TUCCA)][1] RNA-Seq
workflow specifically designed for cellular agriculture projects.

This workflow was developed using the [Snakemake][1.5] workflow management
system. `tucca-rna-seq` is a standardized usage Snakemake workflow that follows
the [best practices][2.5] laid out in the Snakemake documentation (as of
Snakemake v9.3) and can be found in the [Snakemake Workflow Catalog][3].

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


> [!WARNING]
> This workflow is still under construction. [Release v0.9.0][8] marks our first
> public release. v0.9.0 contains all logic to process raw paired-end RNA-Seq
> reads through differential expression. Currently, the workflow can generate a
> large number of DESeq2 result files, especially for experiments with multiple
> conditions and contrasts. While this is thorough, we recognize that navigating
> dozens of individual result files can be challenging. The centerpiece of the
> v1.0.0 release will be an interactive analysis toolkit that allows you to
> dynamically explore and visualize your results. This will include a suite of
> Shiny applications leveraging powerful packages like [`pcaExplorer`][9],
> [`ideal`][10], and [`GeneTonic`][11] to bring your data to life, as well as
> custom scripting to generate our favorite [`clusterProfiler`][12] figures. We
> encourage users to test this v0.9.0 release and provide feedback. Users should
> expect our documentation to be incomplete and continue to have major reworks
> until v1.0.0 is released. Please [open an issue][5] to report any bugs or
> suggest improvements. Additionally, feel free to [contact us][contact] with
> any questions.

## Documentation

The usage of this workflow is described in our documentation at
[tucca-cellag.github.io][2]. If you've found a bug or there is a feature that
we're missing (in the workflow or in our documentation) please
[open an issue][5] to let us know.

## Workflow Overview

<div align="center">
  <img alt="tucca-rna-seq workflow map" src="images/tucca-rna-seq-workflow-no-logo.png" width="700">
  <p>Created in <a href="https://BioRender.com">https://BioRender.com</a></p>
</div>

## Rulegraph

<div align="center">
  <img alt="tucca-rna-seq workflow map" src="images/rulegraph.png" width="700">
  <p>Created via `snakemake --rulegraph`</a></p>
</div>

## How do I get help?

For questions or suggestions regarding the workflow, first, check out our
detailed documentation at [tucca-cellag.github.io][2]. If you can't find the
answer to your question in our documentation you can try checking if someone
has [previously opened an issue][4] answering your question. If you still have
a question please [open an issue][5] so we can help! For any other inquiries,
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
please consult the [contributing guidelines][6]. For all interactions within
the `tucca-cellag` community, we ask that you observe our [code of conduct][7].

&copy; 2025 [Tufts University Center for Cellular Agriculture][1]

[1]: https://cellularagriculture.tufts.edu/
[1.5]: https://snakemake.readthedocs.io/en/stable/
[2]: https://tucca-cellag.github.io/tucca-rna-seq/introduction
[2.5]: https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html
[3]: https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/tucca-cellag%20tucca-rna-seq.html
[4]: https://github.com/tucca-cellag/tucca-rna-seq/issues
[5]: https://github.com/tucca-cellag/tucca-rna-seq/issues/new/choose
[contact]: <mailto:benjamin.bromberg@tufts.edu>
[6]: .github/CONTRIBUTING.md
[7]: CODE_OF_CONDUCT.md
[8]: https://github.com/tucca-cellag/tucca-rna-seq/releases/tag/v0.9.0
[9]: https://bioconductor.org/packages/release/bioc/html/pcaExplorer.html
[10]: https://bioconductor.org/packages/release/bioc/html/ideal.html
[11]: https://bioconductor.org/packages/release/bioc/html/GeneTonic.html
[12]: https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html