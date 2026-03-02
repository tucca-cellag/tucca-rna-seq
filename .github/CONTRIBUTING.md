# Contributing to tucca-rna-seq

Thank you for your interest in contributing to the TUCCA RNA-Seq pipeline!
This document provides guidelines and information to help you get started.

## Code of Conduct

By participating in this project, you agree to abide by our
[Code of Conduct](CODE_OF_CONDUCT.md). Please read it before contributing.

## How to Contribute

### Reporting Bugs

If you find a bug, please
[open an issue](https://github.com/tucca-cellag/tucca-rna-seq/issues/new?template=bug_report.yml)
using the bug report template. Include:

- A clear description of the problem
- The command you ran and the terminal output
- Relevant log files (especially `.snakemake/log/`)
- Your system information (Snakemake version, OS, executor, container engine)

### Suggesting Features

Feature requests are welcome! Please
[open an issue](https://github.com/tucca-cellag/tucca-rna-seq/issues/new?template=feature_request.yml)
using the feature request template. Describe your use case, the proposed
solution, and any alternatives you've considered.

### Submitting Changes

1. **Fork** the repository
2. **Create a branch** from `main` using a descriptive name:
   - `feat/add-kallisto-support`
   - `fix/salmon-index-path`
   - `docs/update-configuration-guide`
3. **Make your changes** (see development guidelines below)
4. **Test your changes** (see testing section)
5. **Submit a pull request** against `main`

## Development Guidelines

### Branch Naming

Use prefixes that describe the type of change:

| Prefix    | Purpose                          |
|-----------|----------------------------------|
| `feat/`   | New features or enhancements     |
| `fix/`    | Bug fixes                        |
| `docs/`   | Documentation changes            |
| `refactor/` | Code refactoring               |
| `test/`   | Adding or updating tests         |
| `chore/`  | Maintenance tasks                |

### Commit Messages

This project uses
[Conventional Commits](https://www.conventionalcommits.org/). PR titles are
automatically linted to enforce this format. Examples:

```
feat(salmon): add support for selective alignment mode
fix(checksum): handle spaces in file paths
docs: update installation instructions
refactor(deseq2): simplify contrast generation logic
test: add integration test for Novogene MD5.txt format
chore: update conda environment versions
```

### Code Style

#### Snakemake files (`.smk`)

- Format with [snakefmt](https://github.com/snakemake/snakefmt) before
  committing
- Follow Snakemake
  [best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html)
- Pass `snakemake --lint` without warnings

#### Python

- Follow [PEP 8](https://peps.python.org/pep-0008/) conventions
- Code is linted with pylint via Super-Linter in CI

#### R scripts

- Follow [tidyverse style guide](https://style.tidyverse.org/) conventions
- Use explicit package namespaces (e.g., `dplyr::filter()`) in workflow
  scripts to avoid conflicts
- R environment is managed with [renv](https://rstudio.github.io/renv/) —
  update `renv.lock` if you add or change R dependencies

#### General

- Use the project's `.editorconfig` settings (UTF-8, LF line endings,
  spaces for indentation)
- Keep lines to 80 characters where practical

### Adding New Rules

When adding a new Snakemake rule:

1. Create a new `.smk` file in `workflow/rules/`
2. Add a corresponding conda environment in `workflow/envs/` if new software
   is required
3. Include the rule file in `workflow/Snakefile`
4. Add helper functions to `workflow/rules/common.smk` as needed
5. Update `config/config.yaml` with any new configuration parameters and
   document them in the schema (`workflow/schemas/`)
6. Add the tool's citation to `CITATIONS.md`

### Adding R Packages

1. Install the package into the renv environment
2. Run `renv::snapshot()` to update `renv.lock`
3. If the package is used in the Snakemake workflow (not just analysis
   notebooks), also add it to `workflow/envs/r_env.yaml`

## Testing

### Running Tests Locally

The test suite uses Snakemake's built-in dry-run and lint capabilities:

```bash
# Lint the workflow
snakemake --lint --configfile .test/local_reads/ensembl/config_complex/config.yaml

# Dry run with a test config
snakemake all --dry-run --configfile .test/local_reads/ensembl/config_complex/config.yaml
```

### Test Data

Test data is organized under `.test/`:

- `.test/data/` — Small FASTQ files and checksum files for local testing
- `.test/local_reads/` — Test configurations for local read scenarios
  (Ensembl and RefSeq genome sources, basic and complex configs)
- `.test/sra_reads/` — Test configurations for SRA download scenarios
  (yeast and mouse)

### CI Pipeline

All pull requests are automatically tested with:

- **Formatting**: Super-Linter (snakefmt, YAML, R, Python/pylint)
- **Linting**: `snakemake --lint` against stable and latest Snakemake versions
- **Integration tests**: Full workflow runs using Conda, Singularity, and
  Apptainer across multiple test scenarios
- **R notebook tests**: Rendering of analysis RMarkdown notebooks

CI must pass before a PR can be merged.

### Adding Test Cases

When contributing new features that affect the workflow:

1. Add test data to `.test/data/` if needed (keep files small)
2. Create or update test configurations in `.test/local_reads/` or
   `.test/sra_reads/`
3. If adding a new test scenario, add it to the CI matrix in
   `.github/workflows/main.yml`

## Release Process

Releases are managed automatically by
[Release Please](https://github.com/googleapis/release-please). When commits
following the Conventional Commits format are merged to `main`, Release Please
will automatically:

1. Create or update a release PR with the changelog
2. Upon merge of the release PR, create a GitHub release with the appropriate
   version bump

## Getting Help

- **Documentation**: [tucca-cellag.github.io/tucca-rna-seq](https://tucca-cellag.github.io/tucca-rna-seq/introduction)
- **Issues**: [GitHub Issues](https://github.com/tucca-cellag/tucca-rna-seq/issues)
- **Discussions**: Open an issue for questions or design discussions

## License

By contributing, you agree that your contributions will be licensed under the
project's [MIT License](../LICENSE).