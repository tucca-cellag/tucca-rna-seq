# GMT Files Directory

This directory contains custom gene set files in GMT format for use with
clusterProfiler's universal enrichment analysis functions.

## File Format

GMT files should follow the standard format:

- Each line represents one gene set
- Fields are tab-separated
- First field: gene set name/ID
- Second field: gene set description (optional, can be empty)
- Remaining fields: gene identifiers (one per field)

Example:

```tsv
GENE_SET_1  Description of gene set 1 GENE1 GENE2 GENE3
GENE_SET_2  Description of gene set 2 GENE4 GENE5 GENE6 GENE7
```

## Available Files

- `muscle_genesets.gmt` - Muscle Gene Sets v3 (March 2019) containing 1,517
  gene sets, added to this workflow via download from [https://www.sys-myo.com/muscle_gene_sets/](https://www.sys-myo.com/muscle_gene_sets/) on 08/05/25

### Muscle Gene Sets Details

The Muscle Gene Sets (version 3, released March 2019) is comprised of 1,517
gene sets:

- **1,156 gene sets** derived from analysis of 302 studies of muscle physiology
  and disease published from 2005-present
- **122 gene sets** derived from published in vitro muscle microarray studies
  carried out from 2005-present
- **185 gene sets** derived from a previous meta-analysis carried out by Jelier
  et al. 2005
- **54 gene sets** from muscle-related gene ontology terms and muscle-relevant
  entries in the MSigDB database

The second column of the GMT file contains descriptions of each gene set. For
muscle gene sets derived from transcriptomic studies, PubMed references or Gene
Expression Omnibus series numbers are provided.

### Citation

If you used the Muscle Gene Sets in your research, please cite:

**Malatras A, Duguez S, Duddy W: Muscle Gene Sets: a versatile methodological**
**aid to functional genomics in the neuromuscular field. Skeletal Muscle 2019,**
**9(1):10.**

[Link to the publication](https://skeletalmusclejournal.biomedcentral.com/articles/10.1186/s13395-019-0196-z)

## Usage in Configuration

To use custom GMT files in your analysis, add them to the `msigdb` section of
your `config.yaml`:

```yaml
enrichment:
  msigdb:
    enabled: true
    collections:
      - "H"  # Hallmark gene sets
      - "C2" # Curated gene sets
    custom_gmt_files:
      - "resources/gmt/muscle_genesets.gmt"
```

## Adding Your Own Files

1. Place your GMT file in this directory
2. Reference it in your `config.yaml` under `enrichment.msigdb.custom_gmt_files`
3. The workflow will automatically load and use your custom gene sets

## Supported Gene ID Types

The workflow supports the following gene ID types for custom GMT files:

- Entrez IDs (recommended)
- Gene symbols
- Ensembl IDs

Make sure your gene IDs match the format used in your differential expression
analysis.
