# Harmonizome Integration

This directory contains cached Harmonizome gene sets downloaded from the [Harmonizome database](https://maayanlab.cloud/Harmonizome/).

## Overview

The Harmonizome integration allows you to perform enrichment analysis against tissue-specific and cell type-specific gene sets from various expression datasets. This provides a data-driven approach to understanding tissue-specific gene expression patterns in your differential expression results.

## Available Datasets

The workflow supports any dataset available in Harmonizome, including:

### GTEx Tissue Gene Expression Profiles

- **Dataset**: [GTEx Tissue Gene Expression Profiles](https://maayanlab.cloud/Harmonizome/dataset/GTEx+Tissue+Gene+Expression+Profiles)
- **Description**: Tissue-specific gene expression from the GTEx project
- **Available gene sets**: muscle, adipose, heart, brain, liver, lung, skin, kidney, thyroid, etc.

### BioGPS Human Cell Type and Tissue Gene Expression Profiles

- **Dataset**: [BioGPS Human Cell Type and Tissue Gene Expression Profiles](https://maayanlab.cloud/Harmonizome/dataset/BioGPS+Human+Cell+Type+and+Tissue+Gene+Expression+Profiles)
- **Description**: Cell type-specific gene expression profiles
- **Available gene sets**: adipocyte, myocyte, neuron, etc.

## Configuration

To use Harmonizome gene sets, configure them in your `config.yaml`:

```yaml
enrichment:
  harmonizome:
    enabled: true
    datasets:
      # GTEx Tissue Gene Expression Profiles
      - name: "GTEx Tissue Gene Expression Profiles"
        gene_sets:
          - "muscle"      # Skeletal muscle tissue
          - "adipose"     # Adipose tissue
          - "heart"       # Heart tissue
      
      # BioGPS Human Cell Type and Tissue Gene Expression Profiles
      - name: "BioGPS Human Cell Type and Tissue Gene Expression Profiles"
        gene_sets:
          - "adipocyte"   # Adipocyte cell type
          - "myocyte"     # Myocyte cell type
    
    # Parameters for ORA analysis
    ora:
      extra: "pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500"
    
    # Parameters for GSEA analysis
    gsea:
      extra: "pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500"
```

## File Format

Gene sets are cached as RDS files with the following naming convention:
```
{gene_set_name}_{dataset_name}.rds
```

For example:
- `muscle_GTEx_Tissue_Gene_Expression_Profiles.rds`
- `adipocyte_BioGPS_Human_Cell_Type_and_Tissue_Gene_Expression_Profiles.rds`

## Usage Examples

### Testing Muscle-Specific Expression

If you're studying muscle-related phenotypes, you might want to test against muscle-specific gene sets:

```yaml
enrichment:
  harmonizome:
    enabled: true
    datasets:
      - name: "GTEx Tissue Gene Expression Profiles"
        gene_sets:
          - "muscle"      # Skeletal muscle
      - name: "BioGPS Human Cell Type and Tissue Gene Expression Profiles"
        gene_sets:
          - "myocyte"     # Muscle cell type
```

### Testing Adipose-Related Expression

For studies involving adipogenesis or adipose tissue:

```yaml
enrichment:
  harmonizome:
    enabled: true
    datasets:
      - name: "GTEx Tissue Gene Expression Profiles"
        gene_sets:
          - "adipose"     # Adipose tissue
      - name: "BioGPS Human Cell Type and Tissue Gene Expression Profiles"
        gene_sets:
          - "adipocyte"   # Adipocyte cell type
```

## Biological Interpretation

### Tissue-Specific Enrichment

When you find enrichment in tissue-specific gene sets, this suggests:

1. **Tissue Identity**: Your samples may have characteristics of that tissue type
2. **Tissue-Specific Processes**: The biological processes you're studying may be particularly active in that tissue
3. **Cell Type Composition**: Your samples may contain cells of that tissue type

### Comparison with Other Gene Sets

- **vs. Hallmark gene sets**: Tissue-specific sets test for constitutive expression, while Hallmarks test for specific biological processes
- **vs. GO terms**: Tissue-specific sets are data-driven, while GO terms are knowledge-based
- **vs. KEGG pathways**: Tissue-specific sets focus on expression patterns, while KEGG focuses on metabolic/signaling pathways

## Technical Details

### API Access

The workflow uses the official Harmonizome Python API via `reticulate`:
- **API Class**: Downloads and uses the official `harmonizomeapi.py` class
- **Python Integration**: Uses `reticulate` to interface with the Python API
- **Official Methods**: Leverages the official `Harmonizome.get()` method for data access
- **API documentation**: Available at the Harmonizome website

### Caching

Gene sets are automatically cached to avoid repeated downloads:
- **Cache location**: `resources/harmonizome/`
- **Cache format**: RDS files containing gene set data and metadata
- **Cache invalidation**: Manual deletion required to force re-download

### Gene ID Conversion

The workflow automatically converts gene identifiers to match your analysis:
- **Input format**: Gene symbols from Harmonizome
- **Output format**: Entrez IDs for clusterProfiler
- **Conversion method**: Uses OrgDb packages (e.g., `org.Hs.eg.db`)

## Troubleshooting

### Common Issues

1. **Gene set not found**: Check the exact spelling and case of gene set names
2. **Download failures**: Verify internet connectivity and Harmonizome API availability
3. **ID conversion errors**: Ensure your OrgDb package supports the target species

### Debugging

Enable verbose logging by checking the workflow logs:
```bash
snakemake --configfile config/config.yaml --cores 1 --verbose
```

## Citation

When using Harmonizome data in your research, please cite:

**Diamant I, Clarke DJB, Evangelista JE, Lingam N, Ma'ayan A. Harmonizome 3.0: integrated knowledge about genes and proteins from diverse multi-omics resources. Nucleic Acids Res. 2024 Nov 20. pii: 53(1):D1016-D1028.**

**Rouillard AD, Gundersen GW, Fernandez NF, Wang Z, Monteiro CD, McDermott MG, Ma'ayan A. The harmonizome: a collection of processed datasets gathered to serve and mine knowledge about genes and proteins. Database (Oxford). 2016 Jul 3;2016. pii: baw100.**

For specific datasets, also cite the original publications as listed above. 