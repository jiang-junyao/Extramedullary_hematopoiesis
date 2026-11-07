# Extramedullary hematopoiesis

Analysis code and input tables for bulk RNA sequencing, Polylox clonal tracing, human HSPCs, niche single-cell RNA sequencing and spatial transcriptomics.

## Analyses

| Directory | Contents |
|---|---|
| [BulkRNA](BulkRNA/README.md) | Differential expression, enrichment, GSEA and gene-expression heatmaps |
| [Polylox](Polylox/README.md) | Barcode tissue distribution, sharing, sampling controls, lineage output and split-sample analyses |
| [HumanHSPC](HumanHSPC/README.md) | Human bone marrow versus spleen HSPC analysis |
| [Niche_scRNA](Niche_scRNA/README.md) | Niche-cell UMAPs and ligand–receptor communication |
| [SpatialRNA](SpatialRNA/README.md) | Spatial preprocessing, domain annotation and RCTD deconvolution |

## Getting started

Run each module from its own working directory. Shared preparation is run once, followed by the required plotting script.

```r
# Start from the repository root in a fresh R session.
setwd('BulkRNA')
source('00_prepare_shared_data.R', encoding = 'UTF-8')
source('Extended_Data_Figure_8.R', encoding = 'UTF-8')
source('Extended_Data_Figure_10.R', encoding = 'UTF-8')
```

For Polylox, start a separate session from the repository root:

```r
setwd('Polylox')
source('00_prepare_shared_data.R', encoding = 'UTF-8')
source('Figure_1.R', encoding = 'UTF-8')
```

The module READMEs describe additional entry points, dependencies, input formats and methods. Do not source a module's scripts directly from the repository root: relative paths are resolved from the module directory.

## Software and data

BulkRNA uses R packages including DESeq2, clusterProfiler and fgsea; Polylox uses dplyr, tidyr, ggplot2 and scales. Human/niche analyses require Seurat and the listed enrichment/communication packages. Spatial preprocessing uses Python; RCTD uses R. Package requirements and available version information are documented per module.

BulkRNA count/annotation tables and Polylox count/Pgen/metadata tables are included under their data/ directories. Large processed human, niche and spatial objects, and certain script-specific external dependencies, are not bundled; see the corresponding input descriptions before running those modules.

Scripts create prepared/ and results/ within each module. These generated files are not included in the repository. Scope distinctions between count-based analyses, source-data plotting and reconstructed examples are described in the relevant README.

## License

See [LICENSE](LICENSE).
