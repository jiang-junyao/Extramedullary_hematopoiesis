# Niche single-cell RNA sequencing

Niche-cell UMAP plots and ligand–receptor communication for Extended Data Figure 7.

Requires Seurat and ggplot2; communication analysis additionally requires CellChat and future.

## Inputs

- data/junyao_processed_v4.rds, including organ, condition and celltype_general2.
- The celltypes_colors palette and jdb_palette() used by the UMAP script.
- obj1/obj2 and their celltype metadata for the CellChat comparison.

These processed objects and script-specific palette/helper dependencies are not bundled. Load the required objects/functions before running the corresponding script.

## Run

From the Niche_scRNA directory:

```r
source('Extended_Data_Figure_7.R', encoding = 'UTF-8')
source('Extended_Data_Figure_7_CellChat.R', encoding = 'UTF-8')
```

R/run_cellchat.R defines communication inference and is sourced by the entry script. Its database classes, thresholds and numeric source/target indices are specified in the code; ensure the supplied celltype identities match those indices.

Outputs are saved under results/Extended_Data_Figure_7/plots/. This module includes UMAP and communication scripts; it does not include the marker DotPlot implementation.
