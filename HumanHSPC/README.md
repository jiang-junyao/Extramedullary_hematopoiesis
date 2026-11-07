# Human HSPCs

Adult bone marrow versus spleen analysis associated with Figure 3c.

Requires Seurat, BuenColors, org.Hs.eg.db, tidyverse and the enrichment packages used by gene_enrich()/setReadable().

## Input and run

Supply data/all_human_atlas.rds with tissue, stage, Batch and scvi reduction. The processed atlas and gene_enrich() helper are not bundled; load the helper in the R session before running.

From the HumanHSPC directory:

```r
source('Figure_3c.R', encoding = 'UTF-8')
```

The script uses the adult tissue comparison, scVI-based UMAP and FindAllMarkers, followed by avg_log2FC>0.5, adjusted P<0.05 and GO qvalue<=0.05 selection. These are not the mouse BulkRNA thresholds.

UMAPs and GO tables are written to results/Figure_3c/plots/ and tables/. This module contains the upstream analysis, not a separate final GO-barplot implementation.
