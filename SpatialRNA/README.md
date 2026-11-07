# Spatial transcriptomics

Spatial preprocessing, domain annotation and RCTD helper code associated with Extended Data Figure 5.

## Python preprocessing

Open spleen_preprocess.ipynb with the working directory set to SpatialRNA. Requires scanpy, omicverse, pandas, numpy, seaborn, matplotlib, scipy and the SpaceFlow/GPU dependencies.

Supply data/bin50_polygon_select.h5ad; it is not bundled. The notebook includes QC, spatially variable genes, SpaceFlow tuning, region annotation and matrix export. Its first cell downloads a pathway database. Organism settings, coordinate splits, region labels, tuning grid and seed 114514 are specified in the cells.

Generated tuning/transfer data are written under prepared/ and saved plots under results/Extended_Data_Figure_5/plots/.

## RCTD

Requires Seurat and spacexr. From the SpatialRNA directory:

```r
source('Extended_Data_Figure_5.R', encoding = 'UTF-8')
# Supply the reference and spatial Seurat objects before calling:
# rctd_decon(ref, querry)
```

The caller supplies ref/querry, including reference cell identities and spatial counts/coordinates. These objects are not bundled. The helper uses doublet_mode='doublet' and returns the RCTD object. It is not a complete final-panel plotting script.
