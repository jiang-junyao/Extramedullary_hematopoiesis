# Existing RCTD helper, formerly Extended Data Figure 6g (now 5g).
# This file defines the original helper; it does not recreate absent final-panel plotting code.
# ref/querry are caller-supplied objects. doublet_mode remains 'doublet'.
# Only non-path correction: coordinates come from the supplied querry argument, not global lwh.
library(Seurat)
library(spacexr)

rctd_decon <- function(ref,querry,nCore = 4){
  # extract information to pass to the RCTD Reference function
  counts <- ref@assays$RNA@counts
  cluster <- as.factor(ref$ct)
  names(cluster) <- colnames(ref)
  nUMI <- ref$nCount_RNA
  names(nUMI) <- colnames(ref)
  reference <- Reference(counts, cluster, nUMI)

  # set up query with the RCTD function SpatialRNA
  counts <- querry@assays$Spatial@counts
  coords <- GetTissueCoordinates(querry)
  colnames(coords) <- c("x", "y")
  coords[is.na(colnames(coords))] <- NULL
  query <- SpatialRNA(coords, counts, colSums(counts))

  RCTD <- create.RCTD(query, reference, max_cores = nCore)
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  return(RCTD)
}
