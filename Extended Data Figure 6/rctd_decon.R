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
  coords <- GetTissueCoordinates(lwh)
  colnames(coords) <- c("x", "y")
  coords[is.na(colnames(coords))] <- NULL
  query <- SpatialRNA(coords, counts, colSums(counts))

  RCTD <- create.RCTD(query, reference, max_cores = nCore)
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  return(RCTD)
}
