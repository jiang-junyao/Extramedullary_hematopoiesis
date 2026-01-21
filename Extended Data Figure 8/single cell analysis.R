run_cellchat <- function(a1,metadata,group,species,
                         search = c('Cell-Cell Contact','ECM-Receptor',
                                    'Secreted Signaling')){
  library(CellChat)
  cellchat <- createCellChat(object = a1, meta = metadata, group.by = group)
  if (species =='mm') {
    db = CellChatDB.mouse
  }
  if (species =='hs') {
    db = CellChatDB.human
  }
  if (species =='zf') {
    db = CellChatDB.zebrafish
  }
  CellChatDB <- db
  CellChatDB.use <- subsetDB(CellChatDB, search = search)
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  future::plan("multicore", workers = 4)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  #cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  return(cellchat)
}
### condition1
cellchat1 = run_cellchat(obj1,obj1@meta.data,'celltype','mm')
### condition2
cellchat2 = run_cellchat(obj2,obj2@meta.data,'celltype','mm')

preg_list = list(cellchat1,cellchat2)
names(preg_list) = c('sp_preg','bm_preg')
cellchat_preg <- mergeCellChat(preg_list, add.names = names(preg_list))

netVisual_heatmap(cellchat_preg, measure = "weight",font.size = 10)

cellchat_preg <- identifyOverExpressedGenes(cellchat_preg, group.dataset = "datasets",
                                            pos.dataset = 'sp_preg',
                                            features.name = 'sp_preg.merged',
                                            only.pos = FALSE, thresh.pc = 0.1,
                                            thresh.fc = 0.05,thresh.p = 0.01,
                                            group.DE.combined = FALSE)
net <- netMappingDEG(cellchat_preg, features.name = 'sp_preg.merged', variable.all = TRUE)
net.up <- subsetCommunication(cellchat_preg, net = net, datasets = "sp_preg",ligand.logFC = 0.05, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat_preg, net = net, datasets = "bm_preg",ligand.logFC = -0.05, receptor.logFC = NULL)

### sources.use, targets.use改一下 代表不同细胞类型

netVisual_chord_gene(preg_list[[1]], sources.use = c(5,6,8,9), targets.use = c(5,6,8,9),
                     slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 50, title.name = paste0("Up-regulated signaling in ", names(preg_list)[1]))

netVisual_chord_gene(preg_list[[2]], sources.use = c(5,6,8), targets.use = c(5,6,8),
                     slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(preg_list)[1]))

