# Existing CellChat comparison, formerly Extended Data Figure 8e (now 7e).
# obj1, obj2 and their original celltype metadata must already be available.
# Source/target numeric indices and all analytical parameters are deliberately retained.
source('R/run_cellchat.R')
dir.create('results/Extended_Data_Figure_7/plots', recursive = TRUE, showWarnings = FALSE)
pdf('results/Extended_Data_Figure_7/plots/CellChat_original_comparison.pdf', width = 10, height = 8)

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


dev.off()
