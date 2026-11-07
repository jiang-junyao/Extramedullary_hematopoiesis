# Existing niche UMAP plots, formerly Extended Data Figure 8b (now 7b).
# celltype_general2, celltypes_colors and jdb_palette() are original-session dependencies.
# Marker DotPlot code was absent from the repository and is not added.
library(Seurat)
library(ggplot2)
seo <- readRDS('data/junyao_processed_v4.rds')
seo <- subset(seo, organ != 'Negative')
dir.create('results/Extended_Data_Figure_7/plots', recursive = TRUE, showWarnings = FALSE)

UMAPPlot(seo,group.by='celltype_general2')+scale_color_manual(values=celltypes_colors)+
  theme_void()+theme(text = element_text(size=16))
ggsave('results/Extended_Data_Figure_7/plots/niche_celltype.pdf',width = 10,height = 8)

UMAPPlot(seo,group.by='celltype_general2')+scale_color_manual(values=celltypes_colors)+
    theme_void()+theme(text = element_text(size=16))+NoLegend()
ggsave('results/Extended_Data_Figure_7/plots/niche_celltype_nolenged.pdf',width = 10,height = 8)
UMAPPlot(seo,group.by='organ')+scale_color_manual(values=jdb_palette('corona'))+
    theme_void()+theme(text = element_text(size=16))+NoLegend()
ggsave('results/Extended_Data_Figure_7/plots/niche_organ_nolenged.pdf',width = 10,height = 8)
UMAPPlot(seo,group.by='organ')+scale_color_manual(values=jdb_palette('corona'))+
    theme_void()+theme(text = element_text(size=16))
ggsave('results/Extended_Data_Figure_7/plots/niche_organ.pdf',width = 10,height = 8)
UMAPPlot(seo,group.by='condition')+
    theme_void()+theme(text = element_text(size=16))+NoLegend()
ggsave('results/Extended_Data_Figure_7/plots/niche_condition_nolenged.pdf',width = 10,height = 8)
UMAPPlot(seo,group.by='condition')+
    theme_void()+theme(text = element_text(size=16))
ggsave('results/Extended_Data_Figure_7/plots/niche_condition.pdf',width = 10,height = 8)
