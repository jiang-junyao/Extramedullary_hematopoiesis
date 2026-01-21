library(org.Hs.eg.db)
library(BuenColors)
library(Seurat)
library(tidyverse)
obj <- readRDS("~/wenqian_sp/Public_human_scrna/all_human_atlas.rds")
obj = subset(obj,tissue %in% c('BM','Spleen'))
obj = RunUMAP(obj,reduction = 'scvi',dims = 1:10,min.dist = 0.8)


UMAPPlot(obj,group.by='tissue',pt.size=0.6)+theme_void()+
  scale_color_manual(values = c(rgb(142/255,224/255,142/255),
                               rgb(250/255,180/255,128/255),
                               rgb(218/255,162/255,200/255),
                               rgb(255/255,170/255,170/255)))
ggsave('/data/jiangjunyao/wenqian_sp/Public_human_scrna/figure/humam_hsc_tissue.pdf',
       width = 10,height = 8)

UMAPPlot(obj,group.by='Batch',pt.size=0.6)+theme_void()+
  scale_color_manual(values = jdb_palette('corona'))
ggsave('/data/jiangjunyao/wenqian_sp/Public_human_scrna/figure/humam_hsc_Batch.pdf',
       width = 10,height = 8)

UMAPPlot(obj,group.by='stage',pt.size=0.6)+theme_void()+
  scale_color_manual(values = jdb_palette('brewer_spectra')[c(1,7)])
ggsave('/data/jiangjunyao/wenqian_sp/Public_human_scrna/figure/humam_hsc_stage.pdf',
       width = 10,height = 8)

obj = NormalizeData(obj)

adult = subset(obj,stage=='adult')
adult@active.ident = as.factor(adult$tissue)
adult_deg = FindAllMarkers(adult)
adult_deg = adult_deg[adult_deg$avg_log2FC>0,]
adult_deg2 = adult_deg[adult_deg$avg_log2FC>0.5,]
adult_deg2 = adult_deg2[adult_deg2$p_val_adj<0.05,]
for (i in unique(adult_deg2$cluster)) {
  tissue_go = gene_enrich(adult_deg2[adult_deg2$cluster==i,'gene'],
                          org.Hs.eg.db,'GO')
  tissue_go = setReadable(tissue_go,org.Hs.eg.db)
  tissue_go = tissue_go@result
  tissue_go = tissue_go[tissue_go$qvalue<=0.05,]
  write.csv(tissue_go,paste0('/data/jiangjunyao/wenqian_sp/Public_human_scrna/adult_hsc_2tissue/',
                             i,'_GO.csv'))
}

