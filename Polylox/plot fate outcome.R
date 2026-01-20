library(pheatmap)
library(tidyverse)
library(FateMapper)
plot_epiblast_heatmap <- function(ratio,group1,sample_name){
  anno_row = as.data.frame(group1[,2])
  rownames(anno_row) = rownames(group1)
  colnames(anno_row) = 'fatebias'
  heatmap_col = c('#F7FBFF','#DEEBF7','#C6DBEF','#9ECAE1','#6BAED6','#4292C6',
                  '#2171B5','#08519C','#08306B')
  fate_name = unique(group1$fate)
  if (nrow(ratio)==1) {
    print('only1')
    ### heatmap
    ratio[ratio==0] = -0.2
    pheatmap(ratio,cluster_rows = F,
             cluster_cols = F,show_rownames = F,
             annotation_row = anno_row,border_color = NA,
             color = colorRampPalette(heatmap_col)(50),
             filename = paste0('E:\\wenqian polylox\\20250423 new fate bias\\fig_all_barcode/',
                               sample_name,'.pdf'),width = 10,height = 8)
  }else{
    row_gap = c()
    for (i in 1:(length(fate_name)-1)) {
      if (i==1) {
        row_gap = c(row_gap,nrow(group1[group1$fate==fate_name[i],]))
      }else{
        row_gap = c(row_gap,nrow(group1[group1$fate==fate_name[i],])+row_gap[i-1])
      }
    }
    ### heatmap
    ratio[ratio==0] = -0.2
    pheatmap(ratio,cluster_rows = F,
             cluster_cols = F,show_rownames = F,
             gaps_row = row_gap,annotation_row = anno_row,border_color = NA,
             color = colorRampPalette(heatmap_col)(50),
             filename = paste0('E:\\wenqian polylox\\20250423 new fate bias\\fig_all_barcode/',
                               sample_name,'.pdf'),width = 10,height = 8)
  }
}

plot_resi <- function(ratio,group1,sample_name){
  anno_row = as.data.frame(group1[,2])
  rownames(anno_row) = rownames(group1)
  colnames(anno_row) = 'fatebias'
  heatmap_col = c('#F7FBFF','#DEEBF7','#C6DBEF','#9ECAE1','#6BAED6','#4292C6',
                  '#2171B5','#08519C','#08306B')
  fate_name = unique(group1$fate)
  if (nrow(ratio)==1) {
    print('only1')
    ### heatmap
    ratio[ratio==0] = -0.2
    pheatmap(ratio,cluster_rows = F,
             cluster_cols = F,show_rownames = F,
             annotation_row = anno_row,border_color = NA,
             color = colorRampPalette(heatmap_col)(50),
             filename = paste0('E:\\wenqian polylox\\20250423 new fate bias\\fig_tissue_resi/',
                               sample_name,'.pdf'),width = 10,height = 8)
  }else{
    row_gap = c()
    for (i in 1:(length(fate_name)-1)) {
      if (i==1) {
        row_gap = c(row_gap,nrow(group1[group1$fate==fate_name[i],]))
      }else{
        row_gap = c(row_gap,nrow(group1[group1$fate==fate_name[i],])+row_gap[i-1])
      }
    }
    ### heatmap
    ratio[ratio==0] = -0.2
    pheatmap(ratio,cluster_rows = F,
             cluster_cols = F,show_rownames = F,
             gaps_row = row_gap,annotation_row = anno_row,border_color = NA,
             color = colorRampPalette(heatmap_col)(50),
             filename = paste0('E:\\wenqian polylox\\20250423 new fate bias\\fig_tissue_resi/',
                               sample_name,'.pdf'),width = 10,height = 8)
  }
}
plot_prog_resi <- function(ratio,group1,sample_name){
  anno_row = as.data.frame(group1[,2])
  rownames(anno_row) = rownames(group1)
  colnames(anno_row) = 'fatebias'
  heatmap_col = c('#F7FBFF','#DEEBF7','#C6DBEF','#9ECAE1','#6BAED6','#4292C6',
                  '#2171B5','#08519C','#08306B')
  fate_name = unique(group1$fate)
  if (nrow(ratio)==1) {
    print('only1')
    ### heatmap
    ratio[ratio==0] = -0.2
    pheatmap(ratio,cluster_rows = F,
             cluster_cols = F,show_rownames = F,
             border_color = NA,
             color = colorRampPalette(heatmap_col)(50),
             filename = paste0('E:\\wenqian polylox\\20250423 new fate bias\\fig_prog_resi/',
                               sample_name,'.pdf'),width = 10,height = 8)
  }else{
    row_gap = c()
    for (i in 1:(length(fate_name)-1)) {
      if (i==1) {
        row_gap = c(row_gap,nrow(group1[group1$fate==fate_name[i],]))
      }else{
        row_gap = c(row_gap,nrow(group1[group1$fate==fate_name[i],])+row_gap[i-1])
      }
    }
    ### heatmap
    ratio[ratio==0] = -0.2
    pheatmap(ratio,cluster_rows = F,
             cluster_cols = F,show_rownames = F,
             border_color = NA,
             color = colorRampPalette(heatmap_col)(50),
             filename = paste0('E:\\wenqian polylox\\20250423 new fate bias\\fig_prog_resi/',
                               sample_name,'.pdf'),width = 10,height = 8)
  }
}
### plot all barcode
dir_all = dir('E:\\wenqian polylox\\20250423 new fate bias\\data')
#dir_all = dir_all[-grep('csv',dir_all)]
for (i in dir_all) {
  dir1 = paste0('E:\\wenqian polylox\\20250423 new fate bias\\data\\',i,'/')
  group1 = read.csv(paste0(dir1,'group_df2.csv'),row.names = 1)
  tissue_ratio = read.csv(paste0(dir1,'tissue_ratio.csv'),row.names = 1)


  resi = group1[group1$hsc_type!='share',]
  resi_ratio = tissue_ratio[rownames(resi),]

  share = group1[group1$hsc_type=='share',]
  share_ratio = tissue_ratio[rownames(share),]

  plot_epiblast_heatmap(tissue_ratio,group1,sample_name=paste0(i,'_all_barcode'))
  plot_epiblast_heatmap(resi_ratio,resi,sample_name=paste0(i,'_resi_barcode'))
  plot_epiblast_heatmap(share_ratio,share,sample_name=paste0(i,'_share_barcode'))
}


### plot single tissue resi barcode
dir_all = dir('E:\\wenqian polylox\\20250423 new fate bias\\data')
for (i in dir_all) {
  dir1 = paste0('E:\\wenqian polylox\\20250423 new fate bias\\data\\',i,'/')
  group1 = read.csv(paste0(dir1,'group_df2.csv'),row.names = 1)
  tissue_ratio = read.csv(paste0(dir1,'tissue_ratio.csv'),row.names = 1)

  resi_type = unique(group1$hsc_type)
  resi_type = resi_type[!resi_type %in% c('share','no_hsc')]
  for (j in resi_type) {
    resi = group1[group1$hsc_type==j,]
    resi_ratio = tissue_ratio[rownames(resi),]
    plot_resi(resi_ratio,resi,sample_name=paste0(i,'_',j))
  }

}
### plot all tissue progenitor resi/share
for (i in c('N22','N18','N16_1','N16_2','N19')) {
  dir1 = paste0('E:\\wenqian polylox\\20250423 new fate bias\\data\\',i,'/')
  group1 = read.csv(paste0(dir1,'group_df2.csv'),row.names = 1)
  tissue_ratio = read.csv(paste0(dir1,'tissue_count.csv'),row.names = 1)
  tissue_ratio = as.data.frame(t(apply(tissue_ratio, 1, function(x){
    row_sum = sum(x)
    return(x/row_sum)
  })))
  resi = rbind(group1[group1$hsc_type =='share',],
               group1[group1$hsc_type =="Legs.LSK-resi",],
               group1[group1$hsc_type =="Spine.LSK-resi",],
               group1[group1$hsc_type =="Spleen.LSK-resi",],
               group1[group1$hsc_type =="Liver.LSM-resi",])
  resi_ratio = tissue_ratio[rownames(resi),]
  if (i %in% c('N19','N30')) {
    plot_prog_resi(resi_ratio,resi,sample_name=paste0(i,'_only_prog_resi_share'))
  }else{
    plot_prog_resi(resi_ratio[,1:4],resi,sample_name=paste0(i,'_only_prog_resi_share'))
  }

}

### N16_1 special issue
n16 = read.csv('E:\\wenqian polylox\\20250423 new fate bias\\data\\N16_1/tissue_ratio.csv',row.names = 1)
group1 = read.csv('E:\\wenqian polylox\\20250423 new fate bias\\data\\N16_1/group_df2.csv',row.names = 1)
group1 = group1[group1$fate!='inactive',]
group1 = group1[group1$hsc_type!='share',]
col_idx = c(grep('Legs',colnames(n16)),
            grep('Spine',colnames(n16)),
            grep('Spleen',colnames(n16)),
            grep('Liver',colnames(n16)))
n16 = n16[,col_idx]

plot_use = rbind(n16[group1[group1$hsc_type=="Legs.LSK-resi",'barcodes'],],
                 n16[group1[group1$hsc_type=="Spine.LSK-resi",'barcodes'],],
                 n16[group1[group1$hsc_type=="Spleen.LSK-resi",'barcodes'],],
                 n16[group1[group1$hsc_type=="Liver.LSM-resi",'barcodes'],])

colSums(plot_use != 0)
df_list = list()
for (i in c("Legs.LSK-resi","Spine.LSK-resi","Spleen.LSK-resi","Liver.LSM-resi")) {
  barcode_use = group1[group1$hsc_type==i,'barcodes']
  sort_df = n16[barcode_use,grep(strsplit(i,'\\.')[[1]][1],colnames(n16))]
  df_list[[i]] =  as.data.frame(n16[barcode_use,][sort_clone_mt(sort_df),])
}
plot_use = do.call(bind_rows,df_list)
plot_use = plot_use[,-grep('CLP',colnames(plot_use))]
plot_use = plot_use[rowSums(plot_use)>0,]
plot_use[plot_use==0]= -0.6
pheatmap(plot_use[,c(1,4,2,3,5:ncol(plot_use))],cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(heatmap_col)(50),
         show_rownames = F,
         filename = 'E:\\wenqian polylox\\20250423 new fate bias\\fig_n16_resi_fate/n16_1.pdf')
### n16_2
n16 = read.csv('E:\\wenqian polylox\\20250423 new fate bias\\data\\N16_2/tissue_ratio.csv',row.names = 1)
group1 = read.csv('E:\\wenqian polylox\\20250423 new fate bias\\data\\N16_2/group_df2.csv',row.names = 1)
group1 = group1[group1$fate!='inactive',]
group1 = group1[group1$hsc_type!='share',]
group1$barcodes = rownames(group1)


col_idx = c(grep('Legs',colnames(n16)),
            grep('Spine',colnames(n16)),
            grep('Spleen',colnames(n16)),
            grep('Liver',colnames(n16)))
n16 = n16[,col_idx]

plot_use = rbind(n16[group1[group1$hsc_type=="Legs.LSK-resi",'barcodes'],],
                 n16[group1[group1$hsc_type=="Spine.LSK-resi",'barcodes'],],
                 n16[group1[group1$hsc_type=="Spleen.LSK-resi",'barcodes'],],
                 n16[group1[group1$hsc_type=="Liver.LSM-resi",'barcodes'],])
df_list = list()
for (i in c("Legs.LSK-resi","Spine.LSK-resi","Spleen.LSK-resi","Liver.LSM-resi")) {
  barcode_use = group1[group1$hsc_type==i,'barcodes']
  sort_df = n16[barcode_use,grep(strsplit(i,'\\.')[[1]][1],colnames(n16))]
  df_list[[i]] =  as.data.frame(n16[barcode_use,][sort_clone_mt(sort_df),])
}
plot_use = do.call(bind_rows,df_list)

plot_use = plot_use[,-grep('CLP',colnames(plot_use))]
plot_use = plot_use[rowSums(plot_use)>0,]
plot_use[plot_use==0]= -0.6
pheatmap(plot_use[,c(1,4,2,3,5:ncol(plot_use))],cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(heatmap_col)(50),
         show_rownames = F,
         filename = 'E:\\wenqian polylox\\20250423 new fate bias\\fig_n16_resi_fate/n16_2.pdf')
dir_all = dir('E:\\wenqian polylox\\20250423 new fate bias\\data')
#dir_all = dir_all[-grep('csv',dir_all)]
all_nor_pgen_mt = readRDS('E:\\wenqian polylox\\20250825_wenqian\\5.Data_Pgen_list\\qc2_4.rds')
for (i in dir_all) {
  dir1 = paste0('E:\\wenqian polylox\\20250423 new fate bias\\data\\',i,'/')
  group1 = read.csv(paste0(dir1,'group_df2.csv'),row.names = 1)
  tissue_ratio = read.csv(paste0(dir1,'tissue_ratio.csv'),row.names = 1)
  colnames(tissue_ratio) <- gsub("\\.", " ", colnames(tissue_ratio))

  data1 = all_nor_pgen_mt[[i]]
  data1 = data1[,colnames(tissue_ratio)]
  data1 = data1[rowSums(data1)>0,]

  other_barcode = rownames(data1)[!rownames(data1) %in% rownames(tissue_ratio)]
  data2 = data1[other_barcode,]
  data2 = as.data.frame(t(apply(data2, 1, function(x){
    row_sum = sum(x)
    return(x/row_sum)
  })))
  other_df = data.frame(
    rep('111',length(other_barcode)),
    rep('no progenitor',length(other_barcode)),
    rep('a',length(other_barcode)))
  rownames(other_df) = rownames(data2)
  final_ratio = rbind(tissue_ratio,data2)
  colnames(other_df) = colnames(group1)[1:3]
  final_df = rbind(group1[,1:3],other_df)

  plot_epiblast_heatmap(final_ratio,final_df,sample_name=paste0(i,'_all_barcode'))

}
heatmap_col = c('#F7FBFF','#DEEBF7','#C6DBEF','#9ECAE1','#6BAED6','#4292C6',
                '#2171B5','#08519C','#08306B')
for (i in c('N6','N9','N13','N31','N32')) {
  dir1 = paste0('E:\\wenqian polylox\\20250423 new fate bias\\data\\',i,'/')
  tissue_ratio = read.csv(paste0(dir1,'tissue_count.csv'),row.names = 1)
  tissue_ratio = tissue_ratio[,grep('HSC|MPP',colnames(tissue_ratio))]
  colnames(tissue_ratio) = gsub('\\.',' ',colnames(tissue_ratio))
  tissue_ratio = merge_to_progenitor(tissue_ratio)
  tissue_ratio = as.data.frame(t(apply(tissue_ratio, 1, function(x){
    row_sum = sum(x)
    return(x/row_sum)
  })))

  num_nonzero <- apply(tissue_ratio, 1, function(x) sum(x != 0))
  resi_ratio <- tissue_ratio[num_nonzero == 1, ]
  share_ratio <- tissue_ratio[num_nonzero != 1, ]
  resi_ratio = resi_ratio[sort_clone_mt(resi_ratio),]


  pheatmap::pheatmap(rbind(share_ratio,resi_ratio),
                     cluster_rows = F,cluster_cols = F,
                     color = colorRampPalette(heatmap_col)(50),
                     show_rownames = F,height = 8,width = 10,
                     filename = paste0('E:\\wenqian polylox\\20250423 new fate bias\\18. hsc mpp heatmap/',
                                       i,'.pdf'))
}

for (i in c('N6','N9','N13','N31','N32')) {
  dir1 = paste0('E:\\wenqian polylox\\20250423 new fate bias\\data\\',i,'/')
  tissue_ratio = read.csv(paste0(dir1,'tissue_count.csv'),row.names = 1)
  tissue_ratio = tissue_ratio[,grep('HSC|MPP',colnames(tissue_ratio))]
  colnames(tissue_ratio) = gsub('\\.',' ',colnames(tissue_ratio))
  tissue_ratio = merge_to_progenitor(tissue_ratio)
  tissue_ratio = as.data.frame(t(apply(tissue_ratio, 1, function(x){
    row_sum = sum(x)
    return(x/row_sum)
  })))

  num_nonzero <- apply(tissue_ratio, 1, function(x) sum(x != 0))
  resi_ratio <- tissue_ratio[num_nonzero == 1, ]
  share_ratio <- tissue_ratio[num_nonzero != 1, ]
  resi_ratio = resi_ratio[sort_clone_mt(resi_ratio),]


  pheatmap::pheatmap(rbind(share_ratio,resi_ratio),
                     cluster_rows = F,cluster_cols = F,
                     color = colorRampPalette(heatmap_col)(50),
                     show_rownames = F,height = 8,width = 10,
                     filename = paste0('E:\\wenqian polylox\\20250423 new fate bias\\18. hsc mpp heatmap/',
                                       i,'.pdf'))
}

for (i in c('N22','N18','N16_1','N16_2','N19','N30')) {
  dir1 = paste0('E:\\wenqian polylox\\20250423 new fate bias\\data\\',i,'/')
  tissue_ratio = read.csv(paste0(dir1,'tissue_count.csv'),row.names = 1)
  colnames(tissue_ratio) = gsub('\\.',' ',colnames(tissue_ratio))
  if (i %in% c('N19','N30')) {
    print('hhh')
  }else{
    tissue_ratio = tissue_ratio[,1:4]
  }
  tissue_ratio = tissue_ratio[rowSums(tissue_ratio)>0,]
  tissue_ratio = tissue_ratio[!is.na(tissue_ratio[,1]),]
  tissue_ratio = as.data.frame(t(apply(tissue_ratio, 1, function(x){
    row_sum = sum(x)
    return(x/row_sum)
  })))

  num_nonzero <- apply(tissue_ratio, 1, function(x) sum(x != 0))
  resi_ratio <- tissue_ratio[num_nonzero == 1, ]
  share_ratio <- tissue_ratio[num_nonzero != 1, ]
  resi_ratio = resi_ratio[sort_clone_mt(resi_ratio),]

  if (i %in% c('N19','N30')) {
    pheatmap::pheatmap(rbind(share_ratio,resi_ratio),
                       cluster_rows = F,cluster_cols = F,
                       color = colorRampPalette(heatmap_col)(50),
                       show_rownames = F,height = 8,width = 10,
                       filename = paste0('E:\\wenqian polylox\\20250423 new fate bias\\18. hsc mpp heatmap/',
                                         i,'.pdf'))
  }else{
    pheatmap::pheatmap(rbind(share_ratio,resi_ratio)[,1:4],
                       cluster_rows = F,cluster_cols = F,
                       color = colorRampPalette(heatmap_col)(50),
                       show_rownames = F,height = 8,width = 10,
                       filename = paste0('E:\\wenqian polylox\\20250423 new fate bias\\18. hsc mpp heatmap/',
                                         i,'.pdf'))
  }
}
