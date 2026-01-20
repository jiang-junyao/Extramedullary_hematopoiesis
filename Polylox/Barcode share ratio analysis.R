### share barcode/tissue barcode
library(tidyverse)
# dir_all = dir('E:\\wenqian polylox\\20250423 new fate bias\\data')
# list_all = list()
# for (i in dir_all) {
#   dir1 = paste0('E:\\wenqian polylox\\20250423 new fate bias\\data\\',i,'/')
#   group1 = read.csv(paste0(dir1,'group_df2.csv'),row.names = 1)
#
#   all_resi = unique(group1$hsc_type)
#   all_resi = all_resi[!all_resi %in% c('share','no_hsc')]
#   share_num = nrow(group1[group1$hsc_type=='share',])
#
#   if (share_num!=0) {
#     all_num = c()
#     all_tissue = c()
#     for (j in all_resi) {
#       resi_num = nrow(group1[group1$hsc_type==j,])
#       all_num = c(all_num,share_num/(share_num+resi_num))
#       j = gsub('HSC-resi','',j)
#       j = gsub('LSK-resi','',j)
#       j = gsub('LSM-resi','',j)
#       all_tissue = c(all_tissue,j)
#     }
#     df1=data.frame(all_num,all_tissue)
#     colnames(df1) = c('share_ratio','tissue')
#     df1$sample = i
#     list_all[[i]] = df1
#   }
# }
# df_all = do.call(bind_rows,list_all)
# write.csv(df_all,'E:\\wenqian polylox\\20250423 new fate bias\\result/all_sample_share_ratio.csv')


## share barcode/tissue barcode
all_nor_pgen_mt = readRDS('E:\\wenqian polylox\\20251231 barcode\\qc2_4.rds')
library(tidyverse)
dir_all = dir('E:\\wenqian polylox\\20250423 new fate bias\\data')
list_all = list()
for (i in dir_all) {
  dir1 = paste0('E:\\wenqian polylox\\20250423 new fate bias\\data\\',i,'/')
  group1 = read.csv(paste0(dir1,'group_df2.csv'),row.names = 1)
  all_df = all_nor_pgen_mt[[i]]
  all_resi = unique(group1$hsc_type)
  all_resi = all_resi[!all_resi %in% c('share','no_hsc')]

  all_num = c()
  all_count = c()
  all_tissue = c()
  all_barcode_num = c()
  all_resi_num = c()
  for (j in all_resi) {

    resi_num = nrow(group1[group1$hsc_type==j,])
    ###
    resi_barcode = rownames(group1[group1$hsc_type==j,])
    
    
    j = gsub('-resi','',j)
    j = gsub('\\.',' ',j)
    barcode_num = all_df[,j]
    barcode_num = barcode_num[barcode_num>0]
    
    
    all_count = c(all_count,sum(all_df[resi_barcode,j])/sum(all_df[,j]))
    all_num = c(all_num,resi_num/length(barcode_num))
    all_tissue = c(all_tissue,j)
    all_resi_num = c(all_resi_num,resi_num)
    all_barcode_num = c(all_barcode_num,length(barcode_num))
  }
  df1=data.frame(all_num,all_count,all_tissue,all_resi_num,all_barcode_num )
  colnames(df1) = c('resi_ratio','resi_ratio_count','tissue','resi barcode num','all barcode num')
  df1$sample = i
  list_all[[i]] = df1

}
df_all = do.call(bind_rows,list_all)
write.csv(df_all,'E:\\wenqian polylox\\20250423 new fate bias\\result/all_sample_resi_ratio.csv')




### progeny barcode source fraction
all_nor_pgen_mt = readRDS('E:\\wenqian polylox\\20250918_wenqian_barcode\\qc2_4.rds')
for (i in dir_all) {
  dir1 = paste0('E:\\wenqian polylox\\20250423 new fate bias\\data\\',i,'/')
  if (i %in% c('N13', 'N9' ,'N6', 'N31' ,'N32' ,'N1')) {
    group1 = read.csv(paste0(dir1,'group_df3.csv'),row.names = 1)
  }else{
    group1 = read.csv(paste0(dir1,'group_df2.csv'),row.names = 1)
  }

  group1 = group1[group1$hsc_type!='no_hsc',]
  tissue_ratio = read.csv(paste0(dir1,'tissue_ratio.csv'),row.names = 1)
  tissue_ratio = tissue_ratio[rownames(group1),]
  colnames(tissue_ratio) <- gsub("\\.", " ", colnames(tissue_ratio))

  sample_ratio = list()
  for (j in colnames(tissue_ratio)) {
    barcode_use = tissue_ratio[,j]
    names(barcode_use) = rownames(tissue_ratio)
    barcode_use = barcode_use[barcode_use>0]
    group_use = group1[names(barcode_use),]
    freq = as.data.frame(table(group_use$hsc_type))
    if (nrow(freq)>0) {
      freq$tissue = j
      freq$ratio = freq$Freq/sum(freq$Freq)
      sample_ratio[[j]] = freq
    }

  }
  sample_ratio = do.call(bind_rows,sample_ratio)
  write.csv(sample_ratio,paste0('E:\\wenqian polylox\\20250423 new fate bias\\result\\downstream tissue barcode source/',
                   i,'.csv'),row.names = F)
  # ggplot(sample_ratio,aes(y=tissue,x=ratio,fill=Var1))+geom_col()+
  #   theme_minimal()+theme(text = element_text(size=16))

}

## group3 merge mpp hsc (from 3.1)
all_nor_pgen_mt = readRDS('E:\\wenqian polylox\\20250423 new fate bias\\result/merge_hsc_mpp_for3.1_nospine.rds')
library(tidyverse)
dir_all = dir('E:\\wenqian polylox\\20250423 new fate bias\\data')
list_all = list()
for (i in c('N32','N31','N13','N6','N9')) {
  dir1 = paste0('E:\\wenqian polylox\\20250423 new fate bias\\data\\',i,'/')
  group1 = read.csv(paste0(dir1,'group_df4.csv'),row.names = 1)
  all_df = all_nor_pgen_mt[[i]]
  colnames(all_df) = gsub('\\.',' ',colnames(all_df))
  all_resi = unique(group1$hsc_type)
  all_resi = all_resi[!all_resi %in% c('share','no_hsc_mpp','no_hsc')]

  all_num = c()
  all_tissue = c()
  all_barcode_num = c()
  all_resi_num = c()
  for (j in all_resi) {

    resi_num = nrow(group1[group1$hsc_type==j,])
    j = gsub('-resi','',j)
    j = gsub('\\.',' ',j)
    barcode_num = all_df[,j]
    barcode_num = barcode_num[barcode_num>0]

    all_num = c(all_num,resi_num/length(barcode_num))
    all_tissue = c(all_tissue,j)
    all_resi_num = c(all_resi_num,resi_num)
    all_barcode_num = c(all_barcode_num,length(barcode_num))
  }
  df1=data.frame(all_num,all_tissue,all_resi_num,all_barcode_num )
  colnames(df1) = c('resi_ratio','tissue','resi barcode num','all barcode num')
  df1$sample = i
  list_all[[i]] = df1

}
df_all = do.call(bind_rows,list_all)
write.csv(df_all,'E:\\wenqian polylox\\20250423 new fate bias\\result/all_sample_resi_ratio_group4.csv')






