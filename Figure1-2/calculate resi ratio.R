### share barcode/tissue barcode
library(tidyverse)



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






