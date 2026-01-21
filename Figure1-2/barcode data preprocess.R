
library(FateMapper)
library(pheatmap)
library(tidyverse)
all_nor_pgen_mt = readRDS('E:\\wenqian polylox\\barcode mt 20250311\\qc2_4.rds')
sample_index = readRDS('E:\\wenqian polylox\\barcode mt 20250311\\Index.rds')

all_progenitor = unique(sample_index[sample_index$Cellfate2 == "Progenitor",'Celltype'])


mt_ratio = as.data.frame(t(apply(mt_use, 1, function(x){
  row_sum = sum(x)
  return(x/row_sum)
})))
merge_sample = function(mt_use,sample_index){

  sample_index$ct_index = paste(sample_index$Local,sample_index$Celltype)
  merge_list = list()
  for (i in unique(sample_index$Cellfate2)) {
    ct_use = sample_index[sample_index$Cellfate2==i,4]
    mt_ct = str_split(colnames(mt_use),' ',simplify = T)[,2]
    if (length(intersect(mt_ct,ct_use))>0) {
      target_mt = as.data.frame(mt_use[,colnames(mt_use)[mt_ct %in% ct_use]])
      merge_list[[i]] = rowSums(target_mt)
    }
  }
  merge_df = do.call(bind_cols,merge_list)
  colnames(merge_df) = names(merge_list)
  return(merge_df)
}
### N1
### pp mt
n1_mt = all_nor_pgen_mt$N22
n1_mt = n1_mt[,-grep('Lung',colnames(n1_mt))]
n1_mt = n1_mt[,colnames(n1_mt)!="Liver LSK"]
n1_mt = n1_mt[rowSums(n1_mt)>0,]
n1_mt_merge = as.data.frame(merge_sample(n1_mt,sample_index))
n1_mt_merge = n1_mt_merge[,c("Progenitor",'Megakaryocyte',"Erythroid","Myeloid","Lymphoid")]
rownames(n1_mt_merge) = rownames(n1_mt)
n1_mt_ratio = as.data.frame(t(apply(n1_mt, 1, function(x){
  row_sum = sum(x)
  return(x/row_sum)
})))
n1_mt_merge_ratio = as.data.frame(t(apply(n1_mt_merge, 1, function(x){
  row_sum = sum(x)
  return(x/row_sum)
})))
active_barcode = rownames(n1_mt_merge_ratio)[n1_mt_merge_ratio$Progenitor!=1 & n1_mt_merge_ratio$Progenitor>0]
inactive_barcode = rownames(n1_mt_merge_ratio)[n1_mt_merge_ratio$Progenitor==1]
n1_mt_merge_ratio_no_pro = as.data.frame(t(apply(n1_mt_merge[active_barcode,-1], 1, function(x){
  row_sum = sum(x)
  return(x/row_sum)
})))
n1_mt_merge_log10 = log10(n1_mt_merge)
n1_mt_merge_log10[n1_mt_merge_log10==-Inf] = 0


group = kmeans(n1_mt_merge_ratio_no_pro,centers =5)$cluster
names(group) = rownames(n1_mt_merge_ratio_no_pro)
group = sort(group)

pheatmap(n1_mt_merge_ratio_no_pro[names(group),],cluster_rows = F,
         cluster_cols = F)

group_df = data.frame(group)
group_df$fate = 'None'
group_df[group_df$group==4,2] = 'Lymphoid-bias'
group_df[group_df$group==2,2] = 'Erythroid-bias'
group_df[group_df$group==5,2] = 'Multilineage'
group_df[group_df$group==3,2] = 'Myeloid-bias'
group_df[group_df$group==1,2] = 'Megakaryocyte-bias'
group_df = rbind(group_df[group_df$fate=='Megakaryocyte-bias',],
                 group_df[group_df$fate=='Erythroid-bias',],
                 group_df[group_df$fate=='Myeloid-bias',],
                 group_df[group_df$fate=='Lymphoid-bias',],
                 group_df[group_df$fate=='Multilineage',])
### order tissue column
sample_index_use = sample_index[sample_index$uni_sample %in% colnames(n1_mt),]
sample_index_use = sample_index_use[order(sample_index_use$order_3),]

n1_mt_ratio_order = n1_mt_ratio[rownames(group_df),unique(sample_index_use$uni_sample)]


fate_list = list()
for (i in unique(group_df$fate)) {
  fate_use = group_df[group_df$fate==i,]
  pro_ct = unique(sample_index_use[sample_index_use$Cellfate2=="Progenitor",'uni_sample'])
  mt_use = n1_mt_ratio_order[rownames(fate_use),pro_ct]
  mt_use = mt_use[sort_clone_mt(mt_use),]
  specific_clone <- c()
  for (colname in names(mt_use)[1:4]) {
    rows <- which(apply(mt_use, 1, function(x) {
      x[which(names(mt_use) == colname)] != 0 && all(x[-which(names(mt_use) == colname)] == 0)
    }))
    specific_clone <- c(specific_clone, rows) # 合并到结果向�?

  }
  other_clone = rownames(fate_use)[!rownames(fate_use) %in% names(specific_clone)]
  df1 = fate_use[c(names(specific_clone),other_clone),]
  df1$hsc_type = c(rep('resi',length(specific_clone)),rep('share',length(other_clone)))
  fate_list[[i]] = df1
}
group_df2 = as.data.frame(do.call(bind_rows,fate_list))

inactivte_group = data.frame(rep(11,length(inactive_barcode)),
                             rep('inactive',length(inactive_barcode)))
colnames(inactivte_group) = colnames(group_df)
rownames(inactivte_group) = inactive_barcode
inactivte_group$hsc_type = 'other'
final_group_df = rbind(group_df2,inactivte_group)

sample_index_use = sample_index[sample_index$uni_sample %in% colnames(n1_mt),]
sample_index_use = sample_index_use[order(sample_index_use$order_3),]
order_use = colnames(n1_mt_ratio[rownames(final_group_df),unique(sample_index_use$uni_sample)])[c(1:4,13:20,5:12)]
pheatmap(n1_mt_ratio[rownames(final_group_df),order_use],cluster_rows = F,
         cluster_cols = F,show_rownames = F)
pheatmap(n1_mt_merge_ratio[rownames(final_group_df),],cluster_rows = F,
         cluster_cols = F,show_rownames = F)
final_group_df$barcodes = rownames(final_group_df)

### output
write.csv(final_group_df,'E:\\wenqian polylox\\20250423 new fate bias\\data\\N22/group_df.csv')
write.csv(n1_mt_ratio[rownames(final_group_df),order_use],'E:\\wenqian polylox\\20250423 new fate bias\\data\\N22/tissue_ratio.csv')
write.csv(n1_mt[rownames(final_group_df),order_use],'E:\\wenqian polylox\\20250423 new fate bias\\data\\N22/tissue_count.csv')

write.csv(n1_mt_merge_ratio[rownames(final_group_df),],'E:\\wenqian polylox\\20250423 new fate bias\\data\\N22/merge_ratio.csv')
write.csv(n1_mt_merge[rownames(final_group_df),],'E:\\wenqian polylox\\20250423 new fate bias\\data\\N22/merge_count.csv')

a1= n1_mt_ratio[rownames(final_group_df),unique(sample_index_use$uni_sample)]
a1 = cbind(final_group_df[,-4],a1)
write.csv(a1,'E:\\wenqian polylox\\20250423 new fate bias\\data\\N22/N22_group_ratio.csv')
