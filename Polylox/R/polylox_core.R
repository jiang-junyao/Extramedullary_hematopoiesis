# 只共享不依赖具体figure的基本操作；组织/重复合并规则由各模块显式决定。
filter_library_counts <- function(counts, min_count = 2) {
  counts <- as.matrix(counts)
  counts[!is.na(counts) & counts < min_count] <- 0
  counts
}

# 在Pgen筛选之前调用，NA分母或全零文库的frequency保持NA。
library_frequency <- function(counts) {
  denominator <- colSums(counts)
  frequency <- sweep(counts, 2, denominator, '/')
  frequency[, !is.finite(denominator) | denominator <= 0] <- NA_real_
  frequency
}

low_pgen <- function(pgen, cutoff = 1e-5) {
  is.finite(pgen) & pgen >= 0 & pgen < cutoff
}

# 输入列必须已经是本次比较使用的组织；缺测不作为未检出。
classify_tissue_barcodes <- function(counts) {
  detected <- counts > 0
  n_tissues <- rowSums(detected)
  type <- ifelse(is.na(n_tissues), NA_character_, ifelse(n_tissues > 1, 'Shared', ifelse(n_tissues == 1, 'Specific', 'Undetected')))
  data.frame(barcode = rownames(counts), n_tissues = n_tissues, barcode_type = type, stringsAsFactors = FALSE)
}
