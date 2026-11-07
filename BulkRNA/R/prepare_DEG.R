library(data.table)
library(DESeq2)

samples <- fread('data/samples.tsv')
membership <- fread('data/comparison_samples.tsv')
comparisons <- fread('data/comparisons.tsv')
count_tables <- lapply(c('Paper', 'EPO0804', 'PregEMP90807', 'Public'), function(x) fread(paste0('data/counts/', x, '.tsv.gz')))
names(count_tables) <- c('Paper', 'EPO0804', 'PregEMP90807', 'Public')

for (i in seq_len(nrow(comparisons))) {
  comparison <- comparisons[i]
  md <- merge(membership[comparison == comparisons$comparison[i]], samples, by = 'sample_key', sort = FALSE)
  setorder(md, sample_order)
  if (comparison$comparison == 'Pregnancy_SP_vs_BM') setorder(md, dataset, group, sample_id)
  datasets <- unique(md$dataset)
  genes <- Reduce(intersect, lapply(count_tables[datasets], function(x) x$gene))
  if (length(datasets) > 1) genes <- sort(genes)
  matrices <- lapply(datasets, function(d) {
    x <- count_tables[[d]][match(genes, gene)]
    as.matrix(x[, md[dataset == d, sample_key], with = FALSE])
  })
  count_matrix <- do.call(cbind, matrices)[, md$sample_key, drop = FALSE]
  rownames(count_matrix) <- genes
  storage.mode(count_matrix) <- 'integer'
  if (comparison$remove_all_zero) count_matrix <- count_matrix[rowSums(count_matrix) > 0, , drop = FALSE]

  md <- as.data.frame(md)
  rownames(md) <- md$sample_key
  md$group <- factor(md$group, levels = c(comparison$reference, comparison$numerator))
  md$cohort <- factor(md$cohort)
  if (comparison$comparison == 'Pregnancy_SP_vs_BM') md$cohort <- factor(md$cohort, levels = c('Paper', '0807'))
  stopifnot(identical(colnames(count_matrix), rownames(md)), !anyDuplicated(rownames(count_matrix)))
  dds <- DESeqDataSetFromMatrix(count_matrix, md, design = as.formula(comparison$design))
  dds <- DESeq(dds, quiet = TRUE)

  # Preserve the original independent-filtering alpha; final DEG selection is below.
  res <- results(dds, contrast = c('group', comparison$numerator, comparison$reference), alpha = 0.25)
  deg <- data.table(gene = rownames(res), as.data.frame(res))
  deg[, is_DEG := is.finite(pvalue) & pvalue < 0.01 & is.finite(padj) & padj < 0.20 & abs(log2FoldChange) >= log2(1.5)]
  deg[, direction := fcase(is_DEG & log2FoldChange > 0, 'Up', is_DEG & log2FoldChange < 0, 'Down', default = 'Not DEG')]
  fwrite(deg, paste0('prepared/tables/', comparison$comparison, '_all_genes.tsv.gz'), sep = '\t', na = 'NA')
  fwrite(deg[is_DEG == TRUE], paste0('prepared/tables/', comparison$comparison, '_DEG.tsv'), sep = '\t')
}
