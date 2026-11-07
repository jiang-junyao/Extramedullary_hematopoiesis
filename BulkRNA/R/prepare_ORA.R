library(data.table)
library(clusterProfiler)

annotation <- unique(fread('data/annotation/TERM2GENE.tsv'), by = c('term', 'gene'))
term2gene <- annotation[, .(term, gene)]

for (comparison in c('EPO_vs_CTRL_SP', 'EMP9_vs_Pregnancy_SP')) {
  deg <- fread(paste0('prepared/tables/', comparison, '_all_genes.tsv.gz'))
  universe <- intersect(deg[is.finite(pvalue), gene], annotation$gene)
  if (comparison == 'EPO_vs_CTRL_SP') genes <- deg[is_DEG == TRUE & log2FoldChange > 0, gene]
  if (comparison == 'EMP9_vs_Pregnancy_SP') genes <- deg[is_DEG == TRUE & log2FoldChange < 0, gene]
  genes <- intersect(genes, universe)
  ora <- enricher(genes, universe = universe, TERM2GENE = term2gene, pAdjustMethod = 'BH',
    pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500)
  fwrite(as.data.table(as.data.frame(ora)), paste0('prepared/tables/', comparison, '_ORA_all.tsv'), sep = '\t')
  fwrite(data.table(gene = universe), paste0('prepared/tables/', comparison, '_ORA_background.tsv'), sep = '\t')
  fwrite(data.table(gene = genes), paste0('prepared/tables/', comparison, '_ORA_query.tsv'), sep = '\t')
}
