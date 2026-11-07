library(data.table)
library(fgsea)
library(BiocParallel)

analyses <- data.table(
  name = c('Pregnancy_signature_in_Public', 'Pregnancy_signature_in_EPO', 'EPO_signature_in_EMP9'),
  comparison = c('Public_Pregnancy_SP_vs_BM', 'EPO_SP_vs_BM', 'EMP9_vs_Pregnancy_SP'),
  signature = c('Pregnancy_SP_high_Top200', 'Pregnancy_SP_high_Top200', 'EPO_up_signature'),
  left = c('Pregnant\nBM HSC', 'EPO\nBM HSC', 'Pregnancy\nSP HSC'),
  right = c('Pregnant\nSP HSC', 'EPO\nSP HSC', 'Pregnancy + EMP9\nSP HSC'))

for (i in seq_len(nrow(analyses))) {
  a <- analyses[i]
  ranked <- fread(paste0('prepared/tables/', a$comparison, '_all_genes.tsv.gz'))[is.finite(log2FoldChange)]
  setorder(ranked, -log2FoldChange, gene)
  # Preserve deterministic tie handling and the exact reversal used in final figures.
  ranked[, score := log2FoldChange - seq_len(.N) * 1e-12]
  ranked <- ranked[.N:1]
  ranked[, score := -score]
  ranked[, rank := seq_len(.N)]
  genes <- fread(paste0('prepared/tables/', a$signature, '.tsv'))$gene
  stats <- setNames(ranked$score, ranked$gene)
  set.seed(20261003)
  fg <- fgseaMultilevel(pathways = list(signature = genes), stats = stats, minSize = 10, maxSize = 1000,
    eps = 0, scoreType = 'std', gseaParam = 1, sampleSize = 101, nPermSimple = 1000, BPPARAM = SerialParam())
  ranked[, hit := gene %in% genes]
  ranked[, running_ES := cumsum(ifelse(hit, abs(score), 0)) / sum(abs(score[hit])) - cumsum(!hit) / sum(!hit)]
  leading <- fg$leadingEdge[[1]]
  fg[, leadingEdge := paste(leading, collapse = ';')]
  fg[, `:=`(left = a$left, right = a$right, signature_genes = length(genes), matched_genes = sum(ranked$hit))]
  fwrite(fg, paste0('prepared/tables/', a$name, '_GSEA.tsv'), sep = '\t', quote = TRUE)
  fwrite(ranked, paste0('prepared/tables/', a$name, '_curve.tsv.gz'), sep = '\t')
}
