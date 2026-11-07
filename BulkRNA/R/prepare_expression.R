library(data.table)

# Run from the repository root. TPM uses every gene, before downstream filtering.
dir.create('prepared/tables', recursive = TRUE, showWarnings = FALSE)

for (dataset in c('Paper', 'EPO0804')) {
  counts <- fread(paste0('data/counts/', dataset, '.tsv.gz'))
  lengths <- fread(paste0('data/annotation/', dataset, '_gene_lengths.tsv.gz'))
  length_bp <- lengths$Length[match(counts$gene, lengths$gene)]
  stopifnot(!anyNA(length_bp), all(length_bp > 0))
  rpk <- as.matrix(counts[, -1]) / (length_bp / 1000)
  tpm <- sweep(rpk, 2, colSums(rpk), '/') * 1e6
  fwrite(data.table(gene = counts$gene, tpm), paste0('prepared/tables/', dataset, '_TPM.tsv.gz'), sep = '\t')
}

# Original full 0807 TPM matrix, retained because its source gene-length file is unavailable.
fwrite(fread('data/tpm/PregEMP90807.tsv.gz'), 'prepared/tables/PregEMP90807_TPM.tsv.gz', sep = '\t')

