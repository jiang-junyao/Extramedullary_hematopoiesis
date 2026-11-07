library(data.table)

# Optional: after 00_process_reads.sh, assemble counts and gene lengths for each batch.
# Write to a separate folder, so the supplied manuscript input matrices remain intact.
samples <- fread('data/samples.tsv')[dataset != 'Public']
dir.create('results/reprocessed', showWarnings = FALSE)
for (dataset_now in unique(samples$dataset)) {
  md <- samples[dataset == dataset_now]
  fc <- lapply(md$sample_key, function(s) fread(paste0('results/featurecounts/', s, '.txt'), skip = 1))
  genes <- fc[[1]]$Geneid
  stopifnot(all(vapply(fc, function(x) identical(x$Geneid, genes), logical(1))))
  m <- do.call(cbind, lapply(fc, function(x) x[[ncol(x)]]))
  colnames(m) <- md$sample_key
  fwrite(data.table(gene = genes, m), paste0('results/reprocessed/', dataset_now, '_counts.tsv.gz'), sep = '\t')
  fwrite(fc[[1]][, .(gene = Geneid, Length)], paste0('results/reprocessed/', dataset_now, '_gene_lengths.tsv.gz'), sep = '\t')
  rpk <- m / (fc[[1]]$Length / 1000)
  tpm <- sweep(rpk, 2, colSums(rpk), '/') * 1e6
  fwrite(data.table(gene = genes, tpm), paste0('results/reprocessed/', dataset_now, '_TPM.tsv.gz'), sep = '\t')
}
