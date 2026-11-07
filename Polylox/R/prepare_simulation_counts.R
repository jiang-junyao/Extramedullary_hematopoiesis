meta <- read.delim(file.path(data_dir, 'library_metadata.tsv'), check.names = FALSE)
prepared <- list()
input_summary <- list()
for (id in ids) {
  raw <- read.csv(file.path(data_dir, 'raw', paste0(id, '_raw_pgen.csv')), check.names = FALSE, colClasses = c(Code = 'character'))
  valid_code <- !is.na(raw$Code) & nzchar(trimws(raw$Code))
  cols <- meta$sample[meta$mouse == id]
  stopifnot(!anyDuplicated(raw$Code), all(cols %in% names(raw)), all(as.matrix(raw[!valid_code, cols, drop = FALSE]) == 0))
  raw <- raw[valid_code, , drop = FALSE]
  library_counts <- as.matrix(raw[, cols, drop = FALSE])
  stopifnot(all(is.finite(library_counts)), all(library_counts >= 0), all(library_counts == floor(library_counts)))
  filtered <- filter_library_counts(library_counts, min_count = 2)
  rownames(filtered) <- raw$Code
  merged <- matrix(0, nrow(raw), 4, dimnames = list(raw$Code, tissues))
  for (tissue in tissues) {
    label <- if (tissue == 'Leg') 'Legs' else tissue
    use <- meta$sample[meta$mouse == id & meta$tissue == label]
    merged[, tissue] <- rowSums(filtered[, use, drop = FALSE])
  }
  good <- low_pgen(raw$pgen, cutoff = 1e-5)
  selected <- good & rowSums(merged) > 0
  progenitor <- merged[selected, , drop = FALSE]
  prepared[[id]] <- list(Code = raw$Code, pgen = raw$pgen, good = good, library_counts = filtered,
    raw_library_counts = library_counts, merged = merged, progenitor = progenitor)
  reference <- as.matrix(read.csv(file.path(data_dir, 'reference', paste0(id, '_progenitor_matrix.csv')), row.names = 1, check.names = FALSE))
  stopifnot(identical(rownames(progenitor), rownames(reference)), identical(colnames(progenitor), colnames(reference)), all(progenitor == reference))
  input_summary[[id]] <- data.frame(analysis_unit = id, independent_animal_id = 'Unknown', libraries = ncol(filtered),
    source_barcode_rows = nrow(raw), singleton_reads_removed = sum(library_counts == 1), pre_pgen_union = sum(rowSums(merged) > 0),
    post_pgen_union = nrow(progenitor), matches_legacy_matrix = TRUE)
}
saveRDS(prepared, file.path(result_dir, 'prepared_counts.rds'))
write.csv(do.call(rbind, input_summary), file.path(result_dir, 'qc/input_summary.csv'), row.names = FALSE)
