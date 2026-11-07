library(dplyr)
library(tidyr)
source('R/polylox_core.R')

out <- 'prepared/tissue/tables'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
pgen_cutoff <- 1e-5
sample_metadata <- read.delim('data/tissue/metadata/all_sample_metadata.tsv', check.names = FALSE)
records <- read.csv('data/tissue/metadata/record_strategy.csv', stringsAsFactors = FALSE)
sample_metadata$tissue[sample_metadata$tissue == 'Legs'] <- 'Leg'

# Each raw library is independently normalized, not pooled by raw read depth.
normalize_cpm <- function(x) {
  library_frequency(x) * 1e6
}

frequency <- list()
classification <- list()
percentages <- list()
denominators <- list()
input_qc <- list()
sample_map <- list()

# The 18 historical records use their saved population strategies.
# Exp1/Exp2 are a separate verified historical branch, not assumed to be PDF animals.
for (unit in c(records$mouse, 'Exp1', 'Exp2')) {
  if (unit %in% records$mouse) {
    dat <- read.csv(paste0('data/tissue/raw/', unit, '_raw_pgen.csv'), colClasses = c(Code = 'character'), check.names = FALSE)
    strategy <- records$strategy[records$mouse == unit]
    meta <- sample_metadata[sample_metadata$mouse == unit, ]
    if (strategy == 'HSC_MPP') meta <- meta[meta$tissue %in% c('Leg', 'Spine', 'Spleen') & meta$cell_type %in% c('HSC', 'MPP'), ]
    if (strategy == 'LSK_LSM') meta <- meta[(meta$tissue %in% c('Leg', 'Spine', 'Spleen') & meta$cell_type == 'LSK') | (meta$tissue == 'Liver' & meta$cell_type == 'LSM'), ]
    if (strategy == 'LSK_3T') meta <- meta[meta$tissue %in% c('Leg', 'Spine', 'Spleen') & meta$cell_type == 'LSK', ]
    input_rows <- nrow(dat)
    missing_pgen <- sum(is.na(dat$pgen))
    # Preserve the legacy cutoff script's normalization population.
    dat <- dat[!is.na(dat$pgen), ]
    pgen_denominator_rule <- 'Missing Pgen rows removed before normalization (legacy cutoff script)'
  } else {
    dat <- read.delim(paste0('data/tissue/raw/', unit, '.tsv'), colClasses = c(Code = 'character'), check.names = FALSE)
    selected <- c('Leg LSK-1', 'Leg LSK-2', 'Spine LSK-1', 'Spine LSK-2', 'Spleen LSK-1', 'Spleen LSK-2')
    if (unit == 'Exp2') selected <- c(selected, 'Liver LSM')
    meta <- data.frame(mouse = unit, sample = selected, tissue = sub(' .*', '', selected),
      cell_type = ifelse(grepl('LSM', selected), 'LSM', 'LSK'), replicate = NA)
    strategy <- ifelse(unit == 'Exp1', 'LSK_3T', 'LSK_LSM')
    input_rows <- nrow(dat)
    missing_pgen <- sum(is.na(dat$pgen))
    pgen_denominator_rule <- 'All barcode rows retained in normalization (September22 script)'
  }

  stopifnot(!anyDuplicated(dat$Code), all(meta$sample %in% names(dat)))
  counts <- as.matrix(dat[, meta$sample, drop = FALSE])
  storage.mode(counts) <- 'numeric'
  stopifnot(all(is.finite(counts)), all(counts >= 0), all(counts == floor(counts)))
  singleton_n <- colSums(counts == 1)
  counts <- filter_library_counts(counts, min_count = 2)
  library_depth <- colSums(counts)
  normalized <- normalize_cpm(counts)
  denominators[[unit]] <- data.frame(analysis_unit = unit, sample = meta$sample, singleton_n, library_depth, pgen_denominator_rule)
  meta$analysis_unit <- unit
  sample_map[[unit]] <- meta

  # Same tissue and same cell type: sum independently normalized repeats, then renormalize.
  meta$profile <- paste(meta$tissue, meta$cell_type, sep = '_')
  profiles <- unique(meta$profile)
  merged <- sapply(profiles, function(x) rowSums(normalized[, meta$sample[meta$profile == x], drop = FALSE]))
  merged <- normalize_cpm(merged)
  tissues <- c('Leg', 'Spine', 'Spleen', 'Liver')
  tissues <- tissues[tissues %in% meta$tissue]

  # HSC and MPP are combined after their individual profile normalization.
  # This old branch deliberately has sum 2e6 per tissue before Pgen, not 1e6.
  tissue_frequency <- sapply(tissues, function(tissue) rowSums(merged[, startsWith(colnames(merged), paste0(tissue, '_')), drop = FALSE]))
  stopifnot(all(is.na(dat$pgen) | dat$pgen >= 0))
  retained <- low_pgen(dat$pgen, cutoff = pgen_cutoff) & rowSums(tissue_frequency) > 0
  tissue_frequency <- tissue_frequency[retained, , drop = FALSE]
  codes <- dat$Code[retained]
  pgen <- dat$pgen[retained]
  rownames(tissue_frequency) <- codes
  n_tissues <- classify_tissue_barcodes(tissue_frequency)$n_tissues
  group <- ifelse(n_tissues >= 2, 'Shared', paste0(tissues[max.col(tissue_frequency, ties.method = 'first')], '-specific'))

  frequency[[unit]] <- data.frame(analysis_unit = unit, Code = codes, pgen, tissue_frequency, check.names = FALSE)
  classification[[unit]] <- data.frame(analysis_unit = unit, Code = codes, pgen, n_tissues, classification = group, strategy)
  percentages[[unit]] <- data.frame(analysis_unit = unit, tissue = tissues,
    detected_n = colSums(tissue_frequency > 0), specific_n = colSums(tissue_frequency > 0 & n_tissues == 1),
    shared_n = colSums(tissue_frequency > 0 & n_tissues >= 2)) %>%
    mutate(specific_percent = 100 * specific_n / detected_n, shared_percent = 100 * shared_n / detected_n)
  input_qc[[unit]] <- data.frame(analysis_unit = unit, input_rows, missing_pgen, normalized_rows = nrow(dat),
    retained_n = length(codes), shared_n = sum(n_tissues >= 2), specific_n = sum(n_tissues == 1), measured_tissues = paste(tissues, collapse = ';'))
}

write.table(bind_rows(frequency), file.path(out, 'tissue_frequency_cpm.tsv'), sep = '\t', row.names = FALSE, quote = FALSE, na = 'NA')
write.table(bind_rows(classification), file.path(out, 'barcode_classification.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
write.table(bind_rows(percentages), file.path(out, 'within_tissue_percentages.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
write.table(bind_rows(denominators), file.path(out, 'library_denominators.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
write.table(bind_rows(sample_map), file.path(out, 'selected_sample_map.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
write.table(bind_rows(input_qc), file.path(out, 'input_record_summary.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)



# Shared row percentages and barcode ordering are prepared once, without making figures.
prepared <- list(frequency = frequency, classification = classification)
row_percentages <- list()
for (unit in names(prepared$frequency)) {
  dat <- left_join(prepared$frequency[[unit]], prepared$classification[[unit]][, c('analysis_unit', 'Code', 'n_tissues', 'classification')], by = c('analysis_unit', 'Code'))
  tissues <- intersect(c('Leg', 'Spine', 'Spleen', 'Liver'), names(dat))
  dat$total <- rowSums(dat[, tissues, drop = FALSE])
  dat$classification <- factor(dat$classification, levels = c('Shared', paste0(tissues, '-specific')))
  dat <- dat %>% arrange(classification, desc(n_tissues), desc(total), desc(Spine), desc(Leg), desc(Spleen))
  percent <- sweep(as.matrix(dat[, tissues, drop = FALSE]), 1, dat$total, '/') * 100
  plot_source <- data.frame(analysis_unit = unit, Code = dat$Code, classification = dat$classification, row_order = seq_len(nrow(dat)), percent, check.names = FALSE)
  row_percentages[[unit]] <- plot_source
}
write.table(bind_rows(row_percentages), file.path(out, 'heatmap_row_percentage.tsv'), sep = '\t', row.names = FALSE, quote = FALSE, na = 'NA')
formal <- read.delim('data/tissue/formal/formal_specific_percent.tsv', stringsAsFactors = FALSE)
write.table(formal, file.path(out, 'formal_specific_percent_replayed.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
saveRDS(list(frequency = frequency, classification = classification, row_percentages = row_percentages, formal = formal), 'prepared/tissue/prepared_tissue.rds')
