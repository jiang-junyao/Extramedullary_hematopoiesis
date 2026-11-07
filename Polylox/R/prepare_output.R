# Run from the package root. Shared preparation computes data and four historical pooled tests, without plotting.
library(dplyr)
library(tidyr)
source('R/polylox_core.R')

input_dir <- 'data/output'
output_result <- 'prepared/output'
for (folder in c('tables', 'qc')) dir.create(file.path(output_result, folder), recursive = TRUE, showWarnings = FALSE)

mice <- read.csv(file.path(input_dir, 'mice.csv'), fileEncoding = 'UTF-8-BOM')
metadata <- read.delim(file.path(input_dir, 'libraries.tsv'), check.names = FALSE, fileEncoding = 'UTF-8-BOM')
metadata$tissue[metadata$tissue == 'Leg'] <- 'Legs'
metadata$cell_type[metadata$cell_type == 'Gr'] <- 'Neu'
tissues <- c('Legs', 'Spine', 'Spleen')
output_types <- c('EryP', 'Mono', 'Neu', 'B', 'MEP', 'GMP')
cell_types <- c('HSC', 'MPP', output_types)
pgen_cutoff <- 1e-5
stopifnot(nrow(mice) == 6, all(table(mice$condition) == 3), !anyDuplicated(metadata[c('mouse', 'sample')]))

# Compute complete profiles once. The normalization denominator precedes Pgen filtering.
profile_rows <- source_rows <- library_rows <- excluded_rows <- list()
for (mouse_id in mice$mouse) {
  raw <- read.csv(file.path(input_dir, 'raw', paste0(mouse_id, '_raw_pgen.csv')), check.names = FALSE, colClasses = c(Code = 'character'))
  valid_code <- !is.na(raw$Code) & trimws(raw$Code) != '' & trimws(raw$Code) != 'NA'
  excluded_rows[[mouse_id]] <- data.frame(mouse = mouse_id, input_rows = nrow(raw), missing_code_excluded = sum(!valid_code))
  raw <- raw[valid_code, ]
  raw$Code <- trimws(raw$Code)
  meta <- metadata[metadata$mouse == mouse_id, ]
  counts <- as.matrix(raw[, meta$sample, drop = FALSE])
  storage.mode(counts) <- 'double'
  stopifnot(!anyDuplicated(raw$Code), all(is.finite(counts)), all(counts >= 0), all(counts == floor(counts)))
  stopifnot(nrow(unique(meta[c('tissue', 'cell_type')])) == 24, !any(raw$pgen < 0, na.rm = TRUE))
  singleton_n <- colSums(counts == 1)
  counts <- filter_library_counts(counts, min_count = 2)
  totals <- colSums(counts)
  stopifnot(all(totals > 0))
  frequency <- library_frequency(counts)
  library_rows[[mouse_id]] <- data.frame(mouse = mouse_id, meta[, c('sample', 'tissue', 'cell_type', 'replicate')], singleton_n, normalization_denominator = totals)

  profiles <- list()
  for (tissue_id in tissues) for (type_id in cell_types) {
    merged <- rowSums(frequency[, meta$tissue == tissue_id & meta$cell_type == type_id, drop = FALSE])
    profiles[[paste(tissue_id, type_id)]] <- merged / sum(merged)
  }
  hspc <- sapply(tissues, function(tissue_id) {
    merged <- profiles[[paste(tissue_id, 'HSC')]] + profiles[[paste(tissue_id, 'MPP')]]
    merged / sum(merged)
  })

  # Each barcode identity remains mouse + Code. Shared is co-detection, not migration.
  eligible <- low_pgen(raw$pgen, cutoff = pgen_cutoff)
  source <- data.frame(mouse = mouse_id, condition = mice$condition[match(mouse_id, mice$mouse)], Code = raw$Code,
    pgen = raw$pgen, eligible_pgen = eligible, HSPC_Legs = hspc[, 'Legs'], HSPC_Spine = hspc[, 'Spine'], HSPC_Spleen = hspc[, 'Spleen'])
  source$origin <- ifelse(rowSums(hspc > 0) >= 2, 'Shared', ifelse(hspc[, 'Spleen'] > 0, 'Non-shared', 'No-splenic-HSPC'))
  source_rows[[mouse_id]] <- source
  profile_rows[[mouse_id]] <- data.frame(mouse = mouse_id, condition = source$condition, Code = raw$Code, pgen = raw$pgen,
    eligible_pgen = eligible, as.data.frame(profiles, check.names = FALSE), check.names = FALSE)
}
profiles_all <- bind_rows(profile_rows)
hspc_sources <- bind_rows(source_rows)
libraries <- bind_rows(library_rows)
write.csv(profiles_all, file.path(output_result, 'tables/complete_profiles.csv'), row.names = FALSE)
write.csv(hspc_sources, file.path(output_result, 'tables/hspc_sources.csv'), row.names = FALSE)
write.csv(libraries, file.path(output_result, 'tables/library_denominators.csv'), row.names = FALSE)
write.csv(bind_rows(excluded_rows), file.path(output_result, 'tables/excluded_rows.csv'), row.names = FALSE)

# Local B and All3T6 use the same anchors and classification; only output scope differs.
barcode_rows <- list()
for (version_id in c('B', 'All3T6')) {
  selected <- hspc_sources$eligible_pgen & hspc_sources$HSPC_Spleen > 0
  output <- sapply(output_types, function(type_id) {
    if (version_id == 'B') return(profiles_all[[paste('Spleen', type_id)]])
    # Normalize separately within each mouse, not across pooled animals.
    merged <- rowSums(profiles_all[, paste(tissues, type_id), drop = FALSE])
    merged / ave(merged, profiles_all$mouse, FUN = sum)
  })
  data <- hspc_sources[selected, ]
  data$version <- version_id
  data$x <- data$HSPC_Spleen
  data <- cbind(data, output[selected, , drop = FALSE])
  data$output_sum <- rowSums(data[, output_types])
  data$active_type_n <- rowSums(data[, output_types] > 0)
  data$active <- data$active_type_n > 0
  barcode_rows[[version_id]] <- data
}
barcode <- bind_rows(barcode_rows)
stopifnot(!anyDuplicated(barcode[c('version', 'mouse', 'Code')]))
write.csv(barcode, file.path(output_result, 'tables/barcode_output.csv'), row.names = FALSE)

counts <- barcode %>% group_by(version, condition, mouse) %>% summarise(total = n(),
  active_specific = sum(active & origin == 'Non-shared'), active_shared = sum(active & origin == 'Shared'),
  inactive_specific = sum(!active & origin == 'Non-shared'), inactive_shared = sum(!active & origin == 'Shared'), .groups = 'drop')
write.csv(counts, file.path(output_result, 'tables/counts_by_mouse.csv'), row.names = FALSE)
counts_pooled <- counts %>% group_by(version, condition) %>% summarise(mouse_n = n(), across(total:inactive_shared, sum), .groups = 'drop')
write.csv(counts_pooled, file.path(output_result, 'tables/counts_pooled.csv'), row.names = FALSE)

# Preserve the historical, exploratory pooled barcode test. This is not mouse-level inference.
boxes <- barcode[barcode$active, ]
statistics <- boxes %>% group_by(version, condition) %>% group_modify(~ {
  a <- .x$output_sum[.x$origin == 'Non-shared']
  b <- .x$output_sum[.x$origin == 'Shared']
  test <- wilcox.test(a, b, exact = FALSE, correct = TRUE)
  data.frame(specific_n = length(a), shared_n = length(b), W = unname(test$statistic), p_unadjusted = test$p.value)
}) %>% ungroup()
write.csv(boxes, file.path(output_result, 'tables/active_boxplot_observations.csv'), row.names = FALSE)
write.csv(statistics, file.path(output_result, 'tables/wilcoxon.csv'), row.names = FALSE)

# RDS preserves calculation precision; each Figure selects and plots its own condition.
saveRDS(list(barcode = barcode, statistics = statistics, counts_by_mouse = counts, counts_pooled = counts_pooled), file.path(output_result, 'tables/output_data.rds'))
qc <- data.frame(check = c('unique_strategy_mouse_Code', 'six_mice_three_per_condition', 'four_exploratory_tests'),
  passed = c(!anyDuplicated(barcode[c('version', 'mouse', 'Code')]), nrow(mice) == 6 && all(table(mice$condition) == 3), nrow(statistics) == 4))
write.csv(qc, file.path(output_result, 'qc/preparation_checks.csv'), row.names = FALSE)
stopifnot(all(qc$passed))

