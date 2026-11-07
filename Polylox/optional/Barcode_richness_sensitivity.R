if (.Platform$OS.type == 'windows') invisible(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
library(ggplot2)
polylox_root <- normalizePath(getwd(), winslash = '/', mustWork = TRUE)
prepared_dir <- file.path(polylox_root, 'prepared/simulation')
data_dir <- file.path(polylox_root, 'data/simulation')
result_dir <- file.path(polylox_root, 'results/Barcode_richness_sensitivity')
for (folder in c('plots', 'tables', 'qc')) dir.create(file.path(result_dir, folder), recursive = TRUE, showWarnings = FALSE)
tissues <- c('Leg', 'Spine', 'Spleen', 'Liver')
ids <- c('N16_1', 'N16_2', 'N22')
# Optional historical alternatives; not the ED2g Pgen-first branch.
summarize_presence <- function(x) {
  detected <- x > 0
  breadth <- rowSums(detected)
  total <- colSums(detected)
  specific <- colSums(detected & (breadth == 1))
  data.frame(tissue = colnames(x), detected = total, specific = specific, shared = total - specific,
    specific_percent = ifelse(total > 0, 100 * specific / total, NA),
    shared_percent = ifelse(total > 0, 100 * (total - specific) / total, NA), row.names = NULL)
}

prepared <- readRDS(file.path(prepared_dir, 'prepared_counts.rds'))
before <- lapply(prepared, function(z) z$merged[rowSums(z$merged) > 0, , drop = FALSE] > 0)
pass <- lapply(prepared, function(z) z$good[rowSums(z$merged) > 0])
union_n <- vapply(before, nrow, integer(1))
stopifnot(all(union_n == c(N16_1 = 2270, N16_2 = 2319, N22 = 1918)))
targets <- lapply(ids, function(id) setNames(as.integer(round(colSums(before[[id]]) * union_n['N22'] / union_n[id])), tissues))
names(targets) <- ids
write.csv(do.call(rbind, lapply(ids, function(id) data.frame(analysis_unit = id, tissue = tissues, original_union = union_n[id],
  reference_union = union_n['N22'], original_tissue_n = colSums(before[[id]]), target_n = targets[[id]]))), file.path(result_dir, 'tables/richness_per_tissue_targets.csv'), row.names = FALSE)

set.seed(20260930)
results <- list()
for (iteration in seq_len(1000)) for (id in ids) {
  b <- before[[id]]
  if (id != 'N22') {
    b[,] <- FALSE
    for (tissue in tissues) {
      eligible <- which(before[[id]][, tissue])
      chosen <- eligible[sample.int(length(eligible), targets[[id]][tissue], replace = FALSE)]
      b[chosen, tissue] <- TRUE
    }
  }
  stopifnot(all(colSums(b) == targets[[id]]))
  results[[length(results) + 1]] <- cbind(analysis_unit = id, iteration = iteration, summarize_presence(b[pass[[id]], , drop = FALSE]))
}
results <- do.call(rbind, results)
write.csv(results, file.path(result_dir, 'tables/richness_per_tissue_iterations.csv'), row.names = FALSE)
selected <- c(1, 112, 223, 334, 445, 556, 667, 778, 889, 1000)
write.csv(results[results$iteration %in% selected, ], file.path(result_dir, 'tables/richness_per_tissue_display_points.csv'), row.names = FALSE)
observed <- read.csv(file.path(prepared_dir, 'tables/observed_tissue_proportions.csv'))
for (adult in c('N16_1', 'N16_2')) {
  a <- observed[observed$analysis_unit %in% c('N22', adult), ]
  a$series <- ifelse(a$analysis_unit == 'N22', 'Embryonic original', 'Adult original')
  b <- results[results$analysis_unit == adult & results$iteration %in% selected, ]
  b$series <- 'Adult normalized'
  p.data <- rbind(a[, c('tissue', 'series', 'specific_percent')], b[, c('tissue', 'series', 'specific_percent')])
  p.data$tissue <- factor(p.data$tissue, levels = tissues)
  p.data$series <- factor(p.data$series, levels = c('Embryonic original', 'Adult original', 'Adult normalized'))
  bars <- aggregate(specific_percent ~ tissue + series, p.data, mean)
  write.csv(bars, file.path(result_dir, 'tables', paste0('not_in_pdf_richness_then_pgen_', adult, '_bar_values.csv')), row.names = FALSE)
  p <- ggplot(bars, aes(tissue, specific_percent, fill = series)) +
    geom_col(position = position_dodge(.82), width = .64, colour = 'black', linewidth = .4) +
    geom_point(data = p.data, aes(fill = series), position = position_jitterdodge(jitter.width = .08, dodge.width = .82, seed = 20261002), shape = 21, size = 2, stroke = .45) +
    scale_fill_manual(values = c('Embryonic original' = '#CCCCCC', 'Adult original' = '#9999A2', 'Adult normalized' = 'white')) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
    labs(title = paste(adult, 'vs N22; per-tissue richness sensitivity'), x = NULL, y = 'Tissue-specific HSPC clones (%)', fill = NULL) +
    theme_classic(base_family = 'Arial', base_size = 11) + theme(axis.text = element_text(colour = 'black'), legend.position = 'right')
  stem <- paste0('not_in_pdf_richness_then_pgen_', adult)
  ggsave(file.path(result_dir, 'plots', paste0(stem, '.png')), p, width = 5.5, height = 3.8, dpi = 300, bg = 'white')
  ggsave(file.path(result_dir, 'plots', paste0(stem, '.pdf')), p, width = 5.5, height = 3.8, device = cairo_pdf, bg = 'white')
}

set.seed(20260930)
union_results <- list()
for (iteration in seq_len(1000)) for (id in ids) {
  counts <- prepared[[id]]$progenitor > 0
  eligible <- which(rowSums(counts) > 0)
  chosen <- if (length(eligible) == 559) eligible else eligible[sample.int(length(eligible), 559, replace = FALSE)]
  union_results[[length(union_results) + 1]] <- cbind(analysis_unit = id, iteration = iteration, summarize_presence(counts[chosen, , drop = FALSE]), final_union = 559)
}
write.csv(do.call(rbind, union_results), file.path(result_dir, 'tables/richness_mouse_union_iterations.csv'), row.names = FALSE)
