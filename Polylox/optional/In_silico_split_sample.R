if (.Platform$OS.type == 'windows') invisible(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
library(ggplot2)
polylox_root <- normalizePath(getwd(), winslash = '/', mustWork = TRUE)
prepared_dir <- file.path(polylox_root, 'prepared/simulation')
data_dir <- file.path(polylox_root, 'data/simulation')
result_dir <- file.path(polylox_root, 'results/In_silico_split_sample')
for (folder in c('plots', 'tables', 'qc')) dir.create(file.path(result_dir, folder), recursive = TRUE, showWarnings = FALSE)
tissues <- c('Leg', 'Spine', 'Spleen', 'Liver')
ids <- c('N16_1', 'N16_2', 'N22')
# Optional: 1000 shared complementary partitions; not part of the supplied PDF panels.
prepared <- readRDS(file.path(prepared_dir, 'prepared_counts.rds'))
meta <- read.delim(file.path(data_dir, 'library_metadata.tsv'), check.names = FALSE)
n_iterations <- 1000L
base_seed <- 20261027L
library_index <- 0L
results <- mappings <- conservation <- list()

split_half <- function(x) {
  total <- sum(x)
  remaining <- rep(floor(total / 2), n_iterations)
  a <- matrix(0L, length(x), n_iterations)
  for (j in which(x > 0)) {
    a[j, ] <- rhyper(n_iterations, x[j], total - x[j], remaining)
    total <- total - x[j]
    remaining <- remaining - a[j, ]
  }
  stopifnot(all(remaining == 0))
  a
}

for (id in ids) for (tissue in tissues) {
  label <- if (tissue == 'Leg') 'Legs' else tissue
  z <- prepared[[id]]
  samples <- meta$sample[meta$mouse == id & meta$tissue == label]
  A <- B <- A_ge2 <- B_ge2 <- matrix(FALSE, length(z$Code), n_iterations)
  original <- rep(0, length(z$Code))
  first <- list()
  for (sample in samples) {
    library_index <- library_index + 1L
    seed <- base_seed + library_index
    set.seed(seed)
    x <- z$library_counts[, sample]
    a <- split_half(x)
    b <- x - a
    stopifnot(all(a + b == x), all(a >= 0), all(b >= 0), all(colSums(a) == floor(sum(x) / 2)))
    A <- A | (a > 0)
    B <- B | (b > 0)
    A_ge2 <- A_ge2 | (a >= 2)
    B_ge2 <- B_ge2 | (b >= 2)
    original <- original + x
    first[[sample]] <- data.frame(analysis_unit = id, tissue = label, sample = sample, Code = z$Code,
      pgen = z$pgen, original_prefiltered_count = x, A = a[, 1], B = b[, 1])
    mappings[[library_index]] <- data.frame(analysis_unit = id, independent_animal_id = 'Unknown', tissue = label, sample = sample,
      seed = seed, raw_reads = sum(z$raw_library_counts[, sample]), removed_singletons = sum(z$raw_library_counts[, sample] == 1),
      sampling_reads = sum(x), reads_A = floor(sum(x) / 2), reads_B = ceiling(sum(x) / 2))
    conservation[[library_index]] <- data.frame(analysis_unit = id, tissue = label, sample = sample, all_barcode_A_plus_B_equals_original = TRUE,
      fixed_half_depth = TRUE, seed = seed, iterations = n_iterations)
  }
  write.csv(do.call(rbind, first), file.path(result_dir, 'tables', paste0('split_', id, '_', label, '_iteration1_library_counts.csv')), row.names = FALSE)
  n0 <- sum(z$good & original > 0)
  for (version in c('post_gt0', 'post_ge2')) {
    aa <- if (version == 'post_gt0') A else A_ge2
    bb <- if (version == 'post_gt0') B else B_ge2
    aa <- aa & z$good
    bb <- bb & z$good
    na <- colSums(aa)
    nb <- colSums(bb)
    ni <- colSums(aa & bb)
    nu <- colSums(aa | bb)
    if (version == 'post_gt0') stopifnot(all(nu == n0))
    results[[length(results) + 1]] <- data.frame(analysis_unit = id, tissue = label, version = version, iteration = seq_len(n_iterations),
      n_original = n0, n_A = na, n_B = nb, n_intersection = ni, n_union = nu, n_A_only = na - ni, n_B_only = nb - ni, n_both_lost = n0 - nu,
      overlap_A = ni / na, overlap_B = ni / nb, jaccard = ni / nu, sample_specific = (nu - ni) / nu, recovery = nu / n0)
  }
}
d <- do.call(rbind, results)
long <- list()
for (metric in c('Shared', 'Specific')) for (side in c('A', 'B')) {
  denominator <- d[[paste0('n_', side)]]
  numerator <- if (metric == 'Shared') d$n_intersection else d[[paste0('n_', side, '_only')]]
  long[[length(long) + 1]] <- data.frame(d[, c('analysis_unit', 'tissue', 'version', 'iteration')], metric = metric, side = side,
    numerator = numerator, denominator = denominator, percent = ifelse(denominator > 0, 100 * numerator / denominator, NA))
}
long <- do.call(rbind, long)
summary <- do.call(rbind, lapply(split(long, interaction(long$analysis_unit, long$tissue, long$version, long$metric, long$side, drop = TRUE)), function(z) {
  data.frame(z[1, c('analysis_unit', 'tissue', 'version', 'metric', 'side')], n_iterations = nrow(z),
    mean_percent = mean(z$percent), lower_percent = unname(quantile(z$percent, .025, type = 7)), upper_percent = unname(quantile(z$percent, .975, type = 7)))
}))
write.csv(d, file.path(result_dir, 'tables/in_silico_iterations.csv'), row.names = FALSE)
write.csv(long, file.path(result_dir, 'tables/in_silico_directional_iterations.csv'), row.names = FALSE)
write.csv(summary, file.path(result_dir, 'tables/in_silico_directional_summary.csv'), row.names = FALSE)
write.csv(do.call(rbind, mappings), file.path(result_dir, 'tables/in_silico_library_seeds.csv'), row.names = FALSE)
write.csv(do.call(rbind, conservation), file.path(result_dir, 'qc/in_silico_conservation.csv'), row.names = FALSE)
for (id in ids) for (metric in c('Shared', 'Specific')) {
  p.data <- summary[summary$analysis_unit == id & summary$metric == metric & summary$version == 'post_gt0', ]
  p.data$tissue <- factor(p.data$tissue, levels = c('Legs', 'Spine', 'Spleen', 'Liver'))
  p.data$side <- factor(p.data$side, levels = c('A', 'B'))
  ymax <- if (metric == 'Shared') 100 else 25
  p <- ggplot(p.data, aes(tissue, mean_percent, fill = side)) +
    geom_col(position = position_dodge(.64), width = .56, colour = 'black', linewidth = .4) +
    geom_errorbar(aes(ymin = lower_percent, ymax = upper_percent), position = position_dodge(.64), width = .13, linewidth = .35) +
    scale_fill_manual(values = c(A = 'black', B = 'white')) +
    scale_x_discrete(labels = c(Legs = 'Leg', Spine = 'Spine', Spleen = 'Spleen', Liver = 'Liver')) +
    scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, ymax / 5), expand = expansion(mult = c(0, .01))) +
    labs(title = id, x = NULL, y = if (metric == 'Shared') 'Shared barcodes (%)' else 'Sample-specific barcodes (%)', fill = NULL) +
    theme_classic(base_family = 'Arial', base_size = 9) +
    theme(axis.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), legend.position = 'right')
  stem <- paste0('not_in_pdf_', id, '_', metric, '_post_gt0')
  ggsave(file.path(result_dir, 'plots', paste0(stem, '.png')), p, width = 8.5 / 2.54, height = 6 / 2.54, dpi = 300, bg = 'white')
  ggsave(file.path(result_dir, 'plots', paste0(stem, '.pdf')), p, width = 8.5 / 2.54, height = 6 / 2.54, device = cairo_pdf, bg = 'white')
}
capture.output(sessionInfo(), file = file.path(result_dir, 'qc/sessionInfo.txt'))
