# Uniform mixing, read depth and Pgen-first richness; shared by ED2 and Responses 1.2/1.3.
prepared <- readRDS(file.path(result_dir, 'prepared_counts.rds'))
history <- read.csv(file.path(data_dir, 'shuffle_data_all.csv'), check.names = FALSE)
history <- history[history$all_sample %in% ids, ]
history$tissue <- unname(c('Legs LSK' = 'Leg', 'Spine LSK' = 'Spine', 'Spleen LSK' = 'Spleen', 'Liver LSM' = 'Liver')[history$all_ct])
history <- data.frame(analysis_unit = history$all_sample, tissue = history$tissue, iteration = history$all_rep,
  detected = history$all_total, specific = history$all_resi_num, shared = history$all_total - history$all_resi_num,
  specific_percent = 100 * history$all_resi_num / history$all_total,
  shared_percent = 100 * (history$all_total - history$all_resi_num) / history$all_total, source = 'historical_100_replicates')
stopifnot(nrow(history) == 1200, !anyNA(history$tissue))

seed <- 20261102L
set.seed(seed)
new <- real <- list()
for (id in ids) {
  counts <- prepared[[id]]$progenitor
  total <- rowSums(counts)
  real[[id]] <- cbind(analysis_unit = id, summarize_presence(counts), source = 'complete_observed_data')
  for (iteration in seq_len(100)) {
    mixed <- t(vapply(total, function(n) as.vector(rmultinom(1, n, rep(.25, 4))), integer(4)))
    colnames(mixed) <- tissues
    stopifnot(all(rowSums(mixed) == total))
    new[[length(new) + 1]] <- cbind(analysis_unit = id, iteration = iteration, summarize_presence(mixed), source = 'regenerated_model_seed20261102')
  }
}
new <- do.call(rbind, new)
real <- do.call(rbind, real)
summary <- list()
for (d in list(history, new)) {
  for (id in ids) for (tissue in tissues) for (metric in c('shared_percent', 'specific_percent')) {
    z <- d[d$analysis_unit == id & d$tissue == tissue, ]
    summary[[length(summary) + 1]] <- data.frame(analysis_unit = id, tissue = tissue, source = z$source[1], metric = metric,
      n_iterations = nrow(z), mean_percent = mean(z[[metric]]), sd_percent = sd(z[[metric]]),
      q025_percent = unname(quantile(z[[metric]], .025, type = 7)), q975_percent = unname(quantile(z[[metric]], .975, type = 7)))
  }
}
summary <- do.call(rbind, summary)
write.csv(history, file.path(result_dir, 'tables/mixing_historical_iterations.csv'), row.names = FALSE)
write.csv(new, file.path(result_dir, 'tables/mixing_regenerated_iterations.csv'), row.names = FALSE)
write.csv(real, file.path(result_dir, 'tables/observed_tissue_proportions.csv'), row.names = FALSE)
write.csv(summary, file.path(result_dir, 'tables/mixing_summary.csv'), row.names = FALSE)


prepared <- readRDS(file.path(result_dir, 'prepared_counts.rds'))

downsample_reads <- function(x, fraction) {
  x <- as.integer(x)
  total <- sum(x)
  remaining <- floor(total * fraction)
  out <- integer(length(x))
  if (total == 0 || remaining == 0) return(out)
  for (j in seq_len(length(x) - 1)) {
    if (remaining <= 0 || total <= 0) break
    out[j] <- rhyper(1, x[j], total - x[j], remaining)
    remaining <- remaining - out[j]
    total <- total - x[j]
  }
  out[length(x)] <- remaining
  out
}

set.seed(123)
results <- list()
for (id in c('N22', 'N16_1', 'N16_2')) {
  counts <- prepared[[id]]$progenitor
  results[[length(results) + 1]] <- cbind(analysis_unit = id, fraction = 1, iteration = 0, summarize_presence(counts))
  for (fraction in seq(.1, .9, .1)) for (iteration in seq_len(100)) {
    sampled <- apply(counts, 2, downsample_reads, fraction = fraction)
    stopifnot(all(colSums(sampled) == floor(colSums(counts) * fraction)), all(sampled <= counts))
    results[[length(results) + 1]] <- cbind(analysis_unit = id, fraction = fraction, iteration = iteration, summarize_presence(sampled))
  }
}
results <- do.call(rbind, results)
write.csv(results, file.path(result_dir, 'tables/read_downsampling_iterations.csv'), row.names = FALSE)
means <- aggregate(specific_percent ~ analysis_unit + fraction + tissue, results, mean)
write.csv(means, file.path(result_dir, 'tables/read_downsampling_means.csv'), row.names = FALSE)
points <- results[results$iteration %in% 0:10, ]
write.csv(points, file.path(result_dir, 'tables/read_downsampling_display_points.csv'), row.names = FALSE)

prepared <- readRDS(file.path(result_dir, 'prepared_counts.rds'))
counts <- lapply(prepared, function(z) z$progenitor > 0)
original_union <- vapply(counts, nrow, integer(1))
stopifnot(all(original_union == c(N16_1 = 667, N16_2 = 681, N22 = 559)))
targets <- lapply(ids, function(id) setNames(as.integer(round(colSums(counts[[id]]) * 559 / original_union[id])), tissues))
names(targets) <- ids
write.csv(do.call(rbind, lapply(ids, function(id) data.frame(analysis_unit = id, tissue = tissues, original_union = unname(original_union[id]),
  reference_union = 559, target_n = targets[[id]], row.names = NULL))), file.path(result_dir, 'tables/richness_pgen_first_targets.csv'), row.names = FALSE)

set.seed(20260930)
results <- list()
for (iteration in seq_len(1000)) for (id in ids) {
  b <- counts[[id]]
  if (id != 'N22') {
    b[,] <- FALSE
    for (tissue in tissues) {
      eligible <- which(counts[[id]][, tissue])
      chosen <- eligible[sample.int(length(eligible), targets[[id]][tissue], replace = FALSE)]
      b[chosen, tissue] <- TRUE
    }
  }
  stopifnot(all(colSums(b) == targets[[id]]))
  results[[length(results) + 1]] <- cbind(analysis_unit = id, iteration = iteration, summarize_presence(b), final_union = sum(rowSums(b) > 0))
}
results <- do.call(rbind, results)
write.csv(results, file.path(result_dir, 'tables/richness_pgen_first_iterations.csv'), row.names = FALSE)
selected <- c(1, 112, 223, 334, 445, 556, 667, 778, 889, 1000)
points <- results[results$iteration %in% selected, ]
write.csv(points, file.path(result_dir, 'tables/richness_pgen_first_display_points.csv'), row.names = FALSE)
means <- aggregate(specific_percent ~ analysis_unit + tissue, results, mean)
write.csv(means, file.path(result_dir, 'tables/richness_pgen_first_means.csv'), row.names = FALSE)
