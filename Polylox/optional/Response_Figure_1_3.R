if (.Platform$OS.type == 'windows') invisible(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
library(ggplot2)
polylox_root <- normalizePath(getwd(), winslash = '/', mustWork = TRUE)
prepared_dir <- file.path(polylox_root, 'prepared/simulation')
data_dir <- file.path(polylox_root, 'data/simulation')
result_dir <- file.path(polylox_root, 'results/Response_Figure_1_3')
for (folder in c('plots', 'tables', 'qc')) dir.create(file.path(result_dir, folder), recursive = TRUE, showWarnings = FALSE)
tissues <- c('Leg', 'Spine', 'Spleen', 'Liver')

# Only Response Figure 1.3 a-e, sharing prepared numerical results with ED2.
points <- read.csv(file.path(prepared_dir, 'tables/read_downsampling_display_points.csv'))
means <- read.csv(file.path(prepared_dir, 'tables/read_downsampling_means.csv'))
adult <- 'N16_2'
cols <- c(Leg = '#8BC98D', Spine = '#77A8E8', Spleen = '#C39F9B', Liver = '#E7B164')
for (tissue in tissues) {
  p.data <- points[points$analysis_unit %in% c(adult, 'N22') & points$tissue == tissue, ]
  q.data <- means[means$analysis_unit %in% c(adult, 'N22') & means$tissue == tissue, ]
  p.data$analysis_unit <- factor(p.data$analysis_unit, levels = c(adult, 'N22'))
  p <- ggplot() +
    geom_smooth(data = q.data, aes(fraction, specific_percent, group = analysis_unit), method = 'loess', se = TRUE, span = .85,
      colour = cols[tissue], fill = cols[tissue], alpha = .20, linewidth = .75) +
    geom_point(data = p.data, aes(fraction, specific_percent, fill = analysis_unit), shape = 21, size = 2, colour = cols[tissue], stroke = .5) +
    scale_fill_manual(values = setNames(c(cols[tissue], 'white'), c(adult, 'N22'))) +
    scale_x_continuous(breaks = seq(.2, 1, .2), limits = c(.05, 1.05), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
    labs(title = tissue, x = 'Fraction of reads used', y = 'Tissue-specific HSPC clones (%)', fill = NULL) +
    theme_classic(base_family = 'Arial', base_size = 11) +
    theme(axis.text = element_text(colour = 'black'), legend.position = c(.2, .12))
  stem <- paste0('ED2_', letters[match(tissue, tissues) + 2], '_', adult, '_', tissue)
  ggsave(file.path(result_dir, 'plots', paste0(stem, '.png')), p, width = 5.5, height = 3.8, dpi = 300, bg = 'white')
  ggsave(file.path(result_dir, 'plots', paste0(stem, '.pdf')), p, width = 5.5, height = 3.8, device = cairo_pdf, bg = 'white')
  write.csv(p.data, file.path(result_dir, 'tables', paste0(stem, '_points.csv')), row.names = FALSE)
  write.csv(q.data, file.path(result_dir, 'tables', paste0(stem, '_means.csv')), row.names = FALSE)
}

points <- read.csv(file.path(prepared_dir, 'tables/richness_pgen_first_display_points.csv'))
means <- read.csv(file.path(prepared_dir, 'tables/richness_pgen_first_means.csv'))
adult <- 'N16_2'
observed <- read.csv(file.path(prepared_dir, 'tables/observed_tissue_proportions.csv'))
{
  a <- observed[observed$analysis_unit %in% c(adult, 'N22'), c('tissue', 'specific_percent', 'analysis_unit')]
  a$series <- ifelse(a$analysis_unit == 'N22', 'Embryonic original', 'Adult original')
  b <- means[means$analysis_unit == adult, ]
  b$series <- 'Adult normalized'
  bars <- rbind(a[, c('tissue', 'series', 'specific_percent')], b[, c('tissue', 'series', 'specific_percent')])
  pts <- points[points$analysis_unit == adult, ]
  pts$series <- 'Adult normalized'
  pts <- rbind(a[, c('tissue', 'series', 'specific_percent')], pts[, c('tissue', 'series', 'specific_percent')])
  bars$tissue <- factor(bars$tissue, levels = tissues)
  bars$series <- factor(bars$series, levels = c('Embryonic original', 'Adult original', 'Adult normalized'))
  pts$tissue <- factor(pts$tissue, levels = tissues)
  pts$series <- factor(pts$series, levels = levels(bars$series))
  write.csv(bars, file.path(result_dir, 'tables', paste0('ED2g_pgen_first_', adult, '_bar_values.csv')), row.names = FALSE)
  p <- ggplot(bars, aes(tissue, specific_percent, fill = series)) +
    geom_col(position = position_dodge(.82), width = .64, colour = 'black', linewidth = .4) +
    geom_point(data = pts, aes(fill = series), position = position_jitterdodge(jitter.width = .08, dodge.width = .82, seed = 20261002), shape = 21, size = 2, stroke = .45) +
    scale_fill_manual(values = c('Embryonic original' = '#CCCCCC', 'Adult original' = '#9999A2', 'Adult normalized' = 'white')) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
    labs(title = paste(adult, 'vs N22; Pgen-first richness sensitivity'), x = NULL, y = 'Tissue-specific HSPC clones (%)', fill = NULL) +
    theme_classic(base_family = 'Arial', base_size = 11) + theme(axis.text = element_text(colour = 'black'), legend.position = 'right')
  stem <- paste0('ED2_g_pgen_first_candidate_', adult)
  ggsave(file.path(result_dir, 'plots', paste0(stem, '.png')), p, width = 5.5, height = 3.8, dpi = 300, bg = 'white')
  ggsave(file.path(result_dir, 'plots', paste0(stem, '.pdf')), p, width = 5.5, height = 3.8, device = cairo_pdf, bg = 'white')
  write.csv(pts, file.path(result_dir, 'tables', paste0(stem, '_points.csv')), row.names = FALSE)
}
