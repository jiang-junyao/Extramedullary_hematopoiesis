if (.Platform$OS.type == 'windows') invisible(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
library(ggplot2)
polylox_root <- normalizePath(getwd(), winslash = '/', mustWork = TRUE)
prepared_dir <- file.path(polylox_root, 'prepared/simulation')
data_dir <- file.path(polylox_root, 'data/simulation')
result_dir <- file.path(polylox_root, 'results/Response_Figure_1_2')
for (folder in c('plots', 'tables', 'qc')) dir.create(file.path(result_dir, folder), recursive = TRUE, showWarnings = FALSE)
tissues <- c('Leg', 'Spine', 'Spleen', 'Liver')

# Only Response Figure 1.2 a-b; no ED2 or other figure output.
summary <- read.csv(file.path(prepared_dir, 'tables/mixing_summary.csv'))
real <- read.csv(file.path(prepared_dir, 'tables/observed_tissue_proportions.csv'))
id <- 'N16_2'
for (metric in c('shared_percent', 'specific_percent')) {
  z <- summary[summary$analysis_unit == id & summary$source == 'historical_100_replicates' & summary$metric == metric, ]
  z$series <- 'Simulation'
  z$lower <- z$q025_percent
  z$upper <- z$q975_percent
  obs <- real[real$analysis_unit == id, ]
  p.data <- rbind(data.frame(tissue = obs$tissue, series = 'Real data', mean_percent = obs[[metric]], lower = obs[[metric]], upper = obs[[metric]]),
    z[, c('tissue', 'series', 'mean_percent', 'lower', 'upper')])
  p.data$tissue <- factor(p.data$tissue, levels = tissues)
  p.data$series <- factor(p.data$series, levels = c('Real data', 'Simulation'))
  p <- ggplot(p.data, aes(tissue, mean_percent, fill = series)) +
    geom_col(position = position_dodge(.72), width = .58, colour = 'black', linewidth = .45) +
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(.72), width = .13, linewidth = .35) +
    scale_fill_manual(values = c('Real data' = 'black', Simulation = 'white')) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = expansion(mult = c(0, .01))) +
    labs(title = paste(id, 'observed vs uniform mixing'), x = NULL, y = if (metric == 'shared_percent') 'Shared HSPC clones (%)' else 'Tissue-specific HSPC clones (%)', fill = NULL) +
    theme_classic(base_family = 'Arial', base_size = 11) +
    theme(axis.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), legend.position = 'right')
  stem <- paste0('ED2_', if (metric == 'shared_percent') 'a' else 'b', '_', id, '_historical_quantile_interval')
  ggsave(file.path(result_dir, 'plots', paste0(stem, '.png')), p, width = 5.5, height = 3.8, dpi = 300, bg = 'white')
  ggsave(file.path(result_dir, 'plots', paste0(stem, '.pdf')), p, width = 5.5, height = 3.8, device = cairo_pdf, bg = 'white')
  write.csv(p.data, file.path(result_dir, 'tables', paste0(stem, '_plot_data.csv')), row.names = FALSE)
}
