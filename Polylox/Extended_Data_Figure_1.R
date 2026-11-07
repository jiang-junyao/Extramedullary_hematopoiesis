# Extended Data Figure 1: panel g candidate
# Candidate/source-replay status is unchanged; run 00_prepare_shared_data.R first.
library(dplyr)
library(tidyr)
library(ggplot2)

input <- 'prepared/tissue/prepared_tissue.rds'
out <- 'results/Extended_Data_Figure_1'
plot_dir <- file.path(out, 'plots')
table_dir <- file.path(out, 'tables')
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
prepared <- readRDS(input)
tissue_order <- c('Leg', 'Spine', 'Spleen', 'Liver')
tissue_colors <- c(Leg = '#8BC98D', Spine = '#77A8E8', Spleen = '#E9D9D7', Liver = '#E7B164')

# Candidate heatmaps for this figure only; no other record is plotted.
units <- c('N16_1')
for (unit in units) {
  plot_source <- prepared$row_percentages[[unit]]
  dat <- plot_source
  tissues <- intersect(tissue_order, names(plot_source))
  plot_data <- plot_source %>% pivot_longer(all_of(tissues), names_to = 'tissue', values_to = 'percentage') %>%
    mutate(tissue = factor(tissue, levels = tissues), y = nrow(dat) - row_order + 1)
  group_n <- table(dat$classification)
  group_n <- group_n[group_n > 0]
  labels <- data.frame(label = paste0(names(group_n), ' (n=', as.integer(group_n), ')'), y = nrow(dat) - cumsum(as.numeric(group_n)) + as.numeric(group_n) / 2 + 0.5)

  p <- ggplot(plot_data, aes(tissue, y, fill = percentage)) +
    geom_tile(width = 1, height = 1) +
    geom_text(data = labels, aes(x = length(tissues) + 0.6, y = y, label = label), inherit.aes = FALSE, hjust = 0, size = 3.4) +
    scale_fill_gradientn(colors = c('#F3F9FD', '#C7E2EF', '#7EBBD5', '#328BC0', '#176EAE', '#2A3967'),
      values = scales::rescale(c(0, 20, 40, 60, 80, 100)), limits = c(0, 100), name = 'Within-barcode\npercentage (%)') +
    scale_x_discrete(position = 'top', drop = FALSE) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0.5, length(tissues) + 2.2), clip = 'off') +
    labs(title = paste(unit, 'progenitor barcode heatmap'), subtitle = 'Historical-data reconstruction; formal panel assignment unconfirmed', x = NULL, y = NULL) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(size = 18, face = 'bold', hjust = 0.5), plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.text.x = element_text(size = 13, angle = 90, hjust = 0, vjust = 0.5), axis.text.y = element_blank(),
      axis.ticks = element_blank(), axis.line = element_blank(), plot.margin = margin(10, 115, 10, 10))
  ggsave(file.path(plot_dir, paste0(unit, '_distribution_candidate.png')), p, width = 9.2, height = 13.568, dpi = 200)
  ggsave(file.path(plot_dir, paste0(unit, '_distribution_candidate.pdf')), p, width = 9.2, height = 13.568, device = cairo_pdf)
  write.table(plot_source, file.path(table_dir, paste0(unit, '_heatmap_row_percentage.tsv')), sep = '\t', row.names = FALSE, quote = FALSE)
}
