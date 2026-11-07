# Original Response3.2 b-e: Non-preg only, without exporting the other condition.
# Run 00_prepare_shared_data.R once at the package root, then select data, plot and save here.
library(dplyr)
library(ggplot2)

output_result <- 'results/Response_Figure_3_2'
cond <- 'Non-preg'
reference_figure <- 'Response3.2'
for (folder in c('plots', 'tables', 'qc')) dir.create(file.path(output_result, folder), recursive = TRUE, showWarnings = FALSE)
prepared <- readRDS('prepared/output/tables/output_data.rds')
barcode <- prepared$barcode
statistics <- prepared$statistics
# Save this Figure's analysis values; axis limits still use both conditions.
write.csv(barcode[barcode$condition == cond, ], file.path(output_result, 'tables/barcode_output.csv'), row.names = FALSE)
write.csv(statistics[statistics$condition == cond, ], file.path(output_result, 'tables/wilcoxon.csv'), row.names = FALSE)
write.csv(barcode[barcode$condition == cond & barcode$active, ], file.path(output_result, 'tables/active_boxplot_observations.csv'), row.names = FALSE)

colors <- c('Non-shared' = '#BE1D2C', Shared = '#000000')
scientific_labels <- function(x) as.expression(lapply(x, function(v) substitute(10^E, list(E = round(log10(v))))))
scatter_theme <- theme_classic(base_size = 16) +
  theme(axis.line = element_line(linewidth = .9), axis.ticks = element_line(linewidth = .9), axis.text = element_text(color = 'black', size = 15),
    axis.title = element_text(size = 18), plot.title = element_text(size = 19, hjust = .5), legend.position = 'right', plot.margin = margin(10, 12, 10, 10))
box_theme <- theme_classic(base_size = 18) +
  theme(axis.line = element_line(linewidth = .9), axis.ticks = element_line(linewidth = .9), axis.ticks.length = grid::unit(.16, 'cm'),
    axis.text.x = element_text(size = 18, color = 'black', angle = 45, hjust = 1), axis.text.y = element_text(size = 16, color = 'black'),
    axis.title.y = element_text(size = 20, margin = margin(r = 12)), plot.title = element_text(size = 16, hjust = .5), plot.title.position = 'plot', legend.position = 'none')

panel_index <- axes <- box_statistics <- scatter_display <- list()
for (version_id in c('B', 'All3T6')) {
  scatter_panel <- if (version_id == 'B') 'd' else 'b'
  box_panel <- if (version_id == 'B') 'e' else 'c'
  scope_label <- if (version_id == 'B') 'Spleen only' else 'Leg, Spine and Spleen'
  all_data <- barcode[barcode$version == version_id, ]
  # True zeros stay zero in the saved biological table and only move in display_y.
  zero_floor <- min(all_data$output_sum[all_data$active]) / 2
  x_limits <- range(all_data$x) * c(1 / 1.4, 1.4)
  y_limits <- c(zero_floor / 1.4, max(all_data$output_sum) * 1.4)
  box_limits <- range(all_data$output_sum[all_data$active]) * c(1 / 1.5, 1.5)
  d <- all_data[all_data$condition == cond, ]
  d <- d[order(d$origin == 'Non-shared'), ]
  d$display_y <- ifelse(d$active, d$output_sum, zero_floor)
  d$origin <- factor(d$origin, levels = c('Shared', 'Non-shared'))
  d$status <- factor(ifelse(d$active, 'Active', 'No-output'), levels = c('Active', 'No-output'))
  scatter_display[[paste(version_id, cond)]] <- d
  p <- ggplot(d, aes(x, display_y, color = origin, shape = status)) +
    geom_point(size = 3.2, alpha = 1, stroke = .9) +
    scale_color_manual(values = colors, labels = c('Shared', 'Spleen-specific'), drop = FALSE) +
    scale_shape_manual(values = c(Active = 16, 'No-output' = 1), drop = FALSE) +
    scale_x_log10(limits = x_limits, labels = scientific_labels, breaks = 10^(-10:0)) +
    scale_y_log10(limits = y_limits, labels = scientific_labels, breaks = 10^(-10:0)) +
    labs(title = paste('Total splenic HSPC clones', cond, scope_label, sep = '\n'),
      x = 'Spleen HSPC frequency', y = 'Progeny frequency (sum of six profiles)', color = 'Tissue specificity', shape = 'Output activity')
  p <- p + scatter_theme

  b.data <- d[d$active, ]
  b.data$box_x <- ifelse(b.data$origin == 'Non-shared', 1.29, 1.71)
  stat <- statistics[statistics$version == version_id & statistics$condition == cond, ]
  p_label <- if (stat$p_unadjusted < .001) format(stat$p_unadjusted, scientific = TRUE, digits = 2) else formatC(stat$p_unadjusted, digits = 3, format = 'f')
  breaks <- 10^(-10:0)
  if (version_id == 'B') breaks <- breaks[log10(breaks) <= log10(box_limits[2]) - .08 * diff(log10(box_limits))]
  b <- ggplot(b.data, aes(box_x, output_sum, group = origin)) +
    geom_boxplot(aes(fill = origin), color = 'black', alpha = .28, width = .28, linewidth = .9, fatten = 2.625, outlier.shape = NA) +
    scale_fill_manual(values = colors) +
    scale_x_continuous(limits = c(1, 2), breaks = c(1.29, 1.71), labels = c('Spleen-specific', 'Shared'), expand = expansion(mult = 0)) +
    scale_y_log10(limits = box_limits, labels = scientific_labels, breaks = breaks) +
    labs(title = paste0('Active HSPC clones\n', cond, '\nP = ', p_label), x = NULL, y = 'Progeny frequency')
  b <- b + box_theme
  built <- ggplot_build(b)$data[[1]]
  built$outliers <- NULL
  box_statistics[[paste(version_id, cond)]] <- cbind(data.frame(version = version_id, condition = cond), built)
  axes[[paste(version_id, cond)]] <- data.frame(version = version_id, condition = cond, x_min = x_limits[1], x_max = x_limits[2],
    scatter_min = y_limits[1], scatter_max = y_limits[2], box_min = box_limits[1], box_max = box_limits[2], zero_floor)
  for (kind in c('scatter', 'boxplot')) {
    panel <- if (kind == 'scatter') scatter_panel else box_panel
    stem <- paste(reference_figure, panel, version_id, 'Sum', kind, sep = '_')
    plot <- if (kind == 'scatter') p else b
    width <- if (kind == 'scatter') 8.2 else 3.3
    height <- if (kind == 'scatter') 6.5 else 6.2 * 5 / 6
    ggsave(file.path(output_result, 'plots', paste0(stem, '.pdf')), plot, width = width, height = height, device = cairo_pdf, bg = 'white')
    ggsave(file.path(output_result, 'plots', paste0(stem, '.png')), plot, width = width, height = height, dpi = 300, bg = 'white')
    panel_index[[stem]] <- data.frame(reference_figure, panel, version = version_id, condition = cond, kind, stem, width, height)
  }
}
write.csv(bind_rows(panel_index), file.path(output_result, 'tables/panel_index.csv'), row.names = FALSE)
write.csv(bind_rows(axes), file.path(output_result, 'qc/axes.csv'), row.names = FALSE)
write.csv(bind_rows(box_statistics), file.path(output_result, 'qc/box_statistics.csv'), row.names = FALSE)
write.csv(bind_rows(scatter_display), file.path(output_result, 'tables/scatter_display.csv'), row.names = FALSE)

