library(data.table)
library(ggplot2)
library(ggtext)
library(patchwork)
library(pheatmap)
library(grid)
library(ragg)

samples <- fread('data/samples.tsv')
membership <- fread('data/comparison_samples.tsv')
comparisons <- fread('data/comparisons.tsv')

save_plot <- function(p, name, width, height) {
  ggsave(file.path(output_dir, 'plots', paste0(panel_names[[name]], '.png')), p, width = width, height = height, dpi = 300, bg = 'white', device = ragg::agg_png)
  ggsave(file.path(output_dir, 'plots', paste0(panel_names[[name]], '.pdf')), p, width = width, height = height, bg = 'white', device = cairo_pdf, family = 'Arial')
}

# Heatmap expression: each source cohort is standardized separately across its samples.
read_tpm <- function(dataset) {
  x <- fread(paste0('prepared/tables/', dataset, '_TPM.tsv.gz'))
  m <- as.matrix(x[, -1])
  rownames(m) <- x$gene
  m
}
tpm <- lapply(c('Paper', 'EPO0804', 'PregEMP90807'), read_tpm)
names(tpm) <- c('Paper', 'EPO0804', 'PregEMP90807')

heatmap_matrix <- function(genes, md, name) {
  pieces <- lapply(unique(md$dataset), function(source_dataset) {
    log_tpm <- log2(tpm[[source_dataset]][genes, md[dataset == source_dataset, sample_key], drop = FALSE] + 1)
    list(log = log_tpm, z = t(scale(t(log_tpm))))
  })
  log_tpm <- do.call(cbind, lapply(pieces, '[[', 'log'))[, md$sample_key, drop = FALSE]
  z <- do.call(cbind, lapply(pieces, '[[', 'z'))[, md$sample_key, drop = FALSE]
  fwrite(data.table(gene = genes, log_tpm), file.path(output_dir, 'tables', paste0(panel_names[[name]], '_log2TPM_plus1.tsv')), sep = '\t')
  fwrite(data.table(gene = genes, z), file.path(output_dir, 'tables', paste0(panel_names[[name]], '_Zscore.tsv')), sep = '\t')
  z[z < -2] <- -2
  z[z > 2] <- 2
  z
}

draw_heatmap <- function(m, md, name, highlight) {
  width <- if (ncol(m) > 12) 10.5 else 9
  height <- if (ncol(m) > 12) 5.6 else 7.4
  ragg::agg_png(file.path(output_dir, 'plots', paste0(panel_names[[name]], '.png')), width = width, height = height, units = 'in', res = 300)
  annotation <- data.frame(condition = md$display_group, row.names = md$sample_key)
  gaps <- which(head(md$display_group, -1) != tail(md$display_group, -1))
  p <- pheatmap(m, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, show_rownames = TRUE,
    labels_row = ifelse(rownames(m) %in% highlight, rownames(m), ''), annotation_col = annotation,
    annotation_names_col = FALSE, annotation_legend = FALSE,
    annotation_colors = list(condition = c('Pregnancy BM' = '#6DC17E', 'Pregnancy SP' = '#F59699', 'EPO BM' = '#6DC17E',
      'EPO SP' = '#F59699', 'Saline' = '#5B8DB8', 'EPO' = '#FBD6D7', 'Pregnancy' = '#D982A6', 'Pregnancy + EMP9' = '#E3DCC2')),
    color = colorRampPalette(c('#0057B8', '#F7F7F7', '#D7191C'))(100), breaks = seq(-2, 2, length.out = 101),
    border_color = 'white', gaps_col = gaps, fontsize = 12, fontsize_row = 13, silent = TRUE)
  gt <- p$gtable
  gt <- gtable::gtable_add_rows(gt, unit(0.42, 'in'), pos = 0)
  matrix_col <- gt$layout$l[gt$layout$name == 'matrix']
  matrix_row <- gt$layout$t[gt$layout$name == 'matrix']
  gt$heights[matrix_row] <- gt$heights[matrix_row] - unit(0.42, 'in')
  groups <- unique(md$display_group)
  centers <- vapply(groups, function(g) mean(which(md$display_group == g)) / ncol(m), numeric(1))
  group_labels <- sub('Pregnancy ', 'Preg ', groups, fixed = TRUE)
  gt <- gtable::gtable_add_grob(gt, textGrob(group_labels, x = centers, y = 0.35, gp = gpar(fontfamily = 'Arial', fontsize = 10)), t = 1, l = matrix_col, name = 'group_labels')
  idx <- which(gt$layout$name == 'row_names')
  labels <- gt$grobs[[idx]]$label
  selected <- which(labels != '')
  origins <- convertY(gt$grobs[[idx]]$y, 'npc', valueOnly = TRUE)
  new_y <- origins[selected]
  gap <- min(0.043, 0.9 / max(1, length(selected) - 1))
  if (length(new_y) > 1) for (j in 2:length(new_y)) new_y[j] <- min(new_y[j], new_y[j - 1] - gap)
  if (length(new_y) && min(new_y) < 0.02) new_y <- new_y + 0.02 - min(new_y)
  gt$grobs[[idx]] <- grobTree(
    segmentsGrob(x0 = unit(0, 'npc'), x1 = unit(0.06, 'npc'), y0 = unit(origins[selected], 'npc'), y1 = unit(new_y, 'npc'), gp = gpar(col = 'grey45', lwd = 0.6)),
    textGrob(labels[selected], x = unit(0.075, 'npc'), y = unit(new_y, 'npc'), hjust = 0,
      gp = gpar(fontfamily = 'Arial', fontsize = 13, fontface = ifelse(labels[selected] %in% c('Klf1', 'Gata1', 'Epor'), 'bold.italic', 'italic'),
        col = ifelse(labels[selected] %in% c('Klf1', 'Gata1', 'Epor'), '#C51B32', '#222222'))))
  grid.newpage(); grid.draw(gt); dev.off()
  cairo_pdf(file.path(output_dir, 'plots', paste0(panel_names[[name]], '.pdf')), width = width, height = height, family = 'Arial')
  grid.newpage(); grid.draw(gt); dev.off()
}
