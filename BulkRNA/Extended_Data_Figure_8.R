# Run 00_prepare_shared_data.R once before this figure.
source('R/plot_common.R')
output_dir <- 'results/Extended_Data_Figure_8'
dir.create(file.path(output_dir, 'plots'), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, 'tables'), recursive = TRUE, showWarnings = FALSE)
panel_names <- c(Pregnancy_SP_vs_BM_volcano = 'b_Pregnancy_SP_vs_BM_volcano',
  Pregnancy_signature_in_Public_GSEA = 'c_Pregnancy_signature_in_Public_GSEA',
  Pregnancy_signature_in_EPO_GSEA = 'e_Pregnancy_signature_in_EPO_GSEA',
  Pregnancy_EPO_heatmap = 'f_Pregnancy_EPO_heatmap')
selected_comparisons <- c('Pregnancy_SP_vs_BM')
selected_gsea <- c('Pregnancy_signature_in_Public', 'Pregnancy_signature_in_EPO')

# Volcanoes: numerical classification uses original statistics; caps affect display only.
labels <- fread('data/annotation/volcano_labels.tsv')
volcano_settings <- data.table(
  comparison = c('Pregnancy_SP_vs_BM', 'EPO_vs_CTRL_SP', 'EMP9_vs_Pregnancy_SP'),
  title = c('Pregnant SP HSC vs pregnant BM HSC', 'EPO vs saline SP HSC', 'Pregnancy + EMP9 vs Pregnancy SP HSC'),
  y_max = c(20, 15, 43), y_cap = c(20, 14, 40),
  down = c('#9CCBA7', '#A6A6A6', '#CBCCCD'), up = c('#F3A0A3', '#D9BEAA', '#E3DCC2'))
volcano_settings <- volcano_settings[comparison %in% selected_comparisons]
for (i in seq_len(nrow(volcano_settings))) {
  v <- volcano_settings[i]
  d <- fread(paste0('prepared/tables/', v$comparison, '_all_genes.tsv.gz'))[is.finite(log2FoldChange)]
  d[, `:=`(x = pmax(-10, pmin(10, log2FoldChange)), y = pmin(v$y_cap, -log10(pmax(pvalue, .Machine$double.xmin))))]
  lab <- labels[comparison == v$comparison]
  colors <- c(Down = v$down, 'Not DEG' = '#CFCFCF', Up = v$up)
  p <- ggplot(d, aes(x, y, color = direction)) +
    geom_point(size = 1.1, na.rm = TRUE) +
    geom_vline(xintercept = c(-1, 1) * log2(1.5), linetype = 'dotted', linewidth = 0.5, color = 'black') +
    geom_hline(yintercept = -log10(0.01), linetype = 'dotted', linewidth = 0.5, color = 'black') +
    geom_text(data = lab, aes(label_x, label_y, label = gene, fontface = fontface, color = I(label_color)), inherit.aes = FALSE, size = 5, family = 'Arial') +
    scale_color_manual(values = colors, guide = 'none') +
    scale_x_continuous(breaks = seq(-10, 10, 5)) +
    coord_cartesian(xlim = c(-10, 10), ylim = c(0, v$y_max), expand = FALSE) +
    labs(title = v$title, x = expression(Log[2]~'Fold Change'), y = expression('-' * Log[10]~italic(P)~'value')) +
    theme_classic(base_family = 'Arial', base_size = 18) +
    theme(plot.title = element_text(size = 18, hjust = 0.5), axis.text = element_text(color = 'black'), plot.margin = margin(12, 15, 12, 12))
  save_plot(p, paste0(v$comparison, '_volcano'), 8.2, 6.2)
  fwrite(d, file.path(output_dir, 'tables', paste0(panel_names[[paste0(v$comparison, '_volcano')]], '_points.tsv.gz')), sep = '\t')
}

# GSEA: labels always read the newly calculated NES and FDR; never hard-code signs.
for (name in selected_gsea) {
  d <- fread(paste0('prepared/tables/', name, '_curve.tsv.gz'))
  s <- fread(paste0('prepared/tables/', name, '_GSEA.tsv'))
  fwrite(d, file.path(output_dir, 'tables', paste0(panel_names[[paste0(name, '_GSEA')]], '_curve.tsv.gz')), sep = '\t')
  fwrite(s, file.path(output_dir, 'tables', paste0(panel_names[[paste0(name, '_GSEA')]], '_statistics.tsv')), sep = '\t')
  exponent <- floor(log10(s$padj))
  fdr <- paste0('FDR: ', formatC(s$padj / 10^exponent, format = 'f', digits = 2), ' × 10<sup>', exponent, '</sup>')
  nes <- paste0('NES: ', formatC(s$NES, format = 'f', digits = 2))
  yr <- range(c(0, d$running_ES))
  span <- diff(yr)
  n <- nrow(d)
  p_curve <- ggplot(d, aes(rank, running_ES)) +
    geom_hline(yintercept = 0, color = 'grey70', linewidth = 0.3) +
    geom_line(linewidth = 0.7, color = 'black') +
    annotate('text', x = n * 0.045, y = yr[1] + span * 0.30, label = nes, hjust = 0, family = 'Arial', size = 5.5) +
    ggtext::geom_richtext(data = data.frame(x = n * 0.045, y = yr[1] + span * 0.12, label = fdr),
      aes(x, y, label = label), inherit.aes = FALSE, hjust = 0, size = 5.5, family = 'Arial', fill = NA, label.color = NA, label.padding = unit(0, 'pt')) +
    scale_x_continuous(limits = c(1, n), expand = c(0, 0)) +
    scale_y_continuous(limits = c(yr[1] - span * 0.07, yr[2] + span * 0.10), expand = c(0, 0)) +
    theme_void() + theme(panel.border = element_rect(fill = NA, color = 'grey35', linewidth = 0.5), plot.margin = margin(0, 4, 0, 4))
  p_hits <- ggplot(d[hit == TRUE], aes(rank)) +
    geom_segment(aes(xend = rank, y = 0, yend = 1), linewidth = 0.3) +
    scale_x_continuous(limits = c(1, n), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() + theme(panel.border = element_rect(fill = NA, color = 'grey35', linewidth = 0.5), plot.margin = margin(0, 4, 0, 4))
  gradient <- data.frame(x = seq_len(500), value = seq(-1, 1, length.out = 500))
  p_gradient <- ggplot(gradient, aes(x, 1, fill = value)) + geom_raster() +
    scale_fill_gradient2(low = '#3045A5', mid = 'white', high = '#D7191C', guide = 'none') +
    scale_x_continuous(expand = c(0, 0)) + theme_void() + theme(plot.margin = margin(0, 4, 0, 4))
  p_labels <- ggplot() +
    annotate('text', x = 0, y = 0.7, label = s$left, hjust = 0, vjust = 1, color = '#3045A5', size = 5, family = 'Arial') +
    annotate('text', x = 1, y = 0.7, label = s$right, hjust = 1, vjust = 1, color = '#D7191C', size = 5, family = 'Arial') +
    xlim(0, 1) + ylim(0, 1) + theme_void() + theme(plot.margin = margin(0, 4, 0, 4))
  p <- (p_curve / p_hits / p_gradient / p_labels) + plot_layout(heights = c(4.5, 1.25, 0.18, 1.25))
  save_plot(p, paste0(name, '_GSEA'), 5.8, 5.0)
}

preg_genes <- fread('prepared/tables/Pregnancy_EPO_heatmap_genes.tsv')$gene
md <- samples[condition %in% c('Pregnancy', 'EPO') & dataset != 'Public']
md[, display_group := paste(ifelse(condition == 'EPO', 'EPO', 'Pregnancy'), tissue)]
md[, group_order := match(display_group, c('Pregnancy BM', 'Pregnancy SP', 'EPO BM', 'EPO SP'))]
setorder(md, group_order, dataset, sample_id)
m <- heatmap_matrix(preg_genes, md, 'Pregnancy_EPO_heatmap')
highlight <- fread('data/annotation/heatmap_highlights.tsv')$gene
draw_heatmap(m, md, 'Pregnancy_EPO_heatmap', highlight)
