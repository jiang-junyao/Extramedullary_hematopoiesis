# Response Figure 3.1: panel g source replay; caption n conflict retained
# Candidate/source-replay status is unchanged; run 00_prepare_shared_data.R first.
library(dplyr)
library(tidyr)
library(ggplot2)

input <- 'prepared/tissue/prepared_tissue.rds'
out <- 'results/Response_Figure_3_1'
plot_dir <- file.path(out, 'plots')
table_dir <- file.path(out, 'tables')
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
prepared <- readRDS(input)
tissue_order <- c('Leg', 'Spine', 'Spleen', 'Liver')
tissue_colors <- c(Leg = '#8BC98D', Spine = '#77A8E8', Spleen = '#E9D9D7', Liver = '#E7B164')

# Formal Source Data replay: preserve original values, without new tests.
job <- 'ResponseFig3_1g'
panel <- 'Fig2i'
dat <- prepared$formal[prepared$formal$panel == panel, ]
dat$condition <- factor(dat$condition, levels = c('Non-preg', 'Pregnancy', 'EPO'))
dat$tissue <- factor(dat$tissue, levels = tissue_order)
means <- dat %>% group_by(condition) %>% summarise(value = mean(specific_percent), .groups = 'drop')
p <- ggplot(dat, aes(condition, specific_percent)) +
  geom_col(data = means, aes(y = value), fill = '#CCCCCC', width = 0.65) +
  geom_point(aes(fill = tissue), shape = 21, position = position_jitter(width = 0.15, height = 0, seed = 20261102), size = 3, stroke = 0.5) +
  scale_fill_manual(values = tissue_colors, drop = TRUE) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(title = paste(job, 'formal Source Data replay'), subtitle = 'Original values retained; no new inferential statistics', x = NULL, y = 'Tissue-specific HSPC clones (%)', fill = 'Tissue') +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(color = 'black'), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot_dir, paste0(job, '_formal_source_replay.png')), p, width = 7.2, height = 5.1, dpi = 300)
ggsave(file.path(plot_dir, paste0(job, '_formal_source_replay.pdf')), p, width = 7.2, height = 5.1, device = cairo_pdf)
write.table(dat, file.path(table_dir, paste0(job, '_points.tsv')), sep = '\t', row.names = FALSE, quote = FALSE)
write.table(means, file.path(table_dir, paste0(job, '_bar_means.tsv')), sep = '\t', row.names = FALSE, quote = FALSE)
