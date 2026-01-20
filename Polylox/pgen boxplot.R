library(tidyverse)
data_list = readRDS('E:\\wenqian polylox\\20250918_wenqian_barcode\\qc2_5.rds')
for (i in names(data_list)) {
  data_use = data_list[[i]]
  pgen = read.table(paste0('E:\\wenqian polylox\\20250825_wenqian\\QC_2_pgen/',
                           i,'_pgen.txt'),header = T)
  rownames(pgen) = pgen$Code
  barcode_summary = data.frame(rownames(data_use),rowSums(data_use))
  barcode_summary$pgen = pgen[barcode_summary[,1],2]
  colnames(barcode_summary) = c('barcode','row_count','pgen')
  total_count <- sum(barcode_summary$row_count)
  barcode_summary$max_barcode_count_freq <- barcode_summary$row_count / total_count
  barcode_summary$log_max_barcode_count_freq <- log10(barcode_summary$max_barcode_count_freq)

  barcode_summary <- barcode_summary[barcode_summary$pgen > 0, ]
  barcode_summary$log_pgen <- log10(barcode_summary$pgen)

  breaks <- seq(-8, -1, by = 0.5)
  barcode_summary$pgen_bin <- cut(
    barcode_summary$log_pgen,
    breaks = breaks,
    include.lowest = TRUE,
    right = FALSE,
    labels = paste0(breaks[-1] - 0.5)
  )

  barcode_summary$pgen_bin <- factor(barcode_summary$pgen_bin, levels = paste0(-breaks[-1] - 8.5))

  p1=ggplot(barcode_summary, aes(x = pgen_bin, y = log_max_barcode_count_freq)) +
    geom_boxplot() +
    xlab("Log probability of generation") +
    ylab("Log frequency") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )+theme(text = element_text(size=16))
  print(p1)
  ggsave(paste0('E:\\wenqian polylox\\20250423 new fate bias\\fig pgen boxplot/',i,
                '.pdf'),width = 10,height = 8)
}



