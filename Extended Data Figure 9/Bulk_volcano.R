
library(ggplot2)
library(ggrepel)
library(dplyr)

data <- as.data.frame(res)
data$gene <- rownames(data)


data <- data[!is.na(data$pvalue), ]

data$logP <- -log10(data$pvalue)
data$is_capped <- data$logP > 30 
data$logP[data$logP > 20] <- 20
logFC_cutoff <- 1.5    #
pval_cutoff <- 0.01    


data$group <- "Stable"
data$group[data$log2FoldChange >= logFC_cutoff & data$pvalue < pval_cutoff] <- "Up"
data$group[data$log2FoldChange <= -logFC_cutoff & data$pvalue < pval_cutoff] <- "Down"
data$group <- factor(data$group, levels = c("Up", "Stable", "Down"))

target_genes <- c("Gata1", "Klf1")

top_up <- data %>%
  filter(group == "Up") %>%
  arrange(pvalue) %>%
  slice_head(n = 10) %>%
  pull(gene)

top_down <- data %>%
  filter(group == "Down") %>%
  arrange(pvalue) %>%
  slice_head(n = 12) %>%
  pull(gene)


label_list <- unique(c(target_genes, top_up, top_down))


data$label <- ifelse(data$gene %in% label_list, data$gene, NA)


my_palette <- c("Up" = rgb(255/255,160/255,160/255),      # Firebrick (深红)
                "Stable" = "#D3D3D3",  # LightGray (浅灰)
                "Down" = rgb(144/255,191/255,249/255))    # DodgerBlue (深蓝)


p <- ggplot(data, aes(x = log2FoldChange, y = logP, color = group)) +


  geom_point(size = 1.5, alpha = 0.7) +


  scale_color_manual(values = my_palette) +


  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "grey40", size = 0.4) +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "grey40", size = 0.4) +


  geom_text_repel(aes(label = label),
                  max.overlaps = Inf,
                  box.padding = 0.6,
                  point.padding = 0.3,
                  size = 4.5,
                  fontface = "italic",
                  segment.size = 0.3,
                  segment.color = "black",
                  show.legend = FALSE) +


  labs(x = bquote(Log[2] ~ "Fold Change"),
       y = bquote(-Log[10] ~ italic("P") ~ "value"), # Y轴标签改为 P value
       title = "Volcano Plot",
       subtitle = paste0("pvalue < ", pval_cutoff, " & |log2FC| > ", logFC_cutoff)) +


  scale_y_continuous(limits = c(0, 22), expand = c(0, 0)) +


  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey50", size = 12),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )


print(p)


ggsave("/data/jiangjunyao/wenqian_smartseq/preg_Volcano_Plot_Top12_PerGroup_2.pdf", plot = p, width = 8, height = 6)
pdf('/data/jiangjunyao/wenqian_smartseq/preg_heatmap.pdf',width = 3,height = 4)
pheatmap(a1[rownames(a2),],cluster_rows = F,cluster_cols = F,show_rownames = F,
         color = colorRampPalette(heatmap_col)(50))
dev.off()
