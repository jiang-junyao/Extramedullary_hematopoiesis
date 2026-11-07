library(data.table)

preg <- fread('prepared/tables/Pregnancy_SP_vs_BM_all_genes.tsv.gz')
epo_tissue <- fread('prepared/tables/EPO_SP_vs_BM_all_genes.tsv.gz')
epo <- fread('prepared/tables/EPO_vs_CTRL_SP_all_genes.tsv.gz')
emp9 <- fread('prepared/tables/EMP9_vs_Pregnancy_SP_all_genes.tsv.gz')
exclusions <- fread('data/annotation/mature_lineage_exclusions.tsv')

preg_up <- preg[is_DEG == TRUE & log2FoldChange > 0]
setorder(preg_up, -log2FoldChange, padj, pvalue, gene)
top200 <- head(preg_up, 200)
top200[, signature_order := .I]
fwrite(top200, 'prepared/tables/Pregnancy_SP_high_Top200.tsv', sep = '\t')
fwrite(epo[is_DEG == TRUE & log2FoldChange > 0], 'prepared/tables/EPO_up_signature.tsv', sep = '\t')

# Pregnancy heatmap: signature overlap, ordered by Pregnancy fold change.
overlap <- merge(top200, epo_tissue[is_DEG == TRUE & log2FoldChange > 0], by = 'gene', suffixes = c('_Pregnancy', '_EPO'))
setorder(overlap, signature_order)
overlap[, excluded := gene %in% exclusions$gene | grepl('^(Igh|Igk|Igl|Tcr|Trbc|Trdc)', gene)]
fwrite(overlap, 'prepared/tables/Pregnancy_EPO_overlap_candidates.tsv', sep = '\t')
fwrite(overlap[excluded == FALSE], 'prepared/tables/Pregnancy_EPO_heatmap_genes.tsv', sep = '\t')

# Treatment heatmap: strict EPO-up AND EMP9-down intersection, EPO-priority order.
intersection <- merge(epo, emp9, by = 'gene', suffixes = c('_EPO', '_EMP9'))
intersection <- intersection[is_DEG_EPO == TRUE & log2FoldChange_EPO > 0 & is_DEG_EMP9 == TRUE & log2FoldChange_EMP9 < 0]
intersection[, excluded := gene %in% exclusions$gene | grepl('^(Igh|Igk|Igl|Tcr)', gene) | grepl('Rik$', gene) | gene == 'Gm15290']
setorder(intersection, -log2FoldChange_EPO, padj_EPO, padj_EMP9, gene)
fwrite(intersection, 'prepared/tables/EPO_EMP9_overlap_candidates.tsv', sep = '\t')
strict <- intersection[excluded == FALSE]
strict[, display_role := 'Strict EPO-up and EMP9-down DEG']
focus <- merge(epo[gene %in% c('Klf1', 'Epor')], emp9[gene %in% c('Klf1', 'Epor')], by = 'gene', suffixes = c('_EPO', '_EMP9'))
focus[, display_role := 'Added display gene; not a strict dual-cutoff DEG']
display <- rbindlist(list(strict, focus), fill = TRUE)
setorder(display, -log2FoldChange_EPO, padj_EPO, padj_EMP9, gene)
fwrite(strict, 'prepared/tables/EPO_EMP9_heatmap_strict_genes.tsv', sep = '\t')
fwrite(display, 'prepared/tables/EPO_EMP9_heatmap_genes.tsv', sep = '\t')

summary <- rbindlist(lapply(fread('data/comparisons.tsv')$comparison, function(name) {
  d <- fread(paste0('prepared/tables/', name, '_all_genes.tsv.gz'))
  data.table(comparison = name, Up = sum(d$direction == 'Up'), Down = sum(d$direction == 'Down'), Total = sum(d$is_DEG))
}))
fwrite(summary, 'prepared/tables/DEG_summary.tsv', sep = '\t')
