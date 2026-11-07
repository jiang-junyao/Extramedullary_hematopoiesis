# Run from this repository root in RStudio or with Rscript.
dir.create('prepared/tables', recursive = TRUE, showWarnings = FALSE)
source('R/prepare_expression.R')
source('R/prepare_DEG.R')
source('R/prepare_gene_sets.R')
source('R/prepare_ORA.R')
source('R/prepare_GSEA.R')
writeLines(capture.output(sessionInfo()), 'prepared/sessionInfo.txt')
