# Shared numerical preparation only: no plotting or figure output.
if (.Platform$OS.type == 'windows') invisible(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
polylox_root <- normalizePath(getwd(), winslash = '/', mustWork = TRUE)
source(file.path(polylox_root, 'R/polylox_core.R'), local = TRUE, encoding = 'UTF-8')
data_dir <- file.path(polylox_root, 'data/simulation')
result_dir <- file.path(polylox_root, 'prepared/simulation')
for (folder in c('tables', 'qc')) dir.create(file.path(result_dir, folder), recursive = TRUE, showWarnings = FALSE)
ids <- c('N16_1', 'N16_2', 'N22')
tissues <- c('Leg', 'Spine', 'Spleen', 'Liver')
summarize_presence <- function(x) {
  detected <- x > 0
  breadth <- rowSums(detected)
  total <- colSums(detected)
  specific <- colSums(detected & (breadth == 1))
  data.frame(tissue = colnames(x), detected = total, specific = specific, shared = total - specific,
    specific_percent = ifelse(total > 0, 100 * specific / total, NA),
    shared_percent = ifelse(total > 0, 100 * (total - specific) / total, NA), row.names = NULL)
}

source(file.path(polylox_root, 'R/prepare_simulation_counts.R'), local = TRUE, encoding = 'UTF-8')
source(file.path(polylox_root, 'R/prepare_simulation_controls.R'), local = TRUE, encoding = 'UTF-8')
