# Prepare shared data only. Run once from the package root; no figure is drawn.
if (.Platform$OS.type == 'windows') invisible(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
source('R/polylox_core.R', encoding = 'UTF-8')

local(source('R/prepare_tissue.R', local = TRUE, encoding = 'UTF-8'))
local(source('R/prepare_simulation.R', local = TRUE, encoding = 'UTF-8'))
local(source('R/prepare_output.R', local = TRUE, encoding = 'UTF-8'))

dir.create('logs', showWarnings = FALSE)
capture.output(sessionInfo(), file = 'logs/preparation_sessionInfo.txt')
