# Plot scatter plots between signature scores and specific cell cycles signatures
source("requirements.R")

# Input paths
cc_signatures_inputpath <- "extdata/signatures/cell_cycle.rds"
tcga_inputpath <- "extdata/tcga-coad/counts_norm_coad_patients.rds"
coad_signatures_inputpath <- "extdata/signatures/signatures_coad.rds"

# Load data
train <- readRDS(tcga_inputpath)
cc_signatures <- readRDS(cc_signatures_inputpath)
coad_signatures <- readRDS(coad_signatures_inputpath)

# Define signatures
signatures <- c(cc_signatures, coad_signatures)

# Calculate signature´s activity
plagePar <- plageParam(t(train), signatures)
gsva <- gsva(plagePar, verbose = FALSE)
gsva <- as.data.frame(t(gsva))

# Plot
s1 = ggscatter(gsva, x = 'ours', y = 'G1Phase', add = 'reg.line', conf.int = T,
               add.params = list(color = viridis(3)[1], fill = 'lightgray')) +
  stat_cor(method = 'pearson')
s2 = ggscatter(gsva, x = 'ours', y = 'M.prophase', add = 'reg.line', conf.int = T,
               add.params = list(color = viridis(3)[1], fill = 'lightgray')) +
  stat_cor(method = 'pearson')
s3 = ggscatter(gsva, x = 'ours', y = 'M.prometaphase', add = 'reg.line', conf.int = T,
               add.params = list(color = viridis(3)[1], fill = 'lightgray')) +
  stat_cor(method = 'pearson')
s4 = ggscatter(gsva, x = 'ours', y = 'M.metaphase.anaphase', add = 'reg.line', conf.int = T,
               add.params = list(color = viridis(3)[1], fill = 'lightgray')) +
  stat_cor(method = 'pearson')
s5 = ggscatter(gsva, x = 'ours', y = 'M.telophase.cytokinesis', add = 'reg.line', conf.int = T,
                add.params = list(color = viridis(3)[1], fill = 'lightgray')) +
  stat_cor(method = 'pearson')

ggarrange(s1, s2, s3, s4, s5, nrow = 3, ncol = 2)
