# Distribution of the yaccs PLAGE score across TCGA-COAD patients
# Dashed line: mean score. Solid line: low/high cutoff used downstream.
source("requirements.R")

# Parameters
# "cutoff"   -> boundary between the two Mclust components (defines low/high)
# "highmean" -> mean of the "high" component
line_type <- "cutoff"

fill_colour <- "#8B7FB0"
line_colour <- "#5B2C8D"

# Inputpath
tcga_inputpath <- "extdata/tcga-coad/counts_norm_coad_patients.rds"
signatures_inputpath <- "extdata/signatures/signatures_coad.rds"

# Outputpath
outputpath <- "Figures/Figure4/supp_distribution.png"

# Load data
signatures <- readRDS(signatures_inputpath)
train <- readRDS(tcga_inputpath)

# Calculate pathway activity (same call as workflow/08_pathway_enrichment.R)
plagePar <- plageParam(
    exprData = t(train),
    geneSets = list(yaccs = signatures[["ours"]])
)

act <- gsva(plagePar, verbose = FALSE) %>%
t() %>%
as.data.frame()

# Cluster patients based on yaccs score (same call as the workflow)
clust <- Mclust(
    act$yaccs,
    modelNames = "E",
    na.rm = TRUE,
    G = 2
)
act$classif <- ifelse(clust$classification == 2, "high", "low")

# Vertical lines
mean_score <- mean(act$yaccs)
cutoff <- mean(c(
    max(act$yaccs[act$classif == "low"]),
    min(act$yaccs[act$classif == "high"])
))
high_mean <- max(clust$parameters$mean)
vline <- if (line_type == "highmean") high_mean else cutoff

# Plot
p <- ggplot(act, aes(x = yaccs)) +
    geom_density(fill = fill_colour, colour = "black", alpha = 0.85) +
    geom_rug(sides = "b", length = unit(0.04, "npc")) +
    geom_vline(xintercept = mean_score, linetype = "dashed", colour = "black") +
    geom_vline(xintercept = vline, colour = line_colour, linewidth = 1.3) +
    labs(x = NULL, y = "density") +
    theme_classic(base_size = 11)

ggsave(outputpath, p, width = 4.5, height = 4, dpi = 200)

cat("n low / high        :", sum(act$classif == "low"), "/", sum(act$classif == "high"), "\n")
cat("mean (dashed)       :", signif(mean_score, 4), "\n")
cat("cutoff low/high     :", signif(cutoff, 4), "\n")
cat("high component mean :", signif(high_mean, 4), "\n")
cat("solid line drawn at :", signif(vline, 4), paste0("(line_type = '", line_type, "')"), "\n")
