# Scatter plots: drug AUC (y) vs signature score (x)
# for every drug within a given MOA
source("requirements.R")

# Parameters
sig <- "ours"
moa <- "tubulin polymerization inhibitor"

# Inputpath
prism_inputpath  <- "extdata/repurposing/PRISM_19Q4_secondary-screen-dose-response-curve-parameters.csv"
scores_inputpath <- "data/repurposing/signatures_scores.rds"

# Outputpath
outputpath <- "Figures/Figure5/scatter_yaccs_tubulin.png"

# Load data
scores <- as.data.frame(readRDS(scores_inputpath))

prism2 <- data.table::fread(
    prism_inputpath,
    select = c("broad_id", "depmap_id", "auc", "moa", "name", "screen_id"),
    showProgress = FALSE
    )

prism2 <- reshape2::dcast(
    dat = as.data.frame(prism2),
    formula = screen_id + name + broad_id + moa ~ depmap_id,
    fun.aggregate = sum,
    value.var = "auc"
    )

prism2[prism2 == 0] <- NA

# remove duplicates (same as the workflow: one row per broad_id)
prism2 <-
    prism2 %>%
    group_by(broad_id) %>%
    arrange(factor(screen_id, levels = c("MTS010", "MTS006", "MTS005", "HTS002"))) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    as.data.frame()

# Metadata (drug name) and AUC matrix
meta <- prism2[, c("broad_id", "name")]
auc_cols <- setdiff(names(prism2), c("screen_id", "name", "broad_id", "moa"))
prism2_m <- as.matrix(prism2[, auc_cols])
rownames(prism2_m) <- prism2$broad_id
prism2_m <- t(prism2_m)

# Signature score
scores_f <- scores[which(scores$signature == sig), ]
signature_score <- scores_f$score
names(signature_score) <- scores_f$cell_lines

z1 <- intersect(names(signature_score), rownames(prism2_m))

# Drugs in the MOA of interest (same construction as drug_type3)
drug_types <- strsplit(prism2$moa, split = ";")
names(drug_types) <- prism2$broad_id
ext <- data.frame(
    broad_id = names(unlist(drug_types)),
    moa = unlist(drug_types),
    stringsAsFactors = FALSE
    )
drug_types <- strsplit(ext$moa, split = ", ")
names(drug_types) <- ext$broad_id
ext2 <- data.frame(
    broad_id = names(unlist(drug_types)),
    moa = unlist(drug_types),
    stringsAsFactors = FALSE
    )
drug_type3 <- split(ext2$broad_id[which(!is.na(ext2$moa))], ext2$moa[which(!is.na(ext2$moa))])

stopifnot(moa %in% names(drug_type3))
drugs <- intersect(drug_type3[[moa]], colnames(prism2_m))

# Long table + per-drug Spearman correlation
plot_df <- list()
stat_df <- list()
for (d in drugs) {
    x <- signature_score[z1]
    y <- prism2_m[z1, d]
    keep <- complete.cases(x, y)
    if (sum(keep) < 3) next

    # exact = FALSE: AUC values are heavily tied, so the exact test is undefined
    ct <- suppressWarnings(
        cor.test(x[keep], y[keep], method = "spearman", exact = FALSE)
        )
    nm <- meta$name[match(d, meta$broad_id)]
    lab <- paste0(nm, " (", d, ")")

    plot_df[[d]] <- data.frame(
        drug = lab, score = x[keep], auc = y[keep], stringsAsFactors = FALSE
        )
    stat_df[[d]] <- data.frame(
        drug = lab, rho = unname(ct$estimate), pval = ct$p.value,
        n = sum(keep), stringsAsFactors = FALSE
        )
}
plot_df <- do.call(rbind, plot_df)
stat_df <- do.call(rbind, stat_df)

# Grid order: highest to lowest correlation
stat_df <- stat_df[order(-stat_df$rho), ]
lv <- stat_df$drug
plot_df$drug <- factor(plot_df$drug, levels = lv)
stat_df$drug <- factor(stat_df$drug, levels = lv)

stat_df$label <- sprintf("rho = %.2f, p = %.2g", stat_df$rho, stat_df$pval)

# Plot
p <- ggplot(plot_df, aes(x = score, y = auc)) +
    geom_point(colour = "black", size = 1.2, alpha = 0.7) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "blue", fill = "grey70") +
    geom_text(
        data = stat_df,
        aes(x = -Inf, y = Inf, label = label),
        hjust = -0.05, vjust = 1.4, size = 3, inherit.aes = FALSE
        ) +
    facet_wrap(~ drug, scales = "free_y", ncol = 4) +
    labs(
        title = paste0(moa, " - ", sig, " signature"),
        subtitle = "Ordered by decreasing Spearman correlation",
        x = paste0(sig, " signature score"),
        y = "Drug AUC (PRISM)"
        ) +
    theme_bw(base_size = 10) +
    theme(strip.text = element_text(size = 7))

ggsave(outputpath, p, width = 13, height = 16, dpi = 150)
cat("drugs plotted:", nrow(stat_df), "\n")
cat("cell lines per panel (range):", min(stat_df$n), "-", max(stat_df$n), "\n")
print(stat_df[, c("drug", "rho", "pval", "n")], row.names = FALSE)
