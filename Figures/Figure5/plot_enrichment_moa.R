# Enrichment plot for a given MOA within a given signature
# ---------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(reshape2)
  library(fgsea); library(ggplot2)
})

sig <- "ours"
moa <- "tubulin polymerization inhibitor"

prism_inputpath  <- "extdata/repurposing/PRISM_19Q4_secondary-screen-dose-response-curve-parameters.csv"
scores_inputpath <- "data/repurposing/signatures_scores.rds"

scores <- as.data.frame(readRDS(scores_inputpath))

prism2 <- fread(prism_inputpath, showProgress = FALSE,
  select = c("broad_id", "depmap_id", "auc", "moa", "name", "screen_id"))

prism2 <- dcast(as.data.frame(prism2),
  screen_id + name + broad_id + moa ~ depmap_id,
  fun.aggregate = sum, value.var = "auc")

prism2[prism2 == 0] <- NA

prism2 <- prism2 %>%
  group_by(broad_id) %>%
  arrange(factor(screen_id, levels = c("MTS010", "MTS006", "MTS005", "HTS002"))) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(-c(screen_id, name))

prism2_m <- as.matrix(prism2[, -c(1, 2)])
rownames(prism2_m) <- prism2$broad_id
prism2_m <- t(prism2_m)

# stats: Spearman correlation between signature score and each drug AUC
scores_f <- scores[scores$signature == sig, ]
signature_score <- scores_f$score
names(signature_score) <- scores_f$cell_lines

z1 <- intersect(names(signature_score), rownames(prism2_m))
cor_z1 <- apply(prism2_m[z1, ], 2, function(x) {
  cor(x, signature_score[z1], use = "complete.obs", method = "sp")
})

# pathways: same drug_type3 construction as the workflow
drug_types <- strsplit(prism2$moa, split = ";")
names(drug_types) <- prism2$broad_id
ext <- data.frame(broad_id = names(unlist(drug_types)),
                  moa = unlist(drug_types), stringsAsFactors = FALSE)
drug_types <- strsplit(ext$moa, split = ", ")
names(drug_types) <- ext$broad_id
ext2 <- data.frame(broad_id = names(unlist(drug_types)),
                   moa = unlist(drug_types), stringsAsFactors = FALSE)
drug_type3 <- split(ext2$broad_id[!is.na(ext2$moa)], ext2$moa[!is.na(ext2$moa)])

stopifnot(moa %in% names(drug_type3))

# plot
p <- plotEnrichment(pathway = drug_type3[[moa]], stats = sort(cor_z1)) +
  labs(title = paste0(moa, " - ", sig),
       x = "rank",
       y = "enrichment score")
