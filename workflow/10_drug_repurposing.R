# Run Drug enrichment analysis
# ---------
source("requirements.R")

# Inputpath
ccle_inputpath <- "extdata/repurposing/sample_info.csv"
prism_inputpath <- "extdata/repurposing/PRISM_19Q4_secondary-screen-dose-response-curve-parameters.csv"
scores_inputpath <- "data/repurposing/signatures_scores.rds"

# Outputpath
outputpath <- "data/repurposing/fgsea_by_signature.rds"

# Load data
scores <- readRDS(scores_inputpath) %>% as.data.frame()

# ccle_samlpe <- read.table(
#     ccle_inputpath,
#     header = TRUE,
#     sep = ",",
#     quote = '"',
#     fill = TRUE,
#     stringsAsFactors = FALSE,
#     row.names = NULL
#     )

prism2 <- read.table(
    prism_inputpath,
    header = TRUE,
    sep = ",",
    quote = '"',
    fill = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
    )

prism2 <- prism2[, c("broad_id", "depmap_id", "auc", "moa", "name", "screen_id")]

prism2 <- reshape2::dcast(
    dat = prism2,
    formula = screen_id + name + broad_id + moa ~ depmap_id,
    fun.aggregate = sum,
    value.var = "auc"
    )

prism2[prism2 == 0] <- NA

# remove duplicates
prism2 <-
    prism2 %>%
    group_by(broad_id) %>%
    arrange(factor(screen_id, levels = c("MTS010", "MTS006", "MTS005", "HTS002"))) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    dplyr::select(-c(screen_id, name))

prism2_m <- as.matrix(prism2[, -c(1, 2)])
rownames(prism2_m) <- prism2$broad_id
prism2_m <- t(prism2_m)


# Load scores from each signature
# =================================
signature <- unique(scores$signature)
fgseaRes <- list()
# i <- 2
for (i in seq_along(signature)) {

  print(paste0("Running ", signature[i]))
  scores_f <- scores[which(scores$signature == signature[i]), ]
  signature_score <- scores_f$score
  names(signature_score) <- scores_f$cell_lines

  z1 <- intersect(names(signature_score), rownames(prism2_m))
  cor_z1 <- apply(prism2_m[z1, ], 2, function(x) {
    cor(x, signature_score[z1], use = "complete.obs", method = "sp")
  })

  # Establishment of drug type list
  # ===========================
  drug_type1 <- split(prism2$broad_id[which(!is.na(prism2$moa))],prism2$moa[which(!is.na(prism2$moa))])
  drug_types <- strsplit(prism2$moa, split = ";")
  names(drug_types) <- prism2$broad_id
  ctrp_cmp_ext = data.frame(broad_id = names(unlist(drug_types)), moa=unlist(drug_types),stringsAsFactors = F)
  drug_type2 = split(ctrp_cmp_ext$broad_id[which(!is.na(ctrp_cmp_ext$moa))],
                     ctrp_cmp_ext$moa[which(!is.na(ctrp_cmp_ext$moa))])
  drug_types = strsplit(ctrp_cmp_ext$moa, split = ", ")
  names(drug_types) = ctrp_cmp_ext$broad_id
  ctrp_cmp_ext2 = data.frame(broad_id = names(unlist(drug_types)), moa=unlist(drug_types),stringsAsFactors = F)
  drug_type3 = split(ctrp_cmp_ext2$broad_id[which(!is.na(ctrp_cmp_ext2$moa))],
                     ctrp_cmp_ext2$moa[which(!is.na(ctrp_cmp_ext2$moa))])
  
  
  # Run GSEA analysis
  # ===========================
  set.seed(1993)
  fgseaRes[[i]] <- fgsea::fgsea(
    pathways = drug_type3,
    stats = sort(cor_z1), sampleSize = 5,
    minSize = 1,
    maxSize = Inf, 
    nproc = 3
    )

  fgseaRes[[i]]$log10adjpval <- -log10(fgseaRes[[i]]$padj)
}
names(fgseaRes) <- signature

saveRDS(fgseaRes, file = outputpath)


fgsea::plotEnrichment(fgseaRes$ours)
