# Caculate pathway activity scores for TCGA-COAD patients and merge with clinical metadata
# =====================
library(TCGAbiolinks)
library(maftools)
require(dplyr)
require(tidyr)
library(GSVA)
require(CMSclassifier)
require(org.Hs.eg.db)

# signatures
tcga_inputpath <- "extdata/tcga-coad/counts_norm_coad_patients.rds"
signatures_inputpath <- "extdata/signatures/signatures_coad.rds"
metadata_inputpath <- "extdata/tcga-coad/metadata_coad_patients.rds"

# Load
gex <- readRDS(tcga_inputpath)
signatures <- readRDS(signatures_inputpath)
metadata <- readRDS(metadata_inputpath)

# Colon Cancer Cluster Subtypes (Guinney, Nat Med)
# =====================
m = mapIds(org.Hs.eg.db,
  keys = colnames(gex),
  column = 'ENTREZID',
  keytype = 'ENSEMBL')
colnames(gex) = m
gex_s = scale(gex)
clusters = CMSclassifier::classifyCMS(t(gex_s), method = 'SSP')[[3]]
metadata$clusters = clusters$SSP.nearestCMS

# Biomarkers Status
# =====================
query <- GDCquery(
  project      = "TCGA-COAD",
  data.category= "Simple Nucleotide Variation",
  data.type    = "Masked Somatic Mutation",
  access       = "open"
  # optionally also specify workflow.type (e.g., MuTect2/MuSE/etc) depending on what you want
)

# GDCdownload(query)
maf <- GDCprepare(query)

genes <- c("BRAF","KRAS","TP53","APC","SMAD4")
mutMat <- maf %>%
  transmute(
    sample = Tumor_Sample_Barcode,
    patient = substr(Tumor_Sample_Barcode, 1, 12),
    gene = Hugo_Symbol
  ) %>%
  filter(gene %in% genes) %>%
  count(patient, gene, name = "n_mut") %>%          # counts per patient+gene
  pivot_wider(names_from = gene, values_from = n_mut, values_fill = 0) %>%
  as.data.frame()

metadata <- metadata %>% left_join(mutMat, by = c("patient" = "patient"))
metadata$BRAF = ifelse(metadata$BRAF > 0, T, F)
metadata$KRAS = ifelse(metadata$KRAS > 0, T, F)
metadata$TP53 = ifelse(metadata$TP53 > 0, T, F)
metadata$APC = ifelse(metadata$APC > 0, T, F)
metadata$SMAD4 = ifelse(metadata$SMAD4 > 0, T, F)

# Calculate pathway activity
method = 'plage'
gsvaPar <- plageParam(gex %>% t(), signatures)
gsva.es <- gsva(gsvaPar, verbose=FALSE)
gsva.es %>% t() %>% as.data.frame() -> gsva

# Convert GSVA to data frame and add sample ID column
gsva_df <- gsva %>%
  as.data.frame() %>%
  tibble::rownames_to_column("barcode")

# Merge metadata with GSVA
merged_data <- metadata %>%
  left_join(gsva_df, by = "barcode") %>%
  rename(
    YACCS = ours,
    kras = KRAS,
    braf = BRAF,
    tp53 = TP53,
    apc = APC,
    smad4 = SMAD4
  ) %>%
  mutate(   # change directionality where needed
    coloGuideEx = coloGuideEx * -1,
    coloGuidePro = coloGuidePro * -1,
    mda114 = mda114 * -1
  )

# Melt by GSVA columns (pathways) and retain metadata variables
gsva_columns <- colnames(gsva)

melted_data <- merged_data %>%
  pivot_longer(
    cols = c("coloGuideEx", "coloGuidePro", "coloPrint", "mda114", "YACCS", "zhang"),
    names_to = "signature",
    values_to = "value"
  ) %>%
  mutate(
    MSI.Status = case_when(
        `MSI Status_MSI-L/MSS` == 1 ~ "MSI-L/MSS",
        `MSI Status_MSI-L/MSS` == 0 ~ "MSI-H",
        TRUE ~ NA_character_
    ),
    gender = ifelse(gender_male == 1, "male", "female"),
    Molecular_Subtype = ifelse(Molecular_Subtype_noCIN == 1, "noCIN", "CIN"),
    Stage = case_when(
      Stage_II == 1 ~ "II",
      Stage_III == 1 ~ "III",
      Stage_IV == 1 ~ "IV",
      TRUE ~ "I"
    )
  ) %>%
  dplyr::select(
    c(barcode, patient, vital_status, gender, age_at_diagnosis, Molecular_Subtype,
        Stage, MSI.Status, SurvTime, clusters, braf, kras, tp53, apc, smad4, signature, value
      )
  )

saveRDS(melted_data, file = "data/signatures/TCGA-COAD_clinvars_sigscores.rds")

