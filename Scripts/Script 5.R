##### 01) Setup & Packages #####
Sys.setenv(LANG = "EN")
options(stringsAsFactors = FALSE, encoding = "UTF-8")
library(tidyverse)
library(data.table)
library(dplyr)
library(tibble)
library(stringr)
library(survival)
library(limma)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(gghalves)
library(patchwork)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(IOBR)
library(TMEscore)
library(ImmuneSubtypeClassifier)
library(Mime1)          # ML.Dev.Prog.Sig and downstream helpers
source("./R/standarize_fun.R")
source("./R/unicox.R")

##### 02) Constants (paths, palette, labels) #####
INPUT_LUAD_DIR      <- "input/tcga_luad"
INPUT_GEO_DIR       <- "input/geo/curated_data"
OUTPUT_FIG_DIR      <- "output/figures/fig3"
OUTPUT_RDS_DIR      <- "output/rds"

PALETTE_CLUSTER3    <- c("#2EC4B6","#BDD5EA","#FFA5AB")  # WoundHealing, IFNGDominant, Inflammatory
dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

clusterDisplayLabels <- c(
  WoundHealing = "Wound Healing",
  IFNGDominant = "IFN-γ Dominant",
  Inflammatory = "Inflammatory"
)

##### 03) Data Loading #####
# Expression (TPM -> log1p)
load(file.path(INPUT_LUAD_DIR, "TCGA-LUAD_mrna_expr_tpm.rdata"))  # -> mrna_expr_tpm
luadTpm <- {
  x <- mrna_expr_tpm[, substr(colnames(mrna_expr_tpm), 14, 15) != "11", drop = FALSE]
  colnames(x) <- substr(colnames(x), 1, 16)
  x <- x[, !duplicated(colnames(x)), drop = FALSE]
  log(x + 1)
}

# Survival & clinical
luadSurvival <- data.table::fread(file.path(INPUT_LUAD_DIR, "TCGA-LUAD.survival.tsv"))
luadClinical <- data.table::fread(file.path(INPUT_LUAD_DIR, "TCGA-LUAD.clinical.tsv")) |>
  dplyr::filter(tissue_type.samples == "Tumor")

# Immune subtypes (from Script 1 standard artifact)
immuneSubtypeAssignments <- readRDS(file.path(OUTPUT_RDS_DIR, "immune_subtypes.rds")) |>
  dplyr::mutate(
    cluster      = factor(as.character(cluster), levels = c("WoundHealing","IFNGDominant","Inflammatory")),
    clusterLabel = clusterDisplayLabels[as.character(cluster)]
  )

##### 04) Differential Expression: Inflammatory vs IFNGDominant #####
immPair <- immuneSubtypeAssignments |>
  dplyr::filter(clusterLabel %in% c("IFN-γ Dominant","Inflammatory")) |>
  dplyr::transmute(Cluster = factor(clusterLabel, levels = c("IFN-γ Dominant","Inflammatory")),
                   ID)

exprPair <- luadTpm[, intersect(colnames(luadTpm), immPair$ID), drop = FALSE]
stopifnot(identical(colnames(exprPair), immPair$ID))

design <- model.matrix(~ 0 + immPair$Cluster)
colnames(design) <- c("IFNGDominant","Inflammatory")
rownames(design) <- colnames(exprPair)

fit <- lmFit(exprPair, design)
fit <- contrasts.fit(fit, makeContrasts(Inflammatory - IFNGDominant, levels = design))
fit <- eBayes(fit)

degInflammatoryVsIFNG <- topTable(fit, number = Inf) |>
  tidyr::drop_na() |>
  dplyr::filter(abs(logFC) > 1.5, adj.P.Val < 0.01)

##### 05) Univariate Cox (filter features by survival) #####
# Align survival rows to expression columns
survAligned <- luadSurvival |> dplyr::filter(sample %in% colnames(exprPair))
exprAligned <- exprPair[, survAligned$sample, drop = FALSE]                   # reorder columns
stopifnot(identical(colnames(exprAligned), survAligned$sample))

# Build univariate Cox input frame: samples x genes (+ OS, OS.time)
uniInput <- exprAligned[rownames(exprAligned) %in% rownames(degInflammatoryVsIFNG), , drop = FALSE] |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("sample") |>
  dplyr::inner_join(survAligned[, c("sample","OS","OS.time")], by = "sample")

# Sanitize column names (coxph-friendly)
names(uniInput) <- gsub("-", "_", names(uniInput), fixed = TRUE)

# Run univariate Cox
geneCols <- setdiff(colnames(uniInput), c("sample","OS","OS.time"))
uniCoxList <- lapply(geneCols, function(g) unicox(uniInput, g))
uniCox <- do.call(rbind, uniCoxList) |> na.omit()
coxSelectedGenes <- uniCox[uniCox$pvalue <= 0.05, "group"]

##### 06) Build Train/Validation Datasets #####
# Training: TCGA
tcgaTrain <- uniInput |>
  dplyr::select(sample, OS, OS.time, dplyr::everything()) |>
  dplyr::rename(ID = sample)

# Helper to make GEO validation cohorts with same columns/order
make_validation <- function(expMat, clinDf, geneSet) {
  expMat |>
    dplyr::filter(rownames(.) %in% geneSet) |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("ID") |>
    dplyr::left_join(
      clinDf |>
        tibble::rownames_to_column("ID") |>
        dplyr::select(ID, OS, OS.time),
      by = "ID"
    ) |>
    dplyr::select(ID, OS.time, OS, dplyr::everything())
}

# Load curated GEO cohorts
load(file.path(INPUT_GEO_DIR, "GSE13213_Curated.Rdata"))   # -> GSE13213_exp, GSE13213_Clinical
load(file.path(INPUT_GEO_DIR, "GSE50081_Curated.Rdata"))   # -> GSE50081_exp, GSE50081_Clinical
load(file.path(INPUT_GEO_DIR, "GSE31210_Curated.Rdata"))   # -> GSE31210_exp, GSE31210_Clinical
load(file.path(INPUT_GEO_DIR, "GSE42127_Curated.Rdata"))   # -> GSE42127_exp, GSE42127_Clinical
load(file.path(INPUT_GEO_DIR, "GSE30219_Curated.Rdata"))   # -> GSE30219_exp, GSE30219_Clinical

val13213 <- make_validation(GSE13213_exp, GSE13213_Clinical, coxSelectedGenes)
val50081 <- make_validation(GSE50081_exp, GSE50081_Clinical, coxSelectedGenes)
val31210 <- make_validation(GSE31210_exp, GSE31210_Clinical, coxSelectedGenes)
val42127 <- make_validation(GSE42127_exp, GSE42127_Clinical, coxSelectedGenes)
val30219 <- make_validation(GSE30219_exp, GSE30219_Clinical, coxSelectedGenes)

# Align common columns across all datasets (order identical)
commonCols <- Reduce(intersect, list(
  colnames(tcgaTrain),
  colnames(val13213), colnames(val50081),
  colnames(val30219), colnames(val31210), colnames(val42127)
))

tcgaTrain <- dplyr::select(tcgaTrain, all_of(commonCols))
val13213  <- dplyr::select(val13213,  all_of(commonCols))
val50081  <- dplyr::select(val50081,  all_of(commonCols))
val30219  <- dplyr::select(val30219,  all_of(commonCols))
val31210  <- dplyr::select(val31210,  all_of(commonCols))
val42127  <- dplyr::select(val42127,  all_of(commonCols))

TrainValidationDatasets <- list(
  Dataset1 = tcgaTrain,
  Dataset2 = val13213,
  Dataset3 = val50081,
  Dataset4 = val30219,
  Dataset5 = val31210,
  Dataset6 = val42127
)

candidateGenes <- intersect(coxSelectedGenes, commonCols)

# saveRDS(candidateGenes, file.path(OUTPUT_RDS_DIR, "ml_candidate_genes.rds"))
# saveRDS(coxSelectedGenes, file.path(OUTPUT_RDS_DIR, "cox_selected_genes.rds"))

##### 07) Train ML Prognostic Signature (Mime1) #####
mlRes <- ML.Dev.Prog.Sig(
  train_data             = TrainValidationDatasets$Dataset1,
  list_train_vali_Data   = TrainValidationDatasets,
  unicox.filter.for.candi= FALSE,
  unicox_p_cutoff        = 0.05,
  candidate_genes        = candidateGenes,
  mode                   = "all",
  nodesize               = 5,
  seed                   = 917
)
# saveRDS(mlRes, file.path(OUTPUT_RDS_DIR, "ml_dev_prog_sig.rds"))

##### 08) C-index distribution & selected model survival plots #####
source("./R/cindex_dis_all.R")
# pdf(file.path(OUTPUT_FIG_DIR, "fig3a_cindex_all.pdf"), width = 7, height = 12)
cindex_dis_all(
  mlRes,
  validate_set = names(TrainValidationDatasets)[-1],
  order        = names(TrainValidationDatasets),
  width        = 0.15,
  n_models     = 100,
  pick_top     = TRUE
)
# dev.off()

cindex_dis_select(
  mlRes,
  model = "StepCox[forward] + Enet[α=0.1]",
  order = names(TrainValidationDatasets)
)

# Per-dataset KM curves for the chosen model
source("./R/rs_sur.R")
survPlots <- vector("list", 6)
for (i in seq_along(TrainValidationDatasets)) {
  survPlots[[i]] <- rs_sur(
    mlRes,
    model_name = "StepCox[forward] + Enet[α=0.1]",
    dataset    = names(TrainValidationDatasets)[i],
    median.line= "hv",
    cutoff     = 0.5,
    conf.int   = TRUE,
    xlab       = "Day",
    pval.coord = c(1000, 0.9)
  )
}
# pdf(file.path(OUTPUT_FIG_DIR, "fig3b_survival_panels.pdf"), width = 13, height = 9)
aplot::plot_list(gglist = survPlots, ncol = 3)
# dev.off()

##### 09) Core feature screening (optional heavy step) #####
# NOTE: The original run is extremely slow; keep commented or use cached results.
# res.feature.all <- ML.Corefeature.Prog.Screen(
#   InputMatrix    = TrainValidationDatasets$Dataset1,
#   candidate_genes= candidateGenes,
#   mode           = "all",
#   nodesize       = 5,
#   seed           = 1314
# )
# saveRDS(res.feature.all, file.path(OUTPUT_RDS_DIR, "ml_core_features_all.rds"))

# If you have res.feature.all from cache:
# source("./R/core_feature_rank.R")
# pdf(file.path(OUTPUT_FIG_DIR, "fig3k_core_feature_rank.pdf"), width = 4, height = 6)
# core_feature_rank(res.feature.all, top = 10)
# dev.off()
# core_feature_select(res.feature.all)

##### 10) AUC @ 1y/3y/5y & selection display #####
allAuc_1y <- cal_AUC_ml_res(
  res.by.ML.Dev.Prog.Sig = mlRes,
  train_data             = TrainValidationDatasets[["Dataset1"]],
  inputmatrix.list       = TrainValidationDatasets,
  mode                   = "all",
  AUC_time               = 1,
  auc_cal_method         = "KM"
)
allAuc_3y <- cal_AUC_ml_res(
  res.by.ML.Dev.Prog.Sig = mlRes,
  train_data             = TrainValidationDatasets[["Dataset1"]],
  inputmatrix.list       = TrainValidationDatasets,
  mode                   = "all",
  AUC_time               = 3,
  auc_cal_method         = "KM"
)
allAuc_5y <- cal_AUC_ml_res(
  res.by.ML.Dev.Prog.Sig = mlRes,
  train_data             = TrainValidationDatasets[["Dataset1"]],
  inputmatrix.list       = TrainValidationDatasets,
  mode                   = "all",
  AUC_time               = 5,
  auc_cal_method         = "KM"
)

auc_plot_1y <- auc_dis_all(
  allAuc_1y,
  dataset      = names(TrainValidationDatasets),
  validate_set = names(TrainValidationDatasets)[-1],
  order        = names(TrainValidationDatasets),
  width        = 0.35,
  year         = 1
)
auc_plot_3y <- auc_dis_all(
  allAuc_3y,
  dataset      = names(TrainValidationDatasets),
  validate_set = names(TrainValidationDatasets)[-1],
  order        = names(TrainValidationDatasets),
  width        = 0.35,
  year         = 3
)
auc_plot_5y <- auc_dis_all(
  allAuc_5y,
  dataset      = names(TrainValidationDatasets),
  validate_set = names(TrainValidationDatasets)[-1],
  order        = names(TrainValidationDatasets),
  width        = 0.35,
  year         = 5
)

source("./R/auc_dis_select.R")
# pdf(file.path(OUTPUT_FIG_DIR, "fig3i_auc_select.pdf"), width = 7, height = 5)
auc_dis_select(
  list(allAuc_1y, allAuc_3y, allAuc_5y),
  model_name = "StepCox[forward] + Enet[α=0.1]",
  dataset    = names(TrainValidationDatasets),
  order      = c("Dataset6","Dataset5","Dataset4","Dataset3","Dataset2","Dataset1"),
  year       = c(1, 3, 5)
)
# dev.off()

##### 11) Meta-analysis of univariate HR across datasets (selected model) #####
uniCoxRs <- cal_unicox_ml_res(
  res.by.ML.Dev.Prog.Sig = mlRes,
  optimal.model          = "StepCox[forward] + Enet[α=0.1]",
  type                   = "categorical"
)
metaModel <- cal_unicox_meta_ml_res(input = uniCoxRs)

source("./R/meta_unicox_vis.R")
pdf(file.path(OUTPUT_FIG_DIR, "fig3j_meta_unicox.pdf"), width = 9, height = 4)
meta_unicox_vis(metaModel, dataset = names(TrainValidationDatasets))
dev.off()
