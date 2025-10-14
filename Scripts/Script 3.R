##### 01) Setup & Packages #####
Sys.setenv(LANG = "EN")
options(stringsAsFactors = FALSE, encoding = "UTF-8")
library(tidyverse)
library(data.table)
library(dplyr)
library(tibble)
library(stringr)
library(RTCGAToolbox)
library(limma)
library(GSVA)
library(clusterProfiler)
library(maftools)
library(ggplot2)
library(ggsci)
library(cowplot)
library(patchwork)
library(ggsignif)
library(gghalves)
library(ComplexHeatmap)
library(pheatmap)
library(circlize)
library(RColorBrewer)
library(viridis)
library(gplots)
library(cowplot)
library(TCellSI)
library(MOVICS)
library(BiocOncoTK)
source("./R/annTrackScale.R")
source("./R/standarize_fun.R")

##### 02) Constants (paths, labels, palette) #####
INPUT_LUAD_DIR      <- "input/tcga_luad"
INPUT_ANALYSIS_DIR  <- "input/analysis_data"
OUTPUT_DATA_DIR     <- "output/data/script2"
OUTPUT_FIG_DIR      <- "output/figures/fig2"
OUTPUT_RDS_DIR      <- "output/rds"
PALETTE_CLUSTER3    <- c("#2EC4B6","#BDD5EA","#FFA5AB")  # WoundHealing, IFNGDominant, Inflammatory
dir.create(OUTPUT_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_RDS_DIR, recursive = TRUE, showWarnings = FALSE)
clusterDisplayLabels <- c(
  WoundHealing = "Wound Healing",
  IFNGDominant = "IFN-γ Dominant",
  Inflammatory = "Inflammatory"
)

##### 03) Load Data #####
load(file.path(INPUT_LUAD_DIR, "TCGA-LUAD_mrna_expr_tpm.rdata"))  # -> mrna_expr_tpm

luadTpm <- {
  x <- mrna_expr_tpm[, substr(colnames(mrna_expr_tpm), 14, 15) != "11", drop = FALSE] # 排除正常样本
  colnames(x) <- substr(colnames(x), 1, 16)
  x <- x[, !duplicated(colnames(x)), drop = FALSE]
  log(x + 1)
}

luadSurvival <- data.table::fread(file.path(INPUT_LUAD_DIR, "TCGA-LUAD.survival.tsv"))
luadClinical <- data.table::fread(file.path(INPUT_LUAD_DIR, "TCGA-LUAD.clinical.tsv")) |>
  dplyr::filter(tissue_type.samples == "Tumor")
immuneSubtypeAssignments <- readRDS(file.path(OUTPUT_RDS_DIR, "immune_subtypes.rds")) |>
  dplyr::mutate(
    cluster      = factor(as.character(cluster), levels = c("WoundHealing","IFNGDominant","Inflammatory")),
    clusterLabel = clusterDisplayLabels[as.character(cluster)]
  )

##### 04) Figure 2A — Immune Landscape (ICI & TIME Heatmaps) #####
immuneSignature <- read.table(file.path(INPUT_ANALYSIS_DIR, "Curated_Immune_Cell_Signature.txt"),
                              sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
cellTypes <- unique(immuneSignature$CellType)
immuneSigList <- lapply(setNames(nm = cellTypes), function(ct) immuneSignature[immuneSignature$CellType == ct, "Symbol"])
iciTargets <- c("CD274","PDCD1","CD247","PDCD1LG2","CTLA4","TNFRSF9","TNFRSF4","TLR9")
immuneSigOrder <- c(
  "T.cells.CD8","T.cells.regulatory..Tregs.","T.cells.CD4.naive","T.cells.follicular.helper",
  "B.cells.naive","B.cells.memory","T.cells.gamma.delta","Dendritic.cells.activated",
  "Macrophages.M1","NK.cells.activated","Plasma.cells","T.cells.CD4.memory.resting",
  "T.cells.CD4.memory.activated","Mast.cells.activated","NK.cells.resting","Macrophages.M0",
  "Macrophages.M2","Eosinophils","Monocytes","Dendritic.cells.resting","Mast.cells.resting",
  "Neutrophils","Endothelial cells","Fibroblasts"
)
write.table(luadTpm,
            file = file.path(OUTPUT_DATA_DIR, "tcga_log2tpm_hugo.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
filterCommonGenes(input.f = file.path(OUTPUT_DATA_DIR, "tcga_log2tpm_hugo.txt"),
                  output.f = file.path(OUTPUT_DATA_DIR, "tcga_log2tpm_hugo_estimate.txt"),
                  id = "GeneSymbol")
estimateScore(file.path(OUTPUT_DATA_DIR, "tcga_log2tpm_hugo_estimate.txt"),
              file.path(OUTPUT_DATA_DIR, "tcga_log2tpm_hugo_estimate_score.txt"),
              platform = "affymetrix")
estScores <- read.table(file = file.path(OUTPUT_DATA_DIR, "tcga_log2tpm_hugo_estimate_score.txt"),
                        header = TRUE, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE, sep = "\t")
rownames(estScores) <- estScores[,2]
colnames(estScores) <- estScores[1,]
estScores <- estScores[-1, c(-1, -2)]
estScores <- sapply(estScores, as.numeric)
rownames(estScores) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity")
estScaled <- annTrackScale(indata = estScores, halfwidth = 2, poolsd = FALSE)
estScaled <- as.data.frame(t(estScaled))  # 行：样本；列：分数
rownames(estScaled) <- colnames(luadTpm)
immuneGsva <- gsva(as.matrix(luadTpm), immuneSigList, method = "gsva")
colnames(immuneGsva) <- substr(colnames(immuneGsva), 1, 16)
df2A <- luadClinical |>
  dplyr::filter(tissue_type.samples == "Tumor") |>
  dplyr::select(sample,
                status = vital_status.demographic,
                stage  = ajcc_pathologic_stage.diagnoses,
                tStage = ajcc_pathologic_t.diagnoses,
                nStage = ajcc_pathologic_n.diagnoses,
                mStage = ajcc_pathologic_m.diagnoses) |>
  dplyr::rename(ID = sample) |>
  dplyr::mutate(
    stage = dplyr::case_when(
      stage %in% c("Stage I","Stage IA","Stage IB") ~ "Stage I",
      stage %in% c("Stage II","Stage IIA","Stage IIB") ~ "Stage II",
      stage %in% c("Stage IIIA","Stage IIIB") ~ "Stage III",
      stage %in% c("Stage IV") ~ "Stage IV",
      .default = "Unknown"
    ),
    tStage = dplyr::case_when(
      tStage %in% c("T1","T1a","T1b") ~ "T1",
      tStage %in% c("T2","T2a","T2b") ~ "T2",
      tStage %in% c("T3","T3a","T3b","T3c") ~ "T3",
      tStage %in% c("T4") ~ "T4",
      .default = "Unknown"
    ),
    mStage = dplyr::case_when(
      mStage %in% c("M0") ~ "M0",
      mStage %in% c("M1","M1a","M1b") ~ "M1",
      mStage %in% c("MX") ~ "MX",
      .default = "Unknown"
    )
  ) |>
  dplyr::right_join(
    immuneSubtypeAssignments |> dplyr::select(ID, cluster, clusterLabel),
    by = "ID"
  ) |>
  tibble::column_to_rownames("ID") |>
  dplyr::mutate(across(everything(), ~ replace_na(as.character(.x), "Unknown"))) |>
  dplyr::mutate(
    clusterLabel = factor(clusterLabel, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory"))
  ) |>
  dplyr::arrange(clusterLabel)

splitByCluster <- factor(df2A$clusterLabel, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory"))
annColors <- list(
  Cluster      = c("Wound Healing" = PALETTE_CLUSTER3[1],
                   "IFN-γ Dominant" = PALETTE_CLUSTER3[2],
                   "Inflammatory" = PALETTE_CLUSTER3[3]),
  ImmuneScore  = inferno(64),
  StromalScore = viridis(64),
  Stage        = c("Stage I"="#8ab1d2","Stage II"="#E58579","Stage III"="#D9BDD8","Stage IV"="#9180AC","Unknown"="#999999"),
  T_Stage      = c("T1"="#FF9F1C","T2"="#FFA5AB","T3"="#023E8A","T4"="#9D4EDD","Unknown"="#999999"),
  N_Stage      = c("N0"="#E64B35FF","N1"="#4DBBD5FF","N2"="#00A087FF","N3"="#8491B4FF","NX"="#3C5488FF","Unknown"="#999999"),
  M_Stage      = c("M0"="#8491B4FF","M1"="#91D1C2FF","MX"="#DC0000FF","Unknown"="#999999"),
  Status       = c("Alive"="#A8817A","Dead"="#E8BE74","Unknown"="#999999")
)
iciExpr <- luadTpm[intersect(rownames(luadTpm), iciTargets), , drop = FALSE]
iciExpr <- iciExpr[iciTargets, rownames(df2A), drop = FALSE]
fig2A_hm1 <- pheatmap(
  standarize.fun(iciExpr, halfwidth = 2),
  border_color = NA,
  annotation_col = df2A[, c("clusterLabel","StromalScore","ImmuneScore")] |>
    dplyr::rename(Cluster = clusterLabel),
  annotation_colors = annColors[c("Cluster","StromalScore","ImmuneScore")],
  color = NMF:::ccRamp(x = c("#2dabb9","white","#FF9F1C"), n = 64),
  column_split = splitByCluster,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cellheight = 12,
  cellwidth  = 0.6,
  name = "ICI",
  cluster_rows = FALSE,
  cluster_cols = FALSE
)

immuneGsvaSub <- immuneGsva[immuneSigOrder, rownames(df2A), drop = FALSE]
fig2A_hm2 <- pheatmap(
  standarize.fun(immuneGsvaSub, halfwidth = 1),
  border_color = NA,
  color = NMF:::ccRamp(x = c("#2dabb9","white","#FF9F1C"), n = 64),
  column_split = splitByCluster,
  gaps_row = c(14, 22),
  show_rownames = TRUE,
  show_colnames = FALSE,
  cellheight = 12,
  cellwidth  = 0.6,
  name = "TIME",
  cluster_rows = FALSE,
  cluster_cols = FALSE
)

pdf(file.path(OUTPUT_FIG_DIR, "fig2a_time_landscape.pdf"), width = 10, height = 10)
draw(fig2A_hm1 %v% fig2A_hm2, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

##### 05) Figure 2B–I — T cell states & Genomic metrics #####
## 2B–I: TCellSI 状态打分
tcellScore <- TCellSI::TCSS_Calculate(luadTpm, ref = TRUE)
tcellDf <- immuneSubtypeAssignments |>
  dplyr::left_join(
    tcellScore |> t() |> as.data.frame() |> tibble::rownames_to_column("ID"),
    by = "ID"
  ) |>
  dplyr::mutate(clusterLabel = factor(clusterLabel, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")))

plotListTcell <- list()
for (state in rownames(tcellScore)) {
  g <- ggplot(tcellDf, aes(clusterLabel, .data[[state]], fill = clusterLabel)) +
    geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.8, color = NA) +
    geom_boxplot(width = 0.4, size = 1.2, outlier.color = NA) +
    geom_jitter(shape = 21, size = 2.5, width = 0.2, alpha = 0.5) +
    geom_signif(
      comparisons = list(c("Wound Healing","IFN-γ Dominant"),
                         c("IFN-γ Dominant","Inflammatory"),
                         c("Wound Healing","Inflammatory")),
      map_signif_level = TRUE, test = t.test, step_increase = 0.08
    ) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(size = 2),
          axis.text.x = element_text(color = "black", size = 13, face = "bold", angle = 45, hjust = 1),
          axis.text.y = element_text(color = "black", size = 13),
          legend.position = "none",
          axis.title.y = element_text(size = 15, face = "bold"),
          axis.ticks = element_line(color = "black", linewidth = 1)) +
    labs(x = NULL, y = paste0(state, " Score")) +
    scale_fill_manual(values = setNames(PALETTE_CLUSTER3, c("Wound Healing","IFN-γ Dominant","Inflammatory")))
  plotListTcell[[state]] <- g
}
pdf(file.path(OUTPUT_FIG_DIR, "fig2e_l_tcell_states.pdf"), width = 12, height = 11)
wrap_plots(plotListTcell, ncol = 4)
dev.off()

## 2B: TMB
load(file.path(INPUT_LUAD_DIR, "TCGA-LUAD_maf.rdata"))  # -> data (MAF-like)
tmbMaf <- read.maf(data)
tmbRes <- tmb(tmbMaf, captureSize = 38, logScale = TRUE)

df2B <- immuneSubtypeAssignments |>
  dplyr::left_join(
    tmbRes |>
      dplyr::rename(ID = Tumor_Sample_Barcode) |>
      dplyr::mutate(ID = substr(ID, 1, 16)),
    by = "ID"
  )
fig2B_TMB <- ggplot(df2B, aes(clusterLabel, total_perMB_log, fill = clusterLabel)) +
  geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.8, color = NA) +
  geom_boxplot(width = 0.4, size = 1.2, outlier.color = NA) +
  geom_jitter(shape = 21, size = 2.5, width = 0.2, alpha = 0.5) +
  geom_signif(comparisons = list(c("Wound Healing","IFN-γ Dominant"),
                                 c("IFN-γ Dominant","Inflammatory"),
                                 c("Wound Healing","Inflammatory")),
              map_signif_level = TRUE, test = t.test, step_increase = 0.08) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.text.x = element_text(color = "black", size = 13, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black", size = 13),
        legend.position = "none",
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.ticks = element_line(color = "black", linewidth = 1)) +
  labs(x = NULL, y = "Total PerMB log") +
  scale_fill_manual(values = setNames(PALETTE_CLUSTER3, c("Wound Healing","IFN-γ Dominant","Inflammatory")))

## 2C: FGA
fgaDat <- data.table::fread(file.path(INPUT_ANALYSIS_DIR, "Fraction_Genome_Altered.txt"))
df2C <- immuneSubtypeAssignments |>
  dplyr::mutate(ID15 = substr(ID, 1, 15)) |>
  dplyr::left_join(fgaDat |> dplyr::rename(ID15 = `Sample ID`), by = "ID15")
fig2C_FGA <- ggplot(df2C, aes(clusterLabel, `Fraction Genome Altered`, fill = clusterLabel)) +
  geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.8, color = NA) +
  geom_boxplot(width = 0.4, size = 1.2, outlier.color = NA) +
  geom_jitter(shape = 21, size = 2.5, width = 0.2, alpha = 0.5) +
  geom_signif(comparisons = list(c("Wound Healing","IFN-γ Dominant"),
                                 c("IFN-γ Dominant","Inflammatory"),
                                 c("Wound Healing","Inflammatory")),
              map_signif_level = TRUE, test = t.test, step_increase = 0.08) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.text.x = element_text(color = "black", size = 13, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black", size = 13),
        legend.position = "none",
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.ticks = element_line(color = "black", linewidth = 1)) +
  labs(x = NULL, y = "Fraction Genome Altered") +
  scale_fill_manual(values = setNames(PALETTE_CLUSTER3, c("Wound Healing","IFN-γ Dominant","Inflammatory")))

## 2D: Mutation Count
mutCountDat <- data.table::fread(file.path(INPUT_ANALYSIS_DIR, "Mutation_Count.txt"))
df2D <- immuneSubtypeAssignments |>
  dplyr::mutate(ID15 = substr(ID, 1, 15)) |>
  dplyr::left_join(mutCountDat |> dplyr::rename(ID15 = `Sample ID`), by = "ID15")
fig2D_MutCount <- ggplot(df2D, aes(clusterLabel, `Mutation Count`, fill = clusterLabel)) +
  geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.8, color = NA) +
  geom_boxplot(width = 0.4, size = 1.2, outlier.color = NA) +
  geom_jitter(shape = 21, size = 2.5, width = 0.2, alpha = 0.5) +
  geom_signif(comparisons = list(c("Wound Healing","IFN-γ Dominant"),
                                 c("IFN-γ Dominant","Inflammatory"),
                                 c("Wound Healing","Inflammatory")),
              map_signif_level = TRUE, test = t.test, step_increase = 0.08) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.text.x = element_text(color = "black", size = 13, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black", size = 13),
        legend.position = "none",
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.ticks = element_line(color = "black", linewidth = 1)) +
  labs(x = NULL, y = "Mutation Count") +
  scale_fill_manual(values = setNames(PALETTE_CLUSTER3, c("Wound Healing","IFN-γ Dominant","Inflammatory")))

## S2D: MSI (MSIsensor.10k 内置数据)
data("MSIsensor.10k")
dfMSI <- immuneSubtypeAssignments |>
  dplyr::mutate(participant_barcode = substr(ID, 1, 12)) |>
  dplyr::left_join(MSIsensor.10k, by = "participant_barcode")
figS2D_MSI <- ggplot(dfMSI, aes(clusterLabel, MSIsensor.score, fill = clusterLabel)) +
  geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.8, color = NA) +
  geom_boxplot(width = 0.4, size = 1.2, outlier.color = NA) +
  geom_jitter(shape = 21, size = 2.5, width = 0.2, alpha = 0.5) +
  geom_signif(comparisons = list(c("Wound Healing","IFN-γ Dominant"),
                                 c("IFN-γ Dominant","Inflammatory"),
                                 c("Wound Healing","Inflammatory")),
              map_signif_level = TRUE, test = t.test, step_increase = 0.08) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.text.x = element_text(color = "black", size = 13, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black", size = 13),
        legend.position = "none",
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.ticks = element_line(color = "black", linewidth = 1)) +
  labs(x = NULL, y = "MSI Score") +
  scale_fill_manual(values = setNames(PALETTE_CLUSTER3, c("Wound Healing","IFN-γ Dominant","Inflammatory")))
pdf(file.path(OUTPUT_FIG_DIR, "fig2b_d_genomic_metrics.pdf"), width = 12, height = 5.5)
print(fig2B_TMB | fig2C_FGA | fig2D_MutCount | figS2D_MSI)
dev.off()

