##### 01) Setup & Packages #####
Sys.setenv(LANG = "EN")
options(stringsAsFactors = FALSE, encoding = "UTF-8")
library(tidyverse)
library(data.table)
library(stringr)
library(dplyr)
library(tibble)
library(limma)
library(GSVA)
library(sva)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
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
library(ggsci)
library(RColorBrewer)
library(ggpp)
library(plotrix)
library(ChAMPdata)
library(IOBR)
library(TMEscore)
library(ImmuneSubtypeClassifier)
library(METAFlux)
library(TCellSI)
source("./R/standarize_fun.R")
source("./R/annTrackScale.R")

##### 02) Constants (paths, palette, labels) #####
INPUT_LUAD_DIR      <- "input/tcga_luad"
INPUT_ANALYSIS_DIR  <- "input/analysis_data"
INPUT_GEO_DIR       <- "input/geo/curated_data"
OUTPUT_DATA_DIR     <- "output/data/script5"
OUTPUT_FIG_DIR      <- "output/figures/fig5"
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

##### 03) Prepare data #####
# Expression
load(file.path(INPUT_LUAD_DIR, "TCGA-LUAD_mrna_expr_tpm.rdata"))  # -> mrna_expr_tpm
luadTpm <- {
  x <- mrna_expr_tpm[, substr(colnames(mrna_expr_tpm), 14, 15) != "11", drop = FALSE]
  colnames(x) <- substr(colnames(x), 1, 16)
  x <- x[, !duplicated(colnames(x)), drop = FALSE]
  log(x + 1)
}

# Survival & clinical (unused below but kept for completeness)
luadSurvival <- data.table::fread(file.path(INPUT_LUAD_DIR, "TCGA-LUAD.survival.tsv"))
luadClinical <- data.table::fread(file.path(INPUT_LUAD_DIR, "TCGA-LUAD.clinical.tsv")) |>
  dplyr::filter(tissue_type.samples == "Tumor")

# Immune subtypes from script 1
immuneSubtypeAssignments <- readRDS(file.path(OUTPUT_RDS_DIR, "immune_subtypes.rds")) |>
  dplyr::mutate(
    cluster      = factor(as.character(cluster), levels = c("WoundHealing","IFNGDominant","Inflammatory")),
    clusterLabel = clusterDisplayLabels[as.character(cluster)]
  )

# Align expression to subtype samples
luadTpm <- luadTpm[, intersect(colnames(luadTpm), immuneSubtypeAssignments$ID), drop = FALSE]

##### 04) A. TIDE — response composition #####
# If you need to export for TIDE webtool:
# tideInput <- t(apply(luadTpm, 1, function(x) x - mean(x)))
# write.table(tideInput, file.path(OUTPUT_DATA_DIR, "tide_input.txt"), sep = "\t", quote = FALSE, col.names = NA)

tideOutput <- data.table::fread(file.path(OUTPUT_DATA_DIR, "tide_output.csv")) |>
  dplyr::left_join(
    immuneSubtypeAssignments[, c("ID","clusterLabel")],
    dplyr::join_by("Patient" == "ID")
  ) |>
  dplyr::rename(Cluster = clusterLabel)

# Composition %
plotData5A <- tideOutput |>
  dplyr::group_by(Cluster, Responder) |>
  dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(Cluster) |>
  dplyr::mutate(Percentage = Count / sum(Count) * 100) |>
  dplyr::mutate(Cluster = factor(Cluster, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")))

xtabTide <- table(tideOutput$Cluster, tideOutput$Responder)
pvalTide <- suppressWarnings(chisq.test(xtabTide)$p.value)

fig5A_TIDE_Stacked <- ggplot(plotData5A, aes(Cluster, Percentage, fill = Responder)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5), size = 4, color = "black") +
  scale_fill_manual(values = c("TRUE"="#E64B35CC","FALSE"="#4DBBD5CC")) +
  labs(x = NULL, y = "Percentage (%)", fill = "Response") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold", size = 13)) +
  annotate("text", x = 2, y = 105,
           label = paste0("italic(P) == ", signif(pvalTide, 3)),
           parse = TRUE, fontface = "italic", size = 5)

##### 05) B. TIDE — score comparisons #####
plotListTide <- list()
tideScoreCols <- colnames(tideOutput)[c(4:9, 11:15)]
tideScoreCols <- intersect(tideScoreCols, colnames(tideOutput))

for (marker in tideScoreCols) {
  g <- ggplot(tideOutput, aes(Cluster, .data[[marker]], fill = Cluster)) +
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
    labs(x = NULL, y = paste0(marker, " Score")) +
    scale_fill_manual(values = setNames(PALETTE_CLUSTER3, c("Wound Healing","IFN-γ Dominant","Inflammatory")))
  plotListTide[[marker]] <- g
}

# Example composite usage:
# ggsave(file.path(OUTPUT_FIG_DIR, "fig5a_d_tide_examples.pdf"),
#        wrap_plots(fig5A_TIDE_Stacked | plotListTide[[1]] | plotListTide[[2]] | plotListTide[[4]]),
#        width = 12, height = 5)

##### 06) E. SubMap (format generation & demo heatmap) #####
generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type") {
  in_gct <- data.frame(GeneID = rownames(in_gct), description = "na", in_gct, check.names = FALSE)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct), "\t", ncol(in_gct) - 2, "\n", file = gct_file, append = TRUE)
  cat(paste(colnames(in_gct), collapse = "\t"), "\n", file = gct_file, append = TRUE)
  apply(in_gct, 1, function(row) cat(paste(row, collapse = "\t"), "\n", file = gct_file, append = TRUE))

  cat(nrow(sam_info), length(levels(factor(sam_info$rank))), 1, "\n", file = cls_file)
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " "), "\n", file = cls_file, sep = "", append = TRUE)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = TRUE)
}

# Public immunotherapy cohort for SubMap reference
skcm.logNC  <- read.table(file.path(INPUT_ANALYSIS_DIR, "skcm.immunotherapy.47samples.log2CountsNorm.txt"),
                          sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
rownames(skcm.logNC) <- toupper(rownames(skcm.logNC))
skcm.info   <- read.table(file.path(INPUT_ANALYSIS_DIR, "skcm.immunotherapy.47sampleInfo.txt"),
                          sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)
skcm.info   <- skcm.info[order(skcm.info$label), ]
skcm.info$rank <- rep(c(1,2,3,4), times = as.character(table(skcm.info$label))) # 1: CTLA4_noR 2: CTLA4_R 3: PD1_noR 4: PD1_R

tmpExpr  <- luadTpm
GENELIST <- intersect(rownames(tmpExpr), rownames(skcm.logNC))
in_gct   <- skcm.logNC[GENELIST, rownames(skcm.info)]

gctRefFile <- file.path(OUTPUT_DATA_DIR, "skcm.immunotherapy.for.submap.gct")
clsRefFile <- file.path(OUTPUT_DATA_DIR, "skcm.immunotherapy.for.submap.cls")
generateInputFileForSubMap(in_gct = in_gct, gct_file = gctRefFile, cls_file = clsRefFile, sam_info = skcm.info, type_name = "rank")

# Our cohort by immune cluster
samplesWH  <- immuneSubtypeAssignments |> dplyr::filter(cluster == "WoundHealing")  |> dplyr::pull(ID)
samplesIFN <- immuneSubtypeAssignments |> dplyr::filter(cluster == "IFNGDominant") |> dplyr::pull(ID)
samplesINF <- immuneSubtypeAssignments |> dplyr::filter(cluster == "Inflammatory") |> dplyr::pull(ID)

sam_info <- data.frame(ImmClust = c(samplesWH, samplesIFN, samplesINF),
                       row.names = c(samplesWH, samplesIFN, samplesINF))
sam_info$rank <- rep(c(1,2,3), times = c(length(samplesWH), length(samplesIFN), length(samplesINF)))

gctCohortFile <- file.path(OUTPUT_DATA_DIR, "immune3.for.submap.gct")
clsCohortFile <- file.path(OUTPUT_DATA_DIR, "immune3.for.submap.cls")
in_gct2 <- log2(tmpExpr[GENELIST, rownames(sam_info), drop = FALSE] + 1)
generateInputFileForSubMap(in_gct = in_gct2, gct_file = gctCohortFile, cls_file = clsCohortFile, sam_info = sam_info, type_name = "rank")

# Demo heatmap matrix (as in original script)
col_fun_submap <- colorRamp2(c(0, 0.5, 1), c("#FFA5AB","white","#BDD5EA"))
# If you want to draw the example:
# pdf(file.path(OUTPUT_FIG_DIR, "fig5e_submap_demo.pdf"), width = 6, height = 6)
# pheatmap(tmpMatrix, cellwidth = 30, cellheight = 30, cluster_rows = FALSE, cluster_cols = FALSE, color = col_fun_submap)
# dev.off()

##### 07) Cohort validations #####

## 7.1 IMvigor210 ----
load(file.path(INPUT_GEO_DIR, "IMvigor210_Cruated.Rdata")) # -> IMvigor210_Clinical
plotData5FG <- IMvigor210_Clinical |>
  tibble::rownames_to_column("ID") |>
  dplyr::filter(`Best Confirmed Overall Response` != "NE") |>
  dplyr::mutate(
    Response_group = dplyr::case_when(
      `Best Confirmed Overall Response` %in% c("CR","PR") ~ "Responder",
      `Best Confirmed Overall Response` %in% c("SD","PD","NE") ~ "Non-responder",
      TRUE ~ NA_character_
    ),
    Cluster = dplyr::case_when(
      `Immune phenotype` == "desert"   ~ "Wound Healing",
      `Immune phenotype` == "excluded" ~ "IFN-γ Dominant",
      `Immune phenotype` == "inflamed" ~ "Inflammatory",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(Cluster)) |>
  dplyr::mutate(
    Cluster = factor(Cluster, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")),
    os = os * 30
  )

plotData5F <- plotData5FG |>
  dplyr::group_by(Cluster, Response_group) |>
  dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(Cluster) |>
  dplyr::mutate(Percentage = Count / sum(Count) * 100)

xtabIMv <- table(plotData5FG$Cluster, plotData5FG$Response_group)
pvalIMv <- suppressWarnings(fisher.test(xtabIMv)$p.value)

fig5F_IMvigor <- ggplot(plotData5F, aes(Cluster, Percentage, fill = Response_group)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5), size = 4, color = "black") +
  scale_fill_manual(values = c("Responder"="#E64B35CC","Non-responder"="#4DBBD5CC")) +
  labs(x = NULL, y = "Percentage (%)", fill = "Response") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold", size = 13),
        legend.position = "none") +
  annotate("text", x = 2, y = 105,
           label = paste0("italic(P) == ", signif(pvalIMv, 3)),
           parse = TRUE, fontface = "italic", size = 5)

# Survival
survObjIMv <- Surv(time = plotData5FG$os, event = plotData5FG$censOS)
fitIMv <- survfit(survObjIMv ~ Cluster, data = plotData5FG)
survPlotIMv <- ggsurvplot(
  fitIMv, data = plotData5FG, pval = TRUE, risk.table = FALSE,
  palette = PALETTE_CLUSTER3, legend.title = "Immune Cluster",
  legend.labs = c("Wound Healing","IFN-γ Dominant","Inflammatory"),
  xlab = "Time (Days)", ylab = "Overall Survival Probability",
  ggtheme = theme_bw() +
    theme(panel.border = element_rect(colour = "black", size = 1.5),
          axis.text = element_text(size = 12, color = "black", face = "bold"),
          axis.title = element_text(size = 12, color = "black", face = "bold"),
          axis.ticks = element_line(size = 1, color = "black"),
          legend.text = element_text(size = 12, color = "black", face = "bold"),
          legend.title = element_text(size = 14, color = "black", face = "bold"))
)

# Save IMvigor figs
pdf(file.path(OUTPUT_FIG_DIR, "fig5f_imvigor_response.pdf"), width = 3, height = 5); print(fig5F_IMvigor); dev.off()
pdf(file.path(OUTPUT_FIG_DIR, "fig5g_imvigor_survival.pdf"), width = 5, height = 5); print(survPlotIMv$plot); dev.off()

## 7.2 GSE135222 ----
load(file.path(INPUT_GEO_DIR, "GSE135222_Curated.Rdata"))
immuneSubtype_GSE135222 <- tme_classifier(eset = GSE135222_exp, scale = TRUE)
plotData5HI <- GSE135222_Clinical |>
  dplyr::left_join(immuneSubtype_GSE135222[, c("ID","TMEcluster")], dplyr::join_by("Sample ID" == "ID")) |>
  dplyr::filter(!is.na(`Clinical benefit`)) |>
  dplyr::mutate(
    Response_group = dplyr::case_when(`Clinical benefit` == "DCB" ~ "Responder",
                                      `Clinical benefit` == "NDB" ~ "Non-responder",
                                      TRUE ~ NA_character_),
    Cluster = dplyr::case_when(
      TMEcluster == "IE" ~ "Wound Healing",
      TMEcluster == "IS" ~ "IFN-γ Dominant",
      TMEcluster == "IA" ~ "Inflammatory",
      TRUE ~ NA_character_
    ),
    Cluster = factor(Cluster, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")),
    PFS.time = as.numeric(PFS.time) * 30
  )

plotData5H <- plotData5HI |>
  dplyr::group_by(Cluster, Response_group) |>
  dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(Cluster) |>
  dplyr::mutate(Percentage = Count / sum(Count) * 100)

xtab135222 <- table(plotData5HI$Cluster, plotData5HI$Response_group)
pval135222 <- suppressWarnings(fisher.test(xtab135222)$p.value)

fig5H_GSE135222 <- ggplot(plotData5H, aes(Cluster, Percentage, fill = Response_group)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5), size = 4, color = "black") +
  scale_fill_manual(values = c("Responder"="#E64B35CC","Non-responder"="#4DBBD5CC")) +
  labs(x = NULL, y = "Percentage (%)", fill = "Response") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold", size = 13),
        legend.position = "none") +
  annotate("text", x = 2, y = 105,
           label = paste0("italic(P) == ", signif(pval135222, 3)),
           parse = TRUE, fontface = "italic", size = 5)

# Survival
survObj135222 <- Surv(time = plotData5HI$PFS.time, event = as.numeric(plotData5HI$PFS))
fit135222 <- survfit(survObj135222 ~ Cluster, data = plotData5HI)
survPlot135222 <- ggsurvplot(
  fit135222, data = plotData5HI, pval = TRUE, risk.table = FALSE,
  palette = PALETTE_CLUSTER3, legend.title = "Immune Cluster",
  legend.labs = c("Wound Healing","IFN-γ Dominant","Inflammatory"),
  xlab = "Time (Days)", ylab = "Overall Survival Probability",
  ggtheme = theme_bw() +
    theme(panel.border = element_rect(colour = "black", size = 1.5),
          axis.text = element_text(size = 12, color = "black", face = "bold"),
          axis.title = element_text(size = 12, color = "black", face = "bold"),
          axis.ticks = element_line(size = 1, color = "black"),
          legend.text = element_text(size = 12, color = "black", face = "bold"),
          legend.title = element_text(size = 14, color = "black", face = "bold"))
)

pdf(file.path(OUTPUT_FIG_DIR, "fig5h_gse135222_response.pdf"), width = 3, height = 5); print(fig5H_GSE135222); dev.off()
pdf(file.path(OUTPUT_FIG_DIR, "fig5i_gse135222_survival.pdf"), width = 5, height = 5); print(survPlot135222$plot); dev.off()

## 7.3 GSE207422 ----
load(file.path(INPUT_GEO_DIR, "GSE207422_Curated.Rdata"))
immuneSubtype_GSE207422 <- tme_classifier(eset = GSE207422_exp, scale = TRUE)
plotData5J_raw <- GSE207422_Clinical |>
  dplyr::left_join(immuneSubtype_GSE207422[, c("ID","TMEcluster")], by = dplyr::join_by("Sample" == "ID")) |>
  dplyr::filter(!is.na(RECIST)) |>
  dplyr::mutate(
    Response_group = dplyr::case_when(RECIST %in% c("CR","PR") ~ "Responder",
                                      RECIST %in% c("SD")      ~ "Non-responder",
                                      TRUE ~ NA_character_),
    Cluster = dplyr::case_when(
      TMEcluster == "IE" ~ "Wound Healing",
      TMEcluster == "IS" ~ "IFN-γ Dominant",
      TMEcluster == "IA" ~ "Inflammatory",
      TRUE ~ NA_character_
    ),
    Cluster = factor(Cluster, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory"))
  )

plotData5J <- plotData5J_raw |>
  dplyr::group_by(Cluster, Response_group) |>
  dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(Cluster) |>
  dplyr::mutate(Percentage = Count / sum(Count) * 100)

pval207422 <- suppressWarnings(fisher.test(table(plotData5J_raw$Cluster, plotData5J_raw$Response_group))$p.value)

fig5J_GSE207422 <- ggplot(plotData5J, aes(Cluster, Percentage, fill = Response_group)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5), size = 4, color = "black") +
  scale_fill_manual(values = c("Responder"="#E64B35CC","Non-responder"="#4DBBD5CC")) +
  labs(x = NULL, y = "Percentage (%)", fill = "Response") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold", size = 13),
        legend.position = "none") +
  annotate("text", x = 2, y = 105,
           label = paste0("italic(P) == ", signif(pval207422, 3)),
           parse = TRUE, fontface = "italic", size = 5)

pdf(file.path(OUTPUT_FIG_DIR, "fig5j_gse207422_response.pdf"), width = 3, height = 5); print(fig5J_GSE207422); dev.off()

## 7.4 GSE126044 ----
load(file.path(INPUT_GEO_DIR, "GSE126044_Curated.Rdata"))
immuneSubtype_GSE126044 <- tme_classifier(eset = GSE126044_exp, scale = TRUE)

plotData5K_raw <- GSE126044_Clinical |>
  dplyr::left_join(immuneSubtype_GSE126044[, c("ID","TMEcluster")], by = dplyr::join_by("ID" == "ID")) |>
  dplyr::mutate(
    Cluster = dplyr::case_when(
      TMEcluster == "IE" ~ "Wound Healing",
      TMEcluster == "IS" ~ "IFN-γ Dominant",
      TMEcluster == "IA" ~ "Inflammatory",
      TRUE ~ NA_character_
    ),
    Cluster = factor(Cluster, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory"))
  )

plotData5K <- plotData5K_raw |>
  dplyr::group_by(Cluster, Responsiveness) |>
  dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(Cluster) |>
  dplyr::mutate(Percentage = Count / sum(Count) * 100)

pval126044 <- suppressWarnings(fisher.test(table(plotData5K_raw$Cluster, plotData5K_raw$Responsiveness))$p.value)

fig5K_GSE126044 <- ggplot(plotData5K, aes(Cluster, Percentage, fill = Responsiveness)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5), size = 4, color = "black") +
  scale_fill_manual(values = c("Responder"="#E64B35CC","Non-responder"="#4DBBD5CC")) +
  labs(x = NULL, y = "Percentage (%)", fill = "Response") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold", size = 13),
        legend.position = "none") +
  annotate("text", x = 2, y = 105,
           label = paste0("italic(P) == ", signif(pval126044, 3)),
           parse = TRUE, fontface = "italic", size = 5)

pdf(file.path(OUTPUT_FIG_DIR, "fig5k_gse126044_response.pdf"), width = 3, height = 5); print(fig5K_GSE126044); dev.off()
