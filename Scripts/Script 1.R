##### 01) Setup & Packages #####
Sys.setenv(LANG = "EN")
options(stringsAsFactors = FALSE, encoding = "UTF-8")
library(tidyverse)
library(data.table)
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(gghalves)
library(ggpp)
library(patchwork)
library(plotrix)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(factoextra)
library(ggConvexHull)
library(IOBR)
library(TMEscore)
library(ImmuneSubtypeClassifier)

source("./R/standarize_fun.R")

##### 02) Constants #####
INPUT_DIR          <- "input/tcga_luad"
OUTPUT_FIG_DIR     <- "output/figures/fig1"
OUTPUT_RDS_DIR     <- "output/rds"
PALETTE_CLUSTER3   <- c("#2EC4B6","#BDD5EA","#FFA5AB")
IMMUNE_SUBTYPE_MIN_SCORE <- 0.6

dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_RDS_DIR, recursive = TRUE, showWarnings = FALSE)

clusterDisplayLabels <- c(
  WoundHealing = "Wound Healing",
  IFNGDominant = "IFN-γ Dominant",
  Inflammatory = "Inflammatory"
)

##### 03) Load & Prepare Data #####
load(file.path(INPUT_DIR, "TCGA-LUAD_mrna_expr_tpm.rdata"))   # provides mrna_expr_tpm

luadTpm <- {
  x <- mrna_expr_tpm[, substr(colnames(mrna_expr_tpm), 14, 15) != "11", drop = FALSE]
  colnames(x) <- substr(colnames(x), 1, 16)
  x[, !duplicated(colnames(x)), drop = FALSE]
}
# log1p
luadTpm <- log(luadTpm + 1)

luadSurvival <- data.table::fread(file.path(INPUT_DIR, "TCGA-LUAD.survival.tsv"))
luadClinical <- data.table::fread(file.path(INPUT_DIR, "TCGA-LUAD.clinical.tsv")) |>
  dplyr::filter(tissue_type.samples == "Tumor")

# Immune subtype calling
immuneSubtypeAssignments <- ImmuneSubtypeClassifier::callEnsemble(X = luadTpm, geneids = "symbol")
calls <- geneMatchErrorReport(X = luadTpm, geneid = "symbol")
match_calls <- geneMatch(X = luadTpm, geneid = "symbol")

# 标准化列名与分群取值（数据层：WoundHealing/IFNGDominant/Inflammatory）
colnames(immuneSubtypeAssignments)[1:2] <- c("ID", "clusterRaw")
colnames(immuneSubtypeAssignments)[3:5] <- c("WoundHealing", "IFNGDominant", "Inflammatory")

immuneSubtypeAssignments <- immuneSubtypeAssignments |>
  dplyr::filter(clusterRaw %in% c("1","2","3")) |>
  dplyr::mutate(
    cluster = dplyr::case_when(
      clusterRaw == "1" ~ "WoundHealing",
      clusterRaw == "2" ~ "IFNGDominant",
      clusterRaw == "3" ~ "Inflammatory",
      .default = NA_character_
    )
  ) |>
  dplyr::mutate(
    cluster = factor(cluster, levels = c("WoundHealing","IFNGDominant","Inflammatory"))
  ) |>
  dplyr::mutate(
    maxScore = apply(dplyr::select(., WoundHealing, IFNGDominant, Inflammatory), 1, max, na.rm = TRUE)
  ) |>
  dplyr::filter(maxScore > IMMUNE_SUBTYPE_MIN_SCORE) |>
  dplyr::mutate(clusterLabel = clusterDisplayLabels[as.character(cluster)])

##### 04) Figure 1A: Cluster Composition (Pie) #####
plotData1A <- immuneSubtypeAssignments

fracWH <- round(sum(plotData1A$cluster == "WoundHealing")/nrow(plotData1A), 2)
fracIF <- round(sum(plotData1A$cluster == "IFNGDominant")/nrow(plotData1A), 2)
fracIN <- round(sum(plotData1A$cluster == "Inflammatory")/nrow(plotData1A), 2)

labels1A <- c(
  paste0("Wound Healing: ", round(fracWH*100,2), "%\n(", sum(plotData1A$cluster=="WoundHealing"), ")"),
  paste0("IFN-γ Dominant: ", round(fracIF*100,2), "%\n(", sum(plotData1A$cluster=="IFNGDominant"), ")"),
  paste0("Inflammatory: ", round(fracIN*100,2), "%\n(", sum(plotData1A$cluster=="Inflammatory"), ")")
)

fig1AClusterPie <- pie3D(
  c(fracWH, fracIF, fracIN),
  labels = labels1A, explode = 0.1, col = PALETTE_CLUSTER3, theta = 1.2, radius = 0.9
)

##### 05) Figure 1B: Clinical Heatmap (Annotation) #####
dfClinicalAnno <- luadClinical |>
  dplyr::filter(tissue_type.samples == "Tumor") |>
  dplyr::select(
    sample,
    status = vital_status.demographic,
    stage = ajcc_pathologic_stage.diagnoses,
    tStage = ajcc_pathologic_t.diagnoses,
    nStage = ajcc_pathologic_n.diagnoses,
    mStage = ajcc_pathologic_m.diagnoses
  ) |>
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
    immuneSubtypeAssignments |>
      dplyr::select(ID, cluster, clusterLabel),
    by = "ID"
  ) |>
  tibble::column_to_rownames(var = "ID") |>
  dplyr::mutate(across(everything(), ~ replace_na(as.character(.x), "Unknown")))

# Annotation
columnAnno <- HeatmapAnnotation(
  Cluster  = dfClinicalAnno$clusterLabel,
  Status   = dfClinicalAnno$status,
  Stage    = dfClinicalAnno$stage,
  T_Stage  = dfClinicalAnno$tStage,
  M_Stage  = dfClinicalAnno$mStage,
  N_Stage  = dfClinicalAnno$nStage,
  col = list(
    Cluster = c("Wound Healing" = "#2EC4B6","IFN-γ Dominant"="#BDD5EA","Inflammatory" = "#FFA5AB"),
    Status  = c("Alive"="#A8817A","Dead"="#E8BE74","Unknown"="#999999"),
    Stage   = c("Stage I"="#8ab1d2","Stage II"="#E58579","Stage III"="#D9BDD8","Stage IV"="#9180AC","Unknown"="#999999"),
    T_Stage = c("T1"="#FF9F1C","T2"="#FFA5AB","T3"="#023E8A","T4"="#9D4EDD","Unknown"="#999999"),
    N_Stage = c("N0"="#E64B35FF","N1"="#4DBBD5FF","N2"="#00A087FF","N3"="#8491B4FF","NX"="#3C5488FF","Unknown"="#999999"),
    M_Stage = c("M0"="#8491B4FF","M1"="#91D1C2FF","MX"="#DC0000FF","Unknown"="#999999")
  )
)

nCols <- nrow(dfClinicalAnno)
mat1bPlaceholder <- matrix(0, nrow = 1, ncol = nCols)
rownames(mat1bPlaceholder) <- "gene"
colnames(mat1bPlaceholder) <- rownames(dfClinicalAnno)

colSplitByCluster <- factor(dfClinicalAnno$clusterLabel, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory"))
colOrderByCluster <- c(
  which(dfClinicalAnno$clusterLabel == "Wound Healing"),
  which(dfClinicalAnno$clusterLabel == "IFN-γ Dominant"),
  which(dfClinicalAnno$clusterLabel == "Inflammatory")
)

fig1BClinicalHeatmap <- ComplexHeatmap::Heatmap(
  mat1bPlaceholder,
  col = c("0" = "white"),
  name = NULL,
  na_col = "white",
  show_column_names = FALSE,
  show_row_names = FALSE,
  row_names_side = "left",
  column_split = colSplitByCluster,
  column_order = colOrderByCluster,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  gap = unit(2, "mm"),
  top_annotation = columnAnno,
  height = unit(0.1, "mm"),
  border = FALSE
)

##### 06) Figure 1C: Clinical Composition Comparisons #####
# Status
df1cStatus <- dfClinicalAnno |>
  dplyr::select(status, clusterLabel) |>
  dplyr::rename(Status = status, Cluster = clusterLabel) |>
  dplyr::group_by(Cluster, Status) |>
  dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(Cluster) |>
  dplyr::mutate(Percentage = Count / sum(Count) * 100) |>
  dplyr::mutate(Cluster = factor(Cluster, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")))
xtabStatus <- table(dfClinicalAnno$clusterLabel, dfClinicalAnno$status)
pvalStatus <- suppressWarnings(chisq.test(xtabStatus)$p.value)

fig1CStatus <- ggplot(df1cStatus, aes(x = Cluster, y = Percentage, fill = Status)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            color = "black", size = 4) +
  scale_fill_manual(values = c("Alive"= "#E64B35CC", "Dead"= "#4DBBD5CC", "Unknown"="#999999")) +
  labs(x = NULL, y = "Percentage (%)", fill = "Status") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x  = element_text(face = "bold", angle = 45, hjust = 1),
    axis.text.y  = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold", size = 13)
  ) +
  annotate("text", x = 2, y = 105,
           label = paste0("italic(P) == ", signif(pvalStatus, 3)),
           parse = TRUE, fontface = "italic", size = 5)

# Stage
df1cStage <- dfClinicalAnno |>
  dplyr::select(stage, clusterLabel) |>
  dplyr::filter(stage != "Unknown") |>
  dplyr::rename(Stage = stage, Cluster = clusterLabel) |>
  dplyr::group_by(Cluster, Stage) |>
  dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(Cluster) |>
  dplyr::mutate(Percentage = Count / sum(Count) * 100) |>
  dplyr::mutate(Cluster = factor(Cluster, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")))
xtabStage <- with(dfClinicalAnno |> dplyr::filter(stage != "Unknown"),
                  table(clusterLabel, stage))
pvalStage <- suppressWarnings(chisq.test(xtabStage)$p.value)

fig1CStage <- ggplot(df1cStage, aes(x = Cluster, y = Percentage, fill = Stage)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            color = "black", size = 4) +
  scale_fill_manual(values = c("Stage I"="#E64B35CC","Stage II"="#4DBBD5CC","Stage III"="#00A087CC","Stage IV"="#3C5488CC")) +
  labs(x = NULL, y = "Percentage (%)", fill = "Stage") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x  = element_text(face = "bold", angle = 45, hjust = 1),
    axis.text.y  = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold", size = 13)
  ) +
  annotate("text", x = 2, y = 105,
           label = paste0("italic(P) == ", signif(pvalStage, 3)),
           parse = TRUE, fontface = "italic", size = 5)

# T Stage
df1cTStage <- dfClinicalAnno |>
  dplyr::select(tStage, clusterLabel) |>
  dplyr::filter(tStage != "Unknown") |>
  dplyr::rename(T_Stage = tStage, Cluster = clusterLabel) |>
  dplyr::group_by(Cluster, T_Stage) |>
  dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(Cluster) |>
  dplyr::mutate(Percentage = Count / sum(Count) * 100) |>
  dplyr::mutate(Cluster = factor(Cluster, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")))
xtabTStage <- with(dfClinicalAnno |> dplyr::filter(tStage != "Unknown"),
                   table(clusterLabel, tStage))
pvalTStage <- suppressWarnings(chisq.test(xtabTStage)$p.value)

fig1CTStage <- ggplot(df1cTStage, aes(x = Cluster, y = Percentage, fill = T_Stage)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            color = "black", size = 4) +
  scale_fill_manual(values = c("T1"="#E64B35CC","T2"="#4DBBD5CC","T3"="#00A087CC","T4"="#3C5488CC")) +
  labs(x = NULL, y = "Percentage (%)", fill = "T_Stage") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x  = element_text(face = "bold", angle = 45, hjust = 1),
    axis.text.y  = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold", size = 13)
  ) +
  annotate("text", x = 2, y = 105,
           label = paste0("italic(P) == ", signif(pvalTStage, 3)),
           parse = TRUE, fontface = "italic", size = 5)

# M Stage
df1cMStage <- dfClinicalAnno |>
  dplyr::select(mStage, clusterLabel) |>
  dplyr::filter(mStage != "Unknown") |>
  dplyr::rename(M_Stage = mStage, Cluster = clusterLabel) |>
  dplyr::group_by(Cluster, M_Stage) |>
  dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(Cluster) |>
  dplyr::mutate(Percentage = Count / sum(Count) * 100) |>
  dplyr::mutate(Cluster = factor(Cluster, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")))
xtabMStage <- with(dfClinicalAnno |> dplyr::filter(mStage != "Unknown"),
                   table(clusterLabel, mStage))
pvalMStage <- suppressWarnings(chisq.test(xtabMStage)$p.value)

fig1CMStage <- ggplot(df1cMStage, aes(x = Cluster, y = Percentage, fill = M_Stage)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            color = "black", size = 4) +
  scale_fill_manual(values = c("M0"="#E64B35CC","M1"="#4DBBD5CC","MX"="#00A087CC")) +
  labs(x = NULL, y = "Percentage (%)", fill = "M_Stage") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x  = element_text(face = "bold", angle = 45, hjust = 1),
    axis.text.y  = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold", size = 13)
  ) +
  annotate("text", x = 2, y = 105,
           label = paste0("italic(P) == ", signif(pvalMStage, 3)),
           parse = TRUE, fontface = "italic", size = 5)

# N Stage
df1cNStage <- dfClinicalAnno |>
  dplyr::select(nStage, clusterLabel) |>
  dplyr::filter(nStage != "Unknown") |>
  dplyr::rename(N_Stage = nStage, Cluster = clusterLabel) |>
  dplyr::group_by(Cluster, N_Stage) |>
  dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(Cluster) |>
  dplyr::mutate(Percentage = Count / sum(Count) * 100) |>
  dplyr::mutate(Cluster = factor(Cluster, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")))
xtabNStage <- with(dfClinicalAnno |> dplyr::filter(nStage != "Unknown"),
                   table(clusterLabel, nStage))
pvalNStage <- suppressWarnings(chisq.test(xtabNStage)$p.value)

fig1CNStage <- ggplot(df1cNStage, aes(x = Cluster, y = Percentage, fill = N_Stage)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            color = "black", size = 4) +
  scale_fill_manual(values = c("N0"="#E64B35CC","N1"="#4DBBD5CC","N2"="#00A087CC","N3"="#8491B4FF","NX"="#3C5488CC")) +
  labs(x = NULL, y = "Percentage (%)", fill = "N_Stage") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x  = element_text(face = "bold", angle = 45, hjust = 1),
    axis.text.y  = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold", size = 13)
  ) +
  annotate("text", x = 2, y = 105,
           label = paste0("italic(P) == ", signif(pvalNStage, 3)),
           parse = TRUE, fontface = "italic", size = 5)

layout1C <- "
ABC
DE.
"
fig1CComposite <- (
  fig1CStage +   # A
  fig1CTStage +  # B
  fig1CStatus +  # C
  fig1CMStage +  # D
  fig1CNStage    # E
) + plot_layout(design = layout1C)

ggsave(file.path(OUTPUT_FIG_DIR, "fig1c_composite.pdf"), fig1CComposite, width = 12, height = 9)

##### 07) Figure 1D: TMEscore #####
tmeScore <- tmescore(eset = luadTpm, method  = "PCA", classify = FALSE)
tmeScoreDf <- immuneSubtypeAssignments |>
  dplyr::left_join(tmeScore, by = "ID")

fig1DTmeScore <- ggplot(tmeScoreDf, aes(clusterLabel, TMEscore, fill = clusterLabel)) +
  geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.8, color = NA) +
  geom_boxplot(width = 0.4, size = 1.2, outlier.color = NA) +
  geom_jitter(aes(fill = clusterLabel), shape = 21, size = 2.5, width = 0.2, alpha = 0.5) +
  geom_signif(comparisons = list(
                c("Wound Healing","IFN-γ Dominant"),
                c("IFN-γ Dominant","Inflammatory"),
                c("Wound Healing","Inflammatory")),
              map_signif_level = TRUE, test = t.test, step_increase = 0.08) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.text.x = element_text(color = "black", size = 13),
        axis.text.y = element_text(color = "black", size = 13),
        legend.position = "none",
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.ticks = element_line(color = "black", linewidth = 1)) +
  labs(x = NULL, y = "TMEscore") +
  scale_fill_manual(values = setNames(PALETTE_CLUSTER3, c("Wound Healing","IFN-γ Dominant","Inflammatory"))) +
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_FIG_DIR, "fig1d_tmescore.pdf"), fig1DTmeScore, width = 3.5, height = 6)

##### 08) Figure 1E: Survival Difference #####
survPlotData <- immuneSubtypeAssignments |>
  dplyr::left_join(luadSurvival, dplyr::join_by("ID" == "sample")) |>
  dplyr::filter(OS.time > 30 & OS.time < 3650) |>
  dplyr::mutate(clusterLabel = factor(clusterLabel, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")))

survObject <- Surv(time = survPlotData$OS.time, event = survPlotData$OS)
fit <- survfit(survObject ~ clusterLabel, data = survPlotData)

survPlot <- ggsurvplot(
  fit,
  data = survPlotData,
  pval = TRUE,
  risk.table = FALSE,
  palette = PALETTE_CLUSTER3,
  legend.title = "Immune Cluster",
  legend.labs = c("Wound Healing","IFN-γ Dominant","Inflammatory"),
  xlab = "Time (Days)",
  ylab = "Overall Survival Probability",
  ggtheme = theme_bw() +
    theme(
      panel.border = element_rect(colour = "black", size = 1.5),
      axis.text = element_text(size = 12, color = "black", face = "bold"),
      axis.title = element_text(size = 12, color = "black", face = "bold"),
      axis.ticks = element_line(size = 1, color = "black"),
      legend.text = element_text(size = 12, color = "black", face = "bold"),
      legend.title = element_text(size = 14, color = "black", face = "bold")
    )
)

survDiffPairwise <- pairwise_survdiff(Surv(time = OS.time, event = OS) ~ clusterLabel, data = survPlotData)
survDiffPvalTable <- as.data.frame(round(survDiffPairwise[["p.value"]], 3))
survDiffPvalTable[is.na(survDiffPvalTable)] <- "-"
survDiffPvalTable <- tibble::rownames_to_column(survDiffPvalTable, var = "  ")

pdf(file.path(OUTPUT_FIG_DIR, "fig1e_survival.pdf"), width = 5, height = 5)
print(survPlot)
dev.off()

##### 09) Figure 1F-G-H: ESTIMATE Scores #####
tmeEstimate <- deconvo_tme(eset = luadTpm, method = "estimate")
tmeEstimateDf <- immuneSubtypeAssignments |>
  dplyr::left_join(tmeEstimate, by = "ID") |>
  dplyr::mutate(clusterLabel = factor(clusterLabel, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")))

fig1FImmuneScore <- ggplot(tmeEstimateDf, aes(clusterLabel, ImmuneScore_estimate, fill = clusterLabel)) +
  geom_half_violin(position = position_nudge(x=0.25), side = "r", width=0.8, color=NA) +
  geom_boxplot(width=0.4, size=1.2, outlier.color = NA) +
  geom_jitter(aes(fill= clusterLabel), shape=21, size=2.5, width=0.2, alpha = 0.5) +
  geom_signif(comparisons = list(
                c("Wound Healing","IFN-γ Dominant"),
                c("IFN-γ Dominant","Inflammatory"),
                c("Wound Healing","Inflammatory")),
              map_signif_level = TRUE, test = t.test, step_increase = 0.08) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.text.x = element_text(color = "black", size = 13),
        axis.text.y = element_text(color = "black", size = 13),
        legend.position = "none",
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.ticks = element_line(color="black", linewidth = 1)) +
  labs(x=NULL, y="Immune Score") +
  scale_fill_manual(values = setNames(PALETTE_CLUSTER3, c("Wound Healing","IFN-γ Dominant","Inflammatory"))) +
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1))

fig1GStromalScore <- ggplot(tmeEstimateDf, aes(clusterLabel, StromalScore_estimate, fill = clusterLabel)) +
  geom_half_violin(position = position_nudge(x=0.25), side = "r", width=0.8, color=NA) +
  geom_boxplot(width=0.4, size=1.2, outlier.color = NA) +
  geom_jitter(aes(fill= clusterLabel), shape=21, size=2.5, width=0.2, alpha = 0.5) +
  geom_signif(comparisons = list(
                c("Wound Healing","IFN-γ Dominant"),
                c("IFN-γ Dominant","Inflammatory"),
                c("Wound Healing","Inflammatory")),
              map_signif_level = TRUE, test = t.test, step_increase = 0.08) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.text.x = element_text(color = "black", size = 13),
        axis.text.y = element_text(color = "black", size = 13),
        legend.position = "none",
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.ticks = element_line(color="black", linewidth = 1)) +
  labs(x=NULL, y="Stromal Score") +
  scale_fill_manual(values = setNames(PALETTE_CLUSTER3, c("Wound Healing","IFN-γ Dominant","Inflammatory"))) +
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1))

fig1HTumorPurity <- ggplot(tmeEstimateDf, aes(clusterLabel, TumorPurity_estimate, fill = clusterLabel)) +
  geom_half_violin(position = position_nudge(x=0.25), side = "r", width=0.8, color=NA) +
  geom_boxplot(width=0.4, size=1.2, outlier.color = NA) +
  geom_jitter(aes(fill= clusterLabel), shape=21, size=2.5, width=0.2, alpha = 0.5) +
  geom_signif(comparisons = list(
                c("Wound Healing","IFN-γ Dominant"),
                c("IFN-γ Dominant","Inflammatory"),
                c("Wound Healing","Inflammatory")),
              map_signif_level = TRUE, test = t.test, step_increase = 0.08) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.text.x = element_text(color = "black", size = 13),
        axis.text.y = element_text(color = "black", size = 13),
        legend.position = "none",
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.ticks = element_line(color="black", linewidth = 1)) +
  labs(x=NULL, y="TumorPurity Score") +
  scale_fill_manual(values = setNames(PALETTE_CLUSTER3, c("Wound Healing","IFN-γ Dominant","Inflammatory"))) +
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1))

fig1FGHComposite <- fig1FImmuneScore + fig1GStromalScore + fig1HTumorPurity
ggsave(file.path(OUTPUT_FIG_DIR, "fig1fgh_composite.pdf"), fig1FGHComposite, width = 10.5, height = 6)

##### 10) Save objects #####
# saveRDS(immuneSubtypeAssignments, file.path(OUTPUT_RDS_DIR, "immune_subtypes.rds"))
# save.image(file.path(OUTPUT_RDS_DIR, "script1_image.Rdata"))
