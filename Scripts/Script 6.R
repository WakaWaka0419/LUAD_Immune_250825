##### 01) Setup & Packages #####
Sys.setenv(LANG = "EN")
options(stringsAsFactors = FALSE, encoding = "UTF-8")
library(Seurat)
library(Matrix)
library(SCP)
library(harmony)
library(patchwork)
library(Nebulosa)
library(ggplot2)
library(data.table)
library(survival)
library(cowplot)
library(dplyr)
library(tibble)
library(ggpubr)
source("./R/standarize_fun.R")

##### 02) Constants (paths, palette) #####
INPUT_SEURAT_RDS   <- "~/Working_folder/DataBase/LUAD/scData/refquery_final.rds"       # source Seurat object
OUTPUT_BASE_DIR    <- "output/script6"
OUTPUT_DATA_DIR    <- file.path(OUTPUT_BASE_DIR, "data")
OUTPUT_FIG_DIR     <- file.path(OUTPUT_BASE_DIR, "figures")

dir.create(OUTPUT_BASE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_FIG_DIR,  recursive = TRUE, showWarnings = FALSE)

PALETTE_CLUSTER3   <- c("#2EC4B6", "#BDD5EA", "#FFA5AB")
PALETTE_BIN_2      <- c("#FF7F00", "#1F78B4")  # high / low binary plots

# DPIS panel
DPIS_GENES <- c("TPX2","CYP4B1","SCGB3A1","CACNA2D2","UBE2C","SFTPB","MYBL2","CDC20","BIRC5","SUSD2")

##### 03) Load #####
seu <- readRDS(INPUT_SEURAT_RDS)
meta_df <- seu@meta.data
print(table(meta_df$id))

##### 04) Figure 6A — baseline annotation UMAP #####
Figure_6A <- SCP::CellDimPlot(
  srt       = seu,
  group.by  = "Cell_Cluster_level1",
  theme_use = "theme_blank"
)

##### 05) Supplement S6 — marker features (DPIS panel) #####
Figure_S6 <- FeatureDimPlot(
  srt        = seu,
  feature    = DPIS_GENES,
  reduction  = "umap",
  theme_use  = "theme_blank"
)
# pdf(file.path(OUTPUT_FIG_DIR, "figure_s6_dpis_feature_umap.pdf"), width = 12, height = 9)
# Figure_S6
# dev.off()

##### 06) DPIS module score & grouping #####
# Seurat::AddModuleScore creates columns DPIS_Score1, DPIS_Score2... here only one set -> "...1"
seu <- AddModuleScore(
  object   = seu,
  features = list(DPIS_GENES),
  name     = "DPIS_Score"
)

# Harmonize the score column name (Seurat creates "DPIS_Score1")
seu@meta.data <- seu@meta.data |>
  dplyr::mutate(DPIS_Score = .data[["DPIS_Score1"]])

Figure_6B <- FeatureDimPlot(
  srt        = seu,
  features   = "DPIS_Score1",
  reduction  = "umap",
  theme_use  = "theme_blank"
)

# Binary group by DPIS score > 0
seu@meta.data <- seu@meta.data |>
  dplyr::mutate(
    DPIS_group = ifelse(.data[["DPIS_Score"]] > 0, "DPIS_High", "DPIS_Low")
  )

Figure_6C <- CellDimPlot(
  srt        = seu,
  group.by   = "DPIS_group",
  reduction  = "umap",
  theme_use  = "theme_blank",
  palcolor   = PALETTE_BIN_2
)

Figure_6D <- CellStatPlot(
  srt       = seu,
  stat.by   = "DPIS_group",
  plot_type = "ring"
)

Figure_6E <- CellStatPlot(
  srt      = seu,
  stat.by  = "Cell_Cluster_level1",
  group.by = "DPIS_group",
  label    = TRUE
)

# pdf(file.path(OUTPUT_FIG_DIR, "figure6_abc.pdf"), width = 12, height = 4.5)
# Figure_6A | Figure_6B | Figure_6C
# dev.off()
# pdf(file.path(OUTPUT_FIG_DIR, "figure6_de.pdf"), width = 8, height = 4.5)
# Figure_6D | Figure_6E
# dev.off()

##### 07) Subset: cancer-cell compartment & re-embed #####
CANCER_L2_LEVELS <- c("CDKN2A Cancer","CXCL1 Cancer","LAMC2 Cancer","Proliferating Cancer","SOX2 Cancer")

seu_cancer <- subset(seu, subset = Cell_Cluster_level2 %in% CANCER_L2_LEVELS) |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA()

# Batch correction & UMAP on Harmony space
seu_cancer <- RunHarmony(
  object         = seu_cancer,
  reduction      = "pca",
  group.by.vars  = "orig.ident",
  reduction.save = "harmony"
)
seu_cancer <- RunUMAP(
  object    = seu_cancer,
  reduction = "harmony",
  dims      = 1:30
)

Figure_6F <- SCP::CellDimPlot(
  srt       = seu_cancer,
  group.by  = "Cell_Cluster_level2",
  theme_use = "theme_blank"
)

Figure_6G <- FeatureDimPlot(
  srt        = seu_cancer,
  features   = "DPIS_Score1",
  theme_use  = "theme_blank"
)

Figure_6H <- FeatureStatPlot(
  srt      = seu_cancer,
  group.by = "Cell_Cluster_level2",
  bg.by    = "Cell_Cluster_level2",
  stat.by  = c("DPIS_Score1"),
  add_box  = TRUE
)

# pdf(file.path(OUTPUT_FIG_DIR, "figure6_fgh.pdf"), width = 12, height = 4.5)
# Figure_6F | Figure_6G | Figure_6H
# dev.off()

##### 08) Relabeling: proliferating ∩ DPIS_high #####
# Keep a copy of original level-2 label
seu_cancer@meta.data <- seu_cancer@meta.data |>
  dplyr::mutate(Cell_Cluster_level2_orig = .data[["Cell_Cluster_level2"]])

# Create combined label for DPIS+ proliferating cancer cells
seu_cancer@meta.data <- seu_cancer@meta.data |>
  dplyr::mutate(
    Cell_Cluster_level2 = ifelse(
      .data[["Cell_Cluster_level2"]] == "Proliferating Cancer" & .data[["DPIS_group"]] == "DPIS_High",
      "DPIS+ Proliferating Cancer",
      .data[["Cell_Cluster_level2"]]
    ),
    Proliferating_Flag = ifelse(.data[["Cell_Cluster_level2_orig"]] == "Proliferating Cancer", "Proliferating Cancer", "Others"),
    DPIS_Flag          = ifelse(.data[["DPIS_group"]] == "DPIS_High", "DPIS Cells", "Others"),
    DPIS_Prolif_Flag   = ifelse(.data[["Cell_Cluster_level2"]] == "DPIS+ Proliferating Cancer", "DPIS+ Proliferating Cancer", "Others")
  )

Figure_6I <- SCP::CellDimPlot(
  srt       = seu_cancer,
  group.by  = "Proliferating_Flag",
  theme_use = "theme_blank",
  palcolor  = PALETTE_BIN_2
)

Figure_6J <- SCP::CellDimPlot(
  srt       = seu_cancer,
  group.by  = "DPIS_Flag",
  theme_use = "theme_blank",
  palcolor  = PALETTE_BIN_2
)

Figure_6K <- SCP::CellDimPlot(
  srt       = seu_cancer,
  group.by  = "DPIS_Prolif_Flag",
  theme_use = "theme_blank",
  palcolor  = PALETTE_BIN_2
)

pdf(file.path(OUTPUT_FIG_DIR, "figure6_ijk.pdf"), width = 12, height = 4.5)
(Figure_6I | Figure_6J | Figure_6K)
dev.off()

# saveRDS(seu_cancer, file.path(OUTPUT_DATA_DIR, "seu_SCENIC.rds"))
