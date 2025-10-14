#####01 Loading Packages#####
Sys.setenv(LANG = 'EN')
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
use_color <- c("#2EC4B6","#BDD5EA","#FFA5AB")

#####02 Main Script#####
seu <- readRDS("~/Working_folder/DataBase/LUAD/scData/refquery_final.rds")
meta.data <- seu@meta.data
table(meta.data$id)
####Figure 6A####
#原始注释图#
Figure_6A <- SCP::CellDimPlot(seu,group.by = "Cell_Cluster_level1",theme_use = "theme_blank")
dpis_genes <- c("TPX2","CYP4B1","SCGB3A1","CACNA2D2","UBE2C","SFTPB","MYBL2","CDC20","BIRC5","SUSD2")
Figure_S6 <- FeatureDimPlot(seu,feature = dpis_genes,reduction = "umap",theme_use = "theme_blank")
#pdf("FigureS6.pdf",width = 12,height = 9)
#Figure_S6
#dev.off()
#创建打分后映射#
dpis_genes <- c("TPX2","CYP4B1","SCGB3A1","CACNA2D2","UBE2C","SFTPB","MYBL2","CDC20","BIRC5","SUSD2")
seu <- AddModuleScore(
  seu,
  features = list(dpis_genes),
  name = "DPIS" 
)

Figure_6B <- FeatureDimPlot(seu,features = "DPIS1",reduction = "umap",theme_use = "theme_blank")
#分组#
seu@meta.data$DPIS_group <- ifelse(seu@meta.data$DPIS1 > 0, "DPIS_High","DPIS_LOW")
Figure_6C <- CellDimPlot(seu,group.by = "DPIS_group",reduction = "umap",theme_use = "theme_blank", palcolor = c("#FF7F00","#1F78B4"))
Figure_6D <- CellStatPlot(seu, stat.by = "DPIS_group", plot_type = "ring")
Figure_6E <- CellStatPlot(seu, stat.by = "Cell_Cluster_level1", group.by = "DPIS_group", label = TRUE)
#pdf("./Output/Script_6/Figure/Figure6ABC.pdf",width = 12,height=4.5)
#Figure_6A | Figure_6B | Figure_6C
#dev.off()
#pdf("./Output/Script_6/Figure/Figure6DE.pdf",width = 8,height=4.5)
#Figure_6D | Figure_6E
#dev.off()
#提取Cancer Cell#
seu_Cancer_Cell_remake <- subset(seu,Cell_Cluster_level2 %in% c("CDKN2A Cancer","CXCL1 Cancer","LAMC2 Cancer",
                                                              "Proliferating Cancer","SOX2 Cancer"))
seu_Cancer_Cell_remake <- seu_Cancer_Cell_remake %>%
  NormalizeData(.) %>%
  FindVariableFeatures(.) %>%
  ScaleData(.) %>%
  RunPCA(.) %>%
  RunHarmony(reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony") %>%
  RunUMAP(reduction = "harmony",1:30)
Figure_6F <- SCP::CellDimPlot(seu_Cancer_Cell_remake,group.by = "Cell_Cluster_level2",theme_use = "theme_blank")
Figure_6G <- FeatureDimPlot(seu_Cancer_Cell_remake,features = "DPIS1",theme_use = "theme_blank")
Figure_6H <- FeatureStatPlot(
  srt = seu_Cancer_Cell_remake, group.by = "Cell_Cluster_level2", bg.by = "Cell_Cluster_level2",
  stat.by = c("DPIS1"), add_box = TRUE
)
#pdf("./Output/Script_6/Figure/Figure6FGH.pdf",width = 12,height=4.5)
#Figure_6F | Figure_6G | Figure_6H
#dev.off()
seu_Cancer_Cell_remake@meta.data$Cell_Cluster_level2_orign <- seu_Cancer_Cell_remake@meta.data$Cell_Cluster_level2
seu_Cancer_Cell_remake@meta.data$Cell_Cluster_level2 <- ifelse(
  seu_Cancer_Cell_remake@meta.data$Cell_Cluster_level2 == "Proliferating Cancer" &
    seu_Cancer_Cell_remake@meta.data$DPIS_group == "DPIS_High",
  "DPIS+ Proliferating Cancer",
  seu_Cancer_Cell_remake@meta.data$Cell_Cluster_level2
)
seu_Cancer_Cell_remake@meta.data$Proliferating <- ifelse(
  seu_Cancer_Cell_remake@meta.data$Cell_Cluster_level2_orign == "Proliferating Cancer",
  "Proliferating Cancer","Others"
)
seu_Cancer_Cell_remake@meta.data$DPIS <- ifelse(
  seu_Cancer_Cell_remake@meta.data$DPIS_group == "DPIS_High",
  "DPIS Cells","Others"
)
seu_Cancer_Cell_remake@meta.data$DPIS_Proliferating <- ifelse(
  seu_Cancer_Cell_remake@meta.data$Cell_Cluster_level2 == "DPIS+ Proliferating Cancer",
  "DPIS+ Proliferating Cancer","Others"
)
Figure_6I <- SCP::CellDimPlot(seu_Cancer_Cell_remake,group.by = "Proliferating",theme_use = "theme_blank",palcolor = c("#FF7F00","#1F78B4"))
Figure_6J <- SCP::CellDimPlot(seu_Cancer_Cell_remake,group.by = "DPIS",theme_use = "theme_blank",palcolor = c("#FF7F00","#1F78B4"))
Figure_6K <- SCP::CellDimPlot(seu_Cancer_Cell_remake,group.by = "DPIS_Proliferating",theme_use = "theme_blank",palcolor = c("#FF7F00","#1F78B4"))
pdf("./Output/Script_6/Figure/Figure6IJK.pdf",width = 12,height=4.5)
Figure_6I | Figure_6J | Figure_6K
dev.off()
#saveRDS(seu_Cancer_Cell_remake,"./Output/Script_6/seu_SCENIC.rds")


