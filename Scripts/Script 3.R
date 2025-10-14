#####Loading Packages#####
library(MOVICS)
library(BiocOncoTK)
library(tidyverse)
library(maftools)
library(GSVA)
library(ComplexHeatmap) 
library(circlize)
library(viridis)
library(gplots)
library(data.table)
library(dplyr)
library(RTCGAToolbox)
library(limma)
library(ComplexHeatmap)
library(RColorBrewer)
library(clusterProfiler)
library(tibble)
library(ggplot2)
library(cowplot)
library(ggsci)
library(viridis)
library(patchwork)
library(TCellSI)
source("./R/annTrackScale.R")
source("./R/standarize_fun.R")
use_color <- c("#2EC4B6","#BDD5EA","#FFA5AB")
#####Data Loading#####
load("./Input/TCGA-LUAD/TCGA-LUAD_mrna_expr_tpm.rdata")
LUAD.TPM <- {
  x <- mrna_expr_tpm[, substr(colnames(mrna_expr_tpm),14,15) != "11", drop = FALSE]
  colnames(x) <- substr(colnames(x), 1, 16)
  x[, !duplicated(colnames(x)), drop = FALSE]
}
duplicated(colnames(LUAD.TPM))
LUAD.Survival <- data.table::fread("./Input/TCGA-LUAD/TCGA-LUAD.survival.tsv")
LUAD.Clinical <- data.table::fread("./Input/TCGA-LUAD/TCGA-LUAD.clinical.tsv") %>%
  dplyr::filter(tissue_type.samples == "Tumor")
LUAD.TPM <- log(LUAD.TPM + 1)
ImmuneSubtypeClass <- readRDS("./Output/Script_1/RData/ImmuneSubtypeClass.RDS")

#####Figure 2A ImmuneLandScape3#####
immune.signature <- read.table("Input/Analysis_data/Curated_Immune_Cell_Signature.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
# 构建计算GSVA的列表
cell.type <- unique(immune.signature$CellType)
immune.sig.ccr <- list()
for (i in cell.type) {
  immune.sig.ccr[[i]] <- immune.signature[which(immune.signature$CellType == i),"Symbol"]
}

# 免疫检查点相关基因
imm.targets <- c("CD274","PDCD1","CD247","PDCD1LG2","CTLA4","TNFRSF9","TNFRSF4","TLR9") 

# 免疫细胞的排序
immune.sig.ccr.order <- c("T.cells.CD8", 
                          "T.cells.regulatory..Tregs.",
                          "T.cells.CD4.naive",
                          "T.cells.follicular.helper",
                          "B.cells.naive",
                          "B.cells.memory",
                          "T.cells.gamma.delta",
                          "Dendritic.cells.activated",
                          "Macrophages.M1",
                          "NK.cells.activated",
                          "Plasma.cells",
                          "T.cells.CD4.memory.resting",
                          "T.cells.CD4.memory.activated",
                          "Mast.cells.activated",
                          "NK.cells.resting",
                          "Macrophages.M0",
                          "Macrophages.M2",
                          "Eosinophils",
                          "Monocytes",
                          "Dendritic.cells.resting",
                          "Mast.cells.resting",
                          "Neutrophils",
                          "Endothelial cells",
                          "Fibroblasts")

# 计算immune/stromal得分
indata <- LUAD.TPM
# 保存到文件
write.table(indata,file = "Output/Script_2/Data/TCGA_log2TPM_hugo.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
filterCommonGenes(input.f = "Output/Script_2/Data/TCGA_log2TPM_hugo.txt" , output.f="Output/Script_2/Data/TCGA_log2TPM_hugo_ESTIMATE.txt", id="GeneSymbol")
estimateScore("Output/Script_2/Data/TCGA_log2TPM_hugo_ESTIMATE.txt","Output/Script_2/Data/TCGA_log2TPM_hugo_estimate_score.txt", platform="affymetrix")
est.tcga <- read.table(file = "Output/Script_2/Data/TCGA_log2TPM_hugo_estimate_score.txt",header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est.tcga) <- est.tcga[,2]; colnames(est.tcga) <- est.tcga[1,]; est.tcga <- est.tcga[-1,c(-1,-2)];
est.tcga <- sapply(est.tcga, as.numeric); rownames(est.tcga) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity"); est.tcga.backup = as.data.frame(est.tcga); colnames(est.tcga.backup) <- colnames(expr)

# 对数值进行标准化并对极端值做截断
est.tcga <- annTrackScale(indata = est.tcga, halfwidth = 2, poolsd = F); est.tcga <- as.data.frame(t(est.tcga)) 
rownames(est.tcga) <- colnames(LUAD.TPM)

tcga.immune.gsva <- gsva(as.matrix(LUAD.TPM),
                         immune.sig.ccr,
                         method = "gsva")

#figure
plot.data.2A <- LUAD.Clinical %>%
  dplyr::filter(tissue_type.samples == "Tumor") %>%
  dplyr::select("sample","vital_status.demographic","ajcc_pathologic_stage.diagnoses",
                "ajcc_pathologic_t.diagnoses","ajcc_pathologic_n.diagnoses","ajcc_pathologic_m.diagnoses") %>%
  setNames(c("ID","Status","Stage","T_Stage","N_Stage","M_Stage")) %>%
  dplyr::mutate(Stage = case_when(
    Stage %in% "Stage I" ~ "Stage I",
    Stage %in% "Stage IA" ~ "Stage I",
    Stage %in% "Stage IB" ~ "Stage I",
    Stage %in% "Stage II" ~ "Stage II",
    Stage %in% "Stage IIA " ~ "Stage II",
    Stage %in% "Stage IIB" ~ "Stage II",
    Stage %in% "Stage IIIA" ~ "Stage III",
    Stage %in% "Stage IIIB" ~ "Stage III",
    Stage %in% "Stage IV" ~ "Stage IV",
    is.na(Stage) ~ "not reported")) %>%
  dplyr::mutate(T_Stage = case_when(
    T_Stage %in% "T1" ~ "T1",
    T_Stage %in% "T1a" ~ "T1",
    T_Stage %in% "T1b" ~ "T1",
    T_Stage %in% "T2" ~ "T2",
    T_Stage %in% "T2a" ~ "T2",
    T_Stage %in% "T2b" ~ "T2",
    T_Stage %in% "T3" ~ "T3",
    T_Stage %in% "T3a" ~ "T3",
    T_Stage %in% "T3b" ~ "T3",
    T_Stage %in% "T3c" ~ "T3",
    T_Stage %in% "T4" ~ "T4",
    is.na(T_Stage) ~ "not reported"))  %>%
  dplyr::mutate(M_Stage = case_when(
    M_Stage %in% "M0" ~ "M0",
    M_Stage %in% "M1" ~ "M1",
    M_Stage %in% "M1a" ~ "M1",
    M_Stage %in% "M1b" ~ "M1",
    M_Stage %in% "MX" ~ "MX",
    is.na(M_Stage) ~ "not reported")) %>%
  right_join(
    ImmuneSubtypeClass %>% as.data.frame(),
    by = "ID"
  ) %>%
  column_to_rownames(var = "ID") %>%
  mutate(across(everything(), ~ replace_na(as.character(.x), "Unknown")))
columnAnno <- HeatmapAnnotation(
  Cluster = plot.data.2A$Cluster,
  Status = plot.data.2A$Status,
  Stage = plot.data.2A$Stage,
  T_Stage = plot.data.2A$T_Stage,
  M_Stage = plot.data.2A$M_Stage,
  N_Stage = plot.data.2A$N_Stage,
  col = list(
    Cluster = c("Wound Healing" = "#2EC4B6","IFN-γ Dominant"="#BDD5EA","Inflammatory" = "#FFA5AB"),
    Status = c("Alive"  = "#A8817A",
               "Dead"    = "#E8BE74"),
    Stage  = c("Stage I" = "#8ab1d2",
               "Stage II"   = "#E58579",
               "Stage III"   = "#D9BDD8",
               "Stage IV"   = "#9180AC",
               "Unknown" = "#999999"),
    T_Stage = c("T1" = "#FF9F1C",
                "T2" = "#FFA5AB",
                "T3" = "#023E8A",
                "T4" = "#9D4EDD",
                "Unknown" = "#999999"),
    N_Stage = c("N0" = "#E64B35FF",
                "N1" = "#4DBBD5FF",
                "N2" = "#00A087FF",
                "N3" = "#8491B4FF",
                "NX" = "#3C5488FF",
                "Unknown" = "#999999"),
    M_Stage = c("M0" = "#8491B4FF",
                "M1" = "#91D1C2FF",
                "MX" = "#DC0000FF",
                "Unknown" = "#999999")
  )
  
)
plot.data.2A <- plot.data.2A %>%
  mutate(Cluster = factor(Cluster, 
                          levels = c("Wound Healing", "IFN-γ Dominant", "Inflammatory"))) %>%
  arrange(Cluster)
split <- factor(plot.data.2A$Cluster, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory"))
# 设置颜色
clust.col <- use_color
heatmap.color <- c("#2dabb9", "white","#FF9F1C" )
blue <- "#5bc0eb"
gold <- "#ECE700"
cyan <- "#00B3D0"

plot.data.2A$ImmuneScore <- as.numeric(est.tcga[rownames(plot.data.2A),"ImmuneScore"])
plot.data.2A$StromalScore <- as.numeric(est.tcga[rownames(plot.data.2A),"StromalScore"])
annColors <- list() # 构建热图的图例颜色
annColors[["Cluster"]] <- c("Wound Healing" = clust.col[1],
                            "IFN-γ Dominant" = clust.col[2],
                            "Inflammatory" = clust.col[3]
)
annColors[["ImmuneScore"]] <- inferno(64)
annColors[["StromalScore"]] <- viridis(64)
## 热图1：免疫检查点基因表达
indata <- LUAD.TPM[intersect(rownames(LUAD.TPM),imm.targets),]

hm1 <- pheatmap(standarize.fun(indata[imm.targets,rownames(plot.data.2A)],halfwidth = 2), # 表达谱数据标准化
                border_color = NA, # 热图单元格无边框
                annotation_col = plot.data.2A[,c("Cluster","StromalScore","ImmuneScore")],
                annotation_colors = annColors[c("Cluster","StromalScore","ImmuneScore")],
                color = NMF:::ccRamp(x = heatmap.color,n = 64),
                show_rownames = T, # 显示行名
                column_split = split,
                show_colnames = F, # 不显示列名
                cellheight = 12, # 热图高度固定
                cellwidth = 0.6, # 热图宽度固定
                name = "ICI", # 图例名字
                cluster_rows = F, # 行不聚类
                cluster_cols = F) # 列不聚类

#pdf("CheckPoint_heatmap.pdf",width = 10,height = 10)
hm1
#dev.off()
## 热图2：肿瘤免疫微环境富集得分
colnames(tcga.immune.gsva) <- substr(colnames(tcga.immune.gsva),1,16)
hm2 <- pheatmap(standarize.fun(tcga.immune.gsva[immune.sig.ccr.order,rownames(plot.data.2A)],halfwidth = 1), # 富集得分标准化并排序
                border_color = NA,
                color = NMF:::ccRamp(x = heatmap.color,n = 64),
                column_split = split,
                gaps_row = c(14,22), # 根据不同类别的免疫细胞分割
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 0.6,
                name = "TIME",
                cluster_rows = F,
                cluster_cols = F)

#pdf("TIMEheatmap.pdf",width = 10,height = 10)
hm2
#dev.off()
draw(hm1 %v% hm2 ,
     heatmap_legend_side = "right", 
     annotation_legend_side = "right")
#pdf("Output/Script_2/Figure/Figure_2ATIME.pdf",width = 10,height = 10)
draw(hm1 %v% hm2 , 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right")
#invisible(dev.off())





######Figure 2 B-I#####
TcellSI.res <- TCellSI::TCSS_Calculate(LUAD.TPM, ref = TRUE)
plot.data.Tcell <- ImmuneSubtypeClass %>%
  left_join(
    TcellSI.res %>%
      t() %>% 
      as.data.frame() %>%
      rownames_to_column(var = "ID"),
    by = "ID"
      
  )
#saveRDS(TcellSI.res,"./Output/Script_2/RDATA/TcellSI_res.RDS")
plot.list.TcellSI <- list()
Tcell.state <- rownames(TcellSI.res)
for(state in Tcell.state){
  Figrue.Tcell <- ggplot(plot.data.Tcell,aes(Cluster,.data[[state]],fill= Cluster))+
    geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
    geom_boxplot(width=0.4,size=1.2,outlier.color =NA)+
    geom_jitter(aes(fill= Cluster),shape=21,size=2.5,width=0.2,alpha = 0.5)+
    geom_signif(
      comparisons = list(c("Wound Healing","IFN-γ Dominant"), c("IFN-γ Dominant","Inflammatory"),c("Wound Healing","Inflammatory")),
      map_signif_level = TRUE, test = t.test, step_increase = 0.08
    )+
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.border = element_rect(size = 2),
          axis.text.x = element_text(color = "black", size = 13),
          axis.text.y = element_text(color = "black",size = 13),
          legend.position = "none",
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.ticks = element_line(color="black",linewidth = 1))+
    labs(x=NULL,y= paste0(state,"Score")) +
    scale_fill_manual(values = use_color) +
    theme(axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1))
  plot.list.TcellSI[[state]] <- Figrue.Tcell
}
pdf("./Output/Script_2/Figure/FigureE_L.pdf",width = 12,height = 11)
wrap_plots(plot.list.TcellSI, ncol=4)
dev.off()

######Figure 2 B C D TMB FGA Mutcount MSI#####
###TMB###
load("./Input/TCGA-LUAD/TCGA-LUAD_maf.rdata")
TMB_data <- data
TMB_maf <- read.maf(TMB_data)
LUAD.tmb <- tmb(TMB_maf, captureSize = 38, logScale = T)
plot.data.2B <- ImmuneSubtypeClass %>%
  left_join(LUAD.tmb %>%
              rename("ID" = `Tumor_Sample_Barcode`) %>%
              mutate(ID = substr(ID,1,16)),
            by = "ID")

Figure2B_TMB <- ggplot(plot.data.2B,aes(Cluster,`total_perMB_log`,fill= Cluster))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
  geom_boxplot(width=0.4,size=1.2,outlier.color =NA)+
  geom_jitter(aes(fill= Cluster),shape=21,size=2.5,width=0.2,alpha = 0.5)+
  geom_signif(
    comparisons = list(c("Wound Healing","IFN-γ Dominant"), c("IFN-γ Dominant","Inflammatory"),c("Wound Healing","Inflammatory")),
    map_signif_level = TRUE, test = t.test, step_increase = 0.08
  )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.text.x = element_text(color = "black", size = 13),
        axis.text.y = element_text(color = "black",size = 13),
        legend.position = "none",
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.ticks = element_line(color="black",linewidth = 1))+
  labs(x=NULL,y="Total PerMB log") +
  scale_fill_manual(values = use_color) +
  theme(axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1))

###FGA###
FAG.dat <- data.table::fread("./Input/Analysis_data/Fraction_Genome_Altered.txt")
plot.data.2C <- ImmuneSubtypeClass %>%
  mutate(ID = substr(ID,1,15)) %>%
  left_join(FAG.dat %>%
              rename("ID" = `Sample ID`),
            by = "ID")

Figure2C_FGA <- ggplot(plot.data.2C,aes(Cluster,`Fraction Genome Altered`,fill= Cluster))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
  geom_boxplot(width=0.4,size=1.2,outlier.color =NA)+
  geom_jitter(aes(fill= Cluster),shape=21,size=2.5,width=0.2,alpha = 0.5)+
  geom_signif(
    comparisons = list(c("Wound Healing","IFN-γ Dominant"), c("IFN-γ Dominant","Inflammatory"),c("Wound Healing","Inflammatory")),
    map_signif_level = TRUE, test = t.test, step_increase = 0.08
  )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.text.x = element_text(color = "black", size = 13),
        axis.text.y = element_text(color = "black",size = 13),
        legend.position = "none",
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.ticks = element_line(color="black",linewidth = 1))+
  labs(x=NULL,y="Fraction Genome Altered") +
  scale_fill_manual(values = use_color) +
  theme(axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1))

###Mut Count###
Mutcount.dat <- data.table::fread("./Input/Analysis_data/Mutation_Count.txt")
plot.data.2D <- ImmuneSubtypeClass %>%
  mutate(ID = substr(ID,1,15)) %>%
  left_join(Mutcount.dat %>%
              rename("ID" = `Sample ID`),
            by = "ID")

Figure2D_Mutcount <- ggplot(plot.data.2D,aes(Cluster,`Mutation Count`,fill= Cluster))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
  geom_boxplot(width=0.4,size=1.2,outlier.color =NA)+
  geom_jitter(aes(fill= Cluster),shape=21,size=2.5,width=0.2,alpha = 0.5)+
  geom_signif(
    comparisons = list(c("Wound Healing","IFN-γ Dominant"), c("IFN-γ Dominant","Inflammatory"),c("Wound Healing","Inflammatory")),
    map_signif_level = TRUE, test = t.test, step_increase = 0.08
  )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.text.x = element_text(color = "black", size = 13),
        axis.text.y = element_text(color = "black",size = 13),
        legend.position = "none",
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.ticks = element_line(color="black",linewidth = 1))+
  labs(x=NULL,y="Mutation Count") +
  scale_fill_manual(values = use_color) +
  theme(axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1))

###MSI###
data("MSIsensor.10k")
plot.data.S2D <- ImmuneSubtypeClass %>%
  mutate(participant_barcode = substr(ID,1,12)) %>%
  left_join(MSIsensor.10k,by = "participant_barcode")

FigureS2D_MSI <- ggplot(plot.data.S2D,aes(Cluster,MSIsensor.score,fill= Cluster))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
  geom_boxplot(width=0.4,size=1.2,outlier.color =NA)+
  geom_jitter(aes(fill= Cluster),shape=21,size=2.5,width=0.2,alpha = 0.5)+
  geom_signif(
    comparisons = list(c("Wound Healing","IFN-γ Dominant"), c("IFN-γ Dominant","Inflammatory"),c("Wound Healing","Inflammatory")),
    map_signif_level = TRUE, test = t.test, step_increase = 0.08
  )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.text.x = element_text(color = "black", size = 13),
        axis.text.y = element_text(color = "black",size = 13),
        legend.position = "none",
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.ticks = element_line(color="black",linewidth = 1))+
  labs(x=NULL,y="MSI Score") +
  scale_fill_manual(values = use_color) +
  theme(axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1))
pdf("./Output/Script_2/Figure/FigureB_D.pdf",width = 12,height = 5.5)
Figure2B_TMB | Figure2C_FGA | Figure2D_Mutcount | FigureS2D_MSI
dev.off()
