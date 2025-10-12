#####01 Loading Packages#####
Sys.setenv(LANG = 'EN')
library(tidyverse)
library(tibble)
library(ConsensusClusterPlus)
library(pheatmap)
library(dendsort)
library(ggplot2)
library(ggsci)
library(scales)
library(ComplexHeatmap)
library(survival)
library(survminer)
library(circlize)
library(MOVICS)
library(IOBR)
library(ggplot2)
library(ggsignif) 
library(gghalves) 
library(tidyverse)
library(factoextra)
library(ggConvexHull)
library(TMEclassifier)
library(plotrix)
library(ggpubr)
library(ComplexHeatmap)
library(grid)
source("./R/standarize_fun.R")
use_color <- c("#2EC4B6","#BDD5EA","#FFA5AB")
#####02 Prepare data#####
load("./Input/TCGA-LUAD/TCGA-LUAD_mrna_expr_tpm.rdata")
LUAD.TPM <- {
  x <- mrna_expr_tpm[, substr(colnames(mrna_expr_tpm),14,15) != "11", drop = FALSE]
  colnames(x) <- substr(colnames(x), 1, 16)
  x[, !duplicated(colnames(x)), drop = FALSE]
}
#duplicated(colnames(LUAD.TPM))
LUAD.Survival <- data.table::fread("./Input/TCGA-LUAD/TCGA-LUAD.survival.tsv")
LUAD.Clinical <- data.table::fread("./Input/TCGA-LUAD/TCGA-LUAD.clinical.tsv") %>%
  dplyr::filter(tissue_type.samples == "Tumor")
LUAD.TPM <- log(LUAD.TPM + 1) 

#####03 TMEclassifier#####
##TMEclassifier##
LUAD_TMEclassifier <- tme_classifier(LUAD.TPM, method = "ensemble", scale = T,tme_deconvolution = F)
table(LUAD_TMEclassifier$TMEcluster)
LUAD_TMEclassifier <- LUAD_TMEclassifier %>%
  column_to_rownames(var = "ID")
##提取特征样本##
delta <- 0.3
# 计算每一类相对其它两类的领先幅度
gap_IE <- LUAD_TMEclassifier$IE - pmax(LUAD_TMEclassifier$IS, LUAD_TMEclassifier$IA, na.rm = TRUE)
gap_IS <- LUAD_TMEclassifier$IS - pmax(LUAD_TMEclassifier$IE, LUAD_TMEclassifier$IA, na.rm = TRUE)
gap_IA <- LUAD_TMEclassifier$IA - pmax(LUAD_TMEclassifier$IE, LUAD_TMEclassifier$IS, na.rm = TRUE)

# 保留规则：TMEcluster 指向哪类，就用对应 gap 判断是否 ≥ 0.3
keep_idx <- (LUAD_TMEclassifier$TMEcluster == "IE" & gap_IE >= delta) |
  (LUAD_TMEclassifier$TMEcluster == "IS" & gap_IS >= delta) |
  (LUAD_TMEclassifier$TMEcluster == "IA" & gap_IA >= delta)
LUAD_TMEclassifier_keep <- LUAD_TMEclassifier[keep_idx, , drop = FALSE]
LUAD_TMEclassifier_keep$TMEcluster <- factor(LUAD_TMEclassifier_keep$TMEcluster,levels = c("IE","IS","IA"))
table(LUAD_TMEclassifier_keep$TMEcluster)
LUAD.Survival <- dplyr::filter(LUAD.Survival,sample %in% rownames(LUAD_TMEclassifier_keep))
colnames(LUAD.Survival)[c(1,2,3)]  <- c("ID","time","status")
LUAD.TPM <- dplyr::select(LUAD.TPM,LUAD.Survival$ID)
tmescore <- tmescore(eset   = LUAD.TPM,#expression data
                     pdata  = LUAD.Survival,
                     method  = "PCA",#default
                     classify = T)#if true, survival data must be provided in pdata

LUAD_TMEclassifier_fina <- LUAD_TMEclassifier_keep %>%
  rownames_to_column(var = "ID") %>%
  left_join(tmescore,by ="ID") %>%
  na.omit(.) %>%
  filter(
    (TMEcluster == "IA" & TMEscore_binary == "High") |
      (TMEcluster == "IE" & TMEscore_binary == "Low")  |
      (TMEcluster == "IS" & TMEscore_binary == "Low")
  )
  
######Figure 1A 分类饼状图 云雨图 #####
plot.data.1A_b <- c(round(sum(LUAD_TMEclassifier_fina$TMEcluster == "IA")/nrow(LUAD_TMEclassifier_fina),2),
                    round(sum(LUAD_TMEclassifier_fina$TMEcluster == "IE")/nrow(LUAD_TMEclassifier_fina),2),
                    round(sum(LUAD_TMEclassifier_fina$TMEcluster == "IS")/nrow(LUAD_TMEclassifier_fina),2))
labels <- c(
  paste0(
    "IA: ", round(sum(LUAD_TMEclassifier_fina$TMEcluster == "IA")/nrow(LUAD_TMEclassifier_fina)*100, 2), "%", 
    "\n(", sum(LUAD_TMEclassifier_fina$TMEcluster == "IA"), ")"
  ),
  paste0(
    "IE: ", round(sum(LUAD_TMEclassifier_fina$TMEcluster == "IE")/nrow(LUAD_TMEclassifier_fina)*100, 2), "%", 
    "\n(", sum(LUAD_TMEclassifier_fina$TMEcluster == "IE"), ")"
  ),
  paste0(
    "IS: ", round(sum(LUAD_TMEclassifier_fina$TMEcluster == "IS")/nrow(LUAD_TMEclassifier_fina)*100, 2), "%", 
    "\n(", sum(LUAD_TMEclassifier_fina$TMEcluster == "IS"), ")"
  )
)

Figure1A_b <- pie3D(plot.data.1A_b, labels = labels, explode = 0.05, col = use_color, theta =0.8)
pdf("./Output/Script_1/Figure_1/Figure1A_b.pdf",width = 4,height = 4)
Figure1A_b
dev.off()
plot.data.1A <- LUAD_TMEclassifier_fina
IA <- ggplot(plot.data.1A,aes(TMEcluster,IA,fill=TMEcluster))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
  geom_boxplot(width=0.4,size=1.2,outlier.color =NA)+
  geom_jitter(aes(fill=TMEcluster),shape=21,size=2.5,width=0.2,alpha = 0.5)+
  geom_signif(
    comparisons = list(c("IE","IA"), c("IS","IA")),
    map_signif_level = TRUE, test = t.test, step_increase = 0.08
  )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.text.x = element_text(color = "black", size = 13),
        axis.text.y = element_text(color = "black",size = 13),
        legend.position = "none",
        axis.title.y = element_text(size = 20),
        axis.ticks = element_line(color="black",linewidth = 1))+
  labs(x=NULL,y="Immune Activated Score") +
  scale_fill_manual(values = use_color)
IE <- ggplot(plot.data.1A,aes(TMEcluster,IE,fill=TMEcluster))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
  geom_boxplot(width=0.4,size=1.2,outlier.color =NA)+
  geom_jitter(aes(fill=TMEcluster),shape=21,size=2.5,width=0.2,alpha = 0.5)+
  geom_signif(
    comparisons = list(c("IE","IS"), c("IE","IA")),
    map_signif_level = TRUE, test = t.test, step_increase = 0.08
  )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.text.x = element_text(color = "black", size = 13),
        axis.text.y = element_text(color = "black",size = 13),
        legend.position = "none",
        axis.title.y = element_text(size = 20),
        axis.ticks = element_line(color="black",linewidth = 1))+
  labs(x=NULL,y="Immune Exclusive Score") +
  scale_fill_manual(values = use_color)

IS <- ggplot(plot.data.1A,aes(TMEcluster,IS,fill=TMEcluster))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
  geom_boxplot(width=0.4,size=1.2,outlier.color =NA)+
  geom_jitter(aes(fill=TMEcluster),shape=21,size=2.5,width=0.2,alpha = 0.5)+
  geom_signif(
    comparisons = list(c("IE","IS"), c("IS","IA")),
    map_signif_level = TRUE, test = t.test, step_increase = 0.08
  )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.text.x = element_text(color = "black", size = 13),
        axis.text.y = element_text(color = "black",size = 13),
        legend.position = "none",
        axis.title.y = element_text(size = 20),
        axis.ticks = element_line(color="black",linewidth = 1))+
  labs(x=NULL,y="Immune Suppressive Score") +
  scale_fill_manual(values = use_color)
IA + IE + IS
Figure_1A <- IA + IE + IS
ggsave("./Output/Script_1/Figure_1/Figrue1A.pdf",Figure_1A,width = 9,height = 5)

######Figure 1B 临床特征分布#####
plot.data.1B <- LUAD.Clinical %>%
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
    LUAD_TMEclassifier_fina %>% as.data.frame(),
    by = "ID"
  ) %>%
  dplyr::select(-c("IA","IE","IS")) %>%
  column_to_rownames(var = "ID") %>%
  mutate(across(everything(), ~ replace_na(as.character(.x), "Unknown")))
show_col(pal_npg("nrc", alpha = 1)(10))
columnAnno <- HeatmapAnnotation(
                                TMECluster = plot.data.1B$TMEcluster,
                                Status = plot.data.1B$Status,
                                Stage = plot.data.1B$Stage,
                                T_Stage = plot.data.1B$T_Stage,
                                M_Stage = plot.data.1B$M_Stage,
                                N_Stage = plot.data.1B$N_Stage,
                                col = list(
                                  TMECluster = c("IA" = "#2EC4B6","IE"="#BDD5EA","IS" = "#FFA5AB"),
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

n <- nrow(plot.data.1B)
mat.plot.data.1B <- matrix(0, nrow = 1, ncol = n)
rownames(mat.plot.data.1B) <- "gene"

# 给列名：优先用行名；没有的话用 ID；再没有就用 S1..Sn
if (!is.null(rownames(plot.data.1B))) {
  colnames(mat.plot.data.1B) <- rownames(plot.data.1B)
} else if ("ID" %in% colnames(plot.data.1B)) {
  colnames(mat.plot.data.1B) <- plot.data.1B$ID
} else {
  colnames(mat.plot.data.1B) <- paste0("S", seq_len(n))
}
split <- factor(plot.data.1B$TMEcluster, levels = c("IA","IE","IS"))
ord   <- c(which(plot.data.1B$TMEcluster == "IE"),
           which(plot.data.1B$TMEcluster == "IS"),
           which(plot.data.1B$TMEcluster == "IA"))
Figure_1B <- ComplexHeatmap::Heatmap(
  mat.plot.data.1B,
  col = c("0" = "white"),
  name = NULL,
  na_col = "white",
  show_column_names = FALSE,
  show_row_names = FALSE,
  row_names_side = "left",
  # 按 IE/IS/IA 分块并指定块顺序
  column_split = split,
  column_order = ord,
  cluster_columns = FALSE,     
  cluster_column_slices = FALSE,
  gap = unit(2, "mm"),
  top_annotation = columnAnno,
  height = unit(0.1, "mm"),
  border = FALSE               # 如果想要格子线，改 TRUE
)
pdf("./Output/Script_1/Figure_1/Figrue1B.pdf",width = 8,height = 4)
Figure_1B
dev.off()

######Figure 1C TME_Classifier#####
##Status##
plot.data.1C_Status <- plot.data.1B %>%
  select(Status, TMEcluster) %>%
  group_by(TMEcluster, Status) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(TMEcluster) %>%
  mutate(Percentage = Count / sum(Count) * 100)

tab <- table(plot.data.1B$TMEcluster, plot.data.1B$Status)
chisq.test(tab)   # 卡方检验
# 或者
fisher.test(tab) 
pval <- chisq.test(tab)$p.value


Figure_1C_status <- ggplot(plot.data.1C_Status, aes(x = TMEcluster, y = Percentage, fill = Status)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),  
            color = "black", size = 4) +
  scale_fill_manual(values = c("Alive"= "#E64B35CC", "Dead"= "#4DBBD5CC")) +
  labs(x = "", y = "Percentage (%)", fill = "Status") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1),
    axis.text.y  = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold",size = 13)
  ) +
  annotate("text", x = 2, y = 105, 
           label = paste0("italic(P) == ", signif(pval, 3)),
           parse = TRUE,  
           face = "bold",
           fontface = "italic",size = 5)
##Stage##
plot.data.1C_Stage <- plot.data.1B %>%
  select(Stage, TMEcluster) %>%
  dplyr::filter(Stage != "Unknown") %>%
  group_by(TMEcluster, Stage) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(TMEcluster) %>%
  mutate(Percentage = Count / sum(Count) * 100)
tab.dat <- plot.data.1B %>%
  select(Stage, TMEcluster) %>%
  dplyr::filter(Stage != "Unknown")
tab <- table(tab.dat$TMEcluster, tab.dat$Stage)
chisq.test(tab)   # 卡方检验
# 或者
pval <- chisq.test(tab)$p.value


Figure.1C_Stage  <- ggplot(plot.data.1C_Stage, aes(x = TMEcluster, y = Percentage, fill = Stage)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),  
            color = "black", size = 4) +
  scale_fill_manual(values = c("Stage I"= "#E64B35CC", "Stage II"= "#4DBBD5CC",
                               "Stage III" = "#00A087CC","Stage IV" = "#3C5488CC"
                               )) +
  labs(x = "", y = "Percentage (%)", fill = "Stage") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1),
    axis.text.y  = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold",size = 13)
  ) +
  annotate("text", x = 2, y = 105, 
           label = paste0("italic(P) == ", signif(pval, 3)),
           parse = TRUE,  
           face = "bold",
           fontface = "italic",size = 5)


##T_Stage##
plot.data.1C_T_Stage <- plot.data.1B %>%
  select(T_Stage, TMEcluster) %>%
  dplyr::filter(T_Stage != "Unknown") %>%
  group_by(TMEcluster, T_Stage) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(TMEcluster) %>%
  mutate(Percentage = Count / sum(Count) * 100)
tab.dat <- plot.data.1B %>%
  select(T_Stage, TMEcluster) %>%
  dplyr::filter(T_Stage != "Unknown")
tab <- table(plot.data.1B$TMEcluster, plot.data.1B$T_Stage)
chisq.test(tab)   # 卡方检验
pval <- chisq.test(tab)$p.value


Figure.1C_T_Stage  <- ggplot(plot.data.1C_T_Stage, aes(x = TMEcluster, y = Percentage, fill = T_Stage)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),  
            color = "black", size = 4) +
  scale_fill_manual(values = c("T1"= "#E64B35CC", "T2"= "#4DBBD5CC",
                               "T3" = "#00A087CC","T4" = "#3C5488CC"
  )) +
  labs(x = "", y = "Percentage (%)", fill = "T_Stage") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1),
    axis.text.y  = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold",size = 13)
  ) +
  annotate("text", x = 2, y = 105, 
           label = paste0("italic(P) == ", signif(pval, 3)),
           parse = TRUE,  
           face = "bold",
           fontface = "italic",size = 5)


##M_Stage##
plot.data.1C_M_Stage <- plot.data.1B %>%
  select(M_Stage, TMEcluster) %>%
  dplyr::filter(M_Stage != "Unknown") %>%
  group_by(TMEcluster, M_Stage) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(TMEcluster) %>%
  mutate(Percentage = Count / sum(Count) * 100)
tab.dat <- plot.data.1B %>%
  select(M_Stage, TMEcluster) %>%
  dplyr::filter(M_Stage != "Unknown")
tab <- table(plot.data.1B$TMEcluster, plot.data.1B$M_Stage)
chisq.test(tab)   # 卡方检验
pval <- chisq.test(tab)$p.value


Figure.1C_M_Stage  <- ggplot(plot.data.1C_M_Stage, aes(x = TMEcluster, y = Percentage, fill = M_Stage)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),  
            color = "black", size = 4) +
  scale_fill_manual(values = c("M0"= "#E64B35CC", "M1"= "#4DBBD5CC",
                               "MX" = "#00A087CC"
  )) +
  labs(x = "", y = "Percentage (%)", fill = "M_Stage") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1),
    axis.text.y  = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold",size = 13)
  ) +
  annotate("text", x = 2, y = 105, 
           label = paste0("italic(P) == ", signif(pval, 3)),
           parse = TRUE,  
           face = "bold",
           fontface = "italic",size = 5)

##N_Stage##
plot.data.1C_N_Stage <- plot.data.1B %>%
  select(N_Stage, TMEcluster) %>%
  dplyr::filter(N_Stage != "Unknown") %>%
  group_by(TMEcluster, N_Stage) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(TMEcluster) %>%
  mutate(Percentage = Count / sum(Count) * 100)
tab.dat <- plot.data.1B %>%
  select(N_Stage, TMEcluster) %>%
  dplyr::filter(N_Stage != "Unknown")
tab <- table(plot.data.1B$TMEcluster, plot.data.1B$N_Stage)
chisq.test(tab)   # 卡方检验
pval <- chisq.test(tab)$p.value


Figure.1C_N_Stage  <- ggplot(plot.data.1C_N_Stage, aes(x = TMEcluster, y = Percentage, fill = N_Stage)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),  
            color = "black", size = 4) +
  scale_fill_manual(values = c("N0"= "#E64B35CC", "N1"= "#4DBBD5CC",
                               "N2" = "#00A087CC","N3" = "#8491B4FF","NX" = "#3C5488CC"
  )) +
  labs(x = "", y = "Percentage (%)", fill = "N_Stage") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1),
    axis.text.y  = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold",size = 13)
  ) +
  annotate("text", x = 2, y = 105, 
           label = paste0("italic(P) == ", signif(pval, 3)),
           parse = TRUE,  
           face = "bold",
           fontface = "italic",size = 5)

Figure_1C <- (Figure.1C_Stage + Figure.1C_T_Stage) / (Figure.1C_M_Stage +Figure.1C_N_Stage)
ggsave("./Output/Script_1/Figure_1/Figure1C.pdf",Figure_1C,width = 8,height = 8)
Figure_1C_status

######Figure 1D Survival#####
IOBR_esmitate <- deconvo_estimate(LUAD.TPM)
IOBR_esmitate.dat <- LUAD_TMEclassifier_fina %>%
  left_join(IOBR_esmitate,by ="ID")
IOBR_esmitate.plot <- ggplot(IOBR_esmitate.dat,aes(TMEcluster,StromalScore_estimate,fill=TMEcluster))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
  geom_boxplot(width=0.4,size=1.2,outlier.color =NA)+
  geom_jitter(aes(fill=TMEcluster),shape=21,size=2.5,width=0.2,alpha = 0.5)+
  geom_signif(
    comparisons = list(c("IE","IS"), c("IS","IA"),c("IA","IE")),
    map_signif_level = TRUE, test = t.test, step_increase = 0.08
  )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.text.x = element_text(color = "black", size = 13),
        axis.text.y = element_text(color = "black",size = 13),
        legend.position = "none",
        axis.title.y = element_text(size = 20),
        axis.ticks = element_line(color="black",linewidth = 1))+
  labs(x=NULL,y="StromalScore_estimate") +
  scale_fill_manual(values = use_color)



tmescore <- tmescore(eset   = LUAD.TPM,#expression data
                     method  ="PCA",#default
                     classify =F)#if true, survival data must be provided in pdata
tmescore.dat <- LUAD_TMEclassifier_fina %>%
  rownames_to_column(var = "ID") %>%
  left_join(tmescore,by ="ID")

ggplot(tmescore.dat,aes(TMEcluster,TMEscore,fill=TMEcluster))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
  geom_boxplot(width=0.4,size=1.2,outlier.color =NA)+
  geom_jitter(aes(fill=TMEcluster),shape=21,size=2.5,width=0.2,alpha = 0.5)+
  geom_signif(
    comparisons = list(c("IE","IS"), c("IS","IA"),c("IA","IE")),
    map_signif_level = TRUE, test = t.test, step_increase = 0.08
  )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.text.x = element_text(color = "black", size = 13),
        axis.text.y = element_text(color = "black",size = 13),
        legend.position = "none",
        axis.title.y = element_text(size = 20),
        axis.ticks = element_line(color="black",linewidth = 1))+
  labs(x=NULL,y="ImmuneScore") +
  scale_fill_manual(values = use_color)


plot.data.1D <- LUAD_TMEclassifier_fina %>%
  left_join(LUAD.Survival, by = "ID") %>%
  filter(time.x > 30 & time.x < 3650)
surv_cluster(plot.data.1D,
             target_group = "TMEcluster",
             ID = "ID",
             levels = c("IA","IE","IS"),
             reference_group = "IS",
             time = "time.x",
             status = "status.x",
             time_type = "day")

