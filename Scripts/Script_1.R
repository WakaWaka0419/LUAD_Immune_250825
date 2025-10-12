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
library(ggplot2)
library(ggsignif) 
library(gghalves) 
library(tidyverse)
library(factoextra)
library(ggConvexHull)
library(plotrix)
library(ggpubr)
library(ComplexHeatmap)
library(grid)
library(gridExtra)
library(ggpp)
library(patchwork)
library(IOBR)
library(TMEscore)
library(ImmuneSubtypeClassifier)
source("./R/standarize_fun.R")
use_color <- c("#2EC4B6","#BDD5EA","#FFA5AB")

#####02 Prepare data#####
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
ImmuneSubtypeClass <- ImmuneSubtypeClassifier::callEnsemble(X = LUAD.TPM, geneids = 'symbol')
calls <- geneMatchErrorReport(X =LUAD.TPM, geneid='symbol')
match_calls <- geneMatch(X =LUAD.TPM, geneid='symbol')
colnames(ImmuneSubtypeClass)[1:2] <- c("ID","Cluster")
colnames(ImmuneSubtypeClass)[3:5] <- c("Wound Healing","IFN-γ Dominant","Inflammatory")
ImmuneSubtypeClass <- dplyr::filter(ImmuneSubtypeClass,Cluster %in% c("1","2","3")) %>%
  mutate(Cluster = factor(case_when(
    Cluster == "1" ~ "Wound Healing",
    Cluster == "2" ~ "IFN-γ Dominant",
    Cluster == "3" ~ "Inflammatory",
    is.na(Cluster) ~ "not reported"
  ),levels = c("Wound Healing","IFN-γ Dominant","Inflammatory"))
  ) %>%
  mutate(Max_Score = apply(select_if(., is.numeric), 1, max, na.rm = TRUE)) %>%
  filter(Max_Score > 0.6)

######Figure 1A 分类饼状图 AND PCA#####
plot.data.1A_b <- c(round(sum(ImmuneSubtypeClass$Cluster == "Wound Healing")/nrow(ImmuneSubtypeClass),2),
                    round(sum(ImmuneSubtypeClass$Cluster == "IFN-γ Dominant")/nrow(ImmuneSubtypeClass),2),
                    round(sum(ImmuneSubtypeClass$Cluster == "Inflammatory")/nrow(ImmuneSubtypeClass),2))
labels <- c(
  paste0(
    "Wound Healing: ", round(sum(ImmuneSubtypeClass$Cluster == "Wound Healing")/nrow(ImmuneSubtypeClass)*100, 2), "%", 
    "\n(", sum(ImmuneSubtypeClass$Cluster == "Wound Healing"), ")"
  ),
  paste0(
    "IFN-γ Dominant: ", round(sum(ImmuneSubtypeClass$Cluster == "IFN-γ Dominant")/nrow(ImmuneSubtypeClass)*100, 2), "%", 
    "\n(", sum(ImmuneSubtypeClass$Cluster == "IFN-γ Dominant"), ")"
  ),
  paste0(
    "Inflammatory: ", round(sum(ImmuneSubtypeClass$Cluster == "Inflammatory")/nrow(ImmuneSubtypeClass)*100, 2), "%", 
    "\n(", sum(ImmuneSubtypeClass$Cluster == "Inflammatory"), ")"
  )
)

Figure1A <- pie3D(plot.data.1A_b, labels = labels, explode = 0.1, col = use_color, theta = 1.2,radius=0.9)

#####Figrue 1B Heatmap for Clinical Characters#####
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
    ImmuneSubtypeClass %>% as.data.frame(),
    by = "ID"
  ) %>%
  column_to_rownames(var = "ID") %>%
  mutate(across(everything(), ~ replace_na(as.character(.x), "Unknown")))
columnAnno <- HeatmapAnnotation(
  Cluster = plot.data.1B$Cluster,
  Status = plot.data.1B$Status,
  Stage = plot.data.1B$Stage,
  T_Stage = plot.data.1B$T_Stage,
  M_Stage = plot.data.1B$M_Stage,
  N_Stage = plot.data.1B$N_Stage,
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

n <- nrow(plot.data.1B)
mat.plot.data.1B <- matrix(0, nrow = 1, ncol = n)
rownames(mat.plot.data.1B) <- "gene"

if (!is.null(rownames(plot.data.1B))) {
  colnames(mat.plot.data.1B) <- rownames(plot.data.1B)
} else if ("ID" %in% colnames(plot.data.1B)) {
  colnames(mat.plot.data.1B) <- plot.data.1B$ID
} else {
  colnames(mat.plot.data.1B) <- paste0("S", seq_len(n))
}
split <- factor(plot.data.1B$Cluster, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory"))
ord   <- c(which(plot.data.1B$Cluster == "Wound Healing"),
           which(plot.data.1B$Cluster == "IFN-γ Dominant"),
           which(plot.data.1B$Cluster == "Inflammatory"))
Figure_1B <- ComplexHeatmap::Heatmap(
  mat.plot.data.1B,
  col = c("0" = "white"),
  name = NULL,
  na_col = "white",
  show_column_names = FALSE,
  show_row_names = FALSE,
  row_names_side = "left",
  column_split = split,
  column_order = ord,
  cluster_columns = FALSE,     
  cluster_column_slices = FALSE,
  gap = unit(2, "mm"),
  top_annotation = columnAnno,
  height = unit(0.1, "mm"),
  border = FALSE              
)



###Figure 1C Character Compare###
plot.data.1C_Status <- plot.data.1B %>%
  select(Status, Cluster) %>%
  group_by(Cluster, Status) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  mutate(Cluster = factor(Cluster,levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")))

tab <- table(plot.data.1B$Cluster, plot.data.1B$Status)
chisq.test(tab)  

Figure_1C_status <- ggplot(plot.data.1C_Status, aes(x = Cluster, y = Percentage, fill = Status)) +
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
  select(Stage, Cluster) %>%
  dplyr::filter(Stage != "Unknown") %>%
  group_by(Cluster, Stage) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  mutate(Cluster = factor(Cluster,levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")))
tab.dat <- plot.data.1B %>%
  select(Stage, Cluster) %>%
  dplyr::filter(Stage != "Unknown")
tab <- table(tab.dat$Cluster, tab.dat$Stage)
chisq.test(tab) 



Figure.1C_Stage  <- ggplot(plot.data.1C_Stage, aes(x = Cluster, y = Percentage, fill = Stage)) +
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
  select(T_Stage, Cluster) %>%
  dplyr::filter(T_Stage != "Unknown") %>%
  group_by(Cluster, T_Stage) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  mutate(Cluster = factor(Cluster,levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")))
tab.dat <- plot.data.1B %>%
  select(T_Stage, Cluster) %>%
  dplyr::filter(T_Stage != "Unknown")
tab <- table(plot.data.1B$Cluster, plot.data.1B$T_Stage)
chisq.test(tab)  


Figure.1C_T_Stage  <- ggplot(plot.data.1C_T_Stage, aes(x = Cluster, y = Percentage, fill = T_Stage)) +
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
  select(M_Stage, Cluster) %>%
  dplyr::filter(M_Stage != "Unknown") %>%
  group_by(Cluster, M_Stage) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  mutate(Cluster = factor(Cluster,levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")))
tab.dat <- plot.data.1B %>%
  select(M_Stage, Cluster) %>%
  dplyr::filter(M_Stage != "Unknown")
tab <- table(plot.data.1B$Cluster, plot.data.1B$M_Stage)
chisq.test(tab)  


Figure.1C_M_Stage  <- ggplot(plot.data.1C_M_Stage, aes(x = Cluster, y = Percentage, fill = M_Stage)) +
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
  select(N_Stage, Cluster) %>%
  dplyr::filter(N_Stage != "Unknown") %>%
  group_by(Cluster, N_Stage) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  mutate(Cluster = factor(Cluster,levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")))
tab.dat <- plot.data.1B %>%
  select(N_Stage, Cluster) %>%
  dplyr::filter(N_Stage != "Unknown") 
tab <- table(plot.data.1B$Cluster, plot.data.1B$N_Stage)
chisq.test(tab)  


Figure.1C_N_Stage  <- ggplot(plot.data.1C_N_Stage, aes(x = Cluster, y = Percentage, fill = N_Stage)) +
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
layout <- "
ABC
DE.
"
Figure_1C_combine <- (
  Figure.1C_Stage   +  # A
    Figure.1C_T_Stage +  # B
    Figure_1C_status  +  # C
    Figure.1C_M_Stage +  # D
    Figure.1C_N_Stage    # E
) + plot_layout(design = layout)

ggsave("./Output/Script_1/Figure_1/Figure_1C_combine.pdf",Figure_1C_combine,width = 12,height = 9)
#####Figure 1D TME Score#####
tmescore <- tmescore(eset   = LUAD.TPM,#expression data
                     method  ="PCA",#default
                     classify =F)#if true, survival data must be provided in pdata
tmescore.dat <- ImmuneSubtypeClass %>%
  left_join(tmescore,by ="ID")

Figure1D_TMEscore <- ggplot(tmescore.dat,aes(Cluster,TMEscore,fill= Cluster))+
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
  labs(x=NULL,y="TMEscore") +
  scale_fill_manual(values = use_color) +
  theme(axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1))
ggsave("./Output/Script_1/Figure_1/Figure1D_TMEscore.pdf",Figure1D_TMEscore,width = 3.5,height = 6)
#####Figure 1E Surv Diff#####
Surv_plot.data <- ImmuneSubtypeClass %>%
  left_join(LUAD.Survival, join_by("ID" == "sample")) %>%
  filter(OS.time > 30 & OS.time < 3650)

Surv_object <- Surv(time = Surv_plot.data$OS.time, event = Surv_plot.data$OS)

fit <- survfit(Surv_object ~ Cluster, data = Surv_plot.data)
surv.plot <- ggsurvplot(
  fit, 
  data = Surv_plot.data, 
  pval = TRUE,               
  risk.table = FALSE,         
  palette = use_color,  
  legend.title = "Immune Cluster",  
  legend.labs = c("Wound Healing","IFN-γ Dominant","Inflammatory"),  
  xlab = "Time (Days)",    
  ylab = "Overall Survival Probability", 
  ggtheme = theme_bw()  +
    theme(
      panel.border = element_rect(colour = "black", size = 1.5), 
      axis.text = element_text(size = 12, color = "black", face = "bold"), 
      axis.title = element_text(size = 12, color = "black", face = "bold"), 
      axis.ticks = element_line(size = 1, color = "black"), 
      legend.text = element_text(size = 12, color = "black", face = "bold"), 
      legend.title = element_text(size = 14, color = "black", face = "bold") 
    )
)
restest <- pairwise_survdiff(Surv(time = OS.time, event = OS) ~ Cluster,
                             data = Surv_plot.data)

restest[["p.value"]]
Surv.diff.p <- as.data.frame(restest[["p.value"]])
Surv.diff.p <- round(Surv.diff.p,3)
Surv.diff.p[is.na(Surv.diff.p)] <- '-'
Surv.diff.p <- rownames_to_column(Surv.diff.p,var = '  ')
pdf("./Output/Script_1/Figure_1/Figure1E_Surv.pdf",width = 5,height = 5)
surv.plot
dev.off()

#####Figure 1 FGH Esmatime#####
TME_estimate <- deconvo_tme(eset = LUAD.TPM, method ="estimate")
TME_estimate.dat <- ImmuneSubtypeClass %>%
  left_join(TME_estimate,by ="ID")

Figure1F_ImmuneScore <- ggplot(TME_estimate.dat,aes(Cluster,ImmuneScore_estimate,fill= Cluster))+
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
  labs(x=NULL,y="Immune Score") +
  scale_fill_manual(values = use_color) +
  theme(axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1))

Figure1G_StromalScore <- ggplot(TME_estimate.dat,aes(Cluster,StromalScore_estimate,fill= Cluster))+
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
  labs(x=NULL,y="Stromal Score") +
  scale_fill_manual(values = use_color) +
  theme(axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1))

Figure1H_TumorPurity <- ggplot(TME_estimate.dat,aes(Cluster,TumorPurity_estimate,fill= Cluster))+
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
  labs(x=NULL,y="TumorPurity Score") +
  scale_fill_manual(values = use_color) +
  theme(axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1))

Figure_FGH <- Figure1F_ImmuneScore + Figure1G_StromalScore + Figure1H_TumorPurity

ggsave("./Output/Script_1/Figure_1/Figure_FGH.pdf",Figure_FGH,width = 10.5,height = 6)


#save.image("./Output/Script_1/RData/Script1_image.Rdata")
#saveRDS(ImmuneSubtypeClass,"./Output/Script_1/ImmuneSubtypeClass.RDS")
