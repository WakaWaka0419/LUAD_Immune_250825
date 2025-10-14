#####01Loading Packeage#####
Sys.setenv(LANG = 'EN')
library(tidyverse)
library(clusterProfiler)
library(tibble)
library(ChAMPdata)
library(circlize)
library(ConsensusClusterPlus)
library(org.Hs.eg.db)
library(pheatmap)
library(dendsort)
library(ggplot2)
library(ggsci)
library(sva) 
library(msigdbr)
library(scales)
library(ComplexHeatmap)
library(survival)
library(survminer)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)
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
library(METAFlux)
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
LUAD.TPM <- dplyr::select(LUAD.TPM,ImmuneSubtypeClass$ID)
LUAD.TPM <- log(LUAD.TPM + 1)
ImmuneSubtypeClass <- readRDS("./Output/Script_1/RData/ImmuneSubtypeClass.RDS")



#####03Main Script#####
####A Tide Compare####
Tide_Input <- t(apply(LUAD.TPM,1,function(x)x - mean(x)))
#write.table(Tide_Input,"./Output/Script_5/Data/Tide_Input.txt",sep = "\t",row.names = T)
Tide_Output <- data.table::fread("./Output/Script_5/Data/Tide_Output.csv")  %>%
  left_join(ImmuneSubtypeClass[,c("ID","Cluster")],
            join_by("Patient" == "ID"))
table(Tide_Output$Cluster, Tide_Output$Responder)

plot.data.5A_Tide <- Tide_Output %>%
  group_by(Cluster, Responder) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  mutate(Cluster = factor(Cluster,levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")))

tbl <- table(Tide_Output$Cluster, Tide_Output$Responder)
pval <- chisq.test(tbl)$p.value


Figure_5A_Tide <- ggplot(plot.data.5A_Tide, aes(x = Cluster, y = Percentage, fill = Responder)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),  
            color = "black", size = 4) +
  scale_fill_manual(values = c("TRUE"= "#E64B35CC", "FALSE"= "#4DBBD5CC")) +
  labs(x = "", y = "Percentage (%)", fill = "Response") +
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

####B Tide Score Compare####
plot.list.Tide <- list()
Tide.Score <- colnames(Tide_Output)[c(4:9,11:15)]
Tide.Score <- intersect(Tide.Score, colnames(Tide_Output))
for(marker in Tide.Score){
  Figrue.Tide <- ggplot(Tide_Output,aes(Cluster,.data[[marker]],fill= Cluster))+
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
    labs(x=NULL,y= paste0(marker,"Score")) +
    scale_fill_manual(values = use_color) +
    theme(axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1))
  plot.list.Tide[[marker]] <- Figrue.Tide
}
wrap_plots(plot.list.Tide, ncol=4)

wrap_plots(Figure_5A_Tide | plot.list.Tide[[1]] |  plot.list.Tide[[2]] | plot.list.Tide[[4]])
#ggsave("./Output/Script_5/Figure/Figure_ABCD.pdf",
#       wrap_plots(Figure_5A_Tide | plot.list.Tide[[1]] |  plot.list.Tide[[2]] | plot.list.Tide[[4]]),
#       width = 12,height = 5)

#ggsave("./Output/Script_5/Figure/Figure_S2.pdf",
#       wrap_plots(plot.list.Tide[-c(1, 2, 4)], ncol = 4, nrow = 2),
#       width = 10.5,height = 10)
####Figure 5E Submap####
# 自定义函数用来产生submap需要的数据格式
generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct, 
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
  
  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

# 创建submap需要的数据格式
skcm.immunotherapy.logNC <- read.table("./Input/Analysis_data/skcm.immunotherapy.47samples.log2CountsNorm.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F) 
rownames(skcm.immunotherapy.logNC) <- toupper(rownames(skcm.immunotherapy.logNC)) 
skcm.immunotherapy.info <- read.table("./Input/Analysis_data/skcm.immunotherapy.47sampleInfo.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

skcm.immunotherapy.info <- skcm.immunotherapy.info[order(skcm.immunotherapy.info$label),]
skcm.immunotherapy.info$rank <- rep(c(1,2,3,4),times=as.character(table(skcm.immunotherapy.info$label))) #1: CTLA4_noR 2: CTLA4_R 3:PD1_noR 4:PD1_R

# 创建submap需要的数据格式
tmp <- LUAD.TPM
GENELIST <- intersect(rownames(tmp),rownames(skcm.immunotherapy.logNC)) 

sam_info <- skcm.immunotherapy.info
in_gct <- skcm.immunotherapy.logNC[GENELIST,rownames(skcm.immunotherapy.info)]

# 产生输出数据的文件名
gct_file <- "./Output/Script_5/Data/skcm.immunotherapy.for.SubMap.gct"
cls_file <- "./Output/Script_5/Data/skcm.immunotherapy.for.SubMap.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

# 提出亚型的样本，顺序排列
match(colnames(LUAD.TPM),ImmuneSubtypeClass$ID)
samples.Wound_Healing <- ImmuneSubtypeClass[which(ImmuneSubtypeClass$Cluster == "Wound Healing"),"ID"]
samples.IFN <- ImmuneSubtypeClass[which(ImmuneSubtypeClass$Cluster == "IFN-γ Dominant"),"ID"]
samples.Inflammatory <- ImmuneSubtypeClass[which(ImmuneSubtypeClass$Cluster == "Inflammatory"),"ID"]

sam_info <- data.frame("ImmClust"=c(samples.Wound_Healing,samples.IFN,samples.Inflammatory),row.names = c(samples.Wound_Healing,samples.IFN,samples.Inflammatory))
sam_info$rank <- rep(c(1,2,3),times=c(length(samples.Wound_Healing),length(samples.IFN),length(samples.Inflammatory))) 

# 产生输出数据的文件名
gct_file <- "./Output/Script_5/Data/Immune2.for.SubMap.gct"
cls_file <- "./Output/Script_5/Data/Immune2.for.SubMap.cls"

in_gct <- log2(tmp[GENELIST,rownames(sam_info)] + 1) # 产生和示例数据类似的形式，log2转化的标准化count值
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")
col_fun <- colorRamp2(c(min(0), median(0),max(1)),c("#FFA5AB", "white","#BDD5EA" ))
cherry    <- "#91D1C299"
lightgrey <- "#849AB499"

tmp <- matrix(
  c(
    # FDR 部分 (A1–A3)
    0.9686314, 0.9008135, 0.2217782, 1.0000000,
    0.9980020, 1.0000000, 1.0000000, 0.01198801,
    0.2237762, 0.8151848, 1.0000000, 1.0000000,
    # Bonferroni 部分 (A1–A3)
    1.0000000, 1.0000000, 0.4435564, 1.0000000,
    1.0000000, 1.0000000, 1.0000000, 0.01198801,
    0.6713287, 1.0000000, 1.0000000, 1.0000000
  ),
  nrow = 6, byrow = TRUE,
  dimnames = list(
    c("Wound Healing p","IFN-γ Dominant p","Inflammatory p",
      "Wound Healing b","IFN-γ Dominant b","Inflammatory b"),
    c("CTAL4-noR","CTLA4-R","PD1-noR","PD1-R")
  )
)

#pdf("./Output/Script_5/Figure/Figure5_E.pdf",width = 6,height = 6)
#pheatmap(tmp, cellwidth = 30, cellheight = 30,
#         cluster_rows = F,cluster_cols = F,
#         color = col_fun,
#         gaps_row = 3,
#         annotation_row = data.frame(pvalue=c("Nominal p value","Nominal p value","Nominal p value","Bonferroni corrected","Bonferroni corrected","Bonferroni corrected"),row.names = rownames(tmp)),
#         annotation_colors = list(pvalue=c("Nominal p value"=lightgrey,"Bonferroni corrected"=cherry)),
#         filename = "heatmap_submap.pdf")
#dev.off()

####Figure 5 Chort Validation####
####IMvigor210####
load("./Input/GEO/Curated_data/IMvigor210_Cruated.Rdata")
###Figure 5F Benefit Compare###
plot.data.5FG <- IMvigor210_Clinical %>%
  rownames_to_column(var = "ID")  %>%
  dplyr::filter(`Best Confirmed Overall Response` != "NE") %>%
  mutate(Response_group = case_when(
    `Best Confirmed Overall Response` %in% c("CR", "PR") ~ "Responder",
    `Best Confirmed Overall Response` %in% c("SD", "PD", "NE") ~ "Non-responder",
    TRUE ~ NA_character_
  )) %>%
  mutate(Cluster = case_when(
    `Immune phenotype` %in% "desert" ~ "Wound Healing",
    `Immune phenotype` %in% "excluded" ~ "IFN-γ Dominant",
    `Immune phenotype` %in% "inflamed" ~ "Inflammatory",
    TRUE ~ NA_character_
  ))   %>%
  dplyr::filter(!is.na(Cluster)) %>%
  mutate(Cluster = factor(Cluster, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory"))) %>%
  dplyr::mutate(os = os*30)

  
plot.data.5F <- plot.data.5FG %>%
  group_by(Cluster, Response_group) %>%
  summarise(Count = n(),
            .groups = "drop")  %>%
  group_by(Cluster) %>%
  mutate(Percentage = Count / sum(Count) * 100)

tbl <- table(plot.data.5FG$Cluster, plot.data.5FG$Response_group)
pval <- fisher.test(tbl)$p.value


Figure_5F_IMvigor <- ggplot(plot.data.5F, aes(x = Cluster, y = Percentage, fill = Response_group)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),  
            color = "black", size = 4) +
  scale_fill_manual(values = c("Responder"= "#E64B35CC", "Non-responder"= "#4DBBD5CC")) +
  labs(x = "", y = "Percentage (%)", fill = "Response") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1),
    axis.text.y  = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold",size = 13),
    legend.position = "none"
  ) +
  annotate("text", x = 2, y = 105, 
           label = paste0("italic(P) == ", signif(pval, 3)),
           parse = TRUE,  
           face = "bold",
           fontface = "italic",size = 5)
#####Figure 5G Surv Diff#####

Surv_object <- Surv(time = plot.data.5FG$os, event = plot.data.5FG$censOS)

fit <- survfit(Surv_object ~ Cluster, data = plot.data.5FG)
surv.plot <- ggsurvplot(
  fit, 
  data = plot.data.5FG, 
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
restest <- pairwise_survdiff(Surv(time = os, event = censOS) ~ Cluster,
                             data = plot.data.5FG)

restest[["p.value"]]
Surv.diff.p <- as.data.frame(restest[["p.value"]])
Surv.diff.p <- round(Surv.diff.p,3)
Surv.diff.p[is.na(Surv.diff.p)] <- '-'
Surv.diff.p <- rownames_to_column(Surv.diff.p,var = '  ')
pdf("./Output/Script_5/Figure/Figure5F_IMvigor.pdf",width = 3,height = 5)
Figure_5F_IMvigor
dev.off()
pdf("./Output/Script_5/Figure/Figure5G_IMvigor.pdf",width = 5,height = 5)
surv.plot$plot
dev.off()

####GSE135222####
load("./Input/GEO/Curated_data/GSE135222_Curated.Rdata")
ImmuneSubtypeClass_GSE135222 <- tme_classifier(eset = GSE135222_exp, scale = TRUE)
plot.data.5HI <-  GSE135222_Clinical %>%
  left_join(ImmuneSubtypeClass_GSE135222[,c("ID","TMEcluster")],
            join_by("Sample ID" == "ID")) %>%
  dplyr::filter(!is.na(`Clinical benefit`)) %>%
  mutate(Response_group = case_when(
    `Clinical benefit` %in% c("DCB") ~ "Responder",
    `Clinical benefit` %in% c("NDB") ~ "Non-responder",
    TRUE ~ NA_character_
  )) %>%
  mutate(Cluster = case_when(
    TMEcluster %in% "IE" ~ "Wound Healing",
    TMEcluster %in% "IS" ~ "IFN-γ Dominant",
    TMEcluster %in% "IA" ~ "Inflammatory",
    TRUE ~ NA_character_
  )) %>%
  mutate(Cluster = factor(Cluster, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory"))) %>%
  mutate(PFS.time = as.numeric(PFS.time) *30)
  
plot.data.5H <- plot.data.5HI %>%
  group_by(Cluster, Response_group) %>%
  summarise(Count = n(),
            .groups = "drop")  %>%
  group_by(Cluster) %>%
  mutate(Percentage = Count / sum(Count) * 100)

tbl <- table(plot.data.5HI$Cluster, plot.data.5HI$Response_group)
pval <- fisher.test(tbl)$p.value


Figure_5H_GSE135222 <- ggplot(plot.data.5H, aes(x = Cluster, y = Percentage, fill = Response_group)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),  
            color = "black", size = 4) +
  scale_fill_manual(values = c("Responder"= "#E64B35CC", "Non-responder"= "#4DBBD5CC")) +
  labs(x = "", y = "Percentage (%)", fill = "Response") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1),
    axis.text.y  = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold",size = 13),
    legend.position = "none"
  ) +
  annotate("text", x = 2, y = 105, 
           label = paste0("italic(P) == ", signif(pval, 3)),
           parse = TRUE,  
           face = "bold",
           fontface = "italic",size = 5)

pdf("./Output/Script_5/Figure/Figure5H_GSE135222.pdf",width = 3,height = 5)
Figure_5H_GSE135222
dev.off()
#####Figure 5G Surv Diff#####

Surv_object <- Surv(time = plot.data.5HI$PFS.time, event = as.numeric(plot.data.5HI$PFS))

fit <- survfit(Surv_object ~ Cluster, data = plot.data.5HI)
surv.plot <- ggsurvplot(
  fit, 
  data = plot.data.5HI, 
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



pdf("./Output/Script_5/Figure/Figure5I_GSE135222.pdf",width = 5,height = 5)
surv.plot$plot
dev.off()


###GSE207422###
load("./Input/GEO/Curated_data/GSE207422_Curated.Rdata")
ImmuneSubtypeClass_GSE207422 <- tme_classifier(eset = GSE207422_exp, scale = TRUE)
plot.data.5J <- GSE207422_Clinical %>%
  left_join(ImmuneSubtypeClass_GSE207422[, c("ID", "TMEcluster")],
            by = join_by("Sample" == "ID")) %>%
  dplyr::filter(!is.na(RECIST)) %>%
  mutate(Response_group = case_when(
    RECIST %in% c("CR","PR") ~ "Responder",
    RECIST %in% c("SD") ~ "Non-responder",
    TRUE ~ NA_character_
  )) %>%
  mutate(Cluster = case_when(
    TMEcluster == "IE" ~ "Wound Healing",
    TMEcluster == "IS" ~ "IFN-γ Dominant",
    TMEcluster == "IA" ~ "Inflammatory",
    TRUE ~ NA_character_
  )) %>%
  mutate(Cluster = factor(Cluster,
                          levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")))


plot.data.5J <- plot.data.5J %>%
  group_by(Cluster, Response_group) %>%
  summarise(Count = n(),
            .groups = "drop")  %>%
  group_by(Cluster) %>%
  mutate(Percentage = Count / sum(Count) * 100)

tbl <- table(plot.data.5J$Cluster, plot.data.5J$Response_group)
pval <- fisher.test(tbl)$p.value


Figure_J_GSE207422 <- ggplot(plot.data.5J, aes(x = Cluster, y = Percentage, fill = Response_group)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),  
            color = "black", size = 4) +
  scale_fill_manual(values = c("Responder"= "#E64B35CC", "Non-responder"= "#4DBBD5CC")) +
  labs(x = "", y = "Percentage (%)", fill = "Response") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1),
    axis.text.y  = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold",size = 13),
    legend.position = "none"
  ) +
  annotate("text", x = 2, y = 105, 
           label = paste0("italic(P) == ", signif(pval, 3)),
           parse = TRUE,  
           face = "bold",
           fontface = "italic",size = 5)

pdf("./Output/Script_5/Figure/Figure5J_GSE207422.pdf",width = 3,height = 5)
Figure_J_GSE207422
dev.off()



###GSE126044###
load("./Input/GEO/Curated_data/GSE126044_Curated.Rdata")
ImmuneSubtypeClass_GSE126044 <- tme_classifier(eset = GSE126044_exp, scale = TRUE)
plot.data.5K <- GSE126044_Clinical %>%
  left_join(ImmuneSubtypeClass_GSE126044[, c("ID", "TMEcluster")],
            by = join_by("ID" == "ID")) %>%
  mutate(Cluster = case_when(
    TMEcluster == "IE" ~ "Wound Healing",
    TMEcluster == "IS" ~ "IFN-γ Dominant",
    TMEcluster == "IA" ~ "Inflammatory",
    TRUE ~ NA_character_
  )) %>%
  mutate(Cluster = factor(Cluster,
                          levels = c("Wound Healing","IFN-γ Dominant","Inflammatory")))


plot.data.5J <- plot.data.5K %>%
  group_by(Cluster, Responsiveness) %>%
  summarise(Count = n(),
            .groups = "drop")  %>%
  group_by(Cluster) %>%
  mutate(Percentage = Count / sum(Count) * 100)

tbl <- table(plot.data.5K$Cluster, plot.data.5K$Responsiveness)
pval <- fisher.test(tbl)$p.value


Figure_K_GSE126044 <- ggplot(plot.data.5J, aes(x = Cluster, y = Percentage, fill = Responsiveness)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5),  
            color = "black", size = 4) +
  scale_fill_manual(values = c("Responder"= "#E64B35CC", "Non-responder"= "#4DBBD5CC")) +
  labs(x = "", y = "Percentage (%)", fill = "Response") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x  = element_text(face = "bold",angle = 45,hjust = 1),
    axis.text.y  = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold",size = 13),
    legend.position = "none"
  ) +
  annotate("text", x = 2, y = 105, 
           label = paste0("italic(P) == ", signif(pval, 3)),
           parse = TRUE,  
           face = "bold",
           fontface = "italic",size = 5)

pdf("./Output/Script_5/Figure/Figure_K_GSE126044.pdf",width = 3,height = 5)
Figure_K_GSE126044
dev.off()
