#####01Loading Packeage#####
Sys.setenv(LANG = 'EN')
library(tidyverse)
library(clusterProfiler)
library(tibble)
library(ConsensusClusterPlus)
library(org.Hs.eg.db)
library(pheatmap)
library(dendsort)
library(ggplot2)
library(ggsci)
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
LUAD.TPM <- log(LUAD.TPM + 1)
ImmuneSubtypeClass <- readRDS("./Output/Script_1/RData/ImmuneSubtypeClass.RDS")
LUAD.TPM <- dplyr::select(LUAD.TPM,ImmuneSubtypeClass$ID)

#####03Main Script#####
###A Compare Cluster Function GO###
##DEG##
#IFN-γ Dominant vs Other
ImmuneSubtype_IFN <- ImmuneSubtypeClass %>%
dplyr::mutate(Cluster = ifelse(Cluster == "IFN-γ Dominant","IFN-γ Dominant","Others")) %>%
  dplyr::select(Cluster,ID) %>%
  mutate(Cluster = factor(Cluster,levels = c("IFN-γ Dominant","Others")))
LUAD.TPM <- dplyr::select(LUAD.TPM,ImmuneSubtype_IFN$ID)
match(ImmuneSubtype_IFN$ID,colnames(LUAD.TPM))
design <- model.matrix(~ 0 + ImmuneSubtype_IFN$Cluster)
colnames(design) <- c("IFNγDominant","Others")
rownames(design) <- colnames(LUAD.TPM)
#差异矩阵
contrast.matrix <- makeContrasts(IFNγDominant - Others,levels = design)
fit <- lmFit(LUAD.TPM,design)
fit <- contrasts.fit(fit,contrast.matrix)
fit <- eBayes(fit)
Diff_IFN <- topTable(fit,number = Inf) %>%
  na.omit() %>%
  filter(abs(logFC > 0.5) & adj.P.Val <0.01) %>%
  rownames_to_column(var = "SYMBOL") %>%
  dplyr::select("SYMBOL","logFC") %>%
  dplyr::mutate(group = "IFN")

#Wound Healing vs Other
ImmuneSubtype_Wound_Healing <- ImmuneSubtypeClass %>%
  dplyr::mutate(Cluster = ifelse(Cluster == "Wound Healing","Wound Healing","Others")) %>%
  dplyr::select(Cluster,ID) %>%
  mutate(Cluster = factor(Cluster,levels = c("Wound Healing","Others")))
LUAD.TPM <- dplyr::select(LUAD.TPM,ImmuneSubtype_Wound_Healing$ID)
match(ImmuneSubtype_Wound_Healing$ID,colnames(LUAD.TPM))
design <- model.matrix(~ 0 + ImmuneSubtype_Wound_Healing$Cluster)
colnames(design) <- c("WoundHealing","Others")
rownames(design) <- colnames(LUAD.TPM)
#差异矩阵
contrast.matrix <- makeContrasts(WoundHealing - Others,levels = design)
fit <- lmFit(LUAD.TPM,design)
fit <- contrasts.fit(fit,contrast.matrix)
fit <- eBayes(fit)
Diff_Wound_Healing <- topTable(fit,number = Inf) %>%
  na.omit() %>%
  filter(abs(logFC > 0.5) & adj.P.Val <0.01) %>%
  rownames_to_column(var = "SYMBOL") %>%
  dplyr::select("SYMBOL","logFC") %>%
  dplyr::mutate(group = "Wound_Healing")


#Inflammatory vs Other
ImmuneSubtype_Inflammatory <- ImmuneSubtypeClass %>%
  dplyr::mutate(Cluster = ifelse(Cluster == "Inflammatory","Inflammatory","Others")) %>%
  dplyr::select(Cluster,ID) %>%
  mutate(Cluster = factor(Cluster,levels = c("Inflammatory","Others")))
LUAD.TPM <- dplyr::select(LUAD.TPM,ImmuneSubtype_Inflammatory$ID)
match(ImmuneSubtype_Inflammatory$ID,colnames(LUAD.TPM))
design <- model.matrix(~ 0 + ImmuneSubtype_Inflammatory$Cluster)
colnames(design) <- c("Inflammatory","Others")
rownames(design) <- colnames(LUAD.TPM)
#差异矩阵
contrast.matrix <- makeContrasts(Inflammatory - Others,levels = design)
fit <- lmFit(LUAD.TPM,design)
fit <- contrasts.fit(fit,contrast.matrix)
fit <- eBayes(fit)
Diff_Inflammatory <- topTable(fit,number = Inf) %>%
  na.omit() %>%
  filter(abs(logFC > 0.5) & adj.P.Val <0.01) %>%
  rownames_to_column(var = "SYMBOL") %>%
  dplyr::select("SYMBOL","logFC") %>%
  dplyr::mutate(group = "Inflammatory")

#Combine DEG#
Combine_DEG <- Reduce(rbind,list(Diff_IFN,Diff_Inflammatory,Diff_Wound_Healing))
Combine_DEG <- Combine_DEG %>% 
  left_join(
    bitr(Combine_DEG$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db),
    by = "SYMBOL"
    )

#CompareCluster#
#enrichGO#
CompareCluster.res <- compareCluster(
  ENTREZID ~ group,
  data = Combine_DEG,
  fun = enrichGO,              # no quotes here
  OrgDb = org.Hs.eg.db,        # correct argument name
  ont = "BP",
  pvalueCutoff = 0.05
)

CompareCluster.sim <- simplify(
  CompareCluster.res,
  cutoff=0.7,
  by="p.adjust",
  select_fun=min)


#write.table(CompareCluster.sim@compareClusterResult,"./Output/Script_4/Data/CompareCluster.sim.txt",sep = "\t",row.names = F)

plot.data.4A <- CompareCluster.sim@compareClusterResult %>%
  dplyr::filter(
    (Cluster == "IFN" & ID %in% c(
      "GO:0034341", # response to type II interferon
      "GO:0140888", # interferon-mediated signaling pathway
      "GO:0019221", # cytokine-mediated signaling pathway
      "GO:0019882", # antigen processing and presentation
      "GO:0042267"  # NK cell mediated cytotoxicity
    )) |
      (Cluster == "Inflammatory" & ID %in% c(
        "GO:0006956", # complement activation
        "GO:0006959", # humoral immune response
        "GO:0030595", # leukocyte chemotaxis
        "GO:0050900", # leukocyte migration
        "GO:0006910"  # phagocytosis, recognition
      )) |
      (Cluster == "Wound_Healing" & ID %in% c(
        "GO:0032964", # collagen biosynthetic process
        "GO:0050673", # epithelial cell proliferation
        "GO:0061448", # connective tissue development
        "GO:0030111", # regulation of Wnt signaling pathway
        "GO:0030510"  # regulation of BMP signaling pathway
      ))
  )  %>%
  mutate(
    # 用 strsplit 分割成两部分，取分子和分母
    GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  ) %>%
  group_by(Cluster) %>%
  arrange(GeneRatio_num, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(
    Description = factor(Description, levels = c( "humoral immune response", 
                                                 "leukocyte migration", "leukocyte chemotaxis", "complement activation", 
                                                 "phagocytosis, recognition","cytokine-mediated signaling pathway", "response to type II interferon", 
                                                 "antigen processing and presentation", "natural killer cell mediated cytotoxicity", 
                                                 "interferon-mediated signaling pathway","epithelial cell proliferation", 
                                                 "connective tissue development", "regulation of Wnt signaling pathway", 
                                                 "regulation of BMP signaling pathway", "collagen biosynthetic process"))
  ) %>%
  mutate(Cluster = factor(Cluster,levels = c("Wound_Healing","IFN","Inflammatory"))) %>%
  mutate(logp = -log(p.adjust))
Figrue_4A <- ggplot(plot.data.4A, aes(Cluster, Description)) +  
  geom_point(aes(fill = logp, size = GeneRatio_num), shape = 21) +  
  theme_bw() +
  scale_fill_gradient(low = "#2dabb9", high = "#FF9F1C") +
  labs(x = "", y = "", fill = "-Log10(P.Adjust)",size = "Gene Ratio") + 
  annotate(geom = "rect",fill = "#2EC4B6",alpha = 0.3,xmax = 1.5,xmin = 0,ymax = 15.5,ymin = 0)+
  annotate(geom = "rect",fill = "#BDD5EA",alpha = 0.3,xmax = 2.5,xmin = 1.5,ymax = 15.5,ymin = 0)+
  annotate(geom = "rect",fill = "#FFA5AB",alpha = 0.3,xmax = Inf,xmin = 2.5,ymax = 15.5,ymin = 0)+
  theme(
    axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 10, angle = 45,colour = "black"), 
    axis.text.y = element_text(size = 10, hjust = 1,colour = "black"), 
    legend.position = "right",
    legend.text = element_text(size = 10),
    #legend.title = element_text(size = 12, face = "bold"),
    panel.border = element_rect(size = 1)
  )
#pdf("./Output/Script_4/Figure/Figure_4A.pdf", width = 6.5, height = 6)
Figrue_4A
#dev.off()
###B Compare Cluster Function GSVA###
#处理原则 ：每个Geneset 不应该出现重复的基因 
hallmarks.data <- read.gmt("./Input/Analysis_data/h.all.v2025.1.Hs.symbols.gmt")
KEGG.data <- read.gmt("./Input/Analysis_data/c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt")
GO_BP.data <- read.gmt("./Input/Analysis_data/c5.go.bp.v2025.1.Hs.symbols.gmt")
GO_CC.data <- read.gmt("./Input/Analysis_data/c5.go.cc.v2025.1.Hs.symbols.gmt")
GO_MF.data<- read.gmt("./Input/Analysis_data/c5.go.mf.v2025.1.Hs.symbols.gmt")
Reactome.data <- read.gmt("./Input/Analysis_data/c2.cp.reactome.v2025.1.Hs.symbols.gmt")
gsva.data <- Reduce(rbind,list(hallmarks.data,KEGG.data,GO_BP.data,GO_CC.data,GO_MF.data,Reactome.data))

h <- dplyr::select(gsva.data, term, gene) %>% #或entrez_gene
  as.data.frame %>% 
  split(., .$term) %>% 
  lapply(., function(x)(x$gene)) #或entrez_gene
gs <- lapply(h, unique)
head(gs)
#saveRDS(gs,"./Output/Script_4/Rdata/gs.RDS")
gsva.res <- gsva(as.matrix(LUAD.TPM),gs)
gsva.res <- as.matrix(gsva.res)
#saveRDS(gsva.data,"./Output/Script_4/Rdata/gsvadata.RDS")

#DE GSVA#
#IFN-γ Dominant vs Other
ImmuneSubtype_IFN <- ImmuneSubtypeClass %>%
  dplyr::mutate(Cluster = ifelse(Cluster == "IFN-γ Dominant","IFN-γ Dominant","Others")) %>%
  dplyr::select(Cluster,ID) %>%
  mutate(Cluster = factor(Cluster,levels = c("IFN-γ Dominant","Others")))
match(ImmuneSubtype_IFN$ID,colnames(LUAD.TPM))
design <- model.matrix(~ 0 + ImmuneSubtype_IFN$Cluster)
colnames(design) <- c("IFNγDominant","Others")
rownames(design) <- colnames(LUAD.TPM)
#差异矩阵
contrast.matrix <- makeContrasts(IFNγDominant - Others,levels = design)
fit <- lmFit(gsva.res,design)
fit <- contrasts.fit(fit,contrast.matrix)
fit <- eBayes(fit)
Diff_IFN <- topTable(fit,number = Inf) %>%
  na.omit() %>%
  rownames_to_column(var = "SYMBOL") %>%
  dplyr::select("SYMBOL","logFC") %>%
  dplyr::mutate(group = "IFN") %>%
  arrange(desc(logFC))

#Wound Healing vs Other
ImmuneSubtype_Wound_Healing <- ImmuneSubtypeClass %>%
  dplyr::mutate(Cluster = ifelse(Cluster == "Wound Healing","Wound Healing","Others")) %>%
  dplyr::select(Cluster,ID) %>%
  mutate(Cluster = factor(Cluster,levels = c("Wound Healing","Others")))
match(ImmuneSubtype_Wound_Healing$ID,colnames(LUAD.TPM))
design <- model.matrix(~ 0 + ImmuneSubtype_Wound_Healing$Cluster)
colnames(design) <- c("WoundHealing","Others")
rownames(design) <- colnames(LUAD.TPM)
#差异矩阵
contrast.matrix <- makeContrasts(WoundHealing - Others,levels = design)
fit <- lmFit(gsva.res,design)
fit <- contrasts.fit(fit,contrast.matrix)
fit <- eBayes(fit)
Diff_Wound_Healing <- topTable(fit,number = Inf) %>%
  na.omit() %>%
  #filter(abs(logFC > 0.5) & adj.P.Val <0.01) %>%
  rownames_to_column(var = "SYMBOL") %>%
  dplyr::select("SYMBOL","logFC") %>%
  dplyr::mutate(group = "Wound_Healing") %>%
  arrange(desc(logFC))


#Inflammatory vs Other
ImmuneSubtype_Inflammatory <- ImmuneSubtypeClass %>%
  dplyr::mutate(Cluster = ifelse(Cluster == "Inflammatory","Inflammatory","Others")) %>%
  dplyr::select(Cluster,ID) %>%
  mutate(Cluster = factor(Cluster,levels = c("Inflammatory","Others")))
LUAD.TPM <- dplyr::select(LUAD.TPM,ImmuneSubtype_Inflammatory$ID)
match(ImmuneSubtype_Inflammatory$ID,colnames(LUAD.TPM))
design <- model.matrix(~ 0 + ImmuneSubtype_Inflammatory$Cluster)
colnames(design) <- c("Inflammatory","Others")
rownames(design) <- colnames(LUAD.TPM)
#差异矩阵
contrast.matrix <- makeContrasts(Inflammatory - Others,levels = design)
fit <- lmFit(gsva.res,design)
fit <- contrasts.fit(fit,contrast.matrix)
fit <- eBayes(fit)
Diff_Inflammatory <- topTable(fit,number = Inf) %>%
  na.omit() %>%
  #filter(abs(logFC > 0.5) & adj.P.Val <0.01) %>%
  rownames_to_column(var = "SYMBOL") %>%
  dplyr::select("SYMBOL","logFC") %>%
  dplyr::mutate(group = "Inflammatory") %>%
  arrange(desc(logFC))


#Select Gene Set#
Cruated_Geneset <- c(
  dplyr::slice(Diff_Wound_Healing,1:5)$SYMBOL,
  dplyr::slice(Diff_IFN,1:5)$SYMBOL,
  dplyr::slice(Diff_Inflammatory,1:5)$SYMBOL
)
Cruated_Geneset <- factor(Cruated_Geneset,levels = Cruated_Geneset)
Cruated_Wound_Healing <- dplyr::filter(Diff_Wound_Healing,SYMBOL %in% Cruated_Geneset)
Cruated_IFN <- dplyr::filter(Diff_IFN,SYMBOL %in% Cruated_Geneset)
Cruated_Inflammatory <- dplyr::filter(Diff_Inflammatory,SYMBOL %in% Cruated_Geneset)

plot.data.4B <- Reduce(function(x, y) merge(x, y, by = "SYMBOL"),
                       list(
                         Cruated_Wound_Healing[c("SYMBOL","logFC")],
                         Cruated_IFN[c("SYMBOL","logFC")],
                         Cruated_Inflammatory[c("SYMBOL","logFC")]
                       )) %>%
  column_to_rownames(var = "SYMBOL") %>%
  { colnames(.) <- c("Wound Healing","IFN","Inflammatory"); . } %>%
  t() %>%
  scale() %>%
  t()
col_fun <- colorRamp2(c(min(-1.5), median(0),max(1.5)),c("#2dabb9", "white","#FF9F1C" ))
rowAnno <- rowAnnotation(
  Cluster = factor(
    rep(c("Wound Healing", "IFN", "Inflammatory"), each = 5),
    levels = c("Wound Healing", "IFN", "Inflammatory")
  ),
  col = list(
    Cluster = c(
      "Wound Healing" = "#2EC4B6",
      "IFN"           = "#BDD5EA",
      "Inflammatory"  = "#FFA5AB"
    )
  )
)

Figure_4B <- Heatmap(
  as.matrix(plot.data.4B[Cruated_Geneset, c("Wound Healing","IFN","Inflammatory")]),
  row_split = factor(
    c(rep("Wound Healing", 5),
      rep("IFN", 5),
      rep("Inflammatory", 5)),
    levels = c("Wound Healing", "IFN", "Inflammatory")
  ),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  left_annotation = rowAnno,
  name = "Z-Score",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10,fontface = "bold",angle = 45)
  
)
#pdf("./Output/Script_4/Figure/Figure_4B.pdf", width = 5.5, height = 6)
Figure_4B
#dev.off()

###C Compare Cluster Metabilism###
data("nutrient_lookup_files")
metabolism_scores <- calculate_reaction_score(LUAD.TPM)
flux<-compute_flux(mras=metabolism_scores ,medium=human_blood) 
#saveRDS(flux,"./Output/Script_4/Rdata/metabolism_flux.RDS")
cbrt <- function(x) {
  sign(x) * abs(x)^(1/3)
}

flux=cbrt(flux)
pathway<-unique(unlist(human_gem$SUBSYSTEM))
pathway_score <- list()
for (i in pathway){
  path=i
  activity_score<-c()
  for (d in 1:ncol(flux)){
    activity_score[d]<-mean(abs(flux[which(unlist(human_gem$SUBSYSTEM)==i),d]))
  } 
  pathway_score[[i]]<-activity_score
}
all_pathway_score<-as.data.frame(do.call(rbind,pathway_score))
colnames(all_pathway_score) <- colnames(flux)

#DE Pathway Activity#
#IFN-γ Dominant vs Other
ImmuneSubtype_IFN <- ImmuneSubtypeClass %>%
  dplyr::mutate(Cluster = ifelse(Cluster == "IFN-γ Dominant","IFN-γ Dominant","Others")) %>%
  dplyr::select(Cluster,ID) %>%
  mutate(Cluster = factor(Cluster,levels = c("IFN-γ Dominant","Others")))
match(ImmuneSubtype_IFN$ID,colnames(all_pathway_score))
design <- model.matrix(~ 0 + ImmuneSubtype_IFN$Cluster)
colnames(design) <- c("IFNγDominant","Others")
rownames(design) <- colnames(LUAD.TPM)
#差异矩阵
contrast.matrix <- makeContrasts(IFNγDominant - Others,levels = design)
fit <- lmFit(all_pathway_score,design)
fit <- contrasts.fit(fit,contrast.matrix)
fit <- eBayes(fit)
Diff_IFN <- topTable(fit,number = Inf) %>%
  na.omit() %>%
  #filter(P.Value <0.05 & logFC > 0) %>%
  rownames_to_column(var = "SYMBOL") %>%
  #dplyr::select("SYMBOL","logFC") %>%
  dplyr::mutate(group = "IFN") %>%
  arrange(desc(logFC))

#Wound Healing vs Other
ImmuneSubtype_Wound_Healing <- ImmuneSubtypeClass %>%
  dplyr::mutate(Cluster = ifelse(Cluster == "Wound Healing","Wound Healing","Others")) %>%
  dplyr::select(Cluster,ID) %>%
  mutate(Cluster = factor(Cluster,levels = c("Wound Healing","Others")))
match(ImmuneSubtype_Wound_Healing$ID,colnames(LUAD.TPM))
design <- model.matrix(~ 0 + ImmuneSubtype_Wound_Healing$Cluster)
colnames(design) <- c("WoundHealing","Others")
rownames(design) <- colnames(all_pathway_score)
#差异矩阵
contrast.matrix <- makeContrasts(WoundHealing - Others,levels = design)
fit <- lmFit(all_pathway_score,design)
fit <- contrasts.fit(fit,contrast.matrix)
fit <- eBayes(fit)
Diff_Wound_Healing <- topTable(fit,number = Inf) %>%
  na.omit() %>%
  #filter(P.Value <0.05 & logFC > 0) %>%
  rownames_to_column(var = "SYMBOL") %>%
  #dplyr::select("SYMBOL","logFC") %>%
  dplyr::mutate(group = "Wound_Healing") %>%
  arrange(desc(logFC))


#Inflammatory vs Other
ImmuneSubtype_Inflammatory <- ImmuneSubtypeClass %>%
  dplyr::mutate(Cluster = ifelse(Cluster == "Inflammatory","Inflammatory","Others")) %>%
  dplyr::select(Cluster,ID) %>%
  mutate(Cluster = factor(Cluster,levels = c("Inflammatory","Others")))
LUAD.TPM <- dplyr::select(LUAD.TPM,ImmuneSubtype_Inflammatory$ID)
match(ImmuneSubtype_Inflammatory$ID,colnames(LUAD.TPM))
design <- model.matrix(~ 0 + ImmuneSubtype_Inflammatory$Cluster)
colnames(design) <- c("Inflammatory","Others")
rownames(design) <- colnames(LUAD.TPM)
#差异矩阵
contrast.matrix <- makeContrasts(Inflammatory - Others,levels = design)
fit <- lmFit(all_pathway_score,design)
fit <- contrasts.fit(fit,contrast.matrix)
fit <- eBayes(fit)
Diff_Inflammatory <- topTable(fit,number = Inf) %>%
  na.omit() %>%
  #filter(P.Value <0.05 & logFC > 0) %>%
  rownames_to_column(var = "SYMBOL") %>%
  #dplyr::select("SYMBOL","logFC") %>%
  dplyr::mutate(group = "Inflammatory") %>%
  arrange(desc(logFC))
Diff_IFN %>%
  dplyr::filter(P.Value < 0.05 & logFC > 0) %>%
  dplyr::slice(1:5) %>%
  dplyr::pull(SYMBOL)


#Select Gene Set#
Cruated_metaPathway <- c(
  Diff_Wound_Healing %>%
    dplyr::filter(P.Value < 0.05 & logFC > 0) %>%
    dplyr::slice(1:5) %>%
    dplyr::pull(SYMBOL),
  Diff_IFN %>%
    dplyr::filter(P.Value < 0.05 & logFC > 0) %>%
    dplyr::slice(1:5) %>%
    dplyr::pull(SYMBOL),
  Diff_Inflammatory %>%
    dplyr::filter(P.Value < 0.05 & logFC > 0) %>%
    dplyr::slice(1:5) %>%
    dplyr::pull(SYMBOL)
)
Cruated_metaPathway <- factor(Cruated_metaPathway,levels = Cruated_metaPathway)
Cruated_Wound_Healing <- dplyr::filter(Diff_Wound_Healing,SYMBOL %in% Cruated_metaPathway)
Cruated_IFN <- dplyr::filter(Diff_IFN,SYMBOL %in% Cruated_metaPathway)
Cruated_Inflammatory <- dplyr::filter(Diff_Inflammatory,SYMBOL %in% Cruated_metaPathway)

plot.data.4C <- Reduce(function(x, y) merge(x, y, by = "SYMBOL"),
                       list(
                         Cruated_Wound_Healing[c("SYMBOL","logFC")],
                         Cruated_IFN[c("SYMBOL","logFC")],
                         Cruated_Inflammatory[c("SYMBOL","logFC")]
                       )) %>%
  column_to_rownames(var = "SYMBOL") %>%
  { colnames(.) <- c("Wound Healing","IFN","Inflammatory"); . } %>%
  t() %>%
  scale() %>%
  t()
col_fun <- colorRamp2(c(min(-1.5), median(0),max(1.5)),c("#2dabb9", "white","#FF9F1C" ))
rowAnno <- rowAnnotation(
  Cluster = factor(
    rep(c("Wound Healing", "IFN", "Inflammatory"), each = 5),
    levels = c("Wound Healing", "IFN", "Inflammatory")
  ),
  col = list(
    Cluster = c(
      "Wound Healing" = "#2EC4B6",
      "IFN"           = "#BDD5EA",
      "Inflammatory"  = "#FFA5AB"
    )
  )
)

Figure_4C <- Heatmap(
  as.matrix(plot.data.4C[Cruated_metaPathway, c("Wound Healing","IFN","Inflammatory")]),
  row_split = factor(
    c(rep("Wound Healing", 5),
      rep("IFN", 5),
      rep("Inflammatory", 5)),
    levels = c("Wound Healing", "IFN", "Inflammatory")
  ),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  left_annotation = rowAnno,
  name = "Z-Score",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10,fontface = "bold",angle = 45)
  
)
#pdf("./Output/Script_4/Figure/Figure_4C.pdf", width = 5.5, height = 6)
Figure_4C
#dev.off()

###D Compare Pathway Activity and Transcription factor###
#Pathway Activity#
net <- decoupleR::get_progeny(organism = 'human', 
                              top = 500)
#saveRDS(net,"./Input/Analysis_data/Pathway_net.RDS")
#saveRDS(net,"./Input/Analysis_data/TF_net.RDS")
Pathway_decouple <- decoupleR::run_mlm(mat = LUAD.TPM, 
                                  net = net, 
                                  .source = 'source', 
                                  .target = 'target',
                                  .mor = 'weight', 
                                  minsize = 5)
Pathway_decouple <- Pathway_decouple %>%
  dplyr::filter(p_value < 0.05) %>%
  left_join(
    ImmuneSubtypeClass %>%
      dplyr::select(ID,Cluster) ,
    join_by("condition" == "ID")
  )
table(Pathway_decouple$source)



Pathway_decouple$Cluster <- factor(
  Pathway_decouple$Cluster,
  levels = c("Wound Healing", "IFN-γ Dominant", "Inflammatory")
)

my_comparisons <- list(
  c("Wound Healing","IFN-γ Dominant"),
  c("Wound Healing","Inflammatory"),
  c("IFN-γ Dominant","Inflammatory")
)


comparisons_list <- list(c("Wound Healing","IFN-γ Dominant"),
                         c("IFN-γ Dominant","Inflammatory"),
                         c("Wound Healing","Inflammatory"))

my_label <- function(p) {
  if (is.na(p)) return("NA")
  if (p <= 0.001) return("***")
  if (p <= 0.01)  return("**")
  if (p <= 0.05)  return("*")
  return("NA")  # 不显著时标 NA
}

Figure_4D <- ggplot(Pathway_decouple, aes(Cluster, score, fill = Cluster)) +
  geom_boxplot(outlier.shape = 21, color = "black") +
  scale_fill_manual(values = use_color) +
  facet_wrap(vars(source), nrow = 1) +
  stat_compare_means(
    method = "t.test",
    comparisons = comparisons_list,
    label = "p.signif",
    step.increase = 0.05,
    label.y = 1.1 * max(Pathway_decouple$score, na.rm = TRUE)  # 适当上移
  ) +
  labs(x = NULL, y = "Pathway activity") +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 15, face = "bold"),
    panel.border = element_rect(size = 1)
  )


pdf("./Output/Script_4/Figure/Figure_4D.pdf",width = 12,height = 4)
Figure_4D
dev.off()

#Compare Transcription factor#
net <- decoupleR::get_collectri(organism = 'human', 
                                split_complexes = FALSE)

TF_decouple <- decoupleR::run_ulm(mat = LUAD.TPM, 
                                  net = net, 
                                  .source = 'source', 
                                  .target = 'target',
                                  .mor = 'mor', 
                                  minsize = 5)

n_tfs <- 500

# Transform to wide matrix
TF_acts_mat <- TF_decouple %>%
  tidyr::pivot_wider(id_cols = 'condition', 
                     names_from = 'source',
                     values_from = 'score') %>%
  tibble::column_to_rownames('condition') %>%
  as.matrix()

# Get top tfs with more variable means across clusters
TFs <- TF_decouple %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(std = stats::sd(score)) %>%
  dplyr::arrange(-abs(std)) %>%
  head(n_tfs) %>%
  dplyr::pull(source)
TF_acts_mat <- as.data.frame(t(TF_acts_mat[,TFs]))

#DE TFS#
#IFN-γ Dominant vs Other
ImmuneSubtype_IFN <- ImmuneSubtypeClass %>%
  dplyr::mutate(Cluster = ifelse(Cluster == "IFN-γ Dominant","IFN-γ Dominant","Others")) %>%
  dplyr::select(Cluster,ID) %>%
  mutate(Cluster = factor(Cluster,levels = c("IFN-γ Dominant","Others")))
TF_acts_mat <- dplyr::select(TF_acts_mat,ImmuneSubtype_IFN$ID)
match(ImmuneSubtype_IFN$ID,colnames(TF_acts_mat))
design <- model.matrix(~ 0 + ImmuneSubtype_IFN$Cluster)
colnames(design) <- c("IFNγDominant","Others")
rownames(design) <- colnames(TF_acts_mat)
#差异矩阵
contrast.matrix <- makeContrasts(IFNγDominant - Others,levels = design)
fit <- lmFit(TF_acts_mat,design)
fit <- contrasts.fit(fit,contrast.matrix)
fit <- eBayes(fit)
Diff_IFN <- topTable(fit,number = Inf) %>%
  na.omit() %>%
  rownames_to_column(var = "SYMBOL") %>%
  dplyr::select("SYMBOL","logFC") %>%
  dplyr::mutate(group = "IFN") %>%
  arrange(desc(logFC)) 

#Wound Healing vs Other
ImmuneSubtype_Wound_Healing <- ImmuneSubtypeClass %>%
  dplyr::mutate(Cluster = ifelse(Cluster == "Wound Healing","Wound Healing","Others")) %>%
  dplyr::select(Cluster,ID) %>%
  mutate(Cluster = factor(Cluster,levels = c("Wound Healing","Others")))
TF_acts_mat <- dplyr::select(TF_acts_mat,ImmuneSubtype_Wound_Healing$ID)
match(ImmuneSubtype_Wound_Healing$ID,colnames(TF_acts_mat))
design <- model.matrix(~ 0 + ImmuneSubtype_Wound_Healing$Cluster)
colnames(design) <- c("WoundHealing","Others")
rownames(design) <- colnames(TF_acts_mat)
#差异矩阵
contrast.matrix <- makeContrasts(WoundHealing - Others,levels = design)
fit <- lmFit(TF_acts_mat,design)
fit <- contrasts.fit(fit,contrast.matrix)
fit <- eBayes(fit)
Diff_Wound_Healing <- topTable(fit,number = Inf) %>%
  na.omit() %>%
  rownames_to_column(var = "SYMBOL") %>%
  dplyr::select("SYMBOL","logFC") %>%
  dplyr::mutate(group = "Wound_Healing") %>%
  arrange(desc(logFC)) 


#Inflammatory vs Other
ImmuneSubtype_Inflammatory <- ImmuneSubtypeClass %>%
  dplyr::mutate(Cluster = ifelse(Cluster == "Inflammatory","Inflammatory","Others")) %>%
  dplyr::select(Cluster,ID) %>%
  mutate(Cluster = factor(Cluster,levels = c("Inflammatory","Others")))
TF_acts_mat <- dplyr::select(TF_acts_mat,ImmuneSubtype_Inflammatory$ID)
match(ImmuneSubtype_Inflammatory$ID,colnames(TF_acts_mat))
design <- model.matrix(~ 0 + ImmuneSubtype_Inflammatory$Cluster)
colnames(design) <- c("Inflammatory","Others")
rownames(design) <- colnames(TF_acts_mat)
#差异矩阵
contrast.matrix <- makeContrasts(Inflammatory - Others,levels = design)
fit <- lmFit(TF_acts_mat,design)
fit <- contrasts.fit(fit,contrast.matrix)
fit <- eBayes(fit)
Diff_Inflammatory <- topTable(fit,number = Inf) %>%
  na.omit() %>%
  rownames_to_column(var = "SYMBOL") %>%
  dplyr::select("SYMBOL","logFC") %>%
  dplyr::mutate(group = "Inflammatory") %>%
  arrange(desc(logFC))  

#Combine DEG#
Cruated_DE_TF <- Reduce(c,
                        list(
                          Diff_Wound_Healing %>% dplyr::slice(1:5) %>% pull(SYMBOL),
                          Diff_IFN %>% dplyr::slice(1:5) %>% pull(SYMBOL),
                          Diff_Inflammatory %>% dplyr::slice(1:5) %>% pull(SYMBOL)
                        ))

plot.data.4E <- Reduce(function(x, y) merge(x, y, by = "SYMBOL"),
                       list(
                         Diff_Wound_Healing[c("SYMBOL","logFC")],
                         Diff_IFN[c("SYMBOL","logFC")],
                         Diff_Inflammatory[c("SYMBOL","logFC")]
                       ))  %>% 
  column_to_rownames(var = "SYMBOL") %>%
  { colnames(.) <- c("Wound Healing","IFN","Inflammatory"); . } %>%
  t() %>%
  scale() %>%
  t() %>%
  as.data.frame() %>%
  dplyr::filter(rownames(.) %in% Cruated_DE_TF) 
col_fun <- colorRamp2(c(min(-1.5), median(0),max(1.5)),c("#2dabb9", "white","#FF9F1C" ))
rowAnno <- rowAnnotation(
  Cluster = factor(
    rep(c("Wound Healing", "IFN", "Inflammatory"), each = 5),
    levels = c("Wound Healing", "IFN", "Inflammatory")
  ),
  col = list(
    Cluster = c(
      "Wound Healing" = "#2EC4B6",
      "IFN"           = "#BDD5EA",
      "Inflammatory"  = "#FFA5AB"
    )
  )
)

Figure_4E <- Heatmap(
  as.matrix(plot.data.4E[Cruated_DE_TF, c("Wound Healing","IFN","Inflammatory")]),
  row_split = factor(
    c(rep("Wound Healing", 5),
      rep("IFN", 5),
      rep("Inflammatory", 5)),
    levels = c("Wound Healing", "IFN", "Inflammatory")
  ),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  left_annotation = rowAnno,
  name = "Z-Score",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10,fontface = "bold",angle = 45)
  
)
pdf("./Output/Script_4/Figure/Figure_4E.pdf", width = 4, height = 7)
Figure_4E
dev.off()

