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
library(Mime1)
source("./R/standarize_fun.R")
source("./R/unicox.R")
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

#####Main Script#####
###Limma###
ImmuneSubtype <- ImmuneSubtypeClass %>%
  dplyr::filter(Cluster %in% c("IFN-γ Dominant","Inflammatory")) %>%
  dplyr::select(Cluster,ID) %>%
  mutate(Cluster = factor(Cluster,levels = c("IFN-γ Dominant","Inflammatory")))
LUAD.TPM <- dplyr::select(LUAD.TPM,ImmuneSubtype$ID)
match(ImmuneSubtype$ID,colnames(LUAD.TPM))
design <- model.matrix(~ 0 + ImmuneSubtype$Cluster)
colnames(design) <- c("IFNγDominant","Inflammatory")
rownames(design) <- colnames(LUAD.TPM)
#差异矩阵
contrast.matrix <- makeContrasts(Inflammatory - IFNγDominant,levels = design)
fit <- lmFit(LUAD.TPM,design)
fit <- contrasts.fit(fit,contrast.matrix)
fit <- eBayes(fit)
allDiff <- topTable(fit,number = Inf) %>%
  na.omit() %>%
  filter(abs(logFC) > 1.5 & adj.P.Val <0.01)
###单因素###
LUAD.Survival <- LUAD.Survival %>%
  dplyr::filter(sample %in% colnames(LUAD.TPM)) 
LUAD.TPM <- dplyr::select(LUAD.TPM,LUAD.Survival$sample)
unicox.dat <- LUAD.TPM %>% 
  dplyr::filter(rownames(.) %in% rownames(allDiff)) %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  inner_join(LUAD.Survival[,c("sample","OS","OS.time")],by="sample") %>% 
  dplyr::rename_all(funs(str_replace_all(.,"-","_")))
BaSurv<-Surv(time = unicox.dat$OS.time,event = unicox.dat$OS)
##提取变量名
names(unicox.dat)
varnames<-colnames(unicox.dat)[2:(ncol(unicox.dat)-2)]
varnames
##循环单因素
univar<-lapply(varnames,function(x){unicox(unicox.dat,x)})
univar<-do.call(rbind,univar)
univar<-na.omit(univar)
cox.res <- univar[univar$pvalue <= 0.05,"group"]


#####Before Machine Learning#####
#TCGA
TCGA_Train_Cohort <- unicox.dat %>%
  dplyr::select(sample, OS, OS.time, dplyr::everything())
#GSE13213
load("./Input/GEO/Curated_data/GSE13213_Curated.Rdata")
GSE13213_Validation_Cohort <- GSE13213_exp %>%
  dplyr::filter(rownames(.) %in% cox.res) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  left_join(
    GSE13213_Clinical %>%
      rownames_to_column(var = "sample") %>%
      dplyr::select(sample,OS,OS.time),by = "sample") %>%
  dplyr::select(sample, OS.time, OS, dplyr::everything())
      

#GSE50081
load("./Input/GEO/Curated_data/GSE50081_Curated.Rdata")
GSE50081_Validation_Cohort <- GSE50081_exp %>%
  dplyr::filter(rownames(.) %in% cox.res) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  left_join(
    GSE50081_Clinical %>%
      rownames_to_column(var = "sample") %>%
      dplyr::select("sample","OS","OS.time"),by = "sample") %>%
  dplyr::select(sample, OS.time, OS, dplyr::everything())

#GSE31210
load("./Input/GEO/Curated_data/GSE31210_Curated.Rdata")
GSE31210_Validation_Cohort <- GSE31210_exp %>%
  dplyr::filter(rownames(.) %in% cox.res) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  left_join(
    GSE31210_Clinical %>%
      rownames_to_column(var = "sample") %>%
      dplyr::select("sample","OS","OS.time"),by = "sample") %>%
  dplyr::select(sample, OS.time, OS, dplyr::everything())

#GSE42127
load("./Input/GEO/Curated_data/GSE42127_Curated.Rdata")
GSE42127_Validation_Cohort <- GSE42127_exp %>%
  dplyr::filter(rownames(.) %in% cox.res) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  left_join(
    GSE42127_Clinical %>%
      rownames_to_column(var = "sample") %>%
      dplyr::select("sample","OS","OS.time"),by = "sample") %>%
  dplyr::select(sample, OS.time, OS, dplyr::everything())

#GSE30219
load("./Input/GEO/Curated_data/GSE30219_Curated.Rdata")
GSE30219_Validation_Cohort <- GSE30219_exp %>%
  dplyr::filter(rownames(.) %in% cox.res) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  left_join(
    GSE30219_Clinical %>%
      rownames_to_column(var = "sample") %>%
      dplyr::select("sample","OS","OS.time"),by = "sample") %>%
  dplyr::select(sample, OS.time, OS, dplyr::everything())

####Data Combine####
colnames(TCGA_Train_Cohort)[1] <- "ID"
colnames(GSE13213_Validation_Cohort)[1] <- "ID"
colnames(GSE50081_Validation_Cohort)[1] <- "ID"
colnames(GSE30219_Validation_Cohort)[1] <- "ID"
colnames(GSE31210_Validation_Cohort)[1] <- "ID"
colnames(GSE42127_Validation_Cohort)[1] <- "ID"
common_cols <- Reduce(intersect, list(
  colnames(GSE13213_Validation_Cohort),
  colnames(GSE50081_Validation_Cohort),
  colnames(GSE30219_Validation_Cohort),
  colnames(GSE31210_Validation_Cohort),
  colnames(GSE42127_Validation_Cohort),
  colnames(TCGA_Train_Cohort)
))
TCGA_Train_Cohort <- dplyr::select(TCGA_Train_Cohort,common_cols)
GSE13213_Validation_Cohort <- dplyr::select(GSE13213_Validation_Cohort,common_cols)
GSE50081_Validation_Cohort <- dplyr::select(GSE50081_Validation_Cohort,common_cols)
GSE30219_Validation_Cohort <- dplyr::select(GSE30219_Validation_Cohort,common_cols)
GSE31210_Validation_Cohort <- dplyr::select(GSE31210_Validation_Cohort,common_cols)
GSE42127_Validation_Cohort <- dplyr::select(GSE42127_Validation_Cohort,common_cols)

Tarin_Validation_Datasets <- list(Dataset1 = TCGA_Train_Cohort,
                                  Dataset2 = GSE13213_Validation_Cohort,
                                  Dataset3 = GSE50081_Validation_Cohort,
                                  Dataset4 = GSE30219_Validation_Cohort,
                                  Dataset5 = GSE31210_Validation_Cohort,
                                  Dataset6 = GSE42127_Validation_Cohort)
Gene <- intersect(cox.res,common_cols)
#saveRDS(Gene,"./Output/Script_3/Rdata/Gene.RDS")
#saveRDS(cox.res,"./Output/Script_3/Rdata/cox_res.RDS")
###Machine Learning Prog Sig###
res <- ML.Dev.Prog.Sig(train_data = Tarin_Validation_Datasets$Dataset1,
                       list_train_vali_Data = Tarin_Validation_Datasets,
                       unicox.filter.for.candi = F,
                       unicox_p_cutoff = 0.05,
                       candidate_genes = Gene,
                       mode = 'all',nodesize =5,seed = 0917)
#saveRDS(res,"./Output/Script_3/Rdata/machineres.RDS")
###C_Index###
source("./R/cindex_dis_all.R")
#pdf("./Output/Script_3/Figure/Figure_3A.pdf",width = 7,height = 12)
cindex_dis_all(res,validate_set = names(Tarin_Validation_Datasets)[-1],order =names(Tarin_Validation_Datasets),width = 0.15,n_models = 100, pick_top = TRUE)
#dev.off()
cindex_dis_select(res,
                  model="StepCox[forward] + Enet[α=0.1]",
                  order= names(Tarin_Validation_Datasets))
###Survival Plot###
source("./R/rs_sur.R")
survplot <- vector("list",6) 
for (i in c(1:6)) {
  print(survplot[[i]]<-rs_sur(res, model_name = "StepCox[forward] + Enet[α=0.1]",dataset = names(Tarin_Validation_Datasets)[i],
                              #color=c("blue","green"),
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = T,
                              xlab="Day",pval.coord=c(1000,0.9)))}
#pdf("./Output/Script_3/Figure/Figure_3B.pdf",width = 13,height =9)
aplot::plot_list(gglist=survplot,ncol=3)
#dev.off()

#运行1h 超慢！！！！！！
#res.feature.all <- ML.Corefeature.Prog.Screen(InputMatrix = Tarin_Validation_Datasets$Dataset1,
#                                              candidate_genes = Gene,
#                                              mode = "all",nodesize =5,seed = 1314 )
source("./R/core_feature_rank.R")
pdf("./Output/Script_3/Figure/Figure_3K.pdf",width = 4,height =6)
core_feature_rank(res.feature.all, top=10)
dev.off()
core_feature_select(res.feature.all)
#saveRDS(res.feature.all,"./Output/Script_3/Rdata/res.feature.all.RDS")
#saveRDS(res,"./Output/Script_3/Rdata/res.RDS")
###AUC Plot###
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = Tarin_Validation_Datasets[["Dataset1"]],
                             inputmatrix.list = Tarin_Validation_Datasets,mode = 'all',AUC_time = 1,
                             auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = Tarin_Validation_Datasets[["Dataset1"]],
                             inputmatrix.list = Tarin_Validation_Datasets,mode = 'all',AUC_time = 3,
                             auc_cal_method="KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = Tarin_Validation_Datasets[["Dataset1"]],
                             inputmatrix.list = Tarin_Validation_Datasets,mode = 'all',AUC_time = 5,
                             auc_cal_method="KM")

auc_plot_1y <- auc_dis_all(all.auc.1y,
            dataset = names(Tarin_Validation_Datasets),
            validate_set=names(Tarin_Validation_Datasets)[-1],
            order= names(Tarin_Validation_Datasets),
            width = 0.35,
            year= 1)
auc_plot_3y <- auc_dis_all(all.auc.3y,
            dataset = names(Tarin_Validation_Datasets),
            validate_set=names(Tarin_Validation_Datasets)[-1],
            order= names(Tarin_Validation_Datasets),
            width = 0.35,
            year=3)
auc_plot_5y <- auc_dis_all(all.auc.5y,
            dataset = names(Tarin_Validation_Datasets),
            validate_set=names(Tarin_Validation_Datasets)[-1],
            order= names(Tarin_Validation_Datasets),
            width = 0.35,
            year=5)
source("./R/auc_dis_select.R")
#pdf("./Output/Script_3/Figure/Figure_3I.pdf",width = 7,height =5)
auc_dis_select(list(all.auc.1y,all.auc.3y,all.auc.5y),
               model_name = "StepCox[forward] + Enet[α=0.1]",
               dataset = names(Tarin_Validation_Datasets),
               order= c("Dataset6","Dataset5","Dataset4","Dataset3","Dataset2","Dataset1"),
               year=c(1,3,5))
#dev.off()
###unicox analysis###
unicox.rs.res <- cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res,optimal.model = "StepCox[forward] + Enet[α=0.1]",type ='categorical')
metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)
source("./R/meta_unicox_vis.R")
pdf("./Output/Script_3/Figure/Figure_3J.pdf",width = 9,height =4)
meta_unicox_vis(metamodel,
                dataset = names(Tarin_Validation_Datasets))
dev.off()
