#####Loading Packages#####
library(ggplot2)
library(data.table)
library(survival)
library(ComplexHeatmap)
library(forestplot)
library(survminer)
library(circlize)
library(data.table)
library(GSVA)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggcor)
#####Main Script#####
####Pancancer TPX2####
rawAnno <- read.delim("~/Working_folder/DataBase/Pancanatlas/merged_sample_quality_annotations.tsv",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) # 数据来自http://api.gdc.cancer.gov/data/1a7d7be8-675d-4e60-a105-19d4121bdebf
rawAnno$simple_barcode <- substr(rawAnno$aliquot_barcode,1,15)
samAnno <- rawAnno[!duplicated(rawAnno$simple_barcode),c("cancer type", "simple_barcode")]
samAnno <- samAnno[which(samAnno$`cancer type` != ""),]
write.table(samAnno,"./Output/Script_7/DATA/output_simple_sample_annotation.txt",sep = "\t",row.names = F,col.names = T,quote = F)
surv <- read.delim("~/Working_folder/DataBase/Pancanatlas/Survival_SupplementalTable_S1_20171025_xena_sp.tsv", sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) # 数据来自https://xenabrowser.net/datapages/?dataset=Survival_SupplementalTable_S1_20171025_xena_sp&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
expr <- fread("~/Working_folder/DataBase/Pancanatlas/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",sep = "\t",stringsAsFactors = F,check.names = F,header = T) # 数据来自http://api.gdc.cancer.gov/data/3586c0da-64d0-4b74-a449-5ff4d9136611
expr <- as.data.frame(expr); rownames(expr) <- expr[,1]; expr <- expr[,-1]
gene <- sapply(strsplit(rownames(expr),"|",fixed = T), "[",1)
expr$gene <- gene
expr <- expr[!duplicated(expr$gene),]
rownames(expr) <- expr$gene; expr <- expr[,-ncol(expr)]
expr[expr < 0] <- 0 
colnames(expr) <- substr(colnames(expr),1,15)
gc()
geneOfInterest <- "TPX2"
if(!is.element(geneOfInterest, rownames(expr))) {
  warning("The gene ", geneOfInterest," cannot be found!")
} else {
  message("The gene ", geneOfInterest," can be matched!")
}
expr.sub <- expr[geneOfInterest, ] 
expr.sub <- as.data.frame(t(na.omit(t(expr.sub)))) 
keepSam <- colnames(expr.sub) 
expr <- expr[geneOfInterest,keepSam] 
rm(expr.sub); gc()
sam <- samAnno[which(samAnno$`cancer type` != "LAML"),"simple_barcode"] 
comsam <- intersect(intersect(colnames(expr), sam), rownames(surv))
tumsam <- comsam[substr(comsam,14,14) == "0"] 
tumAnno <- samAnno[which(samAnno$simple_barcode %in% tumsam),] 
tumAnno <- tumAnno[order(tumAnno$`cancer type`),] 
tumors <- unique(tumAnno$`cancer type`) 
exprSurv <- cbind.data.frame(expr = log2(as.numeric(expr[geneOfInterest,comsam]) + 1),
                             surv[comsam,c("OS","OS.time","DSS","DSS.time","DFI","DFI.time","PFI","PFI.time")])
outTab.cox <- NULL
for(i in tumors) {
  sam <- tumAnno[which(tumAnno$`cancer type` == i),"simple_barcode"]
  exprSurvSub <- exprSurv[sam,]
  
  ## OS
  coxres <- summary(coxph(Surv(OS.time, OS) ~ expr, data = exprSurvSub))
  outTab.cox <- rbind.data.frame(outTab.cox,
                                 data.frame(tumor = i, 
                                            event = "OS", 
                                            beta = coxres$coefficients[1,1], #
                                            hr = coxres$coefficients[1,2], 
                                            lower = coxres$conf.int[1,3], 
                                            upper = coxres$conf.int[1,4], 
                                            p = coxres$coefficients[1,5], 
                                            stringsAsFactors = F),
                                 stringsAsFactors = F)
  
  ## DSS
  coxres <- summary(coxph(Surv(DSS.time, DSS) ~ expr, data = exprSurvSub))
  outTab.cox <- rbind.data.frame(outTab.cox,
                                 data.frame(tumor = i,
                                            event = "DSS",
                                            beta = coxres$coefficients[1,1],
                                            hr = coxres$coefficients[1,2],
                                            lower = coxres$conf.int[1,3],
                                            upper = coxres$conf.int[1,4],
                                            p = coxres$coefficients[1,5],
                                            stringsAsFactors = F),
                                 stringsAsFactors = F)
  
  ## DFI
  if(i %in% c("SKCM","THYM","UVM","GBM")) { 
    outTab.cox <- rbind.data.frame(outTab.cox,
                                   data.frame(tumor = i,
                                              event = "DFI",
                                              beta = NA,
                                              hr = NA,
                                              lower = NA,
                                              upper = NA,
                                              p = NA,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
  } else {
    coxres <- summary(coxph(Surv(DFI.time, DFI) ~ expr, data = exprSurvSub))
    outTab.cox <- rbind.data.frame(outTab.cox,
                                   data.frame(tumor = i,
                                              event = "DFI",
                                              beta = coxres$coefficients[1,1],
                                              hr = coxres$coefficients[1,2],
                                              lower = coxres$conf.int[1,3],
                                              upper = coxres$conf.int[1,4],
                                              p = coxres$coefficients[1,5],
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
  }
  
  ## PFI
  coxres <- summary(coxph(Surv(PFI.time, PFI) ~ expr, data = exprSurvSub))
  outTab.cox <- rbind.data.frame(outTab.cox,
                                 data.frame(tumor = i,
                                            event = "PFI",
                                            beta = coxres$coefficients[1,1],
                                            hr = coxres$coefficients[1,2],
                                            lower = coxres$conf.int[1,4],
                                            upper = coxres$conf.int[1,4],
                                            p = coxres$coefficients[1,5],
                                            stringsAsFactors = F),
                                 stringsAsFactors = F)
}
#write.table(outTab.cox, file = "./Output/Script_7/DATA/output_summary of cox result in pancancer.txt",sep = "\r",row.names = F,col.names = T,quote = F)
minprop <- 0.2
outTab.km <- NULL
for(i in tumors) {
  sam <- tumAnno[which(tumAnno$`cancer type` == i),"simple_barcode"]
  exprSurvSub <- exprSurv[sam,]
  
  ## OS
  bestcut <- surv_cutpoint(exprSurvSub, 
                           time = "OS.time", 
                           event = "OS", 
                           variables = "expr", 
                           minprop = minprop) 
  cutoff <- bestcut$cutpoint[1,1]
  exprSurvSub$group <- factor(ifelse(exprSurvSub$expr > cutoff,"High","Low"), levels = c("Low","High"))
  fitd <- survdiff(Surv(OS.time, OS) ~ group, data=exprSurvSub, na.action=na.exclude)
  p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  HR <- (fitd$obs[2]/fitd$exp[2])/(fitd$obs[1]/fitd$exp[1])
  outTab.km <- rbind.data.frame(outTab.km,
                                data.frame(tumor = i, 
                                           event = "OS", 
                                           hr = HR, 
                                           lower = exp(log(HR) - qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1])), # 置信区间下限
                                           upper = exp(log(HR) + qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1])), # 置信区间上限
                                           p = p.val, 
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
  
  ## DSS
  bestcut <- surv_cutpoint(exprSurvSub, 
                           time = "DSS.time", 
                           event = "DSS", 
                           variables = "expr", 
                           minprop = minprop) 
  cutoff <- bestcut$cutpoint[1,1]
  exprSurvSub$group <- factor(ifelse(exprSurvSub$expr > cutoff,"High","Low"), levels = c("High","Low"))
  fitd <- survdiff(Surv(DSS.time, DSS) ~ group, data=exprSurvSub, na.action=na.exclude)
  p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  outTab.km <- rbind.data.frame(outTab.km,
                                data.frame(tumor = i, 
                                           event = "DSS", 
                                           hr = (fitd$obs[2]/fitd$exp[2])/(fitd$obs[1]/fitd$exp[1]), # Hazard ratio
                                           lower = exp(log(HR) - qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1])), # 置信区间下限
                                           upper = exp(log(HR) + qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1])), # 置信区间上限
                                           p = p.val, # p值
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
  
  ## DFI
  if(i %in% c("SKCM","THYM","UVM","GBM")) { 
    outTab.km <- rbind.data.frame(outTab.km,
                                  data.frame(tumor = i,
                                             event = "DFI",
                                             hr = NA,
                                             lower = NA,
                                             upper = NA,
                                             p = NA,
                                             stringsAsFactors = F),
                                  stringsAsFactors = F)
  } else {
    bestcut <- surv_cutpoint(exprSurvSub, 
                             time = "DFI.time", 
                             event = "DFI", 
                             variables = "expr", 
                             minprop = minprop) 
    cutoff <- bestcut$cutpoint[1,1]
    exprSurvSub$group <- factor(ifelse(exprSurvSub$expr > cutoff,"High","Low"), levels = c("High","Low"))
    fitd <- survdiff(Surv(DFI.time, DFI) ~ group, data=exprSurvSub, na.action=na.exclude)
    p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
    outTab.km <- rbind.data.frame(outTab.km,
                                  data.frame(tumor = i, 
                                             event = "DFI", 
                                             hr = (fitd$obs[2]/fitd$exp[2])/(fitd$obs[1]/fitd$exp[1]), # Hazard ratio
                                             lower = exp(log(HR) - qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1])), # 置信区间下限
                                             upper = exp(log(HR) + qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1])), # 置信区间上限
                                             p = p.val, # p值
                                             stringsAsFactors = F),
                                  stringsAsFactors = F)
  }
  
  ## PFI
  bestcut <- surv_cutpoint(exprSurvSub, 
                           time = "PFI.time", 
                           event = "PFI", 
                           variables = "expr", 
                           minprop = minprop) 
  cutoff <- bestcut$cutpoint[1,1]
  exprSurvSub$group <- factor(ifelse(exprSurvSub$expr > cutoff,"High","Low"), levels = c("High","Low"))
  fitd <- survdiff(Surv(PFI.time, PFI) ~ group, data=exprSurvSub, na.action=na.exclude)
  p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  outTab.km <- rbind.data.frame(outTab.km,
                                data.frame(tumor = i, 
                                           event = "PFI", 
                                           hr = (fitd$obs[2]/fitd$exp[2])/(fitd$obs[1]/fitd$exp[1]), # Hazard ratio
                                           lower = exp(log(HR) - qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1])), # 置信区间下限
                                           upper = exp(log(HR) + qnorm(0.975)*sqrt(1/fitd$exp[2]+1/fitd$exp[1])), # 置信区间上限
                                           p = p.val, 
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}
#write.table(outTab.km, file = "./Output/Script_7/DATA/output_summary of km result in pancancer.txt",sep = "\r",row.names = F,col.names = T,quote = F)

# 设置颜色
blue   <- "#A0CEE3"
yellow <- "#EFEFBE"
sea    <- "#37BCDF"
green  <- "#71BC5D"
cherry <- "#E5588C"
red    <- "#E46A6B"
purple <- "#8959A3"

# 制作热图绘制数据
hmInput <- NULL
for (i in tumors) {
  
  cox.res <- outTab.cox[which(outTab.cox$tumor == i),]
  km.res <- outTab.km[which(outTab.km$tumor == i),]
  
  cox.res$dirct <- ifelse(cox.res$hr > 1 & cox.res$p < 0.05, "Risky",
                          ifelse(cox.res$hr < 1 & cox.res$p < 0.05, "Protective","Nonsense"))
  km.res$dirct <- ifelse(km.res$hr > 1 & km.res$p < 0.05, "Risky",
                         ifelse(km.res$hr < 1 & km.res$p < 0.05, "Protective","Nonsense"))
  hmInput <- rbind.data.frame(hmInput,
                              data.frame(OS.cox = cox.res[1,"dirct"],
                                         OS.km = km.res[1,"dirct"],
                                         DSS.cox = cox.res[2,"dirct"],
                                         DSS.km = km.res[2,"dirct"],
                                         DFI.cox = cox.res[3,"dirct"],
                                         DFI.km = km.res[3,"dirct"],
                                         PFI.cox = cox.res[4,"dirct"],
                                         PFI.km = km.res[4,"dirct"],
                                         row.names = i,
                                         stringsAsFactors = F),
                              stringsAsFactors = F)
}
hmInput[is.na(hmInput)] <- "N/A"
write.table(hmInput, file = "./Output/Script_7/DATA/output_summary of gene prognositicationin pancancer.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
# 修改热图数据
indata <- hmInput
indata[indata == "Risky"] <- 0 #
indata[indata == "Protective"] <- 1 
indata[indata == "N/A"] <- 2 
indata[indata == "Nonsense"] <- 3 

# 构建列注释
annCol <- data.frame("Method" = rep(c("Cox","Log-rank"),4),
                     "Survival Type" = rep(c("OS","DSS","DFI","PFI"), each = 2),
                     check.names = F, 
                     row.names = colnames(indata))
annColors <- list("Method" = c("Cox" = blue, "Log-rank" = yellow),
                  "Survival Type" = c("OS" = sea, "DSS" = green, "DFI" = cherry, "PFI" = purple))
col_fun = circlize::colorRamp2(c(0, 1, 2, 3), c(red,green,"grey50","white")) 

# 绘制热图
hm <- Heatmap(matrix = as.matrix(indata),
              border = "black", 
              rect_gp = gpar(col = "black"), 
              name = "Prognostic role", 
              cluster_rows = F, 
              cluster_columns = F, 
              col = c(red,green,"grey50","white"), 
              show_row_names = T, 
              show_column_names = F, 
              row_names_side = "left", 
              top_annotation = HeatmapAnnotation(df = annCol,
                                                 col = annColors,
                                                 gp = gpar(col = "black"),
                                                 border = TRUE),
              width = grid::unit(8, "cm"),
              height = grid::unit(15, "cm"),
              heatmap_legend_param = list(at = c(0, 1, 2, 3), # 将图例的0123改成对应的文字
                                          legend_gp = grid::gpar(fill = col_fun(c(0,1,2,3))),
                                          labels = c("Risky", "Protective","N/A","Nonsense")))
pdf("./Output/Script_7/Figure/Figure_7Aprognostic heatmap.pdf", width = 10,height = 10)
draw(hm)
invisible(dev.off())
fpInput <- outTab.cox[which(outTab.cox$event == "OS"),]
hrtable <- fpInput[,c("tumor","event","beta","hr","lower","upper","p")]

tabletext <- cbind(c("Cancers",hrtable$tumor),
                   c("p value",ifelse(round(as.numeric(hrtable$p),3) < 0.001,"<0.001",round(as.numeric(hrtable$p),3))),
                   c("HR (95L-95H)",paste0(round(as.numeric(hrtable$hr),3), " (",
                                           round(as.numeric(hrtable$lower),3),"-",
                                           round(as.numeric(hrtable$upper),3),")")))

pdf("./Output/Script_7/Figure/Figure_7Aforestplot of os risk table in pancancer.pdf", width = 8, height = 8)
forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(hrtable$hr)),# HR
           lower=c(NA,as.numeric(hrtable$lower)), # 95%置信区间下限
           upper=c(NA,as.numeric(hrtable$upper)),# 95%置信区间上限
           graph.pos = 4,#图在表中的列位置
           graphwidth = unit(.3,"npc"),# 图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI",# box类型选择钻石
           col=fpColors(box="grey50", lines="grey50", zero = "black"),#box颜色
           boxsize=c(NA,ifelse(as.numeric(hrtable$p) < 0.05,0.8,0.4)),#box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=F,#不显示区
           zero=1,# zero线横坐标
           lwd.zero=2,# zero线宽
           xticks = c(0,1,2,4,8,13),# 横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab="Hazard ratio",
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),
                           "34" = gpar(lwd=2, col="black")),#最后一行底部加黑线???""中数字为nrow(tabletext) + 1
           txt_gp=fpTxtGp(label=gpar(cex=1),#各种字体大小设置
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.2)),
           lineheight = unit(.55,"cm"),#固定行高
           colgap = unit(0.4,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
invisible(dev.off())


####TPX2 CorLink####
gmt2list <- function(annofile){
  if (!file.exists(annofile)) {
    stop("There is no such gmt file.")
  }
  
  if (tools::file_ext(annofile) == "xz") {
    annofile <- xzfile(annofile)
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    close(annofile)
  } else if (tools::file_ext(annofile) == "gmt") {
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
  } else {
    stop ("Only gmt and gmt.xz are accepted for gmt2list")
  }
  
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  
  annoList <- lapply(y, `[`, c(-1,-2))
}
load("./Input/TCGA-LUAD/TCGA-LUAD_mrna_expr_tpm.rdata")
LUAD.TPM <- {
  x <- mrna_expr_tpm[, substr(colnames(mrna_expr_tpm),14,15) != "11", drop = FALSE]
  colnames(x) <- substr(colnames(x), 1, 16)
  x[, !duplicated(colnames(x)), drop = FALSE]
}
duplicated(colnames(LUAD.TPM))
TPX2 <- LUAD.TPM["TPX2",]
Cruated_immune <- data.table::fread("./Input/Analysis_data/Curated_Immune_Cell_Signature.txt")
Cruated_immune <- Cruated_immune %>%
  mutate(Symbol = toupper(str_replace_all(Symbol, "\\s+", ""))) %>%
  filter(!is.na(Symbol), nzchar(Symbol), !is.na(CellType), nzchar(CellType)) %>%
  distinct(CellType, Symbol)
celltype_sets <- Cruated_immune %>%
  group_by(CellType) %>%
  summarise(genes = list(unique(Symbol)), .groups = "drop")

immPath.list <- setNames(celltype_sets$genes, celltype_sets$CellType)
gset <- gmt2list("./Input/Analysis_data/h.all.v2025.1.Hs.symbols.gmt") 
gset.score <- gsva(expr = as.matrix(LUAD.TPM),
                   gset,
                   method = "ssgsea")

# 计算immunotherapy-predicted pathways的富集得分
immPath.score <- gsva(expr = as.matrix(LUAD.TPM),
                      immPath.list, 
                      method = "ssgsea")


immPath.score <- rbind.data.frame(immPath.score,
                                  TPX2)
gset.score <- rbind.data.frame(gset.score,
                                   TPX2)
immCorTPX2 <- NULL
for (i in rownames(immPath.score)) {
  cr <- cor.test(as.numeric(immPath.score[i,]),
                 as.numeric(TPX2),
                 method = "pearson")
  immCorTPX2 <- rbind.data.frame(immCorTPX2,
                                     data.frame(gene = "TPX2",
                                                path = i,
                                                r = cr$estimate,
                                                p = cr$p.value,
                                                stringsAsFactors = F),
                                     stringsAsFactors = F)
}
immCorTPX2$sign <- ifelse(immCorTPX2$r > 0,"pos","neg")
immCorTPX2$absR <- abs(immCorTPX2$r)
immCorTPX2$rSeg <- as.character(cut(immCorTPX2$absR,c(0,0.25,0.5,0.75,1),labels = c("0.25","0.50","0.75","1.00"),include.lowest = T))
immCorTPX2$pSeg <- as.character(cut(immCorTPX2$p,c(0,0.001,0.01,0.05,1),labels = c("<0.001","<0.01","<0.05","ns"),include.lowest = T))
immCorTPX2[nrow(immCorTPX2),"pSeg"] <- "Not Applicable"

immCorTPX2$rSeg <- factor(immCorTPX2$rSeg, levels = c("0.25","0.50","0.75","1.00"))
immCorTPX2$pSeg <- factor(immCorTPX2$pSeg, levels = c("<0.001","<0.01","<0.05","Not Applicable","ns"))
immCorTPX2$sign <- factor(immCorTPX2$sign, levels = c("pos","neg"))

p1 <- quickcor(t(immPath.score), 
               type = "lower",
               show.diag = TRUE) + 
  geom_colour() +
  anno_link(                          
    data    = immCorTPX2,
    mapping = aes(colour = pSeg, size = rSeg, linetype = sign),
    spec.key = "gene",                
    env.key  = "path",
    diag.label = FALSE
  ) +
  scale_size_manual(values = c(0.5, 1, 1.5, 2)) +
  scale_color_manual(values = c("#19A078","#DA6003","#7570B4","#E8288E","#65A818")) +
  scale_fill_gradient2(low = "#9483E1",mid = "white",high = "#E11953",midpoint=0) +
  remove_axis("x")
p1

# 循环计算相关性并绘制右上角
hallmarksCorTPX2 <- NULL
for (i in rownames(gset.score)) {
  cr <- cor.test(as.numeric(gset.score[i,]),
                 as.numeric(TPX2),
                 method = "pearson")
  hallmarksCorTPX2 <- rbind.data.frame(hallmarksCorTPX2,
                                      data.frame(gene = "TPX2",
                                                 path = i,
                                                 r = cr$estimate,
                                                 p = cr$p.value,
                                                 stringsAsFactors = F),
                                      stringsAsFactors = F)
}
hallmarksCorTPX2$sign <- ifelse(hallmarksCorTPX2$r > 0,"pos","neg")
hallmarksCorTPX2$absR <- abs(hallmarksCorTPX2$r)
hallmarksCorTPX2 <- hallmarksCorTPX2 %>%
  arrange(desc(absR)) %>%
  dplyr::slice(1:15)
hallmarksCorTPX2$rSeg <- as.character(cut(hallmarksCorTPX2$absR,c(0,0.25,0.5,0.75,1),labels = c("0.25","0.50","0.75","1.00"),include.lowest = T))
hallmarksCorTPX2$pSeg <- as.character(cut(hallmarksCorTPX2$p,c(0,0.001,0.01,0.05,1),labels = c("<0.001","<0.01","<0.05","ns"),include.lowest = T))
hallmarksCorTPX2[nrow(hallmarksCorTPX2),"pSeg"] <- "Not Applicable"

hallmarksCorTPX2$rSeg <- factor(hallmarksCorTPX2$rSeg, levels = c("0.25","0.50","0.75","1.00"))
hallmarksCorTPX2$pSeg <- factor(hallmarksCorTPX2$pSeg, levels = c("<0.001","<0.01","<0.05","Not Applicable","ns"))
hallmarksCorTPX2$sign <- factor(hallmarksCorTPX2$sign, levels = c("pos","neg"))
gset.score <- gset.score %>%
  as.data.frame(.) %>%
  dplyr::filter(rownames(.) %in% c(hallmarksCorTPX2$path))
  
p2 <- quickcor(t(gset.score), 
               type = "upper",
               show.diag = TRUE) + 
  geom_colour() +
  anno_link(data = hallmarksCorTPX2, 
            mapping = aes(colour = pSeg, size = rSeg, linetype = sign),
            spec.key = "gene",                
            env.key  = "path",
            diag.label = FALSE
  ) +
  scale_size_manual(values = c(0.5, 1, 1.5, 2)) +
  scale_color_manual(values = c("#19A078","#DA6003","#7570B4","#E8288E","#65A818")) +
  scale_fill_gradient2(low = "#9483E1",mid = "white",high = "#E11953",midpoint=0) +
  remove_axis("x")
#ggsave(p1,filename = "./Figure7B_A.pdf", width = 10,height = 8)
#ggsave(p2,filename = "Figure7B_B.pdf", width = 10,height = 8)
