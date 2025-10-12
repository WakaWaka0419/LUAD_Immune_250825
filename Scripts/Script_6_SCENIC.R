#####Loading Packages#####
library(Seurat)
library(tidyverse)
library(ggrepel)
source("./R/SCENIC/makeMetaCells.R")
source("R/SCENIC/compute_module_score.R")
source("R/SCENIC/regulon_specificity.R")
seu_Cancer_Cell_remake <- readRDS("Output/Script_6/DATA/seu_SCENIC.rds")
#####Main Scripts#####
####01 Meta Cell####
metacells <- 
  makeMetaCells(
    seu       = seu_Cancer_Cell_remake,
    min.cells = 10,
    reduction = "umap",
    dims      = 1:2,
    k.param   = 10,
    cores     = 10)

## 保存结果(meta cell matrix)
#saveRDS(metacells, "./Output/Script_6/DATA/00-1.mc.mat.rds")

####02 准备pySCENIC的输入文件####
###TF List###
motif2tfs <- data.table::fread("./Input/Analysis_data/cisTarget_db/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
TFs <- sort(unique(motif2tfs$gene_name))
writeLines(TFs, "./Output/Script_6/DATA/hsa_hgnc_tfs.motifs-v10.txt")
mc.mat <- readRDS("Output/Script_6/DATA/00-1.mc.mat.rds")
mc.meta <- mc.mat$metadata
mc.mat <- mc.mat$mat
###过滤低表达基因###
expr.in.cells <- rowSums(mc.mat > 0)
mc.mat <- mc.mat[expr.in.cells >= 5, ]
###过滤不在cisTargetDB中的基因###
cisdb <- arrow::read_feather("./Input/Analysis_data/cisTarget_db/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
genes.use <- intersect(colnames(cisdb), rownames(mc.mat))
length(genes.use)
dim(mc.mat)
mc.mat <- mc.mat[genes.use, ]
dim(mc.mat)

## 保存为loom文件
loom <- SCopeLoomR::build_loom(
  file.name         = "./Output/Script_6/DATA/00-2.mc_mat_for_step1.loom",
  dgem              = mc.mat,
  default.embedding = NULL
)
loom$close()

## 释放loom对象，允许其它文件占用loom文件
rm(loom)
gc()

###检查文件是否存在###
data_path <- file.path(getwd(), "Output", "Script_6", "DATA")
db_path   <- file.path(getwd(), "Input", "Analysis_data", "cisTarget_db")
# 输入文件
file.exists(file.path(data_path, "00-2.mc_mat_for_step1.loom"))
# TF 文件
file.exists(file.path(db_path, "hsa_hgnc_tfs.motifs-v10.txt"))
# Motif 注释文件
file.exists(file.path(db_path, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"))
# 数据库 feather 文件
file.exists(file.path(db_path, "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"))
file.exists(file.path(db_path, "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"))


####03 after pySCENIC件####
seu_Cancer_Cell_remake <- readRDS("./Output/Script_6/DATA/seu_SCENIC.rds")
regulons <- clusterProfiler::read.gmt("./Output/Script_6/DATA/02-Scenic.regulons.gmt")
rg.names <- unique(regulons$term)
regulon.list <- lapply(rg.names, function(rg) {
  subset(regulons, term == rg)$gene
})
names(regulon.list) <- sub("[0-9]+g", "\\+", rg.names)
summary(sapply(regulon.list, length))
print(regulon.list[1])
saveRDS(regulon.list, "Output/Script_6/DATA/03-1.Scenic.regulons.rds")

###用AUCell计算RAS matrix###
## RAS = regulon activity score ##
seu_Cancer_Cell_remake <- ComputeModuleScore(seu_Cancer_Cell_remake, gene.sets = regulon.list, min.size = 10, cores = 10)
DefaultAssay(seu_Cancer_Cell_remake) <- "AUCell"
seu_Cancer_Cell_remake <- RunUMAP(object = seu_Cancer_Cell_remake,
               features = rownames(seu_Cancer_Cell_remake),
               metric = "correlation", # 注意这里用correlation效果最好
               reduction.name = "umapRAS",
               reduction.key = "umapRAS_")
rasMat <- t(seu_Cancer_Cell_remake[["AUCell"]]@data)
dim(rasMat)
ctMat <- calIndMat(seu_Cancer_Cell_remake$Cell_Cluster_level2)
dim(ctMat)
rssMat <- calRSSMat(rasMat, ctMat)
dim(rssMat)
Figure_6L <- PlotRegulonRank(rssMat,"DPIS+ Proliferating Cancer")
DimPlot2(seu_Cancer_Cell_remake, reduction = "UMAP",group.highlight = "DPIS+ Proliferating Cancer")
Figure_6M <- DimPlot2(seu_Cancer_Cell_remake, reduction = "UMAP",regulon = "ZNF443(+)")
#pdf("./Output/Script_6/Figure/Figure6LM.pdf",width = 8,height=4.5)
#Figure_6L + Figure_6M
#dev.off()
