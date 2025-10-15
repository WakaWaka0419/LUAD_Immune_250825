##### 01) Setup & Packages #####
Sys.setenv(LANG = "EN")
options(stringsAsFactors = FALSE, encoding = "UTF-8")
library(Seurat)
library(tidyverse)
library(ggrepel)
source("./R/SCENIC/makeMetaCells.R")
source("./R/SCENIC/compute_module_score.R")
source("./R/SCENIC/regulon_specificity.R")

##### 02) Constants & IO #####
INPUT_ANALYSIS_DIR      <- "input/analysis_data"
INPUT_CISTARGET_DB_DIR  <- file.path(INPUT_ANALYSIS_DIR, "cisTarget_db")

OUTPUT_BASE_DIR         <- "output/script6"
OUTPUT_DATA_DIR         <- file.path(OUTPUT_BASE_DIR, "data")
OUTPUT_FIG_DIR          <- file.path(OUTPUT_BASE_DIR, "figures")

dir.create(OUTPUT_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

SEURAT_OBJ_RDS          <- file.path(OUTPUT_DATA_DIR, "seu_SCENIC.rds")
MC_RDS                  <- file.path(OUTPUT_DATA_DIR, "00-1.mc.mat.rds")
MC_LOOM                 <- file.path(OUTPUT_DATA_DIR, "00-2.mc_mat_for_step1.loom")
TF_LIST_TXT             <- file.path(OUTPUT_DATA_DIR, "hsa_hgnc_tfs.motifs-v10.txt")
REGULONS_GMT            <- file.path(OUTPUT_DATA_DIR, "02-Scenic.regulons.gmt")
REGULONS_LIST_RDS       <- file.path(OUTPUT_DATA_DIR, "03-1.Scenic.regulons.rds")

# pySCENIC DB files
MOTIF_TBL               <- file.path(INPUT_CISTARGET_DB_DIR, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
RANKINGS_500BP_FEATHER  <- file.path(INPUT_CISTARGET_DB_DIR, "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
RANKINGS_10K_FEATHER    <- file.path(INPUT_CISTARGET_DB_DIR, "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")

##### 03) Load Seurat #####
seuCancer <- readRDS(SEURAT_OBJ_RDS)

##### 04) Build Meta-Cells #####
metacells <- makeMetaCells(
  seu       = seuCancer,
  min.cells = 10,
  reduction = "umap",
  dims      = 1:2,
  k.param   = 10,
  cores     = 10
)

# Persist meta-cell object for downstream steps
saveRDS(metacells, MC_RDS)

##### 05) Prepare pySCENIC Inputs #####
# 5.1 TF list (from motif-gene annotation)
motif2tfs <- data.table::fread(MOTIF_TBL)
tfs <- sort(unique(motif2tfs$gene_name))
writeLines(tfs, TF_LIST_TXT)

# 5.2 Meta-cell matrix -> filter -> genes in cisTarget ranking
mc <- readRDS(MC_RDS)
mc.meta <- mc$metadata
mc.mat  <- mc$mat

# low-expression gene filter
expr.in.cells <- rowSums(mc.mat > 0)
mc.mat <- mc.mat[expr.in.cells >= 5, , drop = FALSE]

# keep only genes present in the 10kb ranking DB
cisdb <- arrow::read_feather(RANKINGS_10K_FEATHER)
genes.use <- intersect(colnames(cisdb), rownames(mc.mat))
mc.mat <- mc.mat[genes.use, , drop = FALSE]

# 5.3 Write loom for pySCENIC step 1
loom <- SCopeLoomR::build_loom(
  file.name         = MC_LOOM,
  dgem              = mc.mat,
  default.embedding = NULL
)
loom$close(); rm(loom); gc()

# 5.4 Sanity checks
data_path <- normalizePath(OUTPUT_DATA_DIR, winslash = "/", mustWork = FALSE)
db_path   <- normalizePath(INPUT_CISTARGET_DB_DIR, winslash = "/", mustWork = FALSE)

file.exists(file.path(data_path, basename(MC_LOOM)))             # input loom
file.exists(file.path(data_path, basename(TF_LIST_TXT)))         # TF list
file.exists(file.path(db_path,   basename(MOTIF_TBL)))           # motif table
file.exists(file.path(db_path,   basename(RANKINGS_500BP_FEATHER)))
file.exists(file.path(db_path,   basename(RANKINGS_10K_FEATHER)))

##### 06) After pySCENIC (import regulons & compute RAS/AUCell) #####
# Assumes external pySCENIC pipeline produced REGULONS_GMT in OUTPUT_DATA_DIR
seuCancer <- readRDS(SEURAT_OBJ_RDS)

regulons_gmt <- clusterProfiler::read.gmt(REGULONS_GMT)
rg.names <- unique(regulons_gmt$term)
regulon.list <- lapply(rg.names, function(rg) regulons_gmt$gene[regulons_gmt$term == rg])
names(regulon.list) <- sub("[0-9]+g", "\\+", rg.names)

saveRDS(regulon.list, REGULONS_LIST_RDS)

# AUCell (regulon activity scores)
seuCancer <- ComputeModuleScore(
  seu        = seuCancer,
  gene.sets  = regulon.list,
  min.size   = 10,
  cores      = 10
)

DefaultAssay(seuCancer) <- "AUCell"
seuCancer <- RunUMAP(
  object        = seuCancer,
  features      = rownames(seuCancer),
  metric        = "correlation",
  reduction.name= "umapRAS",
  reduction.key = "umapRAS_"
)

# 6.1 RSS (regulon specificity score) against cell-type labels
rasMat <- t(seuCancer[["AUCell"]]@data)
ctMat  <- calIndMat(seuCancer$Cell_Cluster_level2)
rssMat <- calRSSMat(rasMat, ctMat)

# 6.2 Example plots (keep objects for saving by caller)
Figure_6L <- PlotRegulonRank(rssMat, "DPIS+ Proliferating Cancer")
DimPlot2(seuCancer, reduction = "UMAP", group.highlight = "DPIS+ Proliferating Cancer")
Figure_6M <- DimPlot2(seuCancer, reduction = "UMAP", regulon = "ZNF443(+)")

# Example save (uncomment if needed)
# pdf(file.path(OUTPUT_FIG_DIR, "figure6_lm.pdf"), width = 8, height = 4.5)
# Figure_6L + Figure_6M
# dev.off()
