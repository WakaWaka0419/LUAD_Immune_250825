##### 01) Setup & Packages #####
Sys.setenv(LANG = "EN")
options(stringsAsFactors = FALSE, encoding = "UTF-8")
library(tidyverse)
library(data.table)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(limma)
library(survival)
library(survminer)
library(GSVA)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(gghalves)
library(ggpp)
library(patchwork)
library(plotrix)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(factoextra)
library(ggConvexHull)
library(IOBR)
library(TMEscore)
library(ImmuneSubtypeClassifier)
library(METAFlux)

##### 02) Constants (paths, palette) #####
INPUT_LUAD_DIR      <- "input/tcga_luad"
INPUT_ANALYSIS_DIR  <- "input/analysis_data"
OUTPUT_FIG_DIR      <- "output/figures/fig4"
OUTPUT_RDS_DIR      <- "output/rds"
PALETTE_CLUSTER3    <- c("#2EC4B6","#BDD5EA","#FFA5AB")  # WoundHealing, IFNGDominant, Inflammatory
IMMUNE_SUBTYPE_MIN_SCORE <- 0.6
dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_RDS_DIR, recursive = TRUE, showWarnings = FALSE)
clusterDisplayLabels <- c(
  WoundHealing = "Wound Healing",
  IFNGDominant = "IFN-γ Dominant",
  Inflammatory = "Inflammatory"
)

##### 03) Load & Prepare Data #####
# 表达矩阵（TPM）
load(file.path(INPUT_LUAD_DIR, "TCGA-LUAD_mrna_expr_tpm.rdata"))  # -> mrna_expr_tpm

luadTpm <- {
  x <- mrna_expr_tpm[, substr(colnames(mrna_expr_tpm), 14, 15) != "11", drop = FALSE]  # 去除正常样本
  colnames(x) <- substr(colnames(x), 1, 16)
  x <- x[, !duplicated(colnames(x)), drop = FALSE]
  log(x + 1)
}

luadSurvival <- data.table::fread(file.path(INPUT_LUAD_DIR, "TCGA-LUAD.survival.tsv"))
luadClinical <- data.table::fread(file.path(INPUT_LUAD_DIR, "TCGA-LUAD.clinical.tsv")) |>
  dplyr::filter(tissue_type.samples == "Tumor")

immuneSubtypeAssignments <- readRDS(file.path(OUTPUT_RDS_DIR, "immune_subtypes.rds")) |>
  dplyr::mutate(
    cluster = factor(as.character(cluster), levels = c("WoundHealing","IFNGDominant","Inflammatory")),
    clusterLabel = clusterDisplayLabels[as.character(cluster)]
  )

subsetExprByIds <- function(expr, ids) {
  keep <- intersect(colnames(expr), ids)
  expr[, keep, drop = FALSE]
}

##### 04) A. Compare Cluster Function — GO（DEG -> compareCluster） #####

# 通用 LIMMA 对比函数：groupName vs Others（列名避免特殊字符）
runLimmaVsOthers <- function(expr, idsGroup, idsAll, coefName) {
  m <- subsetExprByIds(expr, idsAll)
  grp <- ifelse(colnames(m) %in% idsGroup, coefName, "Others")
  design <- model.matrix(~ 0 + factor(grp, levels = c(coefName, "Others")))
  colnames(design) <- c(coefName, "Others")
  fit <- lmFit(m, design)
  fit <- contrasts.fit(fit, makeContrasts(contrasts = paste0(coefName, "-Others"), levels = design))
  eBayes(fit)
}

# 1) IFNGDominant vs Others —— 基因层面
idsIFNG <- immuneSubtypeAssignments |> dplyr::filter(cluster == "IFNGDominant") |> dplyr::pull(ID)
idsALL  <- immuneSubtypeAssignments$ID
fitIFNG_gene <- runLimmaVsOthers(luadTpm, idsIFNG, idsALL, "IFNGDominant")
degIFNG_gene <- topTable(fitIFNG_gene, number = Inf) |>
  dplyr::filter(abs(logFC) > 0.5 & adj.P.Val < 0.01) |>
  tibble::rownames_to_column("SYMBOL") |>
  dplyr::select(SYMBOL, logFC) |>
  dplyr::mutate(group = "IFNG")

# 2) WoundHealing vs Others —— 基因层面
idsWH <- immuneSubtypeAssignments |> dplyr::filter(cluster == "WoundHealing") |> dplyr::pull(ID)
fitWH_gene <- runLimmaVsOthers(luadTpm, idsWH, idsALL, "WoundHealing")
degWH_gene <- topTable(fitWH_gene, number = Inf) |>
  dplyr::filter(abs(logFC) > 0.5 & adj.P.Val < 0.01) |>
  tibble::rownames_to_column("SYMBOL") |>
  dplyr::select(SYMBOL, logFC) |>
  dplyr::mutate(group = "WoundHealing")

# 3) Inflammatory vs Others —— 基因层面
idsINF <- immuneSubtypeAssignments |> dplyr::filter(cluster == "Inflammatory") |> dplyr::pull(ID)
fitINF_gene <- runLimmaVsOthers(luadTpm, idsINF, idsALL, "Inflammatory")
degINF_gene <- topTable(fitINF_gene, number = Inf) |>
  dplyr::filter(abs(logFC) > 0.5 & adj.P.Val < 0.01) |>
  tibble::rownames_to_column("SYMBOL") |>
  dplyr::select(SYMBOL, logFC) |>
  dplyr::mutate(group = "Inflammatory")

degCombined_gene <- dplyr::bind_rows(degIFNG_gene, degINF_gene, degWH_gene) |>
  dplyr::left_join(
    bitr(unique(.$SYMBOL), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    by = "SYMBOL"
  ) |>
  dplyr::filter(!is.na(ENTREZID))

# compareCluster（GO:BP）
cmpGO <- compareCluster(
  ENTREZID ~ group,
  data = degCombined_gene,
  fun  = enrichGO,
  OrgDb = org.Hs.eg.db,
  ont  = "BP",
  pvalueCutoff = 0.05
)

cmpGO_simplified <- simplify(cmpGO, cutoff = 0.7, by = "p.adjust", select_fun = min)

plotData4A <- cmpGO_simplified@compareClusterResult |>
  dplyr::filter(
    (Cluster == "IFNG" & ID %in% c(
      "GO:0034341", "GO:0140888", "GO:0019221", "GO:0019882", "GO:0042267"
    )) |
      (Cluster == "Inflammatory" & ID %in% c(
        "GO:0006956", "GO:0006959", "GO:0030595", "GO:0050900", "GO:0006910"
      )) |
      (Cluster == "WoundHealing" & ID %in% c(
        "GO:0032964", "GO:0050673", "GO:0061448", "GO:0030111", "GO:0030510"
      ))
  ) |>
  dplyr::mutate(
    GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    logp = -log(p.adjust),
    Cluster = factor(Cluster, levels = c("WoundHealing","IFNG","Inflammatory")),
    Description = factor(
      Description,
      levels = c("humoral immune response","leukocyte migration","leukocyte chemotaxis",
                 "complement activation","phagocytosis, recognition",
                 "cytokine-mediated signaling pathway","response to type II interferon",
                 "antigen processing and presentation","natural killer cell mediated cytotoxicity",
                 "interferon-mediated signaling pathway",
                 "epithelial cell proliferation","connective tissue development",
                 "regulation of Wnt signaling pathway","regulation of BMP signaling pathway",
                 "collagen biosynthetic process")
    )
  )

fig4A_CompareClusterGO <- ggplot(plotData4A, aes(Cluster, Description)) +
  geom_point(aes(fill = logp, size = GeneRatio_num), shape = 21) +
  theme_bw() +
  scale_fill_gradient(low = "#2dabb9", high = "#FF9F1C") +
  labs(x = NULL, y = NULL, fill = "-Log10(P.Adjust)", size = "Gene Ratio") +
  annotate("rect", fill = "#2EC4B6", alpha = 0.3, xmax = 1.5, xmin = 0,   ymax = 15.5, ymin = 0) +
  annotate("rect", fill = "#BDD5EA", alpha = 0.3, xmax = 2.5, xmin = 1.5, ymax = 15.5, ymin = 0) +
  annotate("rect", fill = "#FFA5AB", alpha = 0.3, xmax = Inf, xmin = 2.5, ymax = 15.5, ymin = 0) +
  theme(
    axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 10, angle = 45, colour = "black"),
    axis.text.y = element_text(size = 10, hjust = 1, colour = "black"),
    legend.position = "right",
    legend.text = element_text(size = 10),
    panel.border = element_rect(size = 1)
  )

# ggsave(file.path(OUTPUT_FIG_DIR, "fig4a_comparecluster_go.pdf"), fig4A_CompareClusterGO, width = 6.5, height = 6)

##### 05) B. Compare Cluster Function — GSVA #####
hallmarks.data <- read.gmt(file.path(INPUT_ANALYSIS_DIR, "h.all.v2025.1.Hs.symbols.gmt"))
KEGG.data      <- read.gmt(file.path(INPUT_ANALYSIS_DIR, "c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt"))
GO_BP.data     <- read.gmt(file.path(INPUT_ANALYSIS_DIR, "c5.go.bp.v2025.1.Hs.symbols.gmt"))
GO_CC.data     <- read.gmt(file.path(INPUT_ANALYSIS_DIR, "c5.go.cc.v2025.1.Hs.symbols.gmt"))
GO_MF.data     <- read.gmt(file.path(INPUT_ANALYSIS_DIR, "c5.go.mf.v2025.1.Hs.symbols.gmt"))
Reactome.data  <- read.gmt(file.path(INPUT_ANALYSIS_DIR, "c2.cp.reactome.v2025.1.Hs.symbols.gmt"))

gsvaCombined <- dplyr::bind_rows(hallmarks.data, KEGG.data, GO_BP.data, GO_CC.data, GO_MF.data, Reactome.data)
geneSets <- gsvaCombined |>
  dplyr::select(term, gene) |>
  as.data.frame() |>
  split(.$term) |>
  lapply(function(x) unique(x$gene))

gsvaRes <- gsva(as.matrix(luadTpm), geneSets)
gsvaRes <- as.matrix(gsvaRes)

runLimmaOnMatrix <- function(mat, idsGroup, idsAll, coefName) {
  m <- mat[, intersect(colnames(mat), idsAll), drop = FALSE]
  grp <- ifelse(colnames(m) %in% idsGroup, coefName, "Others")
  design <- model.matrix(~ 0 + factor(grp, levels = c(coefName, "Others")))
  colnames(design) <- c(coefName, "Others")
  fit <- lmFit(m, design) |> contrasts.fit(makeContrasts(contrasts = paste0(coefName, "-Others"), levels = design)) |> eBayes()
  tt <- topTable(fit, number = Inf) |> tibble::rownames_to_column("SYMBOL") |> dplyr::select(SYMBOL, logFC)
  tt
}

degWH_gsva   <- runLimmaOnMatrix(gsvaRes, idsWH,  idsALL, "WoundHealing")  |> dplyr::mutate(group = "WoundHealing") |> dplyr::arrange(desc(logFC))
degIFNG_gsva <- runLimmaOnMatrix(gsvaRes, idsIFNG,idsALL, "IFNGDominant")  |> dplyr::mutate(group = "IFNG")         |> dplyr::arrange(desc(logFC))
degINF_gsva  <- runLimmaOnMatrix(gsvaRes, idsINF, idsALL, "Inflammatory")  |> dplyr::mutate(group = "Inflammatory") |> dplyr::arrange(desc(logFC))

curatedGeneSet <- c(
  dplyr::slice(degWH_gsva,   1:5)$SYMBOL,
  dplyr::slice(degIFNG_gsva, 1:5)$SYMBOL,
  dplyr::slice(degINF_gsva,  1:5)$SYMBOL
)
curatedGeneSet <- factor(curatedGeneSet, levels = curatedGeneSet)

plotData4B <- Reduce(function(x, y) merge(x, y, by = "SYMBOL"),
                     list(
                       degWH_gsva[c("SYMBOL","logFC")],
                       degIFNG_gsva[c("SYMBOL","logFC")],
                       degINF_gsva[c("SYMBOL","logFC")]
                     )) |>
  tibble::column_to_rownames("SYMBOL") |>
  { colnames(.) <- c("Wound Healing","IFN","Inflammatory"); . } |>
  t() |>
  scale() |>
  t()

col_fun_4B <- colorRamp2(c(min(-1.5), median(0), max(1.5)), c("#2dabb9","white","#FF9F1C"))
rowAnno_4B <- rowAnnotation(
  Cluster = factor(rep(c("Wound Healing","IFN","Inflammatory"), each = 5),
                   levels = c("Wound Healing","IFN","Inflammatory")),
  col = list(Cluster = c("Wound Healing"="#2EC4B6","IFN"="#BDD5EA","Inflammatory"="#FFA5AB"))
)

fig4B_GSVAHeatmap <- Heatmap(
  as.matrix(plotData4B[curatedGeneSet, c("Wound Healing","IFN","Inflammatory")]),
  row_split = factor(c(rep("Wound Healing",5),rep("IFN",5),rep("Inflammatory",5)),
                     levels = c("Wound Healing","IFN","Inflammatory")),
  cluster_rows = FALSE, cluster_columns = FALSE,
  col = col_fun_4B, left_annotation = rowAnno_4B,
  name = "Z-Score",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10, fontface = "bold", angle = 45)
)
# ggsave 需用 gridGraphics 保存 ComplexHeatmap；这里建议在主脚本中用 pdf(); draw(); dev.off()

##### 06) C. Compare Cluster — Metabolism (METAFlux) #####
data("nutrient_lookup_files")

metabolism_scores <- calculate_reaction_score(luadTpm)
flux <- compute_flux(mras = metabolism_scores, medium = human_blood)

cbrt <- function(x) sign(x) * abs(x)^(1/3)
flux <- cbrt(flux)

pathways <- unique(unlist(human_gem$SUBSYSTEM))
pathway_score <- lapply(pathways, function(pw) {
  idx <- which(unlist(human_gem$SUBSYSTEM) == pw)
  colMeans(abs(flux[idx, , drop = FALSE]))
})
names(pathway_score) <- pathways
all_pathway_score <- as.data.frame(do.call(rbind, pathway_score))
colnames(all_pathway_score) <- colnames(flux)

degIFNG_meta <- runLimmaOnMatrix(all_pathway_score, idsIFNG, idsALL, "IFNGDominant") |> dplyr::mutate(group = "IFNG") |> dplyr::arrange(desc(logFC))
degWH_meta   <- runLimmaOnMatrix(all_pathway_score, idsWH,   idsALL, "WoundHealing") |> dplyr::mutate(group = "WoundHealing") |> dplyr::arrange(desc(logFC))
degINF_meta  <- runLimmaOnMatrix(all_pathway_score, idsINF,  idsALL, "Inflammatory") |> dplyr::mutate(group = "Inflammatory") |> dplyr::arrange(desc(logFC))

topWithP <- function(mat, idsGroup, idsAll, coefName) {
  m <- mat[, intersect(colnames(mat), idsAll), drop = FALSE]
  grp <- ifelse(colnames(m) %in% idsGroup, coefName, "Others")
  design <- model.matrix(~ 0 + factor(grp, levels = c(coefName, "Others")))
  colnames(design) <- c(coefName, "Others")
  fit <- lmFit(m, design) |> contrasts.fit(makeContrasts(contrasts = paste0(coefName, "-Others"), levels = design)) |> eBayes()
  tt <- topTable(fit, number = Inf) |> tibble::rownames_to_column("SYMBOL")
  tt
}
ttIFNG_meta <- topWithP(all_pathway_score, idsIFNG, idsALL, "IFNGDominant")
ttWH_meta   <- topWithP(all_pathway_score, idsWH,   idsALL, "WoundHealing")
ttINF_meta  <- topWithP(all_pathway_score, idsINF,  idsALL, "Inflammatory")

curatedMetaPathway <- c(
  ttWH_meta  |> dplyr::filter(P.Value < 0.05 & logFC > 0) |> dplyr::slice(1:5) |> dplyr::pull(SYMBOL),
  ttIFNG_meta|> dplyr::filter(P.Value < 0.05 & logFC > 0) |> dplyr::slice(1:5) |> dplyr::pull(SYMBOL),
  ttINF_meta |> dplyr::filter(P.Value < 0.05 & logFC > 0) |> dplyr::slice(1:5) |> dplyr::pull(SYMBOL)
)
curatedMetaPathway <- factor(curatedMetaPathway, levels = curatedMetaPathway)

plotData4C <- Reduce(function(x, y) merge(x, y, by = "SYMBOL"),
                     list(
                       ttWH_meta[c("SYMBOL","logFC")],
                       ttIFNG_meta[c("SYMBOL","logFC")],
                       ttINF_meta[c("SYMBOL","logFC")]
                     )) |>
  tibble::column_to_rownames("SYMBOL") |>
  { colnames(.) <- c("Wound Healing","IFN","Inflammatory"); . } |>
  t() |>
  scale() |>
  t()

col_fun_4C <- colorRamp2(c(min(-1.5), median(0), max(1.5)), c("#2dabb9","white","#FF9F1C"))
rowAnno_4C <- rowAnnotation(
  Cluster = factor(rep(c("Wound Healing","IFN","Inflammatory"), each = 5),
                   levels = c("Wound Healing","IFN","Inflammatory")),
  col = list(Cluster = c("Wound Healing"="#2EC4B6","IFN"="#BDD5EA","Inflammatory"="#FFA5AB"))
)

fig4C_MetaHeatmap <- Heatmap(
  as.matrix(plotData4C[curatedMetaPathway, c("Wound Healing","IFN","Inflammatory")]),
  row_split = factor(c(rep("Wound Healing",5),rep("IFN",5),rep("Inflammatory",5)),
                     levels = c("Wound Healing","IFN","Inflammatory")),
  cluster_rows = FALSE, cluster_columns = FALSE,
  col = col_fun_4C, left_annotation = rowAnno_4C,
  name = "Z-Score",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10, fontface = "bold", angle = 45)
)

##### 07) D. Pathway activity (PROGENy) & TF activity (CoReg/DoRothEA) #####
# Pathway activity by decoupleR::PROGENy
net_pathway <- decoupleR::get_progeny(organism = "human", top = 500)

pathway_decouple <- decoupleR::run_mlm(
  mat     = luadTpm,
  net     = net_pathway,
  .source = "source",
  .target = "target",
  .mor    = "weight",
  minsize = 5
) |>
  dplyr::filter(p_value < 0.05) |>
  dplyr::left_join(
    immuneSubtypeAssignments |> dplyr::select(ID, clusterLabel),
    dplyr::join_by("condition" == "ID")
  ) |>
  dplyr::mutate(
    clusterLabel = factor(clusterLabel, levels = c("Wound Healing","IFN-γ Dominant","Inflammatory"))
  )

comparisons_list <- list(
  c("Wound Healing","IFN-γ Dominant"),
  c("IFN-γ Dominant","Inflammatory"),
  c("Wound Healing","Inflammatory")
)

fig4D_PathwayBox <- ggplot(pathway_decouple, aes(clusterLabel, score, fill = clusterLabel)) +
  geom_boxplot(outlier.shape = 21, color = "black") +
  scale_fill_manual(values = setNames(PALETTE_CLUSTER3, c("Wound Healing","IFN-γ Dominant","Inflammatory"))) +
  facet_wrap(vars(source), nrow = 1, scales = "free_y") +
  ggpubr::stat_compare_means(
    method = "t.test",
    comparisons = comparisons_list,
    label = "p.signif",
    step.increase = 0.05,
    label.y = 1.1 * max(pathway_decouple$score, na.rm = TRUE)
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

ggsave(file.path(OUTPUT_FIG_DIR, "fig4d_pathway_boxplots.pdf"), fig4D_PathwayBox, width = 12, height = 4)

# TF activity by decoupleR::get_collectri + ULM
net_tf <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)

tf_decouple <- decoupleR::run_ulm(
  mat     = luadTpm,
  net     = net_tf,
  .source = "source",
  .target = "target",
  .mor    = "mor",
  minsize = 5
)

n_tfs <- 500
tf_mat <- tf_decouple |>
  tidyr::pivot_wider(id_cols = "condition", names_from = "source", values_from = "score") |>
  tibble::column_to_rownames("condition") |>
  as.matrix()
tf_var <- tf_decouple |>
  dplyr::group_by(source) |>
  dplyr::summarise(std = stats::sd(score)) |>
  dplyr::arrange(desc(abs(std))) |>
  dplyr::slice(1:n_tfs) |>
  dplyr::pull(source)
tf_acts_mat <- as.data.frame(t(tf_mat[, tf_var, drop = FALSE]))

runTfDiff <- function(mat, idsGroup, idsAll, coefName) {
  m <- mat[, intersect(colnames(mat), idsAll), drop = FALSE]
  grp <- ifelse(colnames(m) %in% idsGroup, coefName, "Others")
  design <- model.matrix(~ 0 + factor(grp, levels = c(coefName, "Others")))
  colnames(design) <- c(coefName, "Others")
  fit <- lmFit(m, design) |> contrasts.fit(makeContrasts(contrasts = paste0(coefName, "-Others"), levels = design)) |> eBayes()
  topTable(fit, number = Inf) |> tibble::rownames_to_column("SYMBOL") |> dplyr::select(SYMBOL, logFC)
}
degWH_tf   <- runTfDiff(tf_acts_mat, idsWH,  idsALL, "WoundHealing")   |> dplyr::mutate(group = "WoundHealing")   |> dplyr::arrange(desc(logFC))
degIFNG_tf <- runTfDiff(tf_acts_mat, idsIFNG,idsALL, "IFNGDominant")   |> dplyr::mutate(group = "IFNG")           |> dplyr::arrange(desc(logFC))
degINF_tf  <- runTfDiff(tf_acts_mat, idsINF, idsALL, "Inflammatory")   |> dplyr::mutate(group = "Inflammatory")   |> dplyr::arrange(desc(logFC))

curatedDE_TF <- c(
  dplyr::slice(degWH_tf,   1:5) |> dplyr::pull(SYMBOL),
  dplyr::slice(degIFNG_tf, 1:5) |> dplyr::pull(SYMBOL),
  dplyr::slice(degINF_tf,  1:5) |> dplyr::pull(SYMBOL)
)

plotData4E <- Reduce(function(x, y) merge(x, y, by = "SYMBOL"),
                     list(
                       degWH_tf[c("SYMBOL","logFC")],
                       degIFNG_tf[c("SYMBOL","logFC")],
                       degINF_tf[c("SYMBOL","logFC")]
                     )) |>
  tibble::column_to_rownames("SYMBOL") |>
  { colnames(.) <- c("Wound Healing","IFN","Inflammatory"); . } |>
  t() |>
  scale() |>
  t() |>
  as.data.frame() |>
  dplyr::filter(rownames(.) %in% curatedDE_TF)

col_fun_4E <- colorRamp2(c(min(-1.5), median(0), max(1.5)), c("#2dabb9","white","#FF9F1C"))
rowAnno_4E <- rowAnnotation(
  Cluster = factor(rep(c("Wound Healing","IFN","Inflammatory"), each = 5),
                   levels = c("Wound Healing","IFN","Inflammatory")),
  col = list(Cluster = c("Wound Healing"="#2EC4B6","IFN"="#BDD5EA","Inflammatory"="#FFA5AB"))
)

fig4E_TFHeatmap <- Heatmap(
  as.matrix(plotData4E[curatedDE_TF, c("Wound Healing","IFN","Inflammatory")]),
  row_split = factor(c(rep("Wound Healing",5),rep("IFN",5),rep("Inflammatory",5)),
                     levels = c("Wound Healing","IFN","Inflammatory")),
  cluster_rows = FALSE, cluster_columns = FALSE,
  col = col_fun_4E, left_annotation = rowAnno_4E,
  name = "Z-Score",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10, fontface = "bold", angle = 45)
)

# 导出 ComplexHeatmap 建议使用 pdf+draw
# pdf(file.path(OUTPUT_FIG_DIR, "fig4b_gsva_heatmap.pdf"), width = 5.5, height = 6); draw(fig4B_GSVAHeatmap); dev.off()
# pdf(file.path(OUTPUT_FIG_DIR, "fig4c_meta_heatmap.pdf"), width = 5.5, height = 6); draw(fig4C_MetaHeatmap); dev.off()
pdf(file.path(OUTPUT_FIG_DIR, "fig4e_tf_heatmap.pdf"), width = 4, height = 7); draw(fig4E_TFHeatmap); dev.off()

