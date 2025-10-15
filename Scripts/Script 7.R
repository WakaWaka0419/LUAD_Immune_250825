##### 01) Setup & Packages #####
Sys.setenv(LANG = "EN")
options(stringsAsFactors = FALSE, encoding = "UTF-8")
set.seed(20240601)

library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(survival)
library(survminer)
library(ComplexHeatmap)
library(circlize)
library(forestplot)
library(GSVA)
library(ggcor)

##### 02) Constants (paths, palette, gene) #####
# Raw data (PanCanAtlas from local paths)
PAN_ANN_PATH   <- "~/Working_folder/DataBase/Pancanatlas/merged_sample_quality_annotations.tsv"
PAN_SURV_PATH  <- "~/Working_folder/DataBase/Pancanatlas/Survival_SupplementalTable_S1_20171025_xena_sp.tsv"
PAN_EXPR_PATH  <- "~/Working_folder/DataBase/Pancanatlas/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"

# Project IO
OUT_DATA_DIR   <- "output/script7/data"
OUT_FIG_DIR    <- "output/script7/figures"
dir.create(OUT_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG_DIR,  recursive = TRUE, showWarnings = FALSE)

OUT_SIMPLE_ANN <- file.path(OUT_DATA_DIR, "simple_sample_annotation.txt")
OUT_COX_SUM    <- file.path(OUT_DATA_DIR, "summary_cox_pancancer.txt")
OUT_KM_SUM     <- file.path(OUT_DATA_DIR, "summary_km_pancancer.txt")
OUT_PROG_HM    <- file.path(OUT_DATA_DIR, "summary_gene_prognostication_pancancer.txt")
FIG_PROG_HM    <- file.path(OUT_FIG_DIR,  "fig7a_prognostic_heatmap.pdf")
FIG_OS_FP      <- file.path(OUT_FIG_DIR,  "fig7a_forestplot_os.pdf")

# LUAD inputs for correlation analysis
LUAD_TPM_RDS   <- "./Input/TCGA-LUAD/TCGA-LUAD_mrna_expr_tpm.rdata"
CURATED_IMMUNE <- "./Input/Analysis_data/Curated_Immune_Cell_Signature.txt"
HALLMARKS_GMT  <- "./Input/Analysis_data/h.all.v2025.1.Hs.symbols.gmt"

# Colors
COL_BLUE   <- "#A0CEE3"
COL_YELLOW <- "#EFEFBE"
COL_SEA    <- "#37BCDF"
COL_GREEN  <- "#71BC5D"
COL_CHERRY <- "#E5588C"
COL_RED    <- "#E46A6B"
COL_PURPLE <- "#8959A3"

# Gene of interest
GOI <- "TPX2"

##### 03) Load Pan-Cancer data & preprocess #####
rawAnno <- read.delim(PAN_ANN_PATH, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
rawAnno$simple_barcode <- substr(rawAnno$aliquot_barcode, 1, 15)
samAnno <- rawAnno[!duplicated(rawAnno$simple_barcode), c("cancer type", "simple_barcode")]
samAnno <- samAnno[samAnno$`cancer type` != "", ]
write.table(samAnno, OUT_SIMPLE_ANN, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

surv <- read.delim(PAN_SURV_PATH, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

expr <- fread(PAN_EXPR_PATH, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE) |> as.data.frame()
rownames(expr) <- expr[, 1]; expr <- expr[, -1, drop = FALSE]
expr$gene <- sapply(strsplit(rownames(expr), "|", fixed = TRUE), `[`, 1)
expr <- expr[!duplicated(expr$gene), ]
rownames(expr) <- expr$gene; expr$gene <- NULL
expr[expr < 0] <- 0
colnames(expr) <- substr(colnames(expr), 1, 15)
gc()

if (!GOI %in% rownames(expr)) {
  warning("The gene ", GOI, " cannot be found!")
} else {
  message("The gene ", GOI, " can be matched!")
}

# Keep valid samples
expr_row <- expr[GOI, , drop = FALSE]
expr_row <- as.data.frame(t(na.omit(t(expr_row))))  # drop all-NA columns safely
keep_samples <- colnames(expr_row)

# Exclude LAML & keep tumor samples (TCGA code "0" at 14th char)
nonLAML_samples <- samAnno[samAnno$`cancer type` != "LAML", "simple_barcode"]
common_samples  <- intersect(intersect(keep_samples, nonLAML_samples), rownames(surv))
tumor_samples   <- common_samples[substr(common_samples, 14, 14) == "0"]

tumAnno <- samAnno[samAnno$simple_barcode %in% tumor_samples, ]
tumAnno <- tumAnno[order(tumAnno$`cancer type`), ]
tumors  <- unique(tumAnno$`cancer type`)

# Build per-sample survival table with GOI expression
exprSurv <- cbind.data.frame(
  expr = log2(as.numeric(expr[GOI, common_samples]) + 1),
  surv[common_samples, c("OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time", "PFI", "PFI.time")]
)
rownames(exprSurv) <- common_samples

##### 04) Per-cancer Cox (OS/DSS/DFI/PFI) #####
out_cox <- NULL
for (tum in tumors) {
  sam <- tumAnno[tumAnno$`cancer type` == tum, "simple_barcode"]
  sub <- exprSurv[sam, , drop = FALSE]

  # OS
  c1 <- summary(coxph(Surv(OS.time, OS) ~ expr, data = sub))
  out_cox <- rbind(out_cox, data.frame(
    tumor = tum, event = "OS",
    beta = c1$coefficients[1, 1],
    hr   = c1$coefficients[1, 2],
    lower= c1$conf.int[1, 3],
    upper= c1$conf.int[1, 4],
    p    = c1$coefficients[1, 5]
  ))

  # DSS
  c2 <- summary(coxph(Surv(DSS.time, DSS) ~ expr, data = sub))
  out_cox <- rbind(out_cox, data.frame(
    tumor = tum, event = "DSS",
    beta = c2$coefficients[1, 1],
    hr   = c2$coefficients[1, 2],
    lower= c2$conf.int[1, 3],
    upper= c2$conf.int[1, 4],
    p    = c2$coefficients[1, 5]
  ))

  # DFI (skip where not applicable)
  if (tum %in% c("SKCM", "THYM", "UVM", "GBM")) {
    out_cox <- rbind(out_cox, data.frame(
      tumor = tum, event = "DFI",
      beta = NA, hr = NA, lower = NA, upper = NA, p = NA
    ))
  } else {
    c3 <- summary(coxph(Surv(DFI.time, DFI) ~ expr, data = sub))
    out_cox <- rbind(out_cox, data.frame(
      tumor = tum, event = "DFI",
      beta = c3$coefficients[1, 1],
      hr   = c3$coefficients[1, 2],
      lower= c3$conf.int[1, 3],
      upper= c3$conf.int[1, 4],
      p    = c3$coefficients[1, 5]
    ))
  }

  # PFI (fix: lower/upper indices)
  c4 <- summary(coxph(Surv(PFI.time, PFI) ~ expr, data = sub))
  out_cox <- rbind(out_cox, data.frame(
    tumor = tum, event = "PFI",
    beta = c4$coefficients[1, 1],
    hr   = c4$coefficients[1, 2],
    lower= c4$conf.int[1, 3],  # FIXED (was 4)
    upper= c4$conf.int[1, 4],
    p    = c4$coefficients[1, 5]
  ))
}
# write.table(out_cox, OUT_COX_SUM, sep = "\t", row.names = FALSE, quote = FALSE)

##### 05) Per-cancer KM split by best cut (OS/DSS/DFI/PFI) #####
minprop <- 0.2
out_km <- NULL
km_block <- function(sub_df, time_col, event_col) {
  bestcut <- surv_cutpoint(sub_df, time = time_col, event = event_col, variables = "expr", minprop = minprop)
  cutoff  <- bestcut$cutpoint[1, 1]
  sub_df$group <- factor(ifelse(sub_df$expr > cutoff, "High", "Low"), levels = c("Low", "High"))
  fitd <- survdiff(Surv(sub_df[[time_col]], sub_df[[event_col]]) ~ group, data = sub_df, na.action = na.exclude)
  p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  HR    <- (fitd$obs[2] / fitd$exp[2]) / (fitd$obs[1] / fitd$exp[1])
  lower <- exp(log(HR) - qnorm(0.975) * sqrt(1 / fitd$exp[2] + 1 / fitd$exp[1]))
  upper <- exp(log(HR) + qnorm(0.975) * sqrt(1 / fitd$exp[2] + 1 / fitd$exp[1]))
  list(HR = HR, lower = lower, upper = upper, p = p.val)
}

for (tum in tumors) {
  sam <- tumAnno[tumAnno$`cancer type` == tum, "simple_barcode"]
  sub <- exprSurv[sam, , drop = FALSE]

  # OS
  os <- km_block(sub, "OS.time", "OS")
  out_km <- rbind(out_km, data.frame(tumor = tum, event = "OS",  hr = os$HR, lower = os$lower, upper = os$upper, p = os$p))

  # DSS
  dss <- km_block(sub, "DSS.time", "DSS")
  out_km <- rbind(out_km, data.frame(tumor = tum, event = "DSS", hr = dss$HR, lower = dss$lower, upper = dss$upper, p = dss$p))

  # DFI
  if (tum %in% c("SKCM", "THYM", "UVM", "GBM")) {
    out_km <- rbind(out_km, data.frame(tumor = tum, event = "DFI", hr = NA, lower = NA, upper = NA, p = NA))
  } else {
    dfi <- km_block(sub, "DFI.time", "DFI")
    out_km <- rbind(out_km, data.frame(tumor = tum, event = "DFI", hr = dfi$HR, lower = dfi$lower, upper = dfi$upper, p = dfi$p))
  }

  # PFI
  pfi <- km_block(sub, "PFI.time", "PFI")
  out_km <- rbind(out_km, data.frame(tumor = tum, event = "PFI", hr = pfi$HR, lower = pfi$lower, upper = pfi$upper, p = pfi$p))
}
# write.table(out_km, OUT_KM_SUM, sep = "\t", row.names = FALSE, quote = FALSE)

##### 06) Build heatmap matrix (Risky / Protective / N/A / Nonsense) #####
hm_input <- NULL
for (tum in tumors) {
  cx  <- out_cox[out_cox$tumor == tum, ]
  km  <- out_km[out_km$tumor == tum, ]

  cx$dir <- ifelse(cx$hr > 1 & cx$p < 0.05, "Risky",
                   ifelse(cx$hr < 1 & cx$p < 0.05, "Protective", "Nonsense"))
  km$dir <- ifelse(km$hr > 1 & km$p < 0.05, "Risky",
                   ifelse(km$hr < 1 & km$p < 0.05, "Protective", "Nonsense"))

  hm_input <- rbind(hm_input, data.frame(
    OS.cox  = cx[match("OS",  cx$event), "dir"],
    OS.km   = km[match("OS",  km$event), "dir"],
    DSS.cox = cx[match("DSS", cx$event), "dir"],
    DSS.km  = km[match("DSS", km$event), "dir"],
    DFI.cox = cx[match("DFI", cx$event), "dir"],
    DFI.km  = km[match("DFI", km$event), "dir"],
    PFI.cox = cx[match("PFI", cx$event), "dir"],
    PFI.km  = km[match("PFI", km$event), "dir"],
    row.names = tum
  ))
}
hm_input[is.na(hm_input)] <- "N/A"
write.table(hm_input, OUT_PROG_HM, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# Encode to numeric for heatmap colors
mat_hm <- hm_input
mat_hm[mat_hm == "Risky"]      <- 0
mat_hm[mat_hm == "Protective"] <- 1
mat_hm[mat_hm == "N/A"]        <- 2
mat_hm[mat_hm == "Nonsense"]   <- 3
mat_hm <- as.matrix(apply(mat_hm, 2, as.numeric))
rownames(mat_hm) <- rownames(hm_input)

annCol <- data.frame(
  Method = rep(c("Cox", "Log-rank"), 4),
  SurvivalType = rep(c("OS", "DSS", "DFI", "PFI"), each = 2),
  row.names = colnames(mat_hm), check.names = FALSE
)
annColors <- list(
  Method       = c(Cox = COL_BLUE, "Log-rank" = COL_YELLOW),
  SurvivalType = c(OS = COL_SEA, DSS = COL_GREEN, DFI = COL_CHERRY, PFI = COL_PURPLE)
)
col_fun <- c(COL_RED, COL_GREEN, "grey50", "white")  # 0/1/2/3

hm <- Heatmap(
  matrix              = mat_hm,
  border              = "black",
  rect_gp             = gpar(col = "black"),
  name                = "Prognostic role",
  cluster_rows        = FALSE,
  cluster_columns     = FALSE,
  col                 = col_fun,
  show_row_names      = TRUE,
  show_column_names   = FALSE,
  row_names_side      = "left",
  top_annotation      = HeatmapAnnotation(
    df    = annCol,
    col   = annColors,
    gp    = gpar(col = "black"),
    border= TRUE
  ),
  width               = grid::unit(8, "cm"),
  height              = grid::unit(15, "cm"),
  heatmap_legend_param = list(
    at      = c(0, 1, 2, 3),
    labels  = c("Risky", "Protective", "N/A", "Nonsense"),
    legend_gp = grid::gpar(fill = col_fun)
  )
)
pdf(FIG_PROG_HM, width = 10, height = 10)
draw(hm)
invisible(dev.off())

##### 07) Forest plot (OS) #####
fp <- out_cox[out_cox$event == "OS", c("tumor", "event", "beta", "hr", "lower", "upper", "p")]
tabletext <- cbind(
  c("Cancers", fp$tumor),
  c("p value", ifelse(round(fp$p, 3) < 0.001, "<0.001", sprintf("%.3f", round(fp$p, 3)))),
  c("HR (95% CI)", paste0(sprintf("%.3f", fp$hr), " (",
                          sprintf("%.3f", fp$lower), "-",
                          sprintf("%.3f", fp$upper), ")"))
)
pdf(FIG_OS_FP, width = 8, height = 8)
forestplot(
  labeltext = tabletext,
  mean      = c(NA, fp$hr),
  lower     = c(NA, fp$lower),
  upper     = c(NA, fp$upper),
  graph.pos = 4,
  graphwidth= unit(.3, "npc"),
  fn.ci_norm= "fpDrawDiamondCI",
  col       = fpColors(box = "grey50", lines = "grey50", zero = "black"),
  boxsize   = c(NA, ifelse(fp$p < 0.05, 0.8, 0.4)),
  lwd.ci    = 1,
  ci.vertices.height = 0.1, ci.vertices = FALSE,
  zero      = 1,
  lwd.zero  = 2,
  xticks    = c(0, 1, 2, 4, 8, 13),
  lwd.xaxis = 2,
  xlab      = "Hazard ratio",
  hrzl_lines= list("1" = gpar(lwd = 2, col = "black"),
                   as.character(nrow(tabletext) + 1) := gpar(lwd = 2, col = "black")),
  txt_gp    = fpTxtGp(label = gpar(cex = 1), ticks = gpar(cex = 0.85), xlab = gpar(cex = 1), title = gpar(cex = 1.2)),
  lineheight= unit(.55, "cm"),
  colgap    = unit(0.4, "cm"),
  mar       = unit(rep(1.5, 4), "cm"),
  new_page  = FALSE
)
invisible(dev.off())

##### 08) Correlation links (TPX2 vs immune cell sets & HALLMARKS) #####
# Helper: read .gmt to list
gmt2list <- function(annofile){
  if (!file.exists(annofile)) stop("GMT file not found: ", annofile)
  if (tools::file_ext(annofile) == "xz") {
    con <- xzfile(annofile); x <- scan(con, what = "", sep = "\n", quiet = TRUE); close(con)
  } else if (tools::file_ext(annofile) == "gmt") {
    x <- scan(annofile, what = "", sep = "\n", quiet = TRUE)
  } else stop("Only .gmt or .gmt.xz supported")
  y <- strsplit(x, "\t"); names(y) <- sapply(y, `[[`, 1)
  lapply(y, `[`, c(-1, -2))
}

# LUAD TPM
load(LUAD_TPM_RDS) # -> mrna_expr_tpm
luad_tpm <- {
  x <- mrna_expr_tpm[, substr(colnames(mrna_expr_tpm), 14, 15) != "11", drop = FALSE]
  colnames(x) <- substr(colnames(x), 1, 16)
  x[, !duplicated(colnames(x)), drop = FALSE]
}
TPX2_vec <- luad_tpm["TPX2", ]

# Curated immune signatures -> list(CellType = c(genes))
curated_imm <- data.table::fread(CURATED_IMMUNE)
curated_imm <- curated_imm |>
  mutate(Symbol = toupper(str_replace_all(Symbol, "\\s+", ""))) |>
  filter(!is.na(Symbol), nzchar(Symbol), !is.na(CellType), nzchar(CellType)) |>
  distinct(CellType, Symbol)

celltype_sets <- curated_imm |>
  group_by(CellType) |>
  summarise(genes = list(unique(Symbol)), .groups = "drop")
imm_sets <- setNames(celltype_sets$genes, celltype_sets$CellType)

# HALLMARKS
hallmark_sets <- gmt2list(HALLMARKS_GMT)

# GSVA ssGSEA scores
gset_score <- gsva(as.matrix(luad_tpm), hallmark_sets, method = "ssgsea")
imm_score  <- gsva(as.matrix(luad_tpm), imm_sets,      method = "ssgsea")

# Bind TPX2 vector as a row for correlation
imm_score2 <- rbind.data.frame(imm_score, TPX2 = TPX2_vec)
gset_score2<- rbind.data.frame(gset_score, TPX2 = TPX2_vec)

# Correlation: TPX2 vs immune sets
immCor <- do.call(rbind, lapply(rownames(imm_score), function(pth){
  cr <- cor.test(as.numeric(imm_score[pth, ]), as.numeric(TPX2_vec), method = "pearson")
  data.frame(gene = "TPX2", path = pth, r = cr$estimate, p = cr$p.value)
}))
immCor$sign <- ifelse(immCor$r > 0, "pos", "neg")
immCor$absR <- abs(immCor$r)
immCor$rSeg <- cut(immCor$absR, c(0, .25, .5, .75, 1), labels = c("0.25","0.50","0.75","1.00"), include.lowest = TRUE)
immCor$pSeg <- cut(immCor$p,    c(0, .001, .01, .05, 1), labels = c("<0.001","<0.01","<0.05","ns"), include.lowest = TRUE)
immCor[nrow(immCor), "pSeg"] <- "Not Applicable"
immCor$rSeg <- factor(immCor$rSeg, levels = c("0.25","0.50","0.75","1.00"))
immCor$pSeg <- factor(immCor$pSeg, levels = c("<0.001","<0.01","<0.05","Not Applicable","ns"))
immCor$sign <- factor(immCor$sign, levels = c("pos","neg"))

p1 <- quickcor(t(imm_score2), type = "lower", show.diag = TRUE) +
  geom_colour() +
  anno_link(
    data    = immCor,
    mapping = aes(colour = pSeg, size = rSeg, linetype = sign),
    spec.key = "gene",
    env.key  = "path",
    diag.label = FALSE
  ) +
  scale_size_manual(values = c(0.5, 1, 1.5, 2)) +
  scale_color_manual(values = c("#19A078","#DA6003","#7570B4","#E8288E","#65A818")) +
  scale_fill_gradient2(low = "#9483E1", mid = "white", high = "#E11953", midpoint = 0) +
  remove_axis("x")
p1

# Correlation: TPX2 vs HALLMARKS (top 15 by |r|)
hallCor <- do.call(rbind, lapply(rownames(gset_score), function(pth){
  cr <- cor.test(as.numeric(gset_score[pth, ]), as.numeric(TPX2_vec), method = "pearson")
  data.frame(gene = "TPX2", path = pth, r = cr$estimate, p = cr$p.value)
}))
hallCor$absR <- abs(hallCor$r)
hallCor <- hallCor |>
  arrange(desc(absR)) |>
  dplyr::slice(1:15) |>
  mutate(
    sign = ifelse(r > 0, "pos", "neg"),
    rSeg = cut(absR, c(0, .25, .5, .75, 1), labels = c("0.25","0.50","0.75","1.00"), include.lowest = TRUE),
    pSeg = cut(p,    c(0, .001, .01, .05, 1), labels = c("<0.001","<0.01","<0.05","ns"), include.lowest = TRUE)
  )
hallCor[nrow(hallCor), "pSeg"] <- "Not Applicable"
hallCor$rSeg <- factor(hallCor$rSeg, levels = c("0.25","0.50","0.75","1.00"))
hallCor$pSeg <- factor(hallCor$pSeg, levels = c("<0.001","<0.01","<0.05","Not Applicable","ns"))
hallCor$sign <- factor(hallCor$sign, levels = c("pos","neg"))

gset_top <- gset_score2[unique(c(rownames(gset_score)[rownames(gset_score) %in% hallCor$path], "TPX2")), , drop = FALSE]

p2 <- quickcor(t(gset_top), type = "upper", show.diag = TRUE) +
  geom_colour() +
  anno_link(
    data    = hallCor,
    mapping = aes(colour = pSeg, size = rSeg, linetype = sign),
    spec.key = "gene",
    env.key  = "path",
    diag.label = FALSE
  ) +
  scale_size_manual(values = c(0.5, 1, 1.5, 2)) +
  scale_color_manual(values = c("#19A078","#DA6003","#7570B4","#E8288E","#65A818")) +
  scale_fill_gradient2(low = "#9483E1", mid = "white", high = "#E11953", midpoint = 0) +
  remove_axis("x")
p2
# ggsave(file.path(OUT_FIG_DIR, "fig7b_imm_corlinks.pdf"), p1, width = 10, height = 8)
# ggsave(file.path(OUT_FIG_DIR, "fig7b_hallmark_corlinks.pdf"), p2, width = 10, height = 8)
