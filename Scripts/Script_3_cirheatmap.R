library(circlize)
res <- readRDS("./Output/Script_3/Rdata/res.RDS")
res_width <- res$Cindex.res %>%
  pivot_wider(
    names_from = ID,
    values_from = Cindex
  ) %>%
  mutate(MeanCindex = rowMeans(across(starts_with("Dataset")), na.rm = TRUE)) %>%
  arrange(desc(MeanCindex)) %>%
  column_to_rownames(var = "Model")
# 定义颜色梯度：256个颜色
spectral_rev_hex <- c(
  "#5E4FA2","#3288BD","#66C2A5","#ABDDA4","#E6F598",
  "#FFFFBF","#FEE08B","#FDAE61","#F46D43","#D53E4F","#9E0142"
)
grad_cols <- grDevices::colorRampPalette(spectral_rev_hex)(256)

# 定义映射函数（这里按你的C-index取值范围，比如0.5-0.8）
col_fun <- colorRamp2(seq(0.5, 0.7, length.out = 256), grad_cols)

# 转换矩阵（前6列Dataset）
mat <- as.matrix(res_width[, 1:6, drop = FALSE])
rownames(mat) <- rownames(res_width)

pdf("plot.pdf", height = 8, width = 8)
circos.clear()
circos.par(gap.degree = 70, start.degree = 20, track.margin = c(0.001, 0.001))

# 绘制热图
circos.heatmap(
  mat,
  col = col_fun,
  rownames.side = "outside",
  cluster = FALSE,
  cell.border = "white",
  track.height = 0.30
)

# 图例
lgd <- Legend(
  title = "C-index",
  border = "black",
  grid_height = unit(3, "mm"),
  legend_width = unit(30, "mm"),
  at = c(0.5, 0.6, 0.7, 0.8),    # 根据范围设置刻度
  col_fun = col_fun,
  direction = "horizontal"
)
draw(lgd, x = unit(0.55, "npc"), y = unit(0.7, "npc"))

dev.off()
circos.clear()
