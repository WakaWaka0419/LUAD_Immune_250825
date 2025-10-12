cindex_dis_all <- function (object, color = NULL, dataset_col = NULL, validate_set = NULL,
                            order = NULL, width = NULL, height = NULL,
                            tile_border_color = "grey40", tile_border_size = 0.4,
                            n_models = NULL, pick_top = TRUE)
{
  library(ggplot2)
  library(aplot)
  
  # --- 参数默认 ---
  if (is.null(width))  width  <- 0.35
  if (is.null(height)) height <- 0.01
  
  # color: 前3个以前用于梯度，这里仅保留第4个用于右侧单柱填充
  if (is.null(color)) {
    color <- c("#3a83a7", "#f2f1b9", "#DC0000FF", "#3C5488FF")
  }
  
  if (is.null(dataset_col)) {
    dataset_col <- c("#E64B35CC","#4DBBD5CC","#00A087CC","#3C5488CC","#F39B7FCC","#FD8D3CFF",
                     "#756BB1FF", "#636363FF", "#6BAED6FF", "#FD8D3CFF",
                     "#74C476FF", "#9E9AC8FF", "#969696FF", "#9ECAE1FF",
                     "#FDAE6BFF", "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF",
                     "#C6DBEFFF", "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF",
                     "#D9D9D9FF")
  }
  
  # --- 连续渐变色板：Spectral 反向（无需安装任何包） ---
  spectral_rev_hex <- c(
    "#5E4FA2","#3288BD","#66C2A5","#ABDDA4","#E6F598",
    "#FFFFBF","#FEE08B","#FDAE61","#F46D43","#D53E4F","#9E0142"
  )
  grad_cols <- grDevices::colorRampPalette(spectral_rev_hex)(20)
  
  # --- 数据整理 ---
  cindex <- object[["Cindex.res"]]
  cindex$Cindex <- as.numeric(sprintf("%.2f", cindex$Cindex))
  
  # 总体（训练+验证）平均 C-index
  mean_cindex <- aggregate(x = cindex$Cindex, by = list(cindex$Model), FUN = mean)
  colnames(mean_cindex) <- c("Model", "mean")
  mean_cindex$mean  <- as.numeric(sprintf("%.3f", mean_cindex$mean))
  mean_cindex$Value <- "Mean C-index (train + validate)"
  
  # 如果指定了要显示的模型数量，则按表现筛选
  if (!is.null(n_models)) {
    n_models <- max(1, min(n_models, nrow(mean_cindex)))
    ord <- order(mean_cindex$mean, decreasing = pick_top)  # TRUE=从高到低
    keep_models <- mean_cindex$Model[ord][seq_len(n_models)]
    mean_cindex <- mean_cindex[mean_cindex$Model %in% keep_models, , drop = FALSE]
    cindex      <- cindex[cindex$Model %in% keep_models, , drop = FALSE]
  }
  
  # 模型排序：从低到高
  model_levels <- mean_cindex[order(mean_cindex$mean, decreasing = FALSE), "Model"]
  
  # Cohort 标签
  labels <- data.frame(Cohort = unique(cindex$ID))
  if (!is.null(order)) {
    labels$Cohort <- factor(labels$Cohort, levels = order)
    cindex$ID     <- factor(cindex$ID,     levels = order)
  }
  
  cindex$Model      <- factor(cindex$Model,      levels = model_levels)
  mean_cindex$Model <- factor(mean_cindex$Model, levels = model_levels)
  
  # --- 作图 ---
  # 主热图：使用 Spectral(反向) 连续渐变
  p1 <- ggplot(cindex, aes(x = ID, y = Model)) +
    geom_tile(aes(fill = Cindex),
              color = tile_border_color, linewidth = tile_border_size) +
    geom_text(aes(label = Cindex), vjust = 0.5, color = "black", size = 2) +
    scale_fill_gradientn(
      colors = grad_cols,
      # 如需固定色条范围/刻度，可取消下一行两行注释自行设定
       limits = c(0.55, 0.75),
       breaks  = seq(0.55, 0.85, 0.1),
      name   = "C-index",
      guide  = guide_colorbar(
        title.position = "top",
        frame.colour = "black",
        ticks = TRUE,
        ticks.colour = "white",
        ticks.linewidth = 1
      )
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y  = element_text(size = 8),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_blank()
    )
  
  # 顶部 Cohort 色带：加边框
  p2 <- ggplot(labels, aes(Cohort, y = 1)) +
    geom_tile(aes(fill = Cohort),
              color = tile_border_color, linewidth = tile_border_size) +
    scale_fill_manual(values = dataset_col, name = "Cohort") +
    theme_void()
  
  # 右侧单柱：总体平均 C-index
  p3 <- ggplot(mean_cindex, aes(x = Model, y = mean, fill = Value)) +
    geom_col() +
    scale_fill_manual(values = color[4], name = "") +
    geom_text(aes(label = mean), vjust = 0.5, hjust = 1.2, size = 2) +
    coord_flip() +
    labs(y = "") +
    theme(
      axis.title  = element_text(size = 8),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x  = element_blank(),
      axis.text.y  = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white")
    )
  
  # 组合
  print(
    p1 %>%
      insert_top(p2, height = height) %>%
      insert_right(p3, width  = width)
  )
}

