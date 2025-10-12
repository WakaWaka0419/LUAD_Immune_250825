meta_unicox_vis <- function(object, dataset_col = NULL, dataset) {
  library(forestploter)
  library(grid)
  
  if (is.null(dataset_col)) {
    dataset_col <- c("#E64B35CC", "#4DBBD5CC", "#00A087CC", 
                     "#3C5488CC", "#F39B7FCC", "#6BAED6FF", "#FD8D3CFF", 
                     "#FD8D3CFF", "#74C476FF", "#9E9AC8FF",
                     "#969696FF", "#9ECAE1FF", "#FDAE6BFF",
                     "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF",
                     "#C6DBEFFF", "#FDD0A2FF", "#C7E9C0FF",
                     "#DADAEBFF", "#D9D9D9FF")
  }
  
  metamodel <- object
  dt <- metamodel[["data"]]
  
  # === 基础检查与补列 ===
  stopifnot(all(c("HR","LCI","HCI") %in% names(dt)))
  if (!"pvalue" %in% names(dt) && "p" %in% names(dt)) dt$pvalue <- dt$p
  if (!"Group" %in% names(dt) && !"Cohorts" %in% names(dt)) {
    stop("需要有 Group 或 Cohorts 列作为队列名。")
  }
  
  # 权重（字符串百分比）
  dt$Weight_random <- paste0(round(100 * metamodel$w.random / sum(metamodel$w.random), 2), "%")
  dt$Weight_fixed  <- paste0(round(100 * metamodel$w.fixed  / sum(metamodel$w.fixed ), 2), "%")
  
  # === 追加随机/固定效应两行 ===
  dt[nrow(dt)+1, ] <- NA
  dt[nrow(dt), if ("Group" %in% names(dt)) "Group" else "Cohorts"] <- "Random effect model"
  dt[nrow(dt), c("HR","LCI","HCI")] <- exp(c(metamodel$TE.random, metamodel$lower.random, metamodel$upper.random))
  dt[nrow(dt), "pvalue"] <- metamodel$pval.random
  dt[nrow(dt), "Weight_random"] <- "100%"
  dt[nrow(dt), "Weight_fixed"]  <- "--"
  
  dt[nrow(dt)+1, ] <- NA
  dt[nrow(dt), if ("Group" %in% names(dt)) "Group" else "Cohorts"] <- "Fixed effect model"
  dt[nrow(dt), c("HR","LCI","HCI")] <- exp(c(metamodel$TE.fixed, metamodel$lower.fixed, metamodel$upper.fixed))
  dt[nrow(dt), "pvalue"] <- metamodel$pval.fixed
  dt[nrow(dt), "Weight_random"] <- "--"
  dt[nrow(dt), "Weight_fixed"]  <- "100%"
  
  # === 计算 SE（点大小用），以及 TE/SE(TE)（显示用）===
  dt$se <- (log(dt$HCI) - log(dt$HR)) / 1.96
  
  # TE/SE(TE) 优先使用已有列，否则从 HR/SE 推导
  if ("TE" %in% names(dt)) {
    dt$TE <- suppressWarnings(as.numeric(dt$TE))
  } else {
    dt$TE <- log(dt$HR)
  }
  if ("SE(TE)" %in% names(dt)) {
    dt[["SE(TE)"]] <- suppressWarnings(as.numeric(dt[["SE(TE)"]]))
  } else if ("seTE" %in% names(dt)) {
    dt[["SE(TE)"]] <- suppressWarnings(as.numeric(dt[["seTE"]]))
  } else {
    dt[["SE(TE)"]] <- dt$se
  }
  
  # 其他展示列
  show_row <- if ("Group" %in% names(dt)) !is.na(dt$Group) else !is.na(dt$Cohorts)
  dt$` ` <- paste(rep(" ", 20), collapse = " ")
  dt$`HR (95% CI)` <- ifelse(show_row, sprintf("%.2f (%.2f - %.2f)", dt$HR, dt$LCI, dt$HCI), "")
  dt$`P value` <- ifelse(dt$pvalue < 0.001, "P<0.001", sprintf("%.3f", dt$pvalue))
  
  # 统一队列名列为 Cohorts
  if ("Group" %in% names(dt)) names(dt)[names(dt) == "Group"] <- "Cohorts"
  names(dt)[names(dt) == "Weight_random"] <- "Weight(random)"
  names(dt)[names(dt) == "Weight_fixed"]  <- "Weight(fixed)"
  
  # 四舍五入（先数值，后可能置空）
  dt$TE <- round(dt$TE, 2)
  dt[["SE(TE)"]] <- round(dt[["SE(TE)"]], 2)
  
  # 汇总两行不展示 TE/SE(TE)
  dt[(nrow(dt)-1):nrow(dt), c("TE","SE(TE)")] <- NA_real_
  
  # 排序与行名
  rownames(dt) <- dt$Cohorts
  dt <- dt[c(dataset, "Random effect model", "Fixed effect model"), ]
  
  # 主题与绘图
  tm <- forest_theme(
    core = list(bg_params = list(fill = c(dataset_col[seq_along(dataset)], "grey", "grey"), alpha = 0.5)),
    base_size = 10, ci_pch = 16, ci_col = "#762a83", ci_lty = 1, ci_lwd = 1.5, ci_Theight = 0.2,
    refline_lwd = 1, refline_lty = "dashed", refline_col = "grey20",
    vertline_lwd = 1, vertline_lty = "dashed", vertline_col = "grey20",
    summary_fill = "#4575b4", summary_col = "#4575b4",
    footnote_cex = 1, footnote_fontface = "italic", footnote_col = "red"
  )
  
  xmax <- ceiling(max(dt$HCI, na.rm = TRUE))
  ticks <- c(0.5, 2^seq(0, max(0, floor(log2(xmax - 1))), by = 1))
  
  cols_to_show <- c("Cohorts","TE","SE(TE)","Weight(random)","Weight(fixed)"," ","HR (95% CI)","P value")
  
  p <- forestploter::forest(
    dt[, cols_to_show],
    est = dt$HR, lower = dt$LCI, upper = dt$HCI, sizes = dt$se,
    is_summary = c(rep(FALSE, nrow(dt) - 2), TRUE, TRUE),
    ci_column = 6, ref_line = 1, arrow_lab = c("Better","Worse"),
    x_trans = "log2", xlim = c(0, xmax), ticks_at = ticks,
    footnote = " Univariate Cox Regression", theme = tm
  )
  
  p <- add_text(p, text = "Meta analysis of univariate Cox regression",
                part = "header", row = 0, col = 4:6,
                just = "center", gp = gpar(fontface = "bold"))
  p <- add_border(p, part = "header", row = c(0, 1), gp = gpar(lwd = 1))
  p <- insert_text(p, text = "Meta analysis", row = length(dataset) + 1,
                   just = "left", gp = gpar(cex = 0.8, col = "blue", fontface = "italic"))
  print(p)
}
