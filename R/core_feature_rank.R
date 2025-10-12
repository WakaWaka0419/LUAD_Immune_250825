core_feature_rank <- function (object, col = NULL, top = NULL) {
  stopifnot(is.data.frame(object))
  if (!"selected.fea" %in% names(object)) stop("object 中缺少列 'selected.fea'")
  if (is.null(col)) col <- c("#B09C85", "#E18727")
  if (is.null(top)) top <- 50
  
  tmp <- as.data.frame(table(object$selected.fea), stringsAsFactors = FALSE)
  names(tmp) <- c("Var1", "Frequence")
  tmp$Var1 <- gsub("\\.", "-", tmp$Var1)
  
  tmp <- tmp[order(tmp$Frequence, decreasing = FALSE), ]
  top <- min(top, nrow(tmp))
  tmp2 <- utils::tail(tmp, top)
  tmp2$Var1 <- factor(tmp2$Var1, levels = tmp2$Var1)
  
  theme_layer <- if (requireNamespace("DOSE", quietly = TRUE)) DOSE::theme_dose(10) else ggplot2::theme_minimal(base_size = 10)
  
  dataset_col <- c("#B09C85CC","#7E6148CC","#DC0000CC","#91D1C2CC","#8491B4CC","#F39B7FCC",
                   "#3C5488CC","#00A087CC","#4DBBD5CC","#E64B35CC")
  
  p1 <- ggplot2::ggplot(tmp2, ggplot2::aes(x = Var1, y = Frequence, fill = Var1)) +
    ggplot2::geom_bar(position = "dodge", stat = "identity", color = "black", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = dataset_col) +
    theme_layer +
    ggplot2::theme(
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.3),
      plot.title  = ggplot2::element_text(hjust = 0.5),
      panel.background = ggplot2::element_rect(fill = "white")
    ) +
    ggplot2::labs(y = "Frequence of screening", x = "", title = paste("Top", top, "selected genes")) +
    ggplot2::coord_flip()
  
  print(p1)
}
