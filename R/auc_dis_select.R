auc_dis_select <- function (object, model_name, dataset_col = NULL, dataset, order = NULL, 
          year) 
{
  library(ggplot2)
  if (is.null(dataset_col) == T) {
    dataset_col <- c("#E64B35CC", "#4DBBD5CC", "#00A087CC", 
                     "#3C5488CC", "#F39B7FCC", "#6BAED6FF", "#FD8D3CFF", 
                     "#74C476FF", "#9E9AC8FF", "#969696FF", "#9ECAE1FF", 
                     "#FDAE6BFF", "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF", 
                     "#C6DBEFFF", "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", 
                     "#D9D9D9FF")
  }
  else {
    dataset_col <- dataset_col
  }
  auc_select <- data.frame()
  for (j in 1:length(object)) {
    object2 <- object[[j]]
    auc_all <- data.frame()
    for (i in names(object2)) {
      tmp <- object2[[i]]
      tmp2 <- data.frame()
      for (n in dataset) {
        tmp2 <- rbind(tmp2, data.frame(AUC = unique(tmp[[n]][["AUC"]]), 
                                       HR = unique(tmp[[n]][["HR"]])))
      }
      tmp2$ID <- dataset
      tmp2$Model <- i
      auc_all <- rbind(auc_all, tmp2)
    }
    auc_all$HR <- ifelse(auc_all$HR > 1, "HR>1", "HR<1")
    auc_all$AUC <- as.numeric(sprintf("%.2f", auc_all$AUC))
    auc_all_select <- auc_all[auc_all$Model == model_name, 
    ]
    auc_all_select$Year <- paste(year[j], "Year", sep = "-")
    auc_select <- rbind(auc_select, auc_all_select)
  }
  if (is.null(order) == T) {
    order <- dataset
    auc_select$ID <- factor(auc_select$ID, levels = order)
  }
  else {
    order <- order
    auc_select$ID <- factor(auc_select$ID, levels = order)
  }
  mytheme <- theme_bw() +
    theme(axis.text.x = element_text(hjust = 0.5,size = 20), 
          axis.text.y = element_text(size = 20), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 20), 
          plot.title = element_text(hjust = 0.5,size =  22),
          legend.position = "top",
          legend.title = element_blank(),
          legend.background = element_rect(fill = 'transparent'),
          legend.text = element_text(size = 20))
  p1 <- ggplot(data = auc_select, aes(x = ID, y = AUC, fill = ID)) + 
    geom_bar(position = "dodge", stat = "identity") + scale_fill_manual(values = dataset_col, 
                                                                        name = "Cohort") + geom_text(aes(label = AUC), position = position_dodge(width = 0.9), 
                                                                                                     vjust = 0.5, hjust = 1.2, size = 3) + theme(panel.grid = element_blank(), 
                                                                                                                                                 panel.border = element_rect(colour = "black", fill = NA, 
                                                                                                                                                                             size = 0.3), legend.position = "", plot.title = element_text(hjust = 0.5), 
                                                                                                                                                 panel.background = element_rect(fill = "white")) + labs(y = "AUC", 
                                                                                                                                                                                                         x = "", title = model_name) + coord_flip() + facet_grid(. ~ 
                                                                                                                                                                                                                                                                   Year) 
  print(p1)
}
