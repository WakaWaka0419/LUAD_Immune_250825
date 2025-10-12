
rs_sur <- function (object, model_name, dataset, cutoff = NULL, conf.int = F, 
          median.line = NULL, color = NULL, xlab = NULL, pval.coord = NULL) 
{
  library(survival)
  library(survminer)
  if (is.null(color) == T) {
    color <- c("#2EC4B6","#FFA5AB")
  }
  else {
    color <- color
  }
  if (is.null(median.line) == T) {
    median.line <- "none"
  }
  else {
    median.line <- median.line
  }
  if (is.null(xlab) == T) {
    xlab <- "Time"
  }
  else {
    xlab <- xlab
  }
  if (is.null(pval.coord) == T) {
    pval.coord <- c(1, 0.25)
  }
  else {
    pval.coord <- pval.coord
  }
  tmp <- object[["riskscore"]][[model_name]][[dataset]]
  mySurv <- Surv(tmp$OS.time, tmp$OS)
  if (is.null(cutoff) == T) {
    cutoff <- 0.5
    value <- quantile(tmp$RS, probs = c(cutoff))
  }
  else {
    if (cutoff == "mean") {
      value <- mean(tmp$RS)
    }
    else {
      cutoff <- cutoff
      value <- quantile(tmp$RS, probs = c(cutoff))
    }
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
  tmp$Group <- ifelse(tmp$RS > value, "High", "Low")
  Group <- tmp$Group
  Group <- factor(Group, levels = c("Low", "High"))
  tmp$Group <- factor(Group, levels = c("Low", "High"))
  fit <- survfit(Surv(OS.time, OS) ~ Group, data = tmp)
  data.survdiff <- survdiff(mySurv ~ Group)
  p.val <- 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 
                        1)
  HR <- (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 <- exp(log(HR) + qnorm(0.975) * sqrt(1/data.survdiff$exp[2] + 
                                              1/data.survdiff$exp[1]))
  low95 <- exp(log(HR) - qnorm(0.975) * sqrt(1/data.survdiff$exp[2] + 
                                               1/data.survdiff$exp[1]))
  HR <- paste("Hazard Ratio = ", round(HR, 2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95, 2), round(up95, 
                                                       2), sep = " - "), sep = "")
  p1 <- ggsurvplot(fit, data = tmp,  censor = F, 
                   palette = color, legend.title = paste0("Riskscore in ", 
                                                          dataset), font.legend = 11, surv.median.line = median.line, 
                   legend.labs = c(paste0("Low", "(n=", fit$n[1], ")"), 
                                   paste0("High", "(n=", fit$n[2], ")")), pval = paste(pval = ifelse(p.val < 
                                                                                                       0.001, "p < 0.001", paste("p = ", round(p.val, 3), 
                                                                                                                                 sep = "")), HR, CI, sep = "\n"), xlab = xlab, title = model_name, 
                   pval.coord = pval.coord)
  p1 <- p1$plot + theme(plot.title = element_text(hjust = 0.5)) + mytheme
  print(p1)
}
