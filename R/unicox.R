unicox<-function(yT,xVar){
  FML<-as.formula(paste0("BaSurv~",xVar))
  Gcox<-coxph(FML,data=yT)
  Gsum<-summary(Gcox)
  HR<-round(Gsum$coefficients[,2],2)
  Pvalue<-round(Gsum$coefficients[,5],3)
  CI<-data.frame(lower=round(Gsum$conf.int[,3],2),upper=round(Gsum$conf.int[,4],2)) %>% 
    dplyr::mutate(ci=str_c(lower,"-",upper)) %>% 
    .[,"ci"]
  lowCI<-round(Gsum$conf.int[,3],2)
  highCI<-round(Gsum$conf.int[,4],2)
  coef<-round(Gsum$coefficients[,1],2)
  lowCoef<-round(log(Gsum$conf.int[,3]),2)
  highCoef<-round(log(Gsum$conf.int[,4]),2)
  unicox<-data.frame("Characteristics"=xVar,
                     "group"=rownames(Gsum$coefficients),
                     "coef"=coef,
                     "lowCoef"=lowCoef,
                     "highCoef"=highCoef,
                     "hr"=HR,
                     "lowCI"=lowCI,
                     "highCI"=highCI,
                     "CI95"=CI,
                     "pvalue"=Pvalue)
  return(unicox)
}