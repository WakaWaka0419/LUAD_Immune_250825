#####Loading Packages#####
library(tinyarray)
library(GEOquery)
library(tibble)
setwd("./Input/GEO/Orign_data")

#####Main Scripts#####
###01GSE13213###
GSE13213 <- geo_download("GSE13213")
GSE13213_Clinical <- GSE13213$pd %>%
  dplyr::select(c("Annotation:ch1","Survival (days):ch1","Smoking (BI):ch1","Age:ch1","Status:ch1",
                  "Site of relapse:ch1","TNM (Pathological):ch1","Stage (Pathological ):ch1",
                  "Evidence of relapse:ch1","p53 Status:ch1","description","Cohort:ch1","EGFR status:ch1","K-ras Status:ch1"
                  ,"Sex:ch1")) %>%
  dplyr::rename("OS.time" = "Survival (days):ch1") %>%
  dplyr::rename("OS" = "Status:ch1") %>%
  dplyr::mutate(OS = case_when(
    OS == "Alive" ~ 0,
    OS == "Dead" ~ 1
  ))
GSE13213_probe <- getGEO(GSE13213$gpl)
GSE13213_probe <- GSE13213_probe@dataTable@table[,c(1,7)]
GSE13213_exp   <- as.data.frame(GSE13213$exp) %>%
  rownames_to_column(var = "ID") %>%
  left_join(GSE13213_probe,by= "ID") %>%
  dplyr::filter(GENE_SYMBOL != "" & !is.na(GENE_SYMBOL)) %>%
  dplyr::select(-ID) %>%
  aggregate(.~GENE_SYMBOL,FUN = mean) %>%
  column_to_rownames(var = "GENE_SYMBOL")
save(GSE13213_exp,GSE13213_Clinical,file = "../Curated_data/GSE13213_Curated.Rdata")

###02GSE26939###
GSE26939 <- geo_download("GSE26939")
GSE26939_Clinical <- GSE26939$pd %>%
  dplyr::select(c("title", "description", "survival_months:ch1","survival_status(0='alive',1='dead'):ch1","age (90=greater than or equal to 90):ch1", "pack_years:ch1",
                  "fibrosis_percent:ch1","tumor_percent:ch1","Stage:ch1","stage:ch1","grade:ch1","egfr (0='wt',1='mutated'):ch1", "kras (0='wt',1='mutated'):ch1", 
                  "sex:ch1", "Sex:ch1", "smoking_status(0=nonsmoker,1=smoker):ch1", 
                  "stk11 (0='wt',1='mutated'):ch1", "subtype:ch1", 
                  "tp53 (0='wt',1='mutated'):ch1", "adenosquamous_content (0='no', 1='yes'):ch1", 
                  "bronchioloalveolar_content (0='no', 1='yes'):ch1", "marked_lymphocytes:ch1")) %>%
  dplyr::rename("OS.time" = "survival_months:ch1") %>%
  dplyr::rename("OS" = "survival_status(0='alive',1='dead'):ch1") %>%
  dplyr::mutate(OS.time = as.integer(as.numeric(OS.time)*30))
GSE26939_probe <- getGEO(GSE26939$gpl)
GSE26939_probe <- GSE26939_probe@dataTable@table[,c(1,3)] %>%
  dplyr::mutate(ID = as.character(ID))
GSE26939_exp   <- as.data.frame(GSE26939$exp) %>%
  rownames_to_column(var = "ID") %>%
  left_join(GSE26939_probe,by= "ID") %>%
  dplyr::filter(ORF != "" & !is.na(ORF)) %>%
  dplyr::select(-ID) %>%
  aggregate(.~ORF,FUN = mean) %>%
  column_to_rownames(var = "ORF")
save(GSE26939_exp,GSE26939_Clinical,file = "../Curated_data/GSE26939_Curated.Rdata")

###03GSE29019###
GSE29016 <- geo_download("GSE29016")
GSE29016_Clinical <- GSE29016$pd %>%
  dplyr::select(c("title","assay:ch1","os:ch1","age:ch1","tnm:ch1", "source_name_ch1",
                  "Stage:ch1","histology:ch1", "characteristics_ch1.2","ckit_ihc:ch1", 
                  "egfr_amp:ch1", "egfr_ihc:ch1", "egfr_mut:ch1", "gexcluster_ac:ch1", 
                  "her2_ihc:ch1", "kras_mut:ch1", "osbin:ch1", "pakt_ihc:ch1", 
                  "pik3ca_mut:ch1", "pten_ihc:ch1", "smoking:ch1", "characteristics_ch1.5", 
                  "Sex:ch1")) %>%
  dplyr::rename("OS.time" = "os:ch1") %>%
  dplyr::rename("OS" = "characteristics_ch1.2") %>%
  dplyr::filter(source_name_ch1 %in% c("Lung tumor, adenocarcinoma, never smoker","Lung tumor, adenocarcinoma","Lung tumor, adenocarcinoma, smoker")) %>%
  dplyr::mutate(OS.time = as.integer(as.numeric(OS.time)*365)) %>%
  dplyr::mutate(OS = substr(OS,8,8))

GSE29016_probe <- getGEO(GSE29016$gpl)
GSE29016_probe <- GSE29016_probe@dataTable@table[,c(1,14)] 

GSE29016_exp   <- as.data.frame(GSE29016$exp) %>%
  rownames_to_column(var = "ID") %>%
  left_join(GSE29016_probe,by= "ID") %>%
  dplyr::filter(Symbol != "" & !is.na(Symbol)) %>%
  dplyr::select(-ID) %>%
  aggregate(.~Symbol,FUN = mean) %>%
  column_to_rownames(var = "Symbol") %>%
  dplyr::select(rownames(GSE29016_Clinical)) %>%
  log1p(.)
save(GSE29016_exp,GSE29016_Clinical,file = "../Curated_data/GSE29016_Curated.Rdata")


###04GSE30219###
GSE30219 <- geo_download("GSE30219")
GSE30219_Clinical <- GSE30219$pd %>%
  dplyr::select(c("title", "description","follow-up time (months):ch1","disease free survival in months:ch1",
                  "age at surgery:ch1","histology:ch1","pn stage:ch1", "pt stage:ch1","pm stage:ch1", "relapse (event=1; no event=0):ch1",
                  "gender:ch1", "status:ch1", "status", 
                  "submission_date", "last_update_date", "source_name_ch1","data_processing", "tissue:ch1")) %>%
  dplyr::rename("OS.time" = "follow-up time (months):ch1") %>%
  dplyr::rename("DFS.time" = "disease free survival in months:ch1") %>%
  dplyr::rename("OS" = "status:ch1") %>%
  dplyr::filter(`histology:ch1` == "ADC") %>%
  dplyr::mutate(OS.time = as.integer(as.numeric(OS.time)*30)) %>%
  dplyr::mutate(OS = case_when(
    OS == "DEAD" ~ 1,
    OS == "ALIVE" ~ 0
  ))

GSE30219_probe <- getGEO(GSE30219$gpl)
GSE30219_probe <- GSE30219_probe@dataTable@table[,c(1,11)] 

GSE30219_exp   <- as.data.frame(GSE30219$exp) %>%
  rownames_to_column(var = "ID") %>%
  left_join(GSE30219_probe,by= "ID") %>%
  dplyr::filter(`Gene Symbol` != "" & !is.na(`Gene Symbol`)) %>%
  dplyr::select(-ID) %>%
  aggregate(.~`Gene Symbol`,FUN = mean) %>%
  column_to_rownames(var = "Gene Symbol") %>%
  dplyr::select(rownames(GSE30219_Clinical)) 
save(GSE30219_exp,GSE30219_Clinical,file = "../Curated_data/GSE30219_Curated.Rdata")


###05GSE31210###
GSE31210 <- geo_download("GSE31210")
GSE31210_Clinical <- GSE31210$pd %>%
  dplyr::select(c("title", "months before relapse/censor:ch1", "days before relapse/censor:ch1", 
                  "days before death/censor:ch1","myc_copy:ch1","age (years):ch1", "source_name_ch1", 
                  "age:ch1","gene alteration status:ch1","myc:ch1", "pathological stage:ch1","description", 
                  "gender:ch1", "smoking status:ch1", "tissue:ch1","cluster:ch1", "death:ch1", "exclude for prognosis analysis due to incomplete resection or adjuvant therapy:ch1", 
                  "group:ch1", "pstage iorii:ch1", "relapse:ch1", "characteristics_ch1")) %>%
  dplyr::rename("OS.time" = "days before death/censor:ch1" ) %>%
  dplyr::rename("OS" = "death:ch1") %>%
  dplyr::filter(description == "Gene expression data from primary lung ADC") %>%
  dplyr::mutate(OS = case_when(
    OS == "dead" ~ 1,
    OS == "alive" ~ 0
  ))

GSE31210_probe <- getGEO(GSE31210$gpl)
GSE31210_probe <- GSE31210_probe@dataTable@table[,c(1,11)] 

GSE31210_exp   <- as.data.frame(GSE31210$exp) %>%
  rownames_to_column(var = "ID") %>%
  left_join(GSE31210_probe,by= "ID") %>%
  dplyr::filter(`Gene Symbol` != "" & !is.na(`Gene Symbol`)) %>%
  dplyr::select(-ID) %>%
  aggregate(.~`Gene Symbol`,FUN = mean) %>%
  column_to_rownames(var = "Gene Symbol") %>%
  dplyr::select(rownames(GSE31210_Clinical)) %>%
  log1p(.)
save(GSE31210_exp,GSE31210_Clinical,file = "../Curated_data/GSE31210_Curated.Rdata")


###06GSE42127###
GSE42127 <- geo_download("GSE42127")
GSE42127_Clinical <- GSE42127$pd %>%
  dplyr::select(c("title", "description","age at surgery:ch1","overall survival months:ch1","final.pat.stage:ch1",
                  "gender:ch1", "had_adjuvant_chemo:ch1", "histology:ch1", "survival status:ch1")) %>%
  dplyr::rename("OS.time" = "overall survival months:ch1" ) %>%
  dplyr::rename("OS" = "survival status:ch1") %>%
  dplyr::mutate(OS.time = as.integer(as.numeric(OS.time)*30)) %>%
  dplyr::filter(`histology:ch1` == "Adenocarcinoma") %>%
  dplyr::mutate(OS = case_when(
    OS == "D" ~ 1,
    OS == "A" ~ 0
  ))

GSE42127_probe <- getGEO(GSE42127$gpl)
GSE42127_probe <- GSE42127_probe@dataTable@table[,c(1,13)] 

GSE42127_exp   <- as.data.frame(GSE42127$exp) %>%
  rownames_to_column(var = "ID") %>%
  left_join(GSE42127_probe,by= "ID") %>%
  dplyr::filter("Symbol" != "" & !is.na(Symbol)) %>%
  dplyr::select(-ID) %>%
  aggregate(.~Symbol,FUN = mean) %>%
  column_to_rownames(var = "Symbol") %>%
  dplyr::select(rownames(GSE42127_Clinical)) 
save(GSE42127_exp,GSE42127_Clinical,file = "../Curated_data/GSE42127_Curated.Rdata")

###07GSE50081###
GSE50081 <- geo_download("GSE50081")
GSE50081_Clinical <- GSE50081$pd %>%
  dplyr::select(c("title",  "disease-free survival time:ch1", 
                   "age:ch1","survival time:ch1",  "histology:ch1", "smoking:ch1", 
                  "Stage:ch1","recurrence:ch1", "t-stage:ch1", "n-stage:ch1", "Sex:ch1", "status:ch1")) %>%
  dplyr::rename("OS.time" = "survival time:ch1" ) %>%
  dplyr::rename("OS" = "status:ch1") %>%
  dplyr::mutate(OS.time = as.integer(as.numeric(OS.time)*365)) %>%
  dplyr::filter(`histology:ch1` == "adenocarcinoma") %>%
  dplyr::mutate(OS = case_when(
    OS == "dead" ~ 1,
    OS == "alive" ~ 0
  ))

GSE50081_probe <- getGEO(GSE50081$gpl)
GSE50081_probe <- GSE50081_probe@dataTable@table[,c(1,11)] 

GSE50081_exp   <- as.data.frame(GSE50081$exp) %>%
  rownames_to_column(var = "ID") %>%
  left_join(GSE50081_probe,by= "ID") %>%
  dplyr::filter("Gene Symbol" != "" & !is.na(`Gene Symbol`)) %>%
  dplyr::select(-ID) %>%
  aggregate(.~`Gene Symbol`,FUN = mean) %>%
  column_to_rownames(var = "Gene Symbol") %>%
  dplyr::select(rownames(GSE50081_Clinical)) 
save(GSE50081_exp,GSE50081_Clinical,file = "../Curated_data/GSE50081_Curated.Rdata")

###05GSE157010###
GSE157010_exp <- gse157010.expr
GSE157010_Clinical <- gse157010.cli
save(GSE157010_exp,GSE157010_Clinical,file = "../GSE157010_Cruated.Rdata")
