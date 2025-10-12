#####Loading Packages#####
library(tinyarray)
library(GEOquery)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
setwd("./Input/GEO/Orign_data")

#####Main Scripts#####
###01GSE78220###
GSE78220 <- geo_download("GSE78220")
GSE78220_Clinical <- GSE78220$pd %>%
  dplyr::select(c("title", "description","patient id:ch1","overall survival (days):ch1","anatomical location:ch1",
                  "age (yrs):ch1","nf1:ch1","braf:ch1", "nras:ch1","disease status:ch1","source_name_ch1","anti-pd-1 response:ch1","biopsy time:ch1", 
                  "gender:ch1", "previous mapki:ch1", "stranded/unstranded rnaseq:ch1", 
                  "study site:ch1", "vital status:ch1")) %>%
  dplyr::rename("OS.time" = "overall survival (days):ch1") %>%
  dplyr::rename("OS" = "vital status:ch1")
GSE78220_exp <- readxl::read_xlsx("./GSE78220_PatientFPKM.xlsx") %>%
  column_to_rownames(var = "Gene") %>%
  log1p(.)
save(GSE78220_exp,GSE78220_Clinical,file = "../Curated_data/GSE78220_Curated.Rdata")

###02GSE126004###
GSE126044_exp <- data.table::fread("./GSE126044_counts.txt.gz") %>% 
  aggregate(.~V1,FUN = mean) %>%
  column_to_rownames(var = "V1") %>%
  log1p()
GSE126044_Clinical <- readxl::read_xlsx("./GSE126044_Clinical.xlsx") %>%
 dplyr::filter(ID %in% colnames(GSE126044_exp))
save(GSE126044_exp,GSE126044_Clinical,file = "../Curated_data/GSE126044_Curated.Rdata")

###03GSE135225###
GSE135222_exp <- data.table::fread("./GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz") %>%
    mutate(gene_id = sub("\\..*", "", gene_id))
gene_df <- bitr(GSE135222_exp$gene_id,
                fromType = "ENSEMBL",
                toType   = "SYMBOL",
                OrgDb    = org.Hs.eg.db)
GSE135222_exp <- GSE135222_exp %>%
  left_join(gene_df, by = c("gene_id" = "ENSEMBL")) %>%
  dplyr::select(-gene_id) %>%
  aggregate(.~SYMBOL,FUN = mean) %>%
  column_to_rownames(var = "SYMBOL") %>%
  log1p()

GSE135222_Clinical <- readxl::read_xls("./GSE135222.Clinical.xls") %>%
  dplyr::mutate(`Sample ID` = paste0("NSCLC",`Sample ID`)) %>%
  dplyr::filter(`Sample ID` %in% colnames(GSE135222_exp)) %>%
  dplyr::rename("PFS.time" = "PFS") %>%
  dplyr::rename("PFS" = "PD_Event.1_Censoring.0")
save(GSE135222_exp,GSE135222_Clinical,file = "../Curated_data/GSE135222_Curated.Rdata")


###04GSE207422###
GSE207422_exp <- data.table::fread("./GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt") %>%
  column_to_rownames(var = "Gene")
GSE207422_Clinical <- openxlsx::read.xlsx("./GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx") 
save(GSE207422_exp,GSE207422_Clinical,file = "../Curated_data/GSE207422_Curated.Rdata")


###05GSE157010###
GSE157010_exp <- gse157010.expr
GSE157010_Clinical <- gse157010.cli
save(GSE157010_exp,GSE157010_Clinical,file = "../GSE157010_Cruated.Rdata")

###IMvigor210###
load("../Curated_data/IMvigor210CoreBiologies.Rdata")
IMvigor210_exp <- expreSet %>%
  rownames_to_column(var = "entrez_id") %>%  
  left_join(annoData[, c("entrez_id", "symbol")],
            by = "entrez_id") %>%
  dplyr::select(-"entrez_id") %>%
  aggregate(.~symbol,FUN = mean) %>%
  column_to_rownames(var = "symbol")
IMvigor210_Clinical <- phenoData
save(IMvigor210_exp,IMvigor210_Clinical,file = "../Curated_data/IMvigor210_Cruated.Rdata")
