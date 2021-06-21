#Files at the cluster are under /dcl01/lieber/ajaffe/lab/smokingMouse_Indirects/
#processed-data/01_SPEAQeasy/pipeline_output/count_objects/

#with R 4.0.x

### libraries
library(SummarizedExperiment)
library(here)
library(sessioninfo)
library(recount)
library(jaffelab)
library(ggplot2)


## Load objects
load(here("processed-data","01_SPEAQeasy", "pipeline_output", "count_objects",
          "rse_gene_smoking_mouse_n208.Rdata"), verbose = TRUE)
load(here("processed-data","01_SPEAQeasy", "pipeline_output","count_objects",
          "rse_exon_smoking_mouse_n208.Rdata"), verbose = TRUE)
load(here("processed-data","01_SPEAQeasy", "pipeline_output", "count_objects",
          "rse_jx_smoking_mouse_n208.Rdata"), verbose = TRUE)
load(here("processed-data","01_SPEAQeasy", "pipeline_output","count_objects",
          "rse_tx_smoking_mouse_n208.Rdata"), verbose = TRUE)

#Load table of phenotypes
pheno <- read.delim(here("raw-data","Maternal_Smoking_pheno.txt"))


#Rename rownames
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FC1P1-blood"] = "Sample_FC1P1_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FC1P2-blood"] = "Sample_FC1P2_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FC1P3-blood"] = "Sample_FC1P3_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FC2P1-blood"] = "Sample_FC2P1_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FC2P2-blood"] = "Sample_FC2P2_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FC2P3-blood"] = "Sample_FC2P3_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FC2P4-blood"] = "Sample_FC2P4_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FC31-blood"] = "Sample_FC31_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FC32-blood"] = "Sample_FC32_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FC41-blood"] = "Sample_FC41_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FC42-blood"] = "Sample_FC42_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FC43-blood"] = "Sample_FC43_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FE1P1-blood"] = "Sample_FE1P1_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FE1P2-blood"] = "Sample_FE1P2_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FE1P3-blood"] = "Sample_FE1P3_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FE2P1-blood"] = "Sample_FE2P1_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FE2P2-blood"] = "Sample_FE2P2_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FE3P1-blood"] = "Sample_FE3P1_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FE3P2-blood"] = "Sample_FE3P2_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FE3P3-blood"] = "Sample_FE3P3_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FE41-blood"] = "Sample_FE41_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FE42-blood"] = "Sample_FE42_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FE43-blood"] = "Sample_FE43_blood"
pheno$SAMPLE_ID[pheno$SAMPLE_ID == "FE44-blood"] = "Sample_FE44_blood"

rownames(pheno) <- pheno$SAMPLE_ID

#Look at differences
setdiff(rownames(colData(rse_gene)), rownames(pheno))
setdiff(rownames(pheno), rownames(colData(rse_gene)))

## Merge colData with phenotypes table
colData(rse_gene) <- merge(colData(rse_gene), pheno, by="SAMPLE_ID")
colData(rse_exon) <- merge(colData(rse_exon), pheno, by="SAMPLE_ID")
colData(rse_jx) <- merge(colData(rse_jx), pheno, by="SAMPLE_ID")
colData(rse_tx) <- merge(colData(rse_tx), pheno, by="SAMPLE_ID")

#set length to 100 Jx
rowData(rse_jx)$Length <- 100

#Build a matrix with the expression values in RPKM format 
exprs <- list(
  "Gene" = recount::getRPKM(rse_gene, length_var="Length"),
  "Exon" = recount::getRPKM(rse_exon, length_var="Length"),
  "Jx" = recount::getRPKM(rse_jx, length_var= "Length"),
  "Tx" = assays(rse_tx)$tpm 
)



## Identify potential cutoffs
seed <- 20191217
seeds <- seed + 0:3
names(seeds) <- names(exprs)
cutoffs <- sapply(names(exprs), function(type) {
  message(type)
  pdf(paste0('suggested_expr_cutoffs_', tolower(type), '.pdf'), width = 12)
  cuts <- jaffelab::expression_cutoff(exprs[[type]], seed = seeds[type])
  message(paste(cuts, collapse = ' '))
  cut <- max(cuts)
  dev.off()
  return(cut)
})

# Gene
# 2021-06-02 16:53:04 the suggested expression cutoff is 0.16
# 0.18 0.14
# Exon
# 2021-06-02 16:53:27 the suggested expression cutoff is 0.2
# 0.23 0.17
# Jx
# 2021-06-02 16:53:54 the suggested expression cutoff is 0.16
# 0.14 0.17
# Tx
# 2021-06-02 16:53:59 the suggested expression cutoff is 0.23
# 0.25 0.22

cutoffs
# Gene Exon   Jx   Tx 
# 0.18 0.23 0.17 0.25 


### Filter RSEs
means <- lapply(exprs, rowMeans)

rowRanges(rse_gene)$meanExprs <- means[['Gene']]
rowRanges(rse_gene)$passExprsCut <- means[['Gene']] > cutoffs['Gene']
rse_gene <- rse_gene[rowRanges(rse_gene)$passExprsCut]
save(rse_gene, file = 'rse_gene.Rdata')

rowRanges(rse_exon)$meanExprs <- means[['Exon']]
rowRanges(rse_exon)$passExprsCut <- means[['Exon']] > cutoffs['Exon']
rse_exon <- rse_exon[rowRanges(rse_exon)$passExprsCut]
save(rse_exon, file = 'rse_exon.Rdata')

rowRanges(rse_jx)$meanExprs <- means[['Jx']]
rowRanges(rse_jx)$passExprsCut <- means[['Jx']] > cutoffs['Jx'] 
save(rse_jx, file = 'rse_jx.Rdata')

rowRanges(rse_tx)$meanExprs <- means[['Tx']]
rowRanges(rse_tx)$passExprsCut <- means[['Tx']] > cutoffs['Tx']
rse_tx <- rse_tx[rowRanges(rse_tx)$passExprsCut]
save(rse_tx, file = 'rse_tx.Rdata')


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# > session_info()
# oconductor  
# matrixStats          * 0.58.0   2021-01-29 [2] CRAN (R 4.0.3)
# RCurl                  1.98-1.3 2021-03-16 [2] CRAN (R 4.0.4)
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.0.3)
# S4Vectors            * 0.28.1   2020-12-09 [2] Bioconductor  
# sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.0.3)
# SummarizedExperiment * 1.20.0   2020-10-27 [2] Bioconductor  
# withr                  2.4.2    2021-04-18 [2] CRAN (R 4.0.4)
# XVector                0.30.0   2020-10-27 [2] Bioconductor  
# zlibbioc               1.36.0   2020-10-27 [2] Bioconductor  
# 
# [1] /users/bpardo/R/4.0.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/library
