#Files at the cluster are under /dcl01/lieber/ajaffe/lab/smokingMouse_Indirects/
#processed-data/01_SPEAQeasy/pipeline_output/count_objects/

## Libraries
library(SummarizedExperiment)
library(recount)
library(edgeR)
library(jaffelab)
library(here)
library(sessioninfo)

# 1. Data load, exploration and preparation
## Load objects
load("rse_exon_smoking_mouse_n208.Rdata")
load("rse_gene_smoking_mouse_n208.Rdata")
load("rse_jx_smoking_mouse_n208.Rdata")
load("rse_tx_smoking_mouse_n208.Rdata")

## Samples' info
pheno<-read.table(here("raw-data","Maternal_Smoking_pheno.txt"))

## Data dimensions
dim(rse_gene)
# 55401   208
dim(rse_exon)
# 447670    208
dim(rse_jx)
# 1436068     208
dim(rse_tx)
# 142604    208

## Verify data integrity
which(rowSums(is.na(assay(rse_gene)) | assay(rse_gene) == "") > 0)
# named integer(0)
which(rowSums(is.na(assay(rse_exon)) | assay(rse_exon) == "") > 0)
# named integer(0)
which(rowSums(is.na(assay(rse_tx)) | assay(rse_tx) == "") > 0)
# named integer(0)
which(rowSums(is.na(assay(rse_jx)) | assay(rse_jx) == "") > 0)
# named integer(0)

## Correct data
## Samples' ID starting with "F"
unmod_rows<-which(substr(pheno[,1],1,1)=="F")
mod_row_names<-sapply(unmod_rows, function(x){pheno[x,1]<-paste("Sample_",pheno[x,1], sep="")
                                gsub("-","_",pheno[x,1])})
pheno[unmod_rows,1]<-mod_row_names
## Row and Col names
colnames<-pheno[1,-1]
rownames<-pheno[-1,1]
pheno<-pheno[-1,-1]
colnames(pheno)<-colnames
rownames(pheno)<-rownames

## Look at differences
setdiff(rownames(colData(rse_gene)), rownames(pheno))
setdiff(rownames(pheno), rownames(colData(rse_gene)))

## All samples' info in colData of RSE objects
pheno$SAMPLE_ID<-rownames
colData(rse_gene)<-merge(colData(rse_gene), pheno, by="SAMPLE_ID", all=TRUE)
colData(rse_exon)<-merge(colData(rse_exon), pheno, by="SAMPLE_ID", all=TRUE)
colData(rse_jx)<-merge(colData(rse_jx), pheno, by="SAMPLE_ID", all=TRUE)
colData(rse_tx)<-merge(colData(rse_tx), pheno, by="SAMPLE_ID", all=TRUE)




# 2. Data transformation
## 2.1 Data normalization

## Proportion of zeros
length(which(assay(rse_gene)==0))*100/(dim(rse_gene)[1]*dim(rse_gene)[2])
# 56.53367
length(which(assay(rse_exon)==0))*100/(dim(rse_exon)[1]*dim(rse_exon)[2])
# 30.17473
length(which(assay(rse_jx)==0))*100/(dim(rse_jx)[1]*dim(rse_jx)[2])
# 82.34045

## Transform read counts to CPM
norm_counts <- list(
  "Gene" = cpm(calcNormFactors(rse_gene, method = "TMM")),
  "Exon" = cpm(calcNormFactors(rse_exon, method = "TMM")),
   ## TMMwsp method for >80% of zeros
  "Jx" = cpm(calcNormFactors(rse_jx, method = "TMMwsp")),
   ## Scale TPM (Transcripts per million)
  "Tx" = log2(assays(rse_tx)$tpm + 0.5)
)

## 2.2 Data filtering
## Filter genes, exons and jxn with k CPM in at least n samples
length(which(filterByExpr(norm_counts[["Gene"]])))
# 8876
gene_norm_counts<-norm_counts[["Gene"]][which(filterByExpr(norm_counts[["Gene"]])),]
save(gene_norm_counts, file = 'norm_gene.Rdata')

length(which(filterByExpr(norm_counts[["Exon"]])))
# 10538
exon_norm_counts<-norm_counts[["Exon"]][which(filterByExpr(norm_counts[["Exon"]])),]
save(exon_norm_counts, file = 'norm_exon.Rdata')

length(which(filterByExpr(norm_counts[["Jx"]])))
# 11535
jx_norm_counts<-norm_counts[["Jx"]][which(filterByExpr(norm_counts[["Jx"]])),]
save(jx_norm_counts, file = 'norm_jx.Rdata')

## Filter Transcipts
## Identify potential cutoffs
seed <- 20191217
expression_cutoff(norm_counts[["Tx"]], seed = seed, k=2)
# 2022-06-22 13:06:47 the suggested expression cutoff is 0.6
# percent_features_cut  samples_nonzero_cut 
#                 0.63                 0.56 
cutoff<-0.6
length(which(rowMeans(norm_counts[["Tx"]]) > cutoff))
# 32514
tx_norm_counts<-norm_counts[["Tx"]][rowMeans(norm_counts[["Tx"]]) > cutoff,]
save(tx_norm_counts, file = 'norm_tx.Rdata')


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

#  version  R version 4.2.0 (2022-04-22 ucrt)
#  os       Windows 10 x64 (build 19044)
#  system   x86_64, mingw32
#  ui       RStudio
#  language (EN)
#  collate  Spanish_Mexico.utf8
#  ctype    Spanish_Mexico.utf8
#  tz       America/Mexico_City
#  date     2022-06-22
#  rstudio  2022.02.3+492 Prairie Trillium (desktop)
#  pandoc   2.17.1.1 @ C:/Program Files/RStudio/bin/quarto/bin/ (via rmarkdown)
