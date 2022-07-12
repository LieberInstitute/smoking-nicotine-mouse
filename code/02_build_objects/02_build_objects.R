#Files at the cluster are under /dcl01/lieber/ajaffe/lab/smokingMouse_Indirects/
#processed-data/01_SPEAQeasy/pipeline_output/count_objects/


# 1. Build objects
## Load required libraries
library(SummarizedExperiment)
library(recount)
library(edgeR)
library(jaffelab)
library(here)
library(WGCNA)
library(scater)
library(sessioninfo)


## 1.1 Load data, data exploration and preparation
## Load RSE objects
load(here("raw-data/rse_exon_smoking_mouse_n208.Rdata"))
load(here("raw-data/rse_gene_smoking_mouse_n208.Rdata"))
load(here("raw-data/rse_jx_smoking_mouse_n208.Rdata"))
load(here("raw-data/rse_tx_smoking_mouse_n208.Rdata"))
## Samples' info
pheno<-read.table(here("raw-data/Maternal_Smoking_pheno.txt"))
## Flowcell info
flowcells<-read.table(here("processed-data/02_build_objects/flowcell_info.tsv"))

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
which(rowSums(is.na(assay(rse_exon)) | assay(rse_exon) == "") > 0)
which(rowSums(is.na(assay(rse_tx)) | assay(rse_tx) == "") > 0)
which(rowSums(is.na(assay(rse_jx)) | assay(rse_jx) == "") > 0)

## Correct data
tail(pheno[,1])
## Samples' ID starting with "F"
unmod_rows<-which(substr(pheno[,1],1,1)=="F")
mod_row_names<-sapply(unmod_rows, function(x){pheno[x,1]<-paste("Sample_",pheno[x,1], sep="")
                                gsub("-","_",pheno[x,1])})
pheno[unmod_rows,1]<-mod_row_names
## Check
tail(pheno[,1])

## Row and Col names
colnames<-pheno[1,]
pheno<-pheno[-1,]
colnames(pheno)<-colnames

## Look at differences
setdiff(pheno$SAMPLE_ID, rse_gene$SAMPLE_ID)
setdiff(rse_gene$SAMPLE_ID, pheno$SAMPLE_ID)

## Pregnancy qualitative variable
pheno$Pregnancy<-multiGSub(unique(sort(pheno$Pregnant)), c("No", "Yes"), x = pheno$Pregnant)
## Add flowcell info
flowcells<-flowcells[-1,]
colnames(flowcells)<-c("SAMPLE_ID", "flowcell")
pheno<-merge(pheno, flowcells, by="SAMPLE_ID")

## All samples' info in colData of RSE objects
colData(rse_gene)<-merge(colData(rse_gene), pheno, by="SAMPLE_ID")
colData(rse_exon)<-merge(colData(rse_exon), pheno, by="SAMPLE_ID")
colData(rse_jx)<-merge(colData(rse_jx), pheno, by="SAMPLE_ID")
colData(rse_tx)<-merge(colData(rse_tx), pheno, by="SAMPLE_ID")





## 1.2 Data transformation
## Data normalization

## Proportion of zeros
length(which(assay(rse_gene)==0))*100/(dim(rse_gene)[1]*dim(rse_gene)[2])
# 56.53367
length(which(assay(rse_exon)==0))*100/(dim(rse_exon)[1]*dim(rse_exon)[2])
# 30.17473
length(which(assay(rse_jx)==0))*100/(dim(rse_jx)[1]*dim(rse_jx)[2])
# 82.34045

## Transform read counts to CPM
## With log2(CPM + 0.5) they are closer to a Norm distribution
## Genes
assays(rse_gene, withDimnames=FALSE)$norm_counts <- edgeR::cpm(calcNormFactors(rse_gene, method = "TMM"), log = TRUE, prior.count = 0.5)

## Exons
assays(rse_exon, withDimnames=FALSE)$norm_counts<- edgeR::cpm(calcNormFactors(rse_exon, method = "TMM"), log = TRUE, prior.count = 0.5)

## Junctions
## TMMwsp method for >80% of zeros
assays(rse_jx, withDimnames=FALSE)$norm_counts<- edgeR::cpm(calcNormFactors(rse_jx, method = "TMMwsp"), log = TRUE, prior.count = 0.5)

## Transcripts
## Scale TPM (Transcripts per million) to log2(TPM + 0.5)
assays(rse_tx)$norm_tpm<-log2(assays(rse_tx)$tpm + 0.5)





## Compute QC metrics for posterior QCA at gene level
## Mt and ribosomal genes 
subsets=list(Mito=which(seqnames(rse_gene)=="chrM"), 
             Ribo=grep("rRNA",rowData(rse_gene)$gene_type))
## Add general qc-stats based on counts
rse_gene <-addPerCellQC(rse_gene, subsets)




## 1.3 Data filtering 
## Feature filtering based on counts

## Add design matrix to account for group differences
group <- factor(paste(rse_gene$Tissue, rse_gene$Age, rse_gene$Expt, rse_gene$Group, sep="."))
design <- model.matrix(~0 + group) 

## Filter genes with k CPM in at least n samples
genes_kept<-rse_gene[which(filterByExpr(assay(rse_gene), design=design)),]
dim(genes_kept)
# 20891 208 

## Filter exons
exons_kept<-rse_exon[which(filterByExpr(assay(rse_exon), design=design)),]
dim(exons_kept)
# 298076    208

## Filter junctions
jx_kept<-rse_jx[which(filterByExpr(assay(rse_jx), design=design)),]
dim(jx_kept)
# 191553   208

## Filter TPM
## Identify potential cutoffs
seed <- 20191217
expression_cutoff(assays(rse_tx)$tpm, seed = seed, k=2)
# 2022-06-25 16:27:50 the suggested expression cutoff is 0.28
# percent_features_cut  samples_nonzero_cut 
#                 0.29                 0.27 
cutoff<-0.28
## Transcripts that pass cutoff 
tx_kept<-rse_tx[rowMeans(assays(rse_tx)$tpm) > cutoff,]
dim(tx_kept)
# 58693   208



## 1.4 Data separation
## Brain data
## Genes
rse_gene_brain<-rse_gene_norm[,(rse_gene_norm$Tissue=="Brain")]
save(rse_gene_brain, file = 'processed-data/02_build_objects/rse_gene_brain.Rdata')
## Exons
rse_exon_brain<-rse_exon_norm[,(rse_exon_norm$Tissue=="Brain")]
save(rse_exon_brain, file = 'processed-data/02_build_objects/rse_exon_brain.Rdata')
## Jxn
rse_jx_brain<-rse_jx_norm[,(rse_jx_norm$Tissue=="Brain")]
save(rse_jx_brain, file = 'processed-data/02_build_objects/rse_jx_brain.Rdata')
## Tx
rse_tx_brain<-rse_tx_norm[,(rse_tx_norm$Tissue=="Brain")]
save(rse_tx_brain, file = 'processed-data/02_build_objects/rse_tx_brain.Rdata')

## Blood data
## Genes
rse_gene_blood<-rse_gene_norm[,(rse_gene_norm$Tissue=="Blood")]
save(rse_gene_blood, file = 'processed-data/02_build_objects/rse_gene_blood.Rdata')





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
