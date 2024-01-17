#Files at the cluster are under /dcl01/lieber/ajaffe/lab/smokingMouse_Indirects/
#processed-data/01_SPEAQeasy/pipeline_output/count_objects/


# 1. Data preparation


## Load required libraries
library(SummarizedExperiment)
library(recount)
library(edgeR)
library(jaffelab)
library(scuttle)
library(here)
library(WGCNA)
library(scater)
library(biomartr)
library(sessioninfo)


## 1.1 Data exploration and correction

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
save(pheno, file="processed-data/02_build_objects/pheno.tsv")

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
assays(rse_gene, withDimnames=FALSE)$logcounts <- edgeR::cpm(calcNormFactors(rse_gene, method = "TMM"), log = TRUE, prior.count = 0.5)

## Exons
assays(rse_exon, withDimnames=FALSE)$logcounts<- edgeR::cpm(calcNormFactors(rse_exon, method = "TMM"), log = TRUE, prior.count = 0.5)

## Junctions
## TMMwsp method for >80% of zeros
assays(rse_jx, withDimnames=FALSE)$logcounts<- edgeR::cpm(calcNormFactors(rse_jx, method = "TMMwsp"), log = TRUE, prior.count = 0.5)

## Transcripts
## Scale TPM (Transcripts per million) to log2(TPM + 0.5)
assays(rse_tx)$logcounts<-log2(assays(rse_tx)$tpm + 0.5)

## Save rse objects with the assays of normalized counts
save(rse_exon, file="processed-data/02_build_objects/rse_exon_logcounts.Rdata")
save(rse_jx, file="processed-data/02_build_objects/rse_jx_logcounts.Rdata")
save(rse_tx, file="processed-data/02_build_objects/rse_tx_logcounts.Rdata")





## Compute QC metrics for posterior QCA at gene level
## Mt and ribosomal genes 
subsets=list(Mito=which(seqnames(rse_gene)=="chrM"), 
             Ribo=grep("rRNA",rowData(rse_gene)$gene_type))
## Add general qc-stats based on counts
rse_gene <-addPerCellQC(rse_gene, subsets)

## Save rse_gene with logcounts and qc-stats
save(rse_gene, file="processed-data/02_build_objects/rse_gene_logcounts.Rdata")




## 1.3 Data filtering 
## Feature filtering based on counts

## Filter genes with k CPM in at least n samples
## Add design matrix to account for group differences
rse_gene_filt<-rse_gene[which(filterByExpr(assay(rse_gene), 
                design=with(colData(rse_gene), model.matrix(~ Tissue + Age + Expt + Group)))),]
dim(rse_gene_filt)
# 19974   208

## Add actual gene symbols instead of MGI symbols
rowData(rse_gene_filt)$MGI_Symbol<-rowData(rse_gene_filt)$Symbol
symbols<-biomart(genes  = rowData(rse_gene_filt)$ensemblID,
                 mart       = "ENSEMBL_MART_ENSEMBL",
                 dataset    = "mmusculus_gene_ensembl",
                 attributes = c("external_gene_name"),
                 filters    = "ensembl_gene_id")

## Add MGI/ensembl ID for genes without symbol
## Genes without symbol found
no_symbol<-rowData(rse_gene_filt)$ensemblID[(! rowData(rse_gene_filt)$ensemblID 
                                             %in% symbols$ensembl_gene_id)]
## Genes with empty symbol or NA
which_na_symbol<-which(is.na(symbols$external_gene_name) | symbols$external_gene_name=="")
na_symbol<-symbols[which_na_symbol,1]
## Add these genes' IDs
no_symbol<-append(no_symbol, na_symbol)
## Remove these genes from symbols
symbols<-symbols[-which_na_symbol,]
## Add MGI/ensembl IDs for them
for (gene in no_symbol){
    MGI_symbol<-rowData(rse_gene_filt)[which(rowData(rse_gene_filt)$ensemblID==gene), "MGI_Symbol"]
    if (! is.na(MGI_symbol)) {
      symbols[nrow(symbols)+1,]<-c(gene, MGI_symbol)
    }
    else {
      symbols[nrow(symbols)+1,]<-c(gene,gene)
    }
}

## Preserve original genes' order
symbols<-symbols[match(rowData(rse_gene_filt)$ensemblID, symbols$ensembl_gene_id), ]
rowData(rse_gene_filt)$Symbol<-symbols$external_gene_name
save(rse_gene_filt, file = 'processed-data/02_build_objects/rse_gene_filt.Rdata')



## Filter exons
rse_exon_filt<-rse_exon[which(filterByExpr(assay(rse_exon), 
                design=with(colData(rse_exon), model.matrix(~ Tissue + Age + Expt + Group)))),]
dim(rse_exon_filt)
# 290800 208 

## Add gene symbol of all genes associated with exons
rowData(rse_exon_filt)$MGI_Symbol<-rowData(rse_exon_filt)$Symbol
symbols<-biomart(genes  = unique(rowData(rse_exon_filt)$ensemblID),
                 mart       = "ENSEMBL_MART_ENSEMBL",
                 dataset    = "mmusculus_gene_ensembl",
                 attributes = c("external_gene_name"),
                 filters    = "ensembl_gene_id")

## Add MGI/ensembl ID for genes without symbol
## Genes without symbol found
no_symbol<-unique(rowData(rse_exon_filt)$ensemblID)[which(! unique(rowData(rse_exon_filt)$ensemblID) 
                                                  %in%  symbols$ensembl_gene_id)]
## Genes with empty symbol or NA
which_na_symbol<-which(is.na(symbols$external_gene_name) | symbols$external_gene_name=="")
na_symbol<-symbols[which_na_symbol,1]
## Add these genes' IDs
no_symbol<-append(no_symbol, na_symbol)
## Remove these genes from symbols
symbols<-symbols[-which_na_symbol,]
## Add MGI/ensembl IDs for them
for (gene in no_symbol){
  MGI_symbol<-unique(rowData(rse_exon_filt)[which(rowData(rse_exon_filt)$ensemblID==gene), "MGI_Symbol"])
  if (! (is.na(MGI_symbol) | length(MGI_symbol)==0)) {
    symbols[nrow(symbols)+1,]<-c(gene, MGI_symbol)
  }
  else {
    symbols[nrow(symbols)+1,]<-c(gene,gene)
  }
}

## Assign gene symbol to all exons
symbols<-symbols[match(rowData(rse_exon_filt)$ensemblID, symbols$ensembl_gene_id), ]
rowData(rse_exon_filt)$Symbol<-symbols$external_gene_name
save(rse_exon_filt, file = 'processed-data/02_build_objects/rse_exon_filt.Rdata')



## Filter junctions
rse_jx_filt<-rse_jx[which(filterByExpr(assay(rse_jx), 
                design=with(colData(rse_jx), model.matrix(~ Tissue + Age + Expt + Group)))),]
save(rse_jx_filt, file = 'processed-data/02_build_objects/rse_jx_filt.Rdata')
dim(rse_jx_filt)
# 176670 208 


## Filter TPM
## Identify potential cutoffs
seed <- 20191217
expression_cutoff(assays(rse_tx)$tpm, seed = seed, k=2)
# 2022-06-25 16:27:50 the suggested expression cutoff is 0.28
# percent_features_cut  samples_nonzero_cut 
#                 0.29                 0.27 
cutoff<-0.28
## Transcripts that pass cutoff 
rse_tx_filt<-rse_tx[rowMeans(assays(rse_tx)$tpm) > cutoff,]
save(rse_tx_filt, file = 'processed-data/02_build_objects/rse_tx_filt.Rdata')
dim(rse_tx_filt)
# 58693   208





## 1.4 Data separation

## Brain data
## Genes
rse_gene_brain<-rse_gene_filt[,(rse_gene_filt$Tissue=="Brain")]
save(rse_gene_brain, file = 'processed-data/02_build_objects/rse_gene_brain.Rdata')

## Blood data
## Genes
rse_gene_blood<-rse_gene_filt[,(rse_gene_filt$Tissue=="Blood")]
save(rse_gene_blood, file = 'processed-data/02_build_objects/rse_gene_blood.Rdata')


## Adult and brain data
## Genes
rse_gene_brain_adults<-rse_gene_brain[,(rse_gene_brain$Age=="Adult")]
save(rse_gene_brain_adults, file = 'processed-data/02_build_objects/rse_gene_brain_adults.Rdata')
## Exons
rse_exon_brain_adults<-rse_exon_filt[,(rse_exon_filt$Tissue=="Brain" & rse_exon_filt$Age=="Adult")]
save(rse_exon_brain_adults, file = 'processed-data/02_build_objects/rse_exon_brain_adults.Rdata')
## Jxn
rse_jx_brain_adults<-rse_jx_filt[,(rse_jx_filt$Tissue=="Brain" & rse_jx_filt$Age=="Adult")]
save(rse_jx_brain_adults, file = 'processed-data/02_build_objects/rse_jx_brain_adults.Rdata')
## Tx
rse_tx_brain_adults<-rse_tx_filt[,(rse_tx_filt$Tissue=="Brain" & rse_tx_filt$Age=="Adult")]
save(rse_tx_brain_adults, file = 'processed-data/02_build_objects/rse_tx_brain_adults.Rdata')


## Pup and brain data
## Genes
rse_gene_brain_pups<-rse_gene_brain[,(rse_gene_brain$Age=="Pup")]
save(rse_gene_brain_pups, file = 'processed-data/02_build_objects/rse_gene_brain_pups.Rdata')
## Exons
rse_exon_brain_pups<-rse_exon_filt[,(rse_exon_filt$Tissue=="Brain" & rse_exon_filt$Age=="Pup")]
save(rse_exon_brain_pups, file = 'processed-data/02_build_objects/rse_exon_brain_pups.Rdata')
## Jxn
rse_jx_brain_pups<-rse_jx_filt[,(rse_jx_filt$Tissue=="Brain" & rse_jx_filt$Age=="Pup")]
save(rse_jx_brain_pups, file = 'processed-data/02_build_objects/rse_jx_brain_pups.Rdata')
## Tx
rse_tx_brain_pups<-rse_tx_filt[,(rse_tx_filt$Tissue=="Brain" & rse_tx_filt$Age=="Pup")]
save(rse_tx_brain_pups, file = 'processed-data/02_build_objects/rse_tx_brain_pups.Rdata')





## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.0 (2023-04-21)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Mexico_City
# date     2023-12-12
# rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
# pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version    date (UTC) lib source
# AnnotationDbi          1.63.2     2023-07-03 [1] Bioconductor
# aod                    1.3.2      2022-04-02 [1] CRAN (R 4.3.0)
# backports              1.4.1      2021-12-13 [1] CRAN (R 4.3.0)
# base64enc              0.1-3      2015-07-28 [1] CRAN (R 4.3.0)
# beachmat               2.16.0     2023-05-08 [1] Bioconductor
# beeswarm               0.4.0      2021-06-01 [1] CRAN (R 4.3.0)
# Biobase              * 2.61.0     2023-06-02 [1] Bioconductor
# BiocFileCache          2.9.1      2023-07-14 [1] Bioconductor
# BiocGenerics         * 0.47.0     2023-06-02 [1] Bioconductor
# BiocIO                 1.11.0     2023-06-02 [1] Bioconductor
# BiocNeighbors          1.18.0     2023-05-08 [1] Bioconductor
# BiocParallel         * 1.35.3     2023-07-07 [1] Bioconductor
# BiocSingular           1.16.0     2023-05-08 [1] Bioconductor
# biomaRt                2.57.1     2023-06-14 [1] Bioconductor
# Biostrings             2.69.2     2023-07-05 [1] Bioconductor
# bit                    4.0.5      2022-11-15 [1] CRAN (R 4.3.0)
# bit64                  4.0.5      2020-08-30 [1] CRAN (R 4.3.0)
# bitops                 1.0-7      2021-04-24 [1] CRAN (R 4.3.0)
# blob                   1.2.4      2023-03-17 [1] CRAN (R 4.3.0)
# boot                   1.3-28.1   2022-11-22 [1] CRAN (R 4.3.0)
# broom                  1.0.5      2023-06-09 [1] CRAN (R 4.3.0)
# BSgenome               1.69.0     2023-07-07 [1] Bioconductor
# bumphunter             1.43.0     2023-06-14 [1] Bioconductor
# cachem                 1.0.8      2023-05-01 [1] CRAN (R 4.3.0)
# caTools                1.18.2     2021-03-28 [1] CRAN (R 4.3.0)
# checkmate              2.2.0      2023-04-27 [1] CRAN (R 4.3.0)
# cli                    3.6.1      2023-03-23 [1] CRAN (R 4.3.0)
# cluster                2.1.4      2022-08-22 [1] CRAN (R 4.3.0)
# codetools              0.2-19     2023-02-01 [1] CRAN (R 4.3.0)
# colorspace             2.1-0      2023-01-23 [1] CRAN (R 4.3.0)
# corpcor                1.6.10     2021-09-16 [1] CRAN (R 4.3.0)
# cowplot              * 1.1.1      2020-12-30 [1] CRAN (R 4.3.0)
# crayon                 1.5.2      2022-09-29 [1] CRAN (R 4.3.0)
# curl                   5.0.1      2023-06-07 [1] CRAN (R 4.3.0)
# data.table             1.14.8     2023-02-17 [1] CRAN (R 4.3.0)
# DBI                    1.1.3      2022-06-18 [1] CRAN (R 4.3.0)
# dbplyr                 2.3.3      2023-07-07 [1] CRAN (R 4.3.0)
# DelayedArray           0.26.6     2023-07-02 [1] Bioconductor
# DelayedMatrixStats     1.23.0     2023-04-25 [1] Bioconductor
# derfinder              1.35.0     2023-07-07 [1] Bioconductor
# derfinderHelper        1.35.0     2023-06-02 [1] Bioconductor
# digest                 0.6.33     2023-07-07 [1] CRAN (R 4.3.0)
# doParallel             1.0.17     2022-02-07 [1] CRAN (R 4.3.0)
# doRNG                  1.8.6      2023-01-16 [1] CRAN (R 4.3.0)
# downloader             0.4        2015-07-09 [1] CRAN (R 4.3.0)
# dplyr                  1.1.2      2023-04-20 [1] CRAN (R 4.3.0)
# dynamicTreeCut       * 1.63-1     2016-03-11 [1] CRAN (R 4.3.0)
# edgeR                * 3.43.7     2023-06-21 [1] Bioconductor
# EnvStats               2.8.0      2023-07-08 [1] CRAN (R 4.3.0)
# evaluate               0.21       2023-05-05 [1] CRAN (R 4.3.0)
# fANCOVA                0.6-1      2020-11-13 [1] CRAN (R 4.3.0)
# fansi                  1.0.5      2023-10-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1      2022-07-06 [1] CRAN (R 4.3.0)
# fastcluster          * 1.2.3      2021-05-24 [1] CRAN (R 4.3.0)
# fastmap                1.1.1      2023-02-24 [1] CRAN (R 4.3.0)
# filelock               1.0.2      2018-10-05 [1] CRAN (R 4.3.0)
# foreach                1.5.2      2022-02-02 [1] CRAN (R 4.3.0)
# foreign                0.8-84     2022-12-06 [1] CRAN (R 4.3.0)
# Formula                1.2-5      2023-02-24 [1] CRAN (R 4.3.0)
# fs                     1.6.3      2023-07-20 [1] CRAN (R 4.3.0)
# gargle                 1.5.2      2023-07-20 [1] CRAN (R 4.3.0)
# generics               0.1.3      2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.37.2     2023-06-21 [1] Bioconductor
# GenomeInfoDbData       1.2.10     2023-05-28 [1] Bioconductor
# GenomicAlignments      1.37.0     2023-07-07 [1] Bioconductor
# GenomicFeatures        1.53.1     2023-06-22 [1] Bioconductor
# GenomicFiles           1.37.0     2023-07-07 [1] Bioconductor
# GenomicRanges        * 1.53.1     2023-06-02 [1] Bioconductor
# GEOquery               2.69.0     2023-06-16 [1] Bioconductor
# ggbeeswarm             0.7.2      2023-04-29 [1] CRAN (R 4.3.0)
# ggnewscale           * 0.4.9      2023-05-25 [1] CRAN (R 4.3.0)
# ggplot2              * 3.4.4      2023-10-12 [1] CRAN (R 4.3.1)
# ggrepel              * 0.9.3      2023-02-03 [1] CRAN (R 4.3.0)
# glue                   1.6.2      2022-02-24 [1] CRAN (R 4.3.0)
# GO.db                  3.17.0     2023-05-28 [1] Bioconductor
# googledrive            2.1.1      2023-06-11 [1] CRAN (R 4.3.0)
# gplots                 3.1.3      2022-04-25 [1] CRAN (R 4.3.0)
# gridExtra            * 2.3        2017-09-09 [1] CRAN (R 4.3.0)
# gtable                 0.3.4      2023-08-21 [1] CRAN (R 4.3.0)
# gtools                 3.9.4      2022-11-27 [1] CRAN (R 4.3.0)
# here                 * 1.0.1      2020-12-13 [1] CRAN (R 4.3.0)
# Hmisc                * 5.1-0      2023-05-08 [1] CRAN (R 4.3.0)
# hms                    1.1.3      2023-03-21 [1] CRAN (R 4.3.0)
# htmlTable              2.4.1      2022-07-07 [1] CRAN (R 4.3.0)
# htmltools              0.5.5      2023-03-23 [1] CRAN (R 4.3.0)
# htmlwidgets            1.6.2      2023-03-17 [1] CRAN (R 4.3.0)
# httr                   1.4.6      2023-05-08 [1] CRAN (R 4.3.0)
# IRanges              * 2.35.2     2023-06-23 [1] Bioconductor
# irlba                  2.3.5.1    2022-10-03 [1] CRAN (R 4.3.0)
# iterators              1.0.14     2022-02-05 [1] CRAN (R 4.3.0)
# jaffelab             * 0.99.32    2023-05-28 [1] Github (LieberInstitute/jaffelab@21e6574)
# jsonlite               1.8.7      2023-06-29 [1] CRAN (R 4.3.0)
# KEGGREST               1.41.0     2023-07-07 [1] Bioconductor
# KernSmooth             2.23-22    2023-07-10 [1] CRAN (R 4.3.0)
# knitr                  1.43       2023-05-25 [1] CRAN (R 4.3.0)
# labeling               0.4.3      2023-08-29 [1] CRAN (R 4.3.0)
# lattice                0.21-8     2023-04-05 [1] CRAN (R 4.3.0)
# lifecycle              1.0.3      2022-10-07 [1] CRAN (R 4.3.0)
# limma                * 3.57.6     2023-06-21 [1] Bioconductor
# lme4                   1.1-34     2023-07-04 [1] CRAN (R 4.3.0)
# lmerTest               3.1-3      2020-10-23 [1] CRAN (R 4.3.0)
# locfit                 1.5-9.8    2023-06-11 [1] CRAN (R 4.3.0)
# magrittr               2.0.3      2022-03-30 [1] CRAN (R 4.3.0)
# MASS                   7.3-60     2023-05-04 [1] CRAN (R 4.3.0)
# Matrix                 1.6-0      2023-07-08 [1] CRAN (R 4.3.0)
# MatrixGenerics       * 1.13.0     2023-05-20 [1] Bioconductor
# matrixStats          * 1.0.0      2023-06-02 [1] CRAN (R 4.3.0)
# memoise                2.0.1      2021-11-26 [1] CRAN (R 4.3.0)
# minqa                  1.2.5      2022-10-19 [1] CRAN (R 4.3.0)
# munsell                0.5.0      2018-06-12 [1] CRAN (R 4.3.0)
# mvtnorm                1.2-2      2023-06-08 [1] CRAN (R 4.3.0)
# nlme                   3.1-162    2023-01-31 [1] CRAN (R 4.3.0)
# nloptr                 2.0.3      2022-05-26 [1] CRAN (R 4.3.0)
# nnet                   7.3-19     2023-05-03 [1] CRAN (R 4.3.0)
# numDeriv               2016.8-1.1 2019-06-06 [1] CRAN (R 4.3.0)
# pbkrtest               0.5.2      2023-01-19 [1] CRAN (R 4.3.0)
# pillar                 1.9.0      2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.3.0)
# plyr                   1.8.8      2022-11-11 [1] CRAN (R 4.3.0)
# png                    0.1-8      2022-11-29 [1] CRAN (R 4.3.0)
# prettyunits            1.1.1      2020-01-24 [1] CRAN (R 4.3.0)
# progress               1.2.2      2019-05-16 [1] CRAN (R 4.3.0)
# purrr                  1.0.1      2023-01-10 [1] CRAN (R 4.3.0)
# qvalue                 2.33.0     2023-05-11 [1] Bioconductor
# R6                     2.5.1      2021-08-19 [1] CRAN (R 4.3.0)
# rafalib              * 1.0.0      2015-08-09 [1] CRAN (R 4.3.0)
# ragg                   1.2.5      2023-01-12 [1] CRAN (R 4.3.0)
# rappdirs               0.3.3      2021-01-31 [1] CRAN (R 4.3.0)
# rbibutils              2.2.13     2023-01-13 [1] CRAN (R 4.3.0)
# RColorBrewer           1.1-3      2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.11     2023-07-06 [1] CRAN (R 4.3.0)
# RCurl                  1.98-1.12  2023-03-27 [1] CRAN (R 4.3.0)
# Rdpack                 2.4        2022-07-20 [1] CRAN (R 4.3.0)
# readr                  2.1.4      2023-02-10 [1] CRAN (R 4.3.0)
# recount              * 1.27.0     2023-04-25 [1] Bioconductor
# remaCor                0.0.16     2023-06-21 [1] CRAN (R 4.3.0)
# rentrez                1.2.3      2020-11-10 [1] CRAN (R 4.3.0)
# reshape2               1.4.4      2020-04-09 [1] CRAN (R 4.3.0)
# restfulr               0.0.15     2022-06-16 [1] CRAN (R 4.3.0)
# RhpcBLASctl            0.23-42    2023-02-11 [1] CRAN (R 4.3.0)
# rjson                  0.2.21     2022-01-09 [1] CRAN (R 4.3.0)
# rlang                * 1.1.1      2023-04-28 [1] CRAN (R 4.3.0)
# rmarkdown              2.23       2023-07-01 [1] CRAN (R 4.3.0)
# rngtools               1.5.2      2021-09-20 [1] CRAN (R 4.3.0)
# rpart                  4.1.19     2022-10-21 [1] CRAN (R 4.3.0)
# rprojroot              2.0.3      2022-04-02 [1] CRAN (R 4.3.0)
# Rsamtools              2.17.0     2023-07-07 [1] Bioconductor
# RSQLite                2.3.1      2023-04-03 [1] CRAN (R 4.3.0)
# rstudioapi             0.15.0     2023-07-07 [1] CRAN (R 4.3.0)
# rsvd                   1.0.5      2021-04-16 [1] CRAN (R 4.3.0)
# rtracklayer            1.61.0     2023-07-07 [1] Bioconductor
# S4Arrays               1.1.4      2023-06-02 [1] Bioconductor
# S4Vectors            * 0.39.1     2023-06-02 [1] Bioconductor
# ScaledMatrix           1.9.1      2023-05-03 [1] Bioconductor
# scales                 1.2.1      2022-08-20 [1] CRAN (R 4.3.0)
# scater               * 1.29.1     2023-06-15 [1] Github (davismcc/scater@eb6b801)
# scuttle              * 1.9.4      2023-01-23 [1] Bioconductor
# segmented              1.6-4      2023-04-13 [1] CRAN (R 4.3.0)
# sessioninfo          * 1.2.2      2021-12-06 [1] CRAN (R 4.3.0)
# SingleCellExperiment * 1.23.0     2023-04-25 [1] Bioconductor
# sparseMatrixStats      1.13.0     2023-05-20 [1] Bioconductor
# stringi                1.7.12     2023-01-11 [1] CRAN (R 4.3.0)
# stringr              * 1.5.0      2022-12-02 [1] CRAN (R 4.3.0)
# SummarizedExperiment * 1.30.2     2023-06-06 [1] Bioconductor
# survival               3.5-5      2023-03-12 [1] CRAN (R 4.3.0)
# systemfonts            1.0.4      2022-02-11 [1] CRAN (R 4.3.0)
# textshaping            0.3.6      2021-10-13 [1] CRAN (R 4.3.0)
# tibble                 3.2.1      2023-03-20 [1] CRAN (R 4.3.0)
# tidyr                  1.3.0      2023-01-24 [1] CRAN (R 4.3.0)
# tidyselect             1.2.0      2022-10-10 [1] CRAN (R 4.3.0)
# tzdb                   0.4.0      2023-05-12 [1] CRAN (R 4.3.0)
# utf8                   1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
# variancePartition    * 1.32.2     2023-11-14 [1] Bioconductor
# VariantAnnotation      1.47.1     2023-07-07 [1] Bioconductor
# vctrs                  0.6.4      2023-10-12 [1] CRAN (R 4.3.1)
# vipor                  0.4.5      2017-03-22 [1] CRAN (R 4.3.0)
# viridis                0.6.3      2023-05-03 [1] CRAN (R 4.3.0)
# viridisLite            0.4.2      2023-05-02 [1] CRAN (R 4.3.0)
# withr                  2.5.1      2023-09-26 [1] CRAN (R 4.3.1)
# xfun                   0.39       2023-04-20 [1] CRAN (R 4.3.0)
# XML                    3.99-0.14  2023-03-19 [1] CRAN (R 4.3.0)
# xml2                   1.3.5      2023-07-06 [1] CRAN (R 4.3.0)
# XVector                0.41.1     2023-06-02 [1] Bioconductor
# yaml                   2.3.7      2023-01-23 [1] CRAN (R 4.3.0)
# zlibbioc               1.47.0     2023-05-20 [1] Bioconductor
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
