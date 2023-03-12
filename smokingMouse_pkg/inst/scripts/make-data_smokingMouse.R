
## Create complete RSE objects at gene, exon, jxn and tx level adding lognorm counts, info of which features and samples 
## passed filtering steps and which features were DE in the different groups of samples


library(SummarizedExperiment)
library(here)
library(sessioninfo)

## Original rse objects with the correct sample info in colData and logcounts as an assay:
## Gene rse contains QC stats in colData as well
load(here("processed-data/02_build_objects/rse_gene_logcounts.Rdata"), verbose = TRUE)
# rse_gene
load(here("processed-data/02_build_objects/rse_exon_logcounts.Rdata"), verbose = TRUE)
# rse_exon
load(here("processed-data/02_build_objects/rse_jx_logcounts.Rdata"), verbose = TRUE)
# rse_jx
load(here("processed-data/02_build_objects/rse_tx_logcounts.Rdata"), verbose = TRUE)
# rse_tx


####################################
##    02_build_objects analyses
####################################
 
## Analyses:

## 1. Feature filtering

## rse objects with retained features after feature filtering
load(here("processed-data/02_build_objects/rse_gene_filt.Rdata"), verbose = TRUE)
# rse_gene_filt
load(here("processed-data/02_build_objects/rse_exon_filt.Rdata"), verbose = TRUE)
# rse_exon_filt
load(here("processed-data/02_build_objects/rse_jx_filt.Rdata"), verbose = TRUE)
# rse_jx_filt
load(here("processed-data/02_build_objects/rse_tx_filt.Rdata"), verbose = TRUE)
# rse_tx_filt

## Add column to rowData (of the original rse objects) with the info of the features that were retained and dropped after feature filtering
rowData(rse_gene)$retained_after_feature_filtering <- unlist(sapply(rowData(rse_gene)$gencodeID, function(x){if(x %in% rowData(rse_gene_filt)$gencodeID){TRUE} else {FALSE}}))
rowData(rse_exon)$retained_after_feature_filtering <- unlist(sapply(rowData(rse_exon)$exon_gencodeID, function(x){if(x %in% rowData(rse_exon_filt)$exon_gencodeID){TRUE} else {FALSE}}))
rowData(rse_jx)$retained_after_feature_filtering <- unlist(sapply(rownames(rowData(rse_jx)), function(x){if(x %in% rownames(rowData(rse_jx_filt))){TRUE} else {FALSE}}))
rowData(rse_tx)$retained_after_feature_filtering <- unlist(sapply(rowData(rse_tx)$transcript_id, function(x){if(x %in% rowData(rse_tx_filt)$transcript_id){TRUE} else {FALSE}}))



#############################
##    03_EDA analyses
#############################

########### 02_QC ###########

## 1. Sample filtering by QC metrics

## rse objects with samples that passed the QC filtering
load(here("processed-data/03_EDA/02_QC/rse_gene_blood_qc.Rdata"), verbose = TRUE)
# rse_gene_blood_qc
load(here("processed-data/03_EDA/02_QC/rse_gene_brain_pups_qc.Rdata"), verbose = TRUE)
# rse_gene_brain_pups_qc
load(here("processed-data/03_EDA/02_QC/rse_gene_brain_adults_qc.Rdata"), verbose = TRUE)
# rse_gene_brain_adults_qc

## Add column to colData (of the original rse objects) with the info of samples retained and dropped after sample filtering by QC
retained_samples <- union(rse_gene_blood_qc$SAMPLE_ID, union(rse_gene_brain_pups_qc$SAMPLE_ID, rse_gene_brain_adults_qc$SAMPLE_ID))
colData(rse_gene)$retained_after_QC_sample_filtering <- unlist(sapply(rse_gene$SAMPLE_ID, function(x){if(x %in% retained_samples){TRUE} else {FALSE}}))
colData(rse_exon)$retained_after_QC_sample_filtering <- unlist(sapply(rse_exon$SAMPLE_ID, function(x){if(x %in% retained_samples){TRUE} else {FALSE}}))
colData(rse_tx)$retained_after_QC_sample_filtering <- unlist(sapply(rse_tx$SAMPLE_ID, function(x){if(x %in% retained_samples){TRUE} else {FALSE}}))
colData(rse_jx)$retained_after_QC_sample_filtering <- unlist(sapply(rse_jx$SAMPLE_ID, function(x){if(x %in% retained_samples){TRUE} else {FALSE}}))



########### 03_PCA_MDS ###########

## 1. Manual sample filtering (removal of rare samples identified in PCA plots)

## rse objects with samples that passed QC and manual filtering
## (No blood sample was manually removed)
load(here("processed-data/03_EDA/03_PCA/rse_gene_brain_adults_qc_afterPCA.Rdata"), verbose = TRUE)
# rse_gene_brain_adults_qc_afterPCA
load(here("processed-data/03_EDA/03_PCA/rse_gene_brain_pups_qc_afterPCA.Rdata"), verbose = TRUE)
# rse_gene_brain_pups_qc_afterPCA

## Add column to colData (of the original rse objects) with the info of samples retained and dropped after manual sample filtering
retained_samples <- union(rse_gene_blood_qc$SAMPLE_ID, union(rse_gene_brain_adults_qc_afterPCA$SAMPLE_ID, rse_gene_brain_pups_qc_afterPCA$SAMPLE_ID))
colData(rse_gene)$retained_after_manual_sample_filtering <- unlist(sapply(rse_gene$SAMPLE_ID, function(x){if(x %in% retained_samples){TRUE} else {FALSE}}))
colData(rse_exon)$retained_after_manual_sample_filtering <- unlist(sapply(rse_exon$SAMPLE_ID, function(x){if(x %in% retained_samples){TRUE} else {FALSE}}))
colData(rse_tx)$retained_after_manual_sample_filtering <- unlist(sapply(rse_tx$SAMPLE_ID, function(x){if(x %in% retained_samples){TRUE} else {FALSE}}))
colData(rse_jx)$retained_after_manual_sample_filtering <- unlist(sapply(rse_jx$SAMPLE_ID, function(x){if(x %in% retained_samples){TRUE} else {FALSE}}))



#############################
##    04_DEA analyses
#############################

## 1. Differential expression analyses (with fitted models only)

########### Gene level ###########

## Data frames with DEGs 
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata"), verbose = TRUE)
# de_genes_pups_nicotine_fitted
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata"), verbose = TRUE)
# de_genes_pups_smoking_fitted


## Add info of DEGs to rowData of rse objects

## No DEGs in adult blood samples
## Smoking
rowData(rse_gene)$DE_in_adult_blood_smoking <- rep(FALSE, dim(rse_gene)[1])

## No DEGs in adult brain samples 
## Nicotine
rowData(rse_gene)$DE_in_adult_brain_nicotine <- rep(FALSE, dim(rse_gene)[1])
## Smoking
rowData(rse_gene)$DE_in_adult_brain_smoking <- rep(FALSE, dim(rse_gene)[1])

## DEGs in pup brain samples
## Nicotine
rowData(rse_gene)$DE_in_pup_brain_nicotine <- unlist(sapply(rowData(rse_gene)$gencodeID, function(x){if(x %in% de_genes_pups_nicotine_fitted$gencodeID){TRUE} else {FALSE}}))
## Smoking
rowData(rse_gene)$DE_in_pup_brain_smoking <- unlist(sapply(rowData(rse_gene)$gencodeID, function(x){if(x %in% de_genes_pups_smoking_fitted$gencodeID){TRUE} else {FALSE}}))


########### Exon level ###########
## (For pup brain samples only)

## Data frames with DE exons 
load(here("processed-data/04_DEA/Exon_analysis/de_exons_nic.Rdata"), verbose = TRUE)
# de_exons_nic
load(here("processed-data/04_DEA/Exon_analysis/de_exons_smo.Rdata"), verbose = TRUE)
# de_exons_smo

## Add info of DE exons to rowData of rse objects
## DE exons in pup brain samples
## Nicotine
rowData(rse_exon)$DE_in_pup_brain_nicotine <- unlist(sapply(rowData(rse_exon)$exon_gencodeID, function(x){if(x %in% de_exons_nic$exon_gencodeID){TRUE} else {FALSE}}))
## Smoking
rowData(rse_exon)$DE_in_pup_brain_smoking <- unlist(sapply(rowData(rse_exon)$exon_gencodeID, function(x){if(x %in% de_exons_smo$exon_gencodeID){TRUE} else {FALSE}}))


########### Tx level ###########
## (For pup brain samples only)

## Data frames with DE txs
load(here("processed-data/04_DEA/Tx_analysis/de_tx_nic.Rdata"), verbose = TRUE)
# de_tx_nic
load(here("processed-data/04_DEA/Tx_analysis/de_tx_smo.Rdata"), verbose = TRUE)
# de_tx_smo

## Add info of DE txs to rowData of rse objects
## DE txs in pup brain samples
## Nicotine
rowData(rse_tx)$DE_in_pup_brain_nicotine <- unlist(sapply(rowData(rse_tx)$transcript_id, function(x){if(x %in% de_tx_nic$transcript_id){TRUE} else {FALSE}}))
## Smoking
rowData(rse_tx)$DE_in_pup_brain_smoking <- unlist(sapply(rowData(rse_tx)$transcript_id, function(x){if(x %in% de_tx_smo$transcript_id){TRUE} else {FALSE}}))


########### Jxn level ###########
## (For pup brain samples only)

## Data frames with DE jxns
load(here("processed-data/04_DEA/Jx_analysis/de_jxns_nic.Rdata"), verbose = TRUE)
# de_jxns_nic
load(here("processed-data/04_DEA/Jx_analysis/de_jxns_smo.Rdata"), verbose = TRUE)
# de_jxns_smo

## Add info of DE jxns to rowData of rse objects
## DE jxns in pup brain samples
## Nicotine
rowData(rse_jx)$DE_in_pup_brain_nicotine <- unlist(sapply(rownames(rowData(rse_jx)), function(x){if(x %in% rownames(de_jxns_nic)){TRUE} else {FALSE}}))
## Smoking
rowData(rse_jx)$DE_in_pup_brain_smoking <- unlist(sapply(rownames(rowData(rse_jx)), function(x){if(x %in% rownames(de_jxns_smo)){TRUE} else {FALSE}}))



## Save complete datasets
save(rse_gene, file="smokingMouse_pkg/inst/extdata/rse_gene_complete.Rdata")
save(rse_exon, file="smokingMouse_pkg/inst/extdata/rse_exon_complete.Rdata")
save(rse_tx, file="smokingMouse_pkg/inst/extdata/rse_tx_complete.Rdata")
save(rse_jx, file="smokingMouse_pkg/inst/extdata/rse_jx_complete.Rdata")







## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# setting  value
# version  R version 4.2.2 (2022-10-31)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Mexico_City
# date     2023-03-07
# rstudio  2022.12.0+353 Elsbeth Geranium (desktop)
# pandoc   NA


