
library(SummarizedExperiment)
library(here)
library(sessioninfo)

load(here("raw-data/rse_exon_smoking_mouse_n208.Rdata"))
load(here("raw-data/rse_gene_smoking_mouse_n208.Rdata"))
load(here("raw-data/rse_jx_smoking_mouse_n208.Rdata"))
load(here("raw-data/rse_tx_smoking_mouse_n208.Rdata"))



## Add gene and sample information in the rse objects for each analysis done

####################################
##    02_build_objects analyses
####################################

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

 
## Analyses:

## 1. Feature filtering

load(here("processed-data/02_build_objects/rse_gene_filt.Rdata"), verbose = TRUE)
# rse_gene_filt
load(here("processed-data/02_build_objects/rse_exon_filt.Rdata"), verbose = TRUE)
# rse_exon_filt
load(here("processed-data/02_build_objects/rse_jx_filt.Rdata"), verbose = TRUE)
# rse_jx_filt
load(here("processed-data/02_build_objects/rse_tx_filt.Rdata"), verbose = TRUE)
# rse_tx_filt

## Add column to rowData with the info of the features that were retained and dropped after feature filtering
rowData(rse_gene)$retained_after_feature_filtering <- unlist(sapply(rowData(rse_gene)$gencodeID, function(x){if(x %in% rowData(rse_gene_filt)$gencodeID){"Yes"} else {"No"}}))
rowData(rse_exon)$retained_after_feature_filtering <- unlist(sapply(rowData(rse_exon)$exon_gencodeID, function(x){if(x %in% rowData(rse_exon_filt)$exon_gencodeID){"Yes"} else {"No"}}))
rowData(rse_jx)$retained_after_feature_filtering <- unlist(sapply(rownames(rowData(rse_jx)), function(x){if(x %in% rownames(rowData(rse_jx_filt))){"Yes"} else {"No"}}))
rowData(rse_tx)$retained_after_feature_filtering <- unlist(sapply(rowData(rse_tx)$transcript_id, function(x){if(x %in% rowData(rse_tx_filt)$transcript_id){"Yes"} else {"No"}}))



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
colData(rse_gene)$retained_after_QC_sample_filtering <- unlist(sapply(rse_gene$SAMPLE_ID, function(x){if(x %in% retained_samples){"Yes"} else {"No"}}))
colData(rse_exon)$retained_after_QC_sample_filtering <- unlist(sapply(rse_exon$SAMPLE_ID, function(x){if(x %in% retained_samples){"Yes"} else {"No"}}))
colData(rse_tx)$retained_after_QC_sample_filtering <- unlist(sapply(rse_tx$SAMPLE_ID, function(x){if(x %in% retained_samples){"Yes"} else {"No"}}))
colData(rse_jx)$retained_after_QC_sample_filtering <- unlist(sapply(rse_jx$SAMPLE_ID, function(x){if(x %in% retained_samples){"Yes"} else {"No"}}))



########### 03_PCA_MDS ###########

## 1. Manual sample filtering (after removal of rare samples identified in PCA plots)

## rse objects with samples that passed QC and manual filtering
## (No blood sample was manually removed)
load(here("processed-data/03_EDA/03_PCA/rse_gene_brain_adults_qc_afterPCA.Rdata"), verbose = TRUE)
# rse_gene_brain_adults_qc_afterPCA
load(here("processed-data/03_EDA/03_PCA/rse_gene_brain_pups_qc_afterPCA.Rdata"), verbose = TRUE)
# rse_gene_brain_pups_qc_afterPCA

## Add column to colData (of the original rse objects) with the info of samples retained and dropped after manual sample filtering
retained_samples <- union(rse_gene_blood_qc$SAMPLE_ID, union(rse_gene_brain_adults_qc_afterPCA$SAMPLE_ID, rse_gene_brain_pups_qc_afterPCA$SAMPLE_ID))
colData(rse_gene)$retained_after_manual_sample_filtering <- unlist(sapply(rse_gene$SAMPLE_ID, function(x){if(x %in% retained_samples){"Yes"} else {"No"}}))
colData(rse_exon)$retained_after_manual_sample_filtering <- unlist(sapply(rse_exon$SAMPLE_ID, function(x){if(x %in% retained_samples){"Yes"} else {"No"}}))
colData(rse_tx)$retained_after_manual_sample_filtering <- unlist(sapply(rse_tx$SAMPLE_ID, function(x){if(x %in% retained_samples){"Yes"} else {"No"}}))
colData(rse_jx)$retained_after_manual_sample_filtering <- unlist(sapply(rse_jx$SAMPLE_ID, function(x){if(x %in% retained_samples){"Yes"} else {"No"}}))



#############################
##    04_DEA analyses
#############################

## 1. Differential expression analyses (with fitted models only)

########### Gene level ###########

## Data frames with DEGs 
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata"))

## Add info of DEGs to rowData of rse objects
## No DEGs in adult blood samples
## Smoking
rowData(rse_gene)$DE_in_adult_blood_smoking <- rep("No", dim(rse_gene)[1])

## No DEGs in adult brain samples 
## Nicotine
rowData(rse_gene)$DE_in_adult_brain_nicotine <- rep("No", dim(rse_gene)[1])
## Smoking
rowData(rse_gene)$DE_in_adult_brain_smoking <- rep("No", dim(rse_gene)[1])

## DEGs in pup brain samples
## Nicotine
rowData(rse_gene)$DE_in_pup_brain_nicotine <- unlist(sapply(rowData(rse_gene)$gencodeID, function(x){if(x %in% de_genes_pups_nicotine_fitted$gencodeID){"Yes"} else {"No"}}))
## Smoking
rowData(rse_gene)$DE_in_pup_brain_smoking <- unlist(sapply(rowData(rse_gene)$gencodeID, function(x){if(x %in% de_genes_pups_smoking_fitted$gencodeID){"Yes"} else {"No"}}))



########### Exon level ###########
## For pup brain samples only

## Data frames with DE exons 
load(here("processed-data/04_DEA/Exon_analysis/de_exons_nic.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/de_exons_smo.Rdata"))

## DE exons in pup brain samples
## Nicotine
rowData(rse_gene)$DE_in_pup_brain_nicotine <- unlist(sapply(rowData(rse_gene)$gencodeID, function(x){if(x %in% de_genes_pups_nicotine_fitted$gencodeID){"Yes"} else {"No"}}))
## Smoking
rowData(rse_gene)$DE_in_pup_brain_smoking <- unlist(sapply(rowData(rse_gene)$gencodeID, function(x){if(x %in% de_genes_pups_smoking_fitted$gencodeID){"Yes"} else {"No"}}))







