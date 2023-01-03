
## 1. Differential Expression Analysis at the transcript level

## Only for brain and pup samples, and fitted models 
## (Based on EDA and DEA at the gene level)

library(here)
library(SummarizedExperiment)
library(stats)
library(edgeR)
library(limma)
library(ggplot2)
library(rlang)
library(cowplot)
library(ggrepel)
library(jaffelab)
library(VennDiagram) 
library(gridExtra)
library(R.utils)
library(sessioninfo)

load(here("raw-data/rse_tx_smoking_mouse_n208.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_tx_brain_pups_qc.Rdata"))

## Separate samples by Expt
rse_tx_brain_pups_nicotine<-rse_tx_brain_pups_qc[,rse_tx_brain_pups_qc$Expt=="Nicotine"]
rse_tx_brain_pups_smoking<-rse_tx_brain_pups_qc[,rse_tx_brain_pups_qc$Expt=="Smoking"]
save(rse_tx_brain_pups_nicotine, file="processed-data/04_DEA/Tx_analysis/rse_tx_brain_pups_nicotine.Rdata")
save(rse_tx_brain_pups_smoking, file="processed-data/04_DEA/Tx_analysis/rse_tx_brain_pups_smoking.Rdata")



## 1.1 Modelling

## Extract previous output from calcNormFactors 
norm_factors<-calcNormFactors(rse_tx, method = "TMM")
samples_factors<-data.frame(SAMPLE_ID=norm_factors$samples$SAMPLE_ID,
                            norm.factors=norm_factors$samples$norm.factors,
                            lib.size=norm_factors$samples$lib.size)

