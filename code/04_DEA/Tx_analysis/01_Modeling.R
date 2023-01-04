
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



## 1.1 Modeling

## DEA for experiment vs ctrls
DEA_expt_vs_ctl<- function(RSE, name){
  
  pdf(file = paste("plots/04_DEA/01_Modeling/Tx_analysis/DEA_tx_plots_", name, ".pdf", sep="" ))
  par(mfrow=c(2,2))
  
  ## Model matrix using formula for the fitted model
  formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
  model=model.matrix(formula, data=colData(RSE))

  # ## Transform tpm to log2(CPM) 
  # ## Estimate mean-variance relationship for each transcript
  # vExon = voom(assays(RSE)$tpm, design=model, plot=TRUE)
  
  ## Fit linear model for each transcript
  fitTx = lmFit(assays(RSE)$logcounts, design = model)
  
  ## Compute moderated F and t-statistics, and log-odds of DE
  eBTx = eBayes(fitTx)
  
  ## Plot average log expression vs logFC
  limma::plotMA(eBTx, coef = "GroupExperimental", xlab = "Mean of normalized counts", 
                ylab="logFC")
  
  ## Plot -log(p-value) vs logFC
  volcanoplot(eBTx, coef = "GroupExperimental")
  
  ## Select top-ranked exons for Group 
  top_exons = topTable(eBExon, coef="GroupExperimental", p.value = 1, number=nrow(RSE), sort.by="none")
  ## Histogram of adjusted p values
  hist(top_exons$adj.P.Val, xlab="FDR", main="")
  
  dev.off()
  
  return(list(top_exons, vExon, eBExon))
  
}  




