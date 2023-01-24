## 1. Differential Expression Analysis at the junction (jx) level

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

load(here("raw-data/rse_jx_smoking_mouse_n208.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_jx_brain_pups_qc.Rdata"))

## Separate samples by Expt
rse_jx_brain_pups_nicotine<-rse_jx_brain_pups_qc[,rse_jx_brain_pups_qc$Expt=="Nicotine"]
rse_jx_brain_pups_smoking<-rse_jx_brain_pups_qc[,rse_jx_brain_pups_qc$Expt=="Smoking"]
save(rse_jx_brain_pups_nicotine, file="processed-data/04_DEA/Jx_analysis/rse_jx_brain_pups_nicotine.Rdata")
save(rse_jx_brain_pups_smoking, file="processed-data/04_DEA/Jx_analysis/rse_jx_brain_pups_smoking.Rdata")



## 1.1 Modeling

## Extract previous output from calcNormFactors 
norm_factors<-calcNormFactors(rse_jx, method = "TMMwsp")
samples_factors<-data.frame(SAMPLE_ID=norm_factors$samples$SAMPLE_ID,
                            norm.factors=norm_factors$samples$norm.factors,
                            lib.size=norm_factors$samples$lib.size)


## DEA for experiment vs ctrls
DEA_expt_vs_ctl<- function(RSE, name){
  
  ## Previous lib sizes of each sample
  match_samples <- match(RSE$SAMPLE_ID, samples_factors$SAMPLE_ID)
  stopifnot(all(!is.na(match_samples)))
  factors<-samples_factors[match_samples, ]
  
  pdf(file = paste("plots/04_DEA/01_Modeling/Exon_analysis/DEA_exon_plots_", name, ".pdf", sep="" ))
  par(mfrow=c(2,2))
  
  ## Model matrix using formula for the fitted model
  formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
  model=model.matrix(formula, data=colData(RSE))
  ## Use previous norm factors
  RSE_scaled = calcNormFactors(RSE)
  RSE_scaled$samples$lib.size<-factors$lib.size
  RSE_scaled$samples$norm.factors<-factors$norm.factors
  
  ## Transform counts to log2(CPM) 
  ## Estimate mean-variance relationship for each exon
  vExon = voom(RSE_scaled, design=model, plot=TRUE)
  
  ## Fit linear model for each exon
  fitExon = lmFit(vExon)
  
  ## Compute moderated F and t-statistics, and log-odds of DE
  eBExon = eBayes(fitExon)
  
  ## Plot average log expression vs logFC
  limma::plotMA(eBExon, coef = "GroupExperimental", xlab = "Mean of normalized counts", 
                ylab="logFC")
  
  ## Plot -log(p-value) vs logFC
  volcanoplot(eBExon, coef = "GroupExperimental")
  
  ## Select top-ranked exons for Group 
  top_exons = topTable(eBExon, coef="GroupExperimental", p.value = 1, number=nrow(RSE), sort.by="none")
  ## Histogram of adjusted p values
  hist(top_exons$adj.P.Val, xlab="FDR", main="")
  
  dev.off()
  
  return(list(top_exons, vExon, eBExon))
  
}  

