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
  
  pdf(file = paste("plots/04_DEA/01_Modeling/Jx_analysis/DEA_jx_plots_", name, ".pdf", sep="" ))
  par(mfrow=c(2,2))
  
  ## Model matrix using formula for the fitted model
  formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
  model=model.matrix(formula, data=colData(RSE))
  ## Use previous norm factors
  RSE_scaled = calcNormFactors(RSE)
  RSE_scaled$samples$lib.size<-factors$lib.size
  RSE_scaled$samples$norm.factors<-factors$norm.factors
  
  ## Transform counts to log2(CPM) 
  ## Estimate mean-variance relationship for each jxn 
  vJxn = voom(RSE_scaled, design=model, plot=TRUE)
  
  ## Fit linear model for each jxn 
  fitJxn = lmFit(vJxn)
  
  ## Compute moderated F and t-statistics, and log-odds of DE
  eBJxn = eBayes(fitJxn)
  
  ## Plot average log expression vs logFC
  limma::plotMA(eBJxn, coef = "GroupExperimental", xlab = "Mean of normalized counts", 
                ylab="logFC")
  
  ## Plot -log(p-value) vs logFC
  volcanoplot(eBJxn, coef = "GroupExperimental")
  
  ## Select top-ranked jxns for Group 
  top_jxns = topTable(eBJxn, coef="GroupExperimental", p.value = 1, number=nrow(RSE), sort.by="none")
  ## Histogram of adjusted p values
  hist(top_jxns$adj.P.Val, xlab="FDR", main="")
  
  dev.off()
  
  return(list(top_jxns, vJxn, eBJxn))
  
} 



## Perform DEA for each group of samples
apply_DEA<-function(RSE, name){
  ## DEA
  results<-DEA_expt_vs_ctl(RSE, name)
  top_jxns<-results[[1]]
  
  ## If there are DE jxns
  if (length(which(top_jxns$adj.P.Val<0.05))>0){
    ## Signif jxns
    de_jxns<-top_jxns[which(top_jxns$adj.P.Val < 0.05),]
    print(paste(dim(de_jxns)[1], "differentially expressed jxns", sep=" "))
    return(list(results, de_jxns))
  }
  else {
    print("No differentially expressed jxns")
    return(results)
  }
}





##################################
#         Nicotine DEA 
#    Nicotine vs ctrls in pups
##################################

RSE<-rse_jx_brain_pups_nicotine
name<-"nicotine"
results_nic<-apply_DEA(RSE, name)
"205 differentially expressed jxns"
top_jxns_nic<-results_nic[[1]][[1]]
de_jxns_nic<-results_nic[[2]]
save(results_nic, file="processed-data/04_DEA/Jx_analysis/results_nic.Rdata")
save(top_jxns_nic, file="processed-data/04_DEA/Jx_analysis/top_jxns_nic.Rdata")
save(de_jxns_nic, file="processed-data/04_DEA/Jx_analysis/de_jxns_nic.Rdata")



##################################
#          Smoking DEA 
#    Smoking vs ctrls in pups
##################################

RSE<-rse_jx_brain_pups_smoking
name<-"smoking"
results_smo<-apply_DEA(RSE, name)
"9515 differentially expressed jxns"
top_jxns_smo<-results_smo[[1]][[1]]
de_jxns_smo<-results_smo[[2]]
save(results_smo, file="processed-data/04_DEA/Jx_analysis/results_smo.Rdata")
save(top_jxns_smo, file="processed-data/04_DEA/Jx_analysis/top_jxns_smo.Rdata")
save(de_jxns_smo, file="processed-data/04_DEA/Jx_analysis/de_jxns_smo.Rdata")








