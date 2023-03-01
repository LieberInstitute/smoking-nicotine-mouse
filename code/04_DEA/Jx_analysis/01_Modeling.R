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
library(biomartr)
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



## Boxplot of a single jxn
DE_one_boxplot <- function (de_jx, lognorm_DE, jx_ID, jx_name){

  ## Real jx ID
  real_jx_ID <- strsplit(jx_name, "  ")[[1]][1]
  ## q-value for the jx
  q_value<-signif(de_jx[which(rownames(de_jx)==real_jx_ID), "adj.P.Val"], digits = 3)
  ## Class of the jx
  class<-de_jx[which(rownames(de_jx)==real_jx_ID), "Class"]
  
  ## Boxplot for each DE jx
  ggplot(data=as.data.frame(lognorm_DE), 
         aes(x=Group,y=eval(parse_expr(jx_ID)))) + 
    ## Hide outliers
    geom_boxplot(outlier.color = "#FFFFFFFF") +
    ## Samples colored by Group + noise
    geom_jitter(aes(colour=Group),shape=16, 
                position=position_jitter(0.2)) +
    theme_classic() +
    labs(x = "Group", y = "norm counts",
         title = jx_name, 
         subtitle = paste(" FDR:", q_value, '\n', "Class:", class)) +
    theme(plot.margin=unit (c (1,1.5,1,1), 'cm'), legend.position = "none",
          plot.title = element_text(hjust=0.5, size=10, face="bold"), 
          plot.subtitle = element_text(size = 9)) 
  
}



## Boxplot of lognorm counts for top 3 DE jxns
## Obtain lognorm counts of DE jxns
DE_boxplots <- function(RSE, de_jx){
  ## Regress out residuals to remove batch effects
  formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
  model=model.matrix(formula, data=colData(RSE))
  logcounts<-cleaningY(assays(RSE)$logcounts, model, P=2)
  ## Order jxns by q-value
  de_jx<-de_jx[order(de_jx$adj.P.Val),]
  ## Lognorm counts of DE jxs
  lognorm_DE<-logcounts[rownames(de_jx),]
  ## Samples as rows and jxns as columns
  lognorm_DE<-t(lognorm_DE)
  ## Add samples' Group information
  lognorm_DE<-data.frame(lognorm_DE, Group=colData(RSE)$Group)

  plots<-list()
  for (i in 1:3){
    symbol<-biomart(genes  = strsplit(de_jx$newGeneID[i], "[.]")[[1]][1],
                     mart       = "ENSEMBL_MART_ENSEMBL",
                     dataset    = "mmusculus_gene_ensembl",
                     attributes = c("external_gene_name"),
                     filters    = "ensembl_gene_id")
    jx_ID<-colnames(lognorm_DE)[i]
    jx_name<-paste(rownames(de_jx)[i], symbol[2], sep="  ")
    p<-DE_one_boxplot(de_jx, lognorm_DE, jx_ID, jx_name)
    plots[[i]]<-print(p)
  }
  plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol = 3)
  ggsave(here(paste("plots/04_DEA/01_Modeling/Jx_analysis/DE_boxplots_jx_", name, ".pdf", sep="")), 
         width = 35, height = 13, units = "cm")
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
    DE_boxplots(RSE, de_jxns)
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
# date     2023-02-27
# rstudio  2022.12.0+353 Elsbeth Geranium (desktop)
# pandoc   NA






