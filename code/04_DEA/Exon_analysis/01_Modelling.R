
library(here)


## 1. Differential Expression Analysis at the exon level

## Only for brain and pups samples


load(here("raw-data/rse_exon_smoking_mouse_n208.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_exon_brain_pups_qc.Rdata"))



## 1.1 Modelling

## Extract previous output from calcNormFactors for all samples
norm_factors<-calcNormFactors(rse_exon, method = "TMM")
samples_factors<-data.frame(SAMPLE_ID=norm_factors$samples$SAMPLE_ID,
                            norm.factors=norm_factors$samples$norm.factors,
                            lib.size=norm_factors$samples$lib.size)


## DEA for experiment mice (nicotine/smoking) vs ctrls
DEA_expt_vs_ctl<- function(RSE, formula, name, coef){
  
  ## Previous lib sizes of each sample
  match_samples <- match(RSE$SAMPLE_ID, samples_factors$SAMPLE_ID)
  stopifnot(all(!is.na(m)))
  factors<-samples_factors[match_samples, ]
  
  pdf(file = paste("plots/04_DEA/01_Modelling/DEA_plots_", name, ".pdf", sep="" ))
  par(mfrow=c(2,2))
  
  ## Model matrix
  model=model.matrix(formula, data=colData(RSE))
  ## Use previous norm factors to scale the raw library sizes
  RSE_scaled = calcNormFactors(RSE)
  RSE_scaled$samples$lib.size<-factors$lib.size
  RSE_scaled$samples$norm.factors<-factors$norm.factors
  
  ## Transform counts to log2(CPM) 
  ## Estimate mean-variance relationship for each gene
  vGene = voom(RSE_scaled, design=model, plot=TRUE)
  
  ## Fit linear model for each gene 
  fitGene = lmFit(vGene)
  
  ## Empirical Bayesian calculation to obtain our significant genes:
  ## compute F and t-statistics, p-value and logFC for DE 
  eBGene = eBayes(fitGene)
  
  ## Plot average log expression vs logFC
  limma::plotMA(eBGene, coef = coef, xlab = "Mean of normalized counts", 
                ylab="logFC")
  ## Plot -log(p-value) vs logFC
  volcanoplot(eBGene, coef = coef)
  
  ## Select top-ranked genes for Group (expt vs ctrl)
  top_genes = topTable(eBGene, coef=coef, p.value = 1, number=nrow(RSE), sort.by="none")
  ## Histogram of adjusted p values
  hist(top_genes$adj.P.Val, xlab="FDR", main="")
  
  dev.off()
  
  return(list(top_genes, vGene, eBGene))
  
}  

