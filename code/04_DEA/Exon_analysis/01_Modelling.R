
## 1. Differential Expression Analysis at the exon level

## Only for brain and pup samples, and fitted models 
## (Based on EDA and DEA at the gene level)

library(here)
library(SummarizedExperiment)
library(edgeR)
library(sessioninfo)
##CHECK MISSING LIBRARIES

load(here("raw-data/rse_exon_smoking_mouse_n208.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_exon_brain_pups_qc.Rdata"))

## Separate samples by Expt
rse_exon_brain_pups_nicotine<-rse_exon_brain_pups_qc[,rse_exon_brain_pups_qc$Expt=="Nicotine"]
rse_exon_brain_pups_smoking<-rse_exon_brain_pups_qc[,rse_exon_brain_pups_qc$Expt=="Smoking"]



## 1.1 Modelling

## Extract previous output from calcNormFactors 
norm_factors<-calcNormFactors(rse_exon, method = "TMM")
samples_factors<-data.frame(SAMPLE_ID=norm_factors$samples$SAMPLE_ID,
                            norm.factors=norm_factors$samples$norm.factors,
                            lib.size=norm_factors$samples$lib.size)


## DEA for experiment vs ctrls
DEA_expt_vs_ctl<- function(RSE, name){
  
  ## Previous lib sizes of each sample
  match_samples <- match(RSE$SAMPLE_ID, samples_factors$SAMPLE_ID)
  stopifnot(all(!is.na(match_samples)))
  factors<-samples_factors[match_samples, ]
  
  pdf(file = paste("plots/04_DEA/01_Modelling/Exon_analysis/DEA_exon_plots_", name, ".pdf", sep="" ))
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
  
  ## compute F and t-statistics, p-value and logFC 
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



## Plots for DE exons
plots_DE<-function(top_exons, vExon, FDR=0.05, name) {
  ## NS/Down/Upregulated exons
  DE<-vector()
  for (i in 1:dim(top_exons)[1]) {
    if (top_exons$adj.P.Val[i]>FDR) {
      DE<-append(DE, "ns")
    }
    else {
      if (top_exons$logFC[i]>0) {
        DE<-append(DE, "Up")
      }
      else {
        DE<-append(DE, "Down")
      }
    }
  }
  top_exons$DE<- DE
  
  ## Gene symbols for DE exons with |logFC|>1 ??????????
  DEG_symbol<-vector()
  for (i in 1:dim(top_exons)[1]) {
    if (top_exons$DE[i]!="ns" & abs(top_exons$logFC[i])>1) {
      DEG_symbol<-append(DEG_symbol, top_exons$Symbol[i])
    }
    else {
      DEG_symbol<-append(DEG_symbol, NA)
    }
  }
  top_exons$DEG_symbol<- DEG_symbol
  
  
  ## MA plot for DE exons
  cols <- c("Up" = "#ffad73", "Down" = "#26b3ff", "ns" = "grey") 
  sizes <- c("Up" = 2, "Down" = 2, "ns" = 1) 
  alphas <- c("Up" = 1, "Down" = 1, "ns" = 0.5)
  top_exons$mean_log_expr<-apply(vExon$E, 1, mean)
  p1<-ggplot(data = top_exons, 
             aes(x = mean_log_expr,y = logFC,
                 fill = DE,    
                 size = DE,
                 alpha = DE)) + 
    geom_point(shape = 21,    
               colour = "black") +
    scale_fill_manual(values = cols) + 
    scale_size_manual(values = sizes) + 
    scale_alpha_manual(values = alphas) +
    labs(x="Mean of normalized counts")
  
  
  ## Volcano plot for DE exons
  p2<-ggplot(data = top_exons, 
             aes(x = logFC,y = -log10(adj.P.Val),
                 fill = DE,    
                 size = DE,
                 alpha = DE,
                 label= DEG_symbol)) +
    geom_point(shape = 21) + 
    geom_hline(yintercept = -log10(FDR),
               linetype = "dashed") + 
    geom_vline(xintercept = c(-1,1),
               linetype = "dashed") +
    geom_label_repel(fill="white", size=2, max.overlaps = Inf,  
                     box.padding = 0.2, 
                     show.legend=FALSE) +
    labs(y="-log10(FDR)")+
    scale_fill_manual(values = cols) + 
    scale_size_manual(values = sizes) + 
    scale_alpha_manual(values = alphas) 
  
  plot_grid(p1, p2, ncol=2)
  ggsave(paste("plots/04_DEA/01_Modelling/Gene_analysis/DEG_exon_plots_", name, ".pdf", sep=""), 
         width = 35, height = 15, units = "cm")
}



##################################
#         Nicotine DEA 
#    Nicotine vs ctrls in pups
##################################

RSE<-rse_exon_brain_pups_nicotine

## Naive model
name<-"nicotine"
results<-DEA_expt_vs_ctl(RSE, name)
top_exons<-results[[1]]
vExon<-results[[2]]
plots_DE(top_exons, vExon, name = name)




##################################
#          Smoking DEA 
#    Smoking vs ctrls in pups
##################################

RSE<-rse_exon_brain_pups_smoking

## Naive model
name<-"smoking"
results<-DEA_expt_vs_ctl(RSE, name)
top_exons<-results[[1]]
vExon<-results[[2]]
plots_DE(top_exons, vExon, name = name)







