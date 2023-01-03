
## 1. Differential Expression Analysis at the exon level

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

load(here("raw-data/rse_exon_smoking_mouse_n208.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_exon_brain_pups_qc.Rdata"))

## Separate samples by Expt
rse_exon_brain_pups_nicotine<-rse_exon_brain_pups_qc[,rse_exon_brain_pups_qc$Expt=="Nicotine"]
rse_exon_brain_pups_smoking<-rse_exon_brain_pups_qc[,rse_exon_brain_pups_qc$Expt=="Smoking"]
save(rse_exon_brain_pups_nicotine, file="processed-data/04_DEA/Exon_analysis/rse_exon_brain_pups_nicotine.Rdata")
save(rse_exon_brain_pups_smoking, file="processed-data/04_DEA/Exon_analysis/rse_exon_brain_pups_smoking.Rdata")



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
  
  ## Compute F and t-statistics, p-value and logFC 
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
plots_DE<-function(top_exons, vExon, FDR=0.05, logFC=0.25, name) {
  ## NS/Down/Upregulated exons
  DE<-vector()
  for (i in 1:dim(top_exons)[1]) {
    ## NS exons are those with p-values>0.05 and |logFC|<0.25
    if (top_exons$adj.P.Val[i]>FDR || abs(top_exons$logFC[i])<logFC) {
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
  
  ## Gene symbols of DE exons with |logFC|>1 
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
  ggsave(paste("plots/04_DEA/01_Modelling/Exon_analysis/DEG_exon_plots_", name, ".pdf", sep=""), 
         width = 35, height = 15, units = "cm")
}



## Boxplot of a single DE exon
DE_one_boxplot <- function (de_exons, lognorm_DE, exon_ID, exon_name){
  
  ## q-value for the exon
  q_value<-signif(de_exons[which(de_exons$exon_libdID==exon_ID), "adj.P.Val"], digits = 3)
  
  ## Boxplot for the DE exon
  ggplot(data=as.data.frame(lognorm_DE), 
         aes(x=Group,y=eval(parse_expr(exon_ID)))) + 
    ## Hide outliers
    geom_boxplot(outlier.color = "#FFFFFFFF") +
    ## Samples colored by Group + noise
    geom_jitter(aes(colour=Group),shape=16, 
                position=position_jitter(0.2)) +
    theme_classic() +
    labs(x = "Group", y = "norm counts",
         title = exon_name,
         subtitle = paste("FDR:", q_value)) +
    theme(plot.margin=unit (c (1,1.5,1,1), 'cm'), legend.position = "none",
          plot.title = element_text(hjust=0.5, size=10, face="bold"), 
          plot.subtitle = element_text(size = 9)) 
  
}



## Boxplot of lognorm counts for top 3 DE exons
## Obtain lognorm counts of DE exons
DE_boxplots <- function(RSE, vExon, de_exons){
  ## Regress out residuals to remove batch effects
  formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
  model=model.matrix(formula, data=colData(RSE))
  vExon$E<-cleaningY(vExon$E, model, P=2)
  ## Order exons by q-value
  de_exons<-de_exons[order(de_exons$adj.P.Val),]
  lognorm_DE<-vExon$E[rownames(de_exons),]
  ## Samples as rows and exons as columns
  lognorm_DE<-t(lognorm_DE)
  ## Exons' IDs as colnames
  colnames(lognorm_DE)<-de_exons$exon_libdID
  ## Add samples' Group information
  lognorm_DE<-data.frame(lognorm_DE, "Group"=colData(RSE)$Group)
  
  
  plots<-list()
  for (i in 1:3){
    exon_ID<-colnames(lognorm_DE)[i]
    exon_name<-paste(de_exons$Symbol[i]," - ", de_exons$seqnames[i], ":", de_exons$start[i], "-", de_exons$end[i], sep="")
    p<-DE_one_boxplot(de_exons, lognorm_DE, exon_ID, exon_name)
    plots[[i]]<-p
  }
  plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol = 3)
  ggsave(here(paste("plots/04_DEA/01_Modelling/Exon_analysis/DE_boxplots_exons_",name, ".pdf", sep="")), 
         width = 25, height = 10, units = "cm")
}



## Perform DEA for each group of samples
apply_DEA<-function(RSE, name){
  ## DEA
  results<-DEA_expt_vs_ctl(RSE, name)
  top_exons<-results[[1]]
  
  ## If there are DE exons
  if (length(which(top_exons$adj.P.Val<0.05))>0){
    vExon<-results[[2]]
    ## Retain only DE exons with |logFC|>0.25
    de_exons<-top_exons[which(top_exons$adj.P.Val < 0.05 & abs(top_exons$logFC)>0.25),]
    ## Plots for DE genes
    plots_DE(top_exons, vExon, 0.05, 0.25, name)
    DE_boxplots(RSE, vExon, de_exons)
    print(paste(dim(de_exons)[1], "differentially expressed exons", sep=" "))
    return(list(results, de_exons))
  }
  else {
    print("No differentially expressed exons")
    return(results)
  }
}





##################################
#         Nicotine DEA 
#    Nicotine vs ctrls in pups
##################################

RSE<-rse_exon_brain_pups_nicotine
name<-"nicotine"
results_nic<-apply_DEA(RSE, name)
"1115 differentially expressed exons"
top_exons_nic<-results_nic[[1]][[1]]
de_exons_nic<-results_nic[[2]]
save(results_nic, file="processed-data/04_DEA/Exon_analysis/results_nic.Rdata")
save(top_exons_nic, file="processed-data/04_DEA/Exon_analysis/top_exons_nic.Rdata")
save(de_exons_nic, file="processed-data/04_DEA/Exon_analysis/de_exons_nic.Rdata")



##################################
#          Smoking DEA 
#    Smoking vs ctrls in pups
##################################

RSE<-rse_exon_brain_pups_smoking
name<-"smoking"
results_smo<-apply_DEA(RSE, name)
"5983 differentially expressed exons"
top_exons_smo<-results_smo[[1]][[1]]
de_exons_smo<-results_smo[[2]]
save(results_smo, file="processed-data/04_DEA/Exon_analysis/results_smo.Rdata")
save(top_exons_smo, file="processed-data/04_DEA/Exon_analysis/top_exons_smo.Rdata")
save(de_exons_smo, file="processed-data/04_DEA/Exon_analysis/de_exons_smo.Rdata")





## Reproducibility information

options(width = 120)
session_info()

# setting  value
# version  R version 4.2.0 (2022-04-22 ucrt)
# os       Windows 10 x64 (build 19044)
# system   x86_64, mingw32
# ui       RStudio
# language (EN)
# collate  Spanish_Mexico.utf8
# ctype    Spanish_Mexico.utf8
# tz       America/Mexico_City
# date     2022-11-18
# rstudio  2022.07.2+576 Spotted Wakerobin (desktop)
# pandoc   2.19.2 @ C:/Program Files/RStudio/bin/quarto/bin/tools/ (via rmarkdown)




