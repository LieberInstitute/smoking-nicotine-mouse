
## 1. Differential Expression Analysis at the transcript (tx) level

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
  
  ## Fit linear model for each transcript
  fitTx = lmFit(assays(RSE)$logcounts, design = model)
  
  ## Compute moderated F and t-statistics, and log-odds of DE
  eBTx = eBayes(fitTx)
  
  ## Plot average log expression vs logFC
  limma::plotMA(eBTx, coef = "GroupExperimental", xlab = "Mean of normalized counts", 
                ylab="logFC")
  
  ## Plot -log(p-value) vs logFC
  volcanoplot(eBTx, coef = "GroupExperimental")
  
  ## Select top-ranked transcripts for Group 
  top_tx = topTable(eBTx, coef="GroupExperimental", p.value = 1, number=nrow(RSE), sort.by="none")
  ## Histogram of adjusted p values
  hist(top_tx$adj.P.Val, xlab="FDR", main="")
  
  dev.off()
  
  ## Add relevant info 
  top_tx$Symbol <- rowData(RSE)$gene_name
  top_tx$ensembl_id <- rowData(RSE)$gene_id
  top_tx$transcript_id<-rowData(RSE)$transcript_id
  top_tx$transcript_name<-rowData(RSE)$transcript_name
  
  return(list(top_tx, eBTx))
  
}  



## Plots for DE transcripts
plots_DE<-function(top_tx, RSE, FDR=0.05, name) {
  
  ## NS/Down/Upregulated transcripts
  DE<-vector()
  for (i in 1:dim(top_tx)[1]) {
    if (top_tx$adj.P.Val[i]>FDR) {
      DE<-append(DE, "ns")
    }
    else {
      if (top_tx$logFC[i]>0) {
        DE<-append(DE, "Up")
      }
      else {
        DE<-append(DE, "Down")
      }
    }
  }
  top_tx$DE<- DE
  
  ## Gene symbols for DE tx with |logFC|>1
  DEtx_symbol<-vector()
  for (i in 1:dim(top_tx)[1]) {
    if (top_tx$DE[i]!="ns" & abs(top_tx$logFC[i])>1) {
      DEtx_symbol<-append(DEtx_symbol, top_tx$Symbol[i])
    }
    else {
      DEtx_symbol<-append(DEtx_symbol, NA)
    }
  }
  top_tx$DEtx_symbol<- DEtx_symbol
  
  
  ## MA plot for DE trascripts
  cols <- c("Up" = "#ffad73", "Down" = "#26b3ff", "ns" = "grey") 
  sizes <- c("Up" = 2, "Down" = 2, "ns" = 1) 
  alphas <- c("Up" = 1, "Down" = 1, "ns" = 0.5)
  top_tx$mean_log_expr<-apply(assays(RSE)$logcounts, 1, mean)
  p1<-ggplot(data = top_tx, 
             aes(x = mean_log_expr, y = logFC,
                 fill = DE,    
                 size = DE,
                 alpha = DE)) + 
    geom_point(shape = 21,    
               colour = "black") +
    scale_fill_manual(values = cols) + 
    scale_size_manual(values = sizes) + 
    scale_alpha_manual(values = alphas) +
    labs(x="Mean of normalized counts")
  
  
  ## Volcano plot for DE tx
  p2<-ggplot(data = top_tx, 
             aes(x = logFC,y = -log10(adj.P.Val),
                 fill = DE,    
                 size = DE,
                 alpha = DE,
                 label= DEtx_symbol)) +
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
  ggsave(paste("plots/04_DEA/01_Modeling/Tx_analysis/DEtx_plots_", name, ".pdf", sep=""), 
         width = 35, height = 15, units = "cm")
}



## Boxplot of a single transcript
DE_one_boxplot <- function (de_tx, lognorm_DE, tx_ID, tx_name){
  
  ## q-value for the tx
  q_value<-signif(de_tx[which(de_tx$transcript_id==tx_ID), "adj.P.Val"], digits = 3)
  
  ## Boxplot for each DE tx
  ggplot(data=as.data.frame(lognorm_DE), 
         aes(x=Group,y=eval(parse_expr(tx_ID)))) + 
    ## Hide outliers
    geom_boxplot(outlier.color = "#FFFFFFFF") +
    ## Samples colored by Group + noise
    geom_jitter(aes(colour=Group),shape=16, 
                position=position_jitter(0.2)) +
    theme_classic() +
    labs(x = "Group", y = "norm counts",
         title = tx_name, 
         subtitle = paste("FDR:", q_value)) +
    theme(plot.margin=unit (c (1,1.5,1,1), 'cm'), legend.position = "none",
          plot.title = element_text(hjust=0.5, size=10, face="bold"), 
          plot.subtitle = element_text(size = 9)) 
  
}



## Boxplot of lognorm counts for top 3 DE tx
## Obtain lognorm counts of DE tx
DE_boxplots <- function(RSE, de_tx){
  ## Regress out residuals to remove batch effects
  formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
  model=model.matrix(formula, data=colData(RSE))
  logcounts<-cleaningY(assays(RSE)$logcounts, model, P=2)
  ## Order tx by q-value
  de_tx<-de_tx[order(de_tx$adj.P.Val),]
  ## Lognorm counts of DE tx
  lognorm_DE<-logcounts[rownames(de_tx),]
  ## Samples as rows and tx as columns
  lognorm_DE<-t(lognorm_DE)
  ## Add samples' Group information
  lognorm_DE<-data.frame(lognorm_DE, "Group"=colData(RSE)$Group)
  
  
  plots<-list()
  for (i in 1:3){
    tx_ID<-colnames(lognorm_DE)[i]
    tx_name<-paste(de_tx$Symbol[i]," - ", de_tx$transcript_id[i], sep="")
    p<-DE_one_boxplot(de_tx, lognorm_DE, tx_ID, tx_name)
    plots[[i]]<-p
  }
  plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol = 3)
  ggsave(here(paste("plots/04_DEA/01_Modeling/Tx_analysis/DE_boxplots_tx_",name, ".pdf", sep="")), 
         width = 25, height = 10, units = "cm")
}



## Perform DEA for each group of samples
apply_DEA<-function(RSE, name){
  ## DEA
  results<-DEA_expt_vs_ctl(RSE, name)
  top_tx<-results[[1]]
  
  ## If there are DE tx
  if (length(which(top_tx$adj.P.Val<0.05))>0){
    model=model.matrix(formula, data=colData(RSE))
    de_tx<-top_tx[top_tx$adj.P.Val < 0.05,]
    ## Plots for DE tx
    plots_DE(top_tx, RSE, 0.05, name)
    DE_boxplots(RSE, de_tx)
    print(paste(length(which(top_tx$adj.P.Val<0.05)), "differentially expressed tx", sep=" "))
    return(list(results, de_tx))
  }
  else {
    print("No differentially expressed tx")
    return(results)
  }
}





##################################
#         Nicotine DEA 
#    Nicotine vs ctrls in pups
##################################

RSE<-rse_tx_brain_pups_nicotine
name<-"nicotine"
results_nic<-apply_DEA(RSE, name)
"232 differentially expressed tx"
top_tx_nic<-results_nic[[1]][[1]]
de_tx_nic<-results_nic[[2]]
save(results_nic, file="processed-data/04_DEA/Tx_analysis/results_nic.Rdata")
save(top_tx_nic, file="processed-data/04_DEA/Tx_analysis/top_tx_nic.Rdata")
save(de_tx_nic, file="processed-data/04_DEA/Tx_analysis/de_tx_nic.Rdata")



##################################
#         Smoking DEA 
#    Smoking vs ctrls in pups
##################################

RSE<-rse_tx_brain_pups_smoking
name<-"smoking"
results_smo<-apply_DEA(RSE, name)
"4059 differentially expressed tx"
top_tx_smo<-results_smo[[1]][[1]]
de_tx_smo<-results_smo[[2]]
save(results_smo, file="processed-data/04_DEA/Tx_analysis/results_smo.Rdata")
save(top_tx_smo, file="processed-data/04_DEA/Tx_analysis/top_tx_smo.Rdata")
save(de_tx_smo, file="processed-data/04_DEA/Tx_analysis/de_tx_smo.Rdata")







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
# date     2023-01-05
# rstudio  2022.07.2+576 Spotted Wakerobin (desktop)
# pandoc   2.19.2 @ C:/Program Files/RStudio/bin/quarto/bin/tools/ (via rmarkdown)


