
## 1. Differential Expression Analysis at the exon level

## Only for pup brain samples, and fitted models 
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
load(here("processed-data/03_EDA/03_PCA/rse_exon_brain_pups_qc_afterPCA.Rdata"))

## Simplify rse name
rse_exon_brain_pups_qc <- rse_exon_brain_pups_qc_afterPCA

## Separate samples by Expt
rse_exon_brain_pups_nicotine<-rse_exon_brain_pups_qc[,rse_exon_brain_pups_qc$Expt=="Nicotine"]
rse_exon_brain_pups_smoking<-rse_exon_brain_pups_qc[,rse_exon_brain_pups_qc$Expt=="Smoking"]
save(rse_exon_brain_pups_nicotine, file="processed-data/04_DEA/Exon_analysis/rse_exon_brain_pups_nicotine.Rdata")
save(rse_exon_brain_pups_smoking, file="processed-data/04_DEA/Exon_analysis/rse_exon_brain_pups_smoking.Rdata")



## 1.1 Modeling

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



## Plots for DE exons
plots_DE<-function(top_exons, vExon, FDR=0.05, logFC=0.25, name) {
  ## NS/Down/Upregulated exons
  DE<-vector()
  for (i in 1:dim(top_exons)[1]) {
    ## NS exons are those with FDR>=0.05 and |logFC|<=0.25
    if (top_exons$adj.P.Val[i]>=FDR || abs(top_exons$logFC[i])<=logFC) {
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
  cols <- c("Up" = "red", "Down" = "#26b3ff", "ns" = "grey") 
  sizes <- c("Up" = 2, "Down" = 2, "ns" = 1) 
  alphas <- c("Up" = 1, "Down" = 1, "ns" = 0.3)
  top_exons$mean_log_expr<-apply(vExon$E, 1, mean)
  
  p1<-ggplot(data = top_exons,
             aes(x = mean_log_expr,y = logFC,
                 fill = DE,
                 size = DE,
                 alpha = DE)) +
    geom_point(shape = 21) +
    theme_bw() +
    scale_fill_manual(values = cols, name=NULL) +
    scale_size_manual(values = sizes, name=NULL) +
    scale_alpha_manual(values = alphas, name=NULL) +
    labs(x="Mean of log-normalized counts", y="Log2FC (Exposed vs Ctrl)") +
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          legend.position = c(0.82, 0.15),
          legend.background = element_rect(fill=NA, color='black'),
          legend.key.height = unit(0.15,"cm"),
          axis.title = element_text(size = (10)),
          legend.text = element_text(size=11))
  
  
  ## Volcano plot for DE exons
  p2<-ggplot(data = top_exons,
             aes(x = logFC,y = -log10(adj.P.Val),
                 fill = DE,
                 size = DE,
                 alpha = DE,
                 label= DEG_symbol)) +
    theme_bw() +
    geom_point(shape =21) +
    geom_hline(yintercept = -log10(FDR),
               linetype = "dashed", color = 'gray65', linewidth=0.5) +
    geom_vline(xintercept = c(-1,1),
               linetype = "dashed", color = 'gray65', linewidth=0.5) +
    geom_label_repel(aes(fontface = 'bold', fill=NULL, alpha=NULL),
                     size=2.3,
                     max.overlaps = Inf,
                     min.segment.length = unit(0.2, "cm"),
                     point.padding = unit(0.1, "cm"),
                     box.padding = 0.2,
                     label.padding = 0.2,
                     label.size = 0.2,
                     show.legend=FALSE) +
    labs(y="-log10(FDR)", x="Log2FC (Exposed vs Ctrl)")+
    scale_fill_manual(values = cols, name=NULL) +
    scale_size_manual(values = sizes, name=NULL) +
    scale_alpha_manual(values = alphas, name=NULL) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = c(0.13, 0.85),
          legend.background = element_rect(fill='white', color='black'),
          legend.key.height = unit(0.15,"cm"),
          legend.text = element_text(size=11),
          legend.title = element_text(size = 12),
          axis.title = element_text(size = (12)),
          axis.text = element_text(size = 10)) 
  
  plot_grid(p1, p2, ncol=2)
  ggsave(paste("plots/04_DEA/01_Modeling/Exon_analysis/DE_exons_plots_", name, ".pdf", sep=""), 
         width = 35, height = 15, units = "cm")
}



## Boxplot of a single DE exon
DE_one_boxplot <- function (de_exons, lognorm_DE, exon_ID, exon_name){
  
  ## q-value for the exon
  q_value<-signif(de_exons[which(de_exons$exon_libdID==exon_ID), "adj.P.Val"], digits = 3)
  
  ## FC
  FC<-signif(2**(de_exons[which(de_exons$exon_libdID==exon_ID), "logFC"]), digits=2)
  
  ## Boxplot for the DE exon
  ggplot(data=as.data.frame(lognorm_DE), 
         aes(x=Group,y=eval(parse_expr(exon_ID)))) + 
    ## Hide outliers
    geom_boxplot(outlier.color = "#FFFFFFFF", width=0.35) +
    ## Samples colored by Group + noise
    geom_jitter(aes(colour=Group),shape=16, 
                position=position_jitter(0.2), size=2.1) +
    theme_bw() +
    scale_color_manual(values=c("Control" = "seashell3", "Experimental" = "orange3")) +
    scale_x_discrete(labels=c("Control"="Ctrl","Experimental"="Expt")) +
    labs(x = "Group", y = "lognorm counts",
         title = exon_name, 
         subtitle = paste("FDR:", q_value, '      ', 'FC:', FC)) +
    theme(plot.margin=unit (c(0.4,0.4,0.4,0.4), 'cm'), 
          legend.position = "none",
          plot.title = element_text(hjust=0.5, size=12, face="bold"), 
          plot.subtitle = element_text(size = 10),
          axis.title = element_text(size = (12)),
          axis.text = element_text(size = 10.5)) 
  
}



## Boxplot of lognorm counts for top 3 DE exons
## Obtain lognorm counts of DE exons
DE_boxplots <- function(RSE, vExon, de_exons){
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
  ggsave(here(paste("plots/04_DEA/01_Modeling/Exon_analysis/DE_boxplots_exons_",name, ".pdf", sep="")), 
         width = 28, height = 10, units = "cm")
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
    ## Plots for DE exons
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

## Create csv with results of DE exons
de_exons_nic <- de_exons_nic[,c('exon_gencodeID', 'seqnames', 'start',  'end', 'strand', 'gencodeID', 'Symbol', 'logFC', "t" , "P.Value", "adj.P.Val")]
colnames(de_exons_nic)[7] <- 'gene_Symbol'
de_exons_nic <- de_exons_nic[order(de_exons_nic$adj.P.Val, decreasing = FALSE),]
write.table(de_exons_nic, file = "processed-data/04_DEA/Exon_analysis/de_exons_brain_pup_nicotine.csv", row.names = FALSE, col.names = TRUE, sep = '\t')



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

## Create csv with results of DE exons
de_exons_smo <- de_exons_smo[,c('exon_gencodeID', 'seqnames', 'start',  'end', 'strand', 'gencodeID', 'Symbol', 'logFC', "t" , "P.Value", "adj.P.Val")]
colnames(de_exons_smo)[7] <- 'gene_Symbol' 
de_exons_smo <- de_exons_smo[order(de_exons_smo$adj.P.Val, decreasing = FALSE),]
write.table(de_exons_smo, file = "processed-data/04_DEA/Exon_analysis/de_exons_brain_pup_smoking.csv", row.names = FALSE, col.names = TRUE, sep = '\t')







## Reproducibility information

options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.0 (2023-04-21)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Mexico_City
# date     2023-12-25
# rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
# pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# AnnotationDbi          1.63.2    2023-07-03 [1] Bioconductor
# backports              1.4.1     2021-12-13 [1] CRAN (R 4.3.0)
# base64enc              0.1-3     2015-07-28 [1] CRAN (R 4.3.0)
# Biobase              * 2.61.0    2023-06-02 [1] Bioconductor
# BiocFileCache          2.9.1     2023-07-14 [1] Bioconductor
# BiocGenerics         * 0.47.0    2023-06-02 [1] Bioconductor
# biomaRt                2.57.1    2023-06-14 [1] Bioconductor
# biomartr             * 1.0.7     2023-12-02 [1] CRAN (R 4.3.1)
# Biostrings             2.69.2    2023-07-05 [1] Bioconductor
# bit                    4.0.5     2022-11-15 [1] CRAN (R 4.3.0)
# bit64                  4.0.5     2020-08-30 [1] CRAN (R 4.3.0)
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
# blob                   1.2.4     2023-03-17 [1] CRAN (R 4.3.0)
# cachem                 1.0.8     2023-05-01 [1] CRAN (R 4.3.0)
# cellranger             1.1.0     2016-07-27 [1] CRAN (R 4.3.0)
# checkmate              2.2.0     2023-04-27 [1] CRAN (R 4.3.0)
# cli                    3.6.1     2023-03-23 [1] CRAN (R 4.3.0)
# cluster                2.1.4     2022-08-22 [1] CRAN (R 4.3.0)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# cowplot              * 1.1.1     2020-12-30 [1] CRAN (R 4.3.0)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# curl                   5.0.1     2023-06-07 [1] CRAN (R 4.3.0)
# data.table             1.14.8    2023-02-17 [1] CRAN (R 4.3.0)
# DBI                    1.1.3     2022-06-18 [1] CRAN (R 4.3.0)
# dbplyr                 2.3.3     2023-07-07 [1] CRAN (R 4.3.0)
# DelayedArray           0.26.6    2023-07-02 [1] Bioconductor
# digest                 0.6.33    2023-07-07 [1] CRAN (R 4.3.0)
# dplyr                  1.1.2     2023-04-20 [1] CRAN (R 4.3.0)
# edgeR                * 3.43.7    2023-06-21 [1] Bioconductor
# evaluate               0.21      2023-05-05 [1] CRAN (R 4.3.0)
# fansi                  1.0.5     2023-10-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
# fastmap                1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
# filelock               1.0.2     2018-10-05 [1] CRAN (R 4.3.0)
# foreign                0.8-84    2022-12-06 [1] CRAN (R 4.3.0)
# Formula                1.2-5     2023-02-24 [1] CRAN (R 4.3.0)
# fs                     1.6.3     2023-07-20 [1] CRAN (R 4.3.0)
# gargle                 1.5.2     2023-07-20 [1] CRAN (R 4.3.0)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.37.2    2023-06-21 [1] Bioconductor
# GenomeInfoDbData       1.2.10    2023-05-28 [1] Bioconductor
# GenomicRanges        * 1.53.1    2023-06-02 [1] Bioconductor
# ggplot2              * 3.4.4     2023-10-12 [1] CRAN (R 4.3.1)
# ggrepel              * 0.9.3     2023-02-03 [1] CRAN (R 4.3.0)
# glue                   1.6.2     2022-02-24 [1] CRAN (R 4.3.0)
# googledrive            2.1.1     2023-06-11 [1] CRAN (R 4.3.0)
# gridExtra            * 2.3       2017-09-09 [1] CRAN (R 4.3.0)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.0)
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
# Hmisc                * 5.1-0     2023-05-08 [1] CRAN (R 4.3.0)
# hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
# htmlTable              2.4.1     2022-07-07 [1] CRAN (R 4.3.0)
# htmltools              0.5.5     2023-03-23 [1] CRAN (R 4.3.0)
# htmlwidgets            1.6.2     2023-03-17 [1] CRAN (R 4.3.0)
# httr                   1.4.6     2023-05-08 [1] CRAN (R 4.3.0)
# IRanges              * 2.35.2    2023-06-23 [1] Bioconductor
# jaffelab             * 0.99.32   2023-05-28 [1] Github (LieberInstitute/jaffelab@21e6574)
# KEGGREST               1.41.0    2023-07-07 [1] Bioconductor
# knitr                  1.43      2023-05-25 [1] CRAN (R 4.3.0)
# labeling               0.4.3     2023-08-29 [1] CRAN (R 4.3.0)
# lattice                0.21-8    2023-04-05 [1] CRAN (R 4.3.0)
# lifecycle              1.0.3     2022-10-07 [1] CRAN (R 4.3.0)
# limma                * 3.57.6    2023-06-21 [1] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [1] CRAN (R 4.3.0)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# MASS                   7.3-60    2023-05-04 [1] CRAN (R 4.3.0)
# Matrix                 1.6-0     2023-07-08 [1] CRAN (R 4.3.0)
# MatrixGenerics       * 1.13.0    2023-05-20 [1] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.3.0)
# memoise                2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# nlme                   3.1-162   2023-01-31 [1] CRAN (R 4.3.0)
# nnet                   7.3-19    2023-05-03 [1] CRAN (R 4.3.0)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# png                    0.1-8     2022-11-29 [1] CRAN (R 4.3.0)
# prettyunits            1.1.1     2020-01-24 [1] CRAN (R 4.3.0)
# progress               1.2.2     2019-05-16 [1] CRAN (R 4.3.0)
# purrr                  1.0.1     2023-01-10 [1] CRAN (R 4.3.0)
# R.methodsS3          * 1.8.2     2022-06-13 [1] CRAN (R 4.3.0)
# R.oo                 * 1.25.0    2022-06-12 [1] CRAN (R 4.3.0)
# R.utils              * 2.12.2    2022-11-11 [1] CRAN (R 4.3.0)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.3.0)
# ragg                   1.2.5     2023-01-12 [1] CRAN (R 4.3.0)
# rappdirs               0.3.3     2021-01-31 [1] CRAN (R 4.3.0)
# RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.11    2023-07-06 [1] CRAN (R 4.3.0)
# RCurl                  1.98-1.12 2023-03-27 [1] CRAN (R 4.3.0)
# readxl               * 1.4.3     2023-07-06 [1] CRAN (R 4.3.0)
# rlang                * 1.1.1     2023-04-28 [1] CRAN (R 4.3.0)
# rmarkdown              2.23      2023-07-01 [1] CRAN (R 4.3.0)
# rpart                  4.1.19    2022-10-21 [1] CRAN (R 4.3.0)
# rprojroot              2.0.3     2022-04-02 [1] CRAN (R 4.3.0)
# RSQLite                2.3.1     2023-04-03 [1] CRAN (R 4.3.0)
# rstudioapi             0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
# S4Arrays               1.1.4     2023-06-02 [1] Bioconductor
# S4Vectors            * 0.39.1    2023-06-02 [1] Bioconductor
# scales                 1.2.1     2022-08-20 [1] CRAN (R 4.3.0)
# segmented              1.6-4     2023-04-13 [1] CRAN (R 4.3.0)
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
# stringi                1.7.12    2023-01-11 [1] CRAN (R 4.3.0)
# stringr                1.5.0     2022-12-02 [1] CRAN (R 4.3.0)
# SummarizedExperiment * 1.30.2    2023-06-06 [1] Bioconductor
# systemfonts            1.0.4     2022-02-11 [1] CRAN (R 4.3.0)
# textshaping            0.3.6     2021-10-13 [1] CRAN (R 4.3.0)
# tibble                 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.4     2023-10-12 [1] CRAN (R 4.3.1)
# withr                  2.5.1     2023-09-26 [1] CRAN (R 4.3.1)
# xfun                   0.39      2023-04-20 [1] CRAN (R 4.3.0)
# XML                    3.99-0.14 2023-03-19 [1] CRAN (R 4.3.0)
# xml2                   1.3.5     2023-07-06 [1] CRAN (R 4.3.0)
# XVector                0.41.1    2023-06-02 [1] Bioconductor
# zlibbioc               1.47.0    2023-05-20 [1] Bioconductor
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────


