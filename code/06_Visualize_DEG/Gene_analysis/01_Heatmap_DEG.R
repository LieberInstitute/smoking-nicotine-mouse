
# 1. Visualize differentially expressed genes

library(here)
library(SummarizedExperiment)
library(stats)
library(pheatmap)
library(rlang)
library(sessioninfo)


load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_nicotine.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_smoking.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/results_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/results_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/05_GO_KEGG/Gene_analysis/all_DEG.Rdata"))
load(here("processed-data/05_GO_KEGG/Gene_analysis/intersections.Rdata"))


## (For fitted models only)


## Manually determine coloring for plot annotation of all heatmaps
ann_colors = list()
ann_colors[["Sex"]]=c("F"="hotpink1", "M"="dodgerblue")
ann_colors[["Group"]]=c("Ctrl"="seashell3", "Expt"="orange3")
ann_colors[["FDR nic"]]=c("0.01-0.05"="#9900CC", "<0.01"="#CC99FF")
ann_colors[["FDR smo"]]=c("0.01-0.05"="#FF33FF", "<0.01"="#FF99FF")
ann_colors[["FC smo"]]=c(">=1"="#00CC66", "<1"="#CCFF99")
ann_colors[["FC nic"]]=c(">=1"="#00CC66", "<1"="#CCFF99")
ann_colors[["Significance"]]=c("Nic only"="navajowhite2", "Smo only"="thistle3", "Both"="deeppink3")
ann_colors[["Expt & Group"]]=c("Nicotine Control"="lightblue3", "Nicotine Experimental"="skyblue3", 
                               "Smoking Control"="lightsalmon1", "Smoking Experimental"="salmon3")




## 1.1 Heatmaps of all DEG in Nicotine or Smoking separately

DEG_heatmaps<- function(rse, results, de_genes, filename){
  
  if (filename=="nicotine"){
    name="nic"
  }
  else {
    name="smo"
  }
  
  ## Extract lognorm counts of all genes
  vGene<-results[[1]][[2]]
  DEG<-de_genes$ensemblID
  rownames(vGene)<-vGene$genes$ensemblID
  colnames(vGene)<-vGene$targets$SAMPLE_ID
  
  ## Retain lognorm counts of DEG only
  vGene_DEG <-vGene$E[which(vGene$genes$ensemblID %in% DEG), ]
  
  ## Center the data to make differences more evident
  vGene_DEG<-(vGene_DEG-rowMeans(vGene_DEG))/rowSds(vGene_DEG)
  
  
  ## Samples' info
  ann_cols <- as.data.frame(vGene$targets[, c("Group", "Sex")])
  rownames(ann_cols) <- vGene$targets$SAMPLE_ID
  colnames(ann_cols) <- c("Group", "Sex")
  ann_cols$Group <- gsub('Experimental', 'Expt', gsub('Control', 'Ctrl', ann_cols$Group))
  
  ## Genes' info
  ## FDRs
  FDRs<-data.frame(signif(de_genes$adj.P.Val, digits = 3))
  rownames(FDRs)<-de_genes$ensemblID
  FDRs<-data.frame("FDR"=apply(FDRs, 1, function(x){if(x>0.01|x==0.01){paste("0.01-0.05")} 
                                                 else {paste("<0.01")}}))
  FDRs[, paste("FDR", name)]<-FDRs$FDR
  FDRs$FDR<-NULL
  
  ## FCs
  FCs<-data.frame(signif(2**(de_genes$logFC), digits = 3))
  rownames(FCs)<-de_genes$ensemblID
  FCs<-data.frame("FC"=apply(FCs, 1, function(x){if(x>1|x==1){paste(">=1")} else {paste("<1")}}))
  FCs[, paste("FC", name)]<-FCs$FC
  FCs$FC<-NULL
  ## Ignore FDRs
  ann_rows<-FCs
  
  ## Join annotation info
  anns<-list("Group"=ann_cols$Group, "Sex"=ann_cols$Sex)
  anns[[paste("FDR", name)]]=ann_rows$FDRs
  anns[[paste("FC", name)]]=ann_rows$FCs 

  
  
  ## Heatmap colors
  break1<-seq(min(vGene_DEG),0,by=0.001)
  break2<-seq(0,max(vGene_DEG),by=0.001)
  my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(break1)),
                  c(colorRampPalette(colors = c("lightsalmon", "darkred"))(n = length(break2)))
  )
  
  
  
  if (filename=="nicotine"){
    width=6.5
    height=5.5
    cellwidth=5
    cellheight=0.25
  }
  else {
    width=8.5
    height=8
    cellwidth=2.5
    cellheight=0.06
  }
  
  
  
  ## Display heatmap
  pheatmap(
    vGene_DEG,
    breaks = c(break1, break2),
    color=my_palette,
    cluster_rows = TRUE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    cluster_cols = TRUE,
    annotation_col = ann_cols,
    annotation_row = ann_rows,
    annotation_colors = ann_colors, 
    cellwidth = cellwidth,
    cellheight = cellheight,
    fontsize=11, 
    width = width,
    height = height,
    filename=paste("plots/06_Heatmap_DEG/Gene_analysis/Heatmap_DEG_", filename, ".pdf", sep="")
    
  )
}


## Nicotine DEG
DEG_heatmaps(rse_gene_brain_pups_nicotine, results_pups_nicotine_fitted, de_genes_pups_nicotine_fitted, "nicotine")
## Smoking DEG
DEG_heatmaps(rse_gene_brain_pups_smoking, results_pups_smoking_fitted, de_genes_pups_smoking_fitted, "smoking")









## 1.2 Heatmaps for Nicotine VS Smoking DEG

nic_vs_smo_heatmaps<- function(DEG_list, option, filename){

  ## Counts for both smo and nic DEG
  vGene_nic<-results_pups_nicotine_fitted[[1]][[2]]
  colnames(vGene_nic)<-vGene_nic$targets$SAMPLE_ID
  vGene_smo<-results_pups_smoking_fitted[[1]][[2]]
  colnames(vGene_smo)<-vGene_smo$targets$SAMPLE_ID

  
  
  ## Heatmaps comparing DEG in nic VS smo 
  if (option=="nic_and_smo"){
  
    ## Extract lognorm counts of the genes
    vGene_DEG_nic <-vGene_nic$E[which(vGene_nic$genes$ensemblID %in% DEG_list$ensemblID), ]
    vGene_DEG_smo <-vGene_smo$E[which(vGene_smo$genes$ensemblID %in% DEG_list$ensemblID), ]
    
    ## Center the data
    vGene_DEG_smo<-(vGene_DEG_smo-rowMeans(vGene_DEG_smo))/rowSds(vGene_DEG_smo)
    vGene_DEG_nic<-(vGene_DEG_nic-rowMeans(vGene_DEG_nic))/rowSds(vGene_DEG_nic)
    
    ## Join data of nic and smo samples 
    vGene_DEG<-cbind(vGene_DEG_nic, vGene_DEG_smo)
    
    
    
    ## Samples' info
    ann_cols_nic <- data.frame("Expt & Group"=paste(vGene_nic$targets$Expt, vGene_nic$targets$Group), 
                               "Sex"=vGene_nic$targets$Sex, check.names=FALSE)
    rownames(ann_cols_nic)<-vGene_nic$targets$SAMPLE_ID
    ann_cols_smo <- data.frame("Expt & Group"=paste(vGene_smo$targets$Expt, vGene_smo$targets$Group), 
                               "Sex"=vGene_smo$targets$Sex, check.names=FALSE)
    rownames(ann_cols_smo)<-vGene_smo$targets$SAMPLE_ID
    ann_cols<-rbind(ann_cols_nic, ann_cols_smo)
    
    ## Genes' info
    ## FDRs in nic and smo
    top_genes_nic<-results_pups_nicotine_fitted[[1]][[1]]
    data_nic<-top_genes_nic[which(top_genes_nic$ensemblID %in% DEG_list$ensemblID), c("ensemblID","adj.P.Val", "logFC")]
    rownames<-rownames(data_nic)
    FDR_nic<-data.frame("FDR nic"=signif(data_nic$adj.P.Val, digits = 3), check.names=FALSE)
    rownames(FDR_nic)<-rownames
    
    top_genes_smo<-results_pups_smoking_fitted[[1]][[1]]
    data_smo<-top_genes_smo[which(top_genes_smo$ensemblID %in% DEG_list$ensemblID), c("ensemblID","adj.P.Val", "logFC")]
    rownames<-rownames(data_smo)
    FDR_smo<-data.frame("FDR smo"=signif(data_smo$adj.P.Val, digits = 3), check.names=FALSE)
    rownames(FDR_smo)<-rownames
    
    FDRs<-cbind(FDR_nic, FDR_smo)
    
    FDRs$Significance<-apply(FDRs, 1, function(x){if (x[1]<=0.05 && x[2]<=0.05){paste("Both")} 
                                                  else if (x[1]<=0.05 && x[2]>0.05){paste("Nic only")}
                                                  else if (x[1]>0.05 && x[2]<=0.05){paste("Smo only")}})
    
    ## FCs in nic and smo 
    FC_nic<-data.frame(signif(2**data_nic$logFC, digits = 3))
    rownames(FC_nic)<-rownames(data_nic)
    FC_nic<-data.frame("FC nic"=apply(FC_nic, 1, function(x){if(x>1|x==1){paste(">=1")} else {paste("<1")}}),
                       check.names=FALSE)
    
    FC_smo<-data.frame(signif(2**data_smo$logFC, digits = 3))
    rownames(FC_smo)<-rownames(data_smo)
    FC_smo<-data.frame("FC smo"=apply(FC_smo, 1, function(x){if(x>1|x==1){paste(">=1")} else {paste("<1")}}),
                       check.names=FALSE)
    
    FCs<-cbind(FC_nic, FC_smo)
    # (Ignore FDRs)
    ann_rows<-FCs
    
    ## Join annotation info 
    anns<-list("Expt & Group"=ann_cols$`Expt & Group`, "Sex"=ann_cols$Sex,
               "FC nic"=ann_rows$`FC nic`, "FC smo"=ann_rows$`FC smo`)
    
    
    
    ## Heatmap colors
    break1<-seq(min(vGene_DEG),0,by=0.001)
    break2<-seq(0,max(vGene_DEG),by=0.001)
    my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(break1)),
                    c(colorRampPalette(colors = c("lightsalmon", "darkred"))(n = length(break2)))
    )
  
  
    
    ## Display heatmap
    pheatmap(
      vGene_DEG,
      breaks = c(break1, break2),
      color=my_palette,
      cluster_rows = TRUE,
      show_rownames = FALSE,
      cluster_cols = TRUE,
      annotation_col = ann_cols,
      annotation_row = ann_rows,
      annotation_colors = ann_colors, 
      fontsize=10, 
      width = 13,
      height = 11,
      filename=paste("plots/06_Heatmap_DEG/Gene_analysis/Heatmap_", filename, "_DEG.pdf", sep="")
    )
  }
  
  
  
  
  ## Heatmaps of DEG list in only nic or smo
  else {
    
    if (option=="nic"){
      expt_name="nicotine"
    }
    else {
      expt_name="smoking"
    }
    
    vGene=eval(parse_expr(paste("vGene", option, sep="_")))
    de_genes=eval(parse_expr(paste("de_genes_pups", expt_name, "fitted", sep="_")))
    top_genes<-eval(parse_expr(paste("results_pups", expt_name, "fitted", sep="_")))[[1]][[1]]
    rse<-eval(parse_expr(paste("rse_gene_brain_pups", expt_name, sep="_")))
    
    
    
    ## Extract lognorm counts of the genes
    vGene_DEG <-vGene$E[which(vGene$genes$ensemblID %in% DEG_list$ensemblID), ]
    
    ## Center the data to make differences more evident
    vGene_DEG<-(vGene_DEG-rowMeans(vGene_DEG))/rowSds(vGene_DEG)
    
    
    
    ## Samples' info
    ann_cols <- as.data.frame(vGene$targets[, c("Group", "Sex")])
    rownames(ann_cols) <- vGene$targets$SAMPLE_ID
    ann_cols$Group <- gsub('Experimental', 'Expt', gsub('Control', 'Ctrl', ann_cols$Group))
    
    ## Genes' info
    ## FDRs
    data<-top_genes[which(top_genes$ensemblID %in% DEG_list$ensemblID), c("ensemblID","adj.P.Val", "logFC")]
    rownames<-rownames(data)
    FDRs<-data.frame("FDR"=apply(as.data.frame(data$adj.P.Val), 1, function(x){if(x>0.01|x==0.01){paste("0.01-0.05")} 
      else {paste("<0.01")}}))
    FDRs[,paste("FDR", option)]<-FDRs$FDR
    FDRs$FDR<-NULL
    rownames(FDRs)<-rownames
    
    ## FCs in nic or smo 
    FCs<-data.frame(signif(2**(data$logFC), digits = 3))
    rownames(FCs)<-rownames
    FCs<-data.frame("FC"=apply(FCs, 1, function(x){if(x>1|x==1){paste(">=1")} else {paste("<1")}}))
    FCs[, paste("FC", option)]<-FCs$FC
    FCs$FC<-NULL
    ann_rows<-cbind(FDRs, FCs)

    ## Join annotation info
    anns<-list("Group"=ann_cols$Group, "Sex"=ann_cols$Sex)
    anns[[paste("FDR", option)]]=ann_rows$FDR
    anns[[paste("FC", option)]]=ann_rows$FC
    
    
  
    ## Heatmap colors
    break1<-seq(min(vGene_DEG),0,by=0.001)
    break2<-seq(0,max(vGene_DEG),by=0.001)
    my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(break1)),
                    c(colorRampPalette(colors = c("lightsalmon", "darkred"))(n = length(break2)))
    )
    
    
    
    ## Display heatmap
    pheatmap(
      vGene_DEG,
      breaks = c(break1, break2),
      color=my_palette,
      cluster_rows = TRUE,
      show_rownames = FALSE,
      cluster_cols = TRUE,
      annotation_col = ann_cols,
      annotation_row = ann_rows,
      annotation_colors = ann_colors, 
      fontsize=8, 
      width = 15,
      height = 13,
      filename=paste("plots/06_Heatmap_DEG/Gene_analysis/Heatmap_", filename, "_DEG_", option, ".pdf", sep="")
    )  
  }
}





## Create heatmaps

### All DEG in either smo or nic

DEG_list<-all
option<-"nic_and_smo"
name<-"all"
nic_vs_smo_heatmaps(DEG_list, option, name)



### Heatmap for DEG Up regulated in nicotine only

DEG_list<-intersections[["only up nic"]]
name<-"only_Up_nic"

option<-"nic"
nic_vs_smo_heatmaps(DEG_list, option, name)

option<-"nic_and_smo"
nic_vs_smo_heatmaps(DEG_list, option, name)



### Heatmap for DEG Up regulated in smoking only

DEG_list<-intersections[["only up smo"]]
name<-"only_Up_smo"

option<-"smo"
nic_vs_smo_heatmaps(DEG_list, option, name)

option<-"nic_and_smo"
nic_vs_smo_heatmaps(DEG_list, option, name)



### Heatmap for DEG Up regulated in both nicotine and smoking

DEG_list<-intersections[["smo Up nic Up"]]
name<-"smoUp_nicUp"

option<-"smo"
nic_vs_smo_heatmaps(DEG_list, option, name)

option<-"nic"
nic_vs_smo_heatmaps(DEG_list, option, name)

option<-"nic_and_smo"
nic_vs_smo_heatmaps(DEG_list, option, name)



### Heatmap for DEG Down regulated in nicotine only

DEG_list<-intersections[["only down nic"]]
name<-"only_Down_nic"

option<-"nic"
nic_vs_smo_heatmaps(DEG_list, option, name)

option<-"nic_and_smo"
nic_vs_smo_heatmaps(DEG_list, option, name)



### Heatmap for DEG Down regulated in smoking only

DEG_list<-intersections[["only down smo"]]
name<-"only_Down_smo"

option<-"smo"
nic_vs_smo_heatmaps(DEG_list, option, name)

option<-"nic_and_smo"
nic_vs_smo_heatmaps(DEG_list, option, name)



### Heatmap for DEG Down regulated in both nicotine and smoking

DEG_list<-intersections[["smo Down nic Down"]]
name<-"smoDown_nicDown"

option<-"smo"
nic_vs_smo_heatmaps(DEG_list, option, name)

option<-"nic"
nic_vs_smo_heatmaps(DEG_list, option, name)

option<-"nic_and_smo"
nic_vs_smo_heatmaps(DEG_list, option, name)



### Hetamap for DEG Up regulated in nicotine and down in smoking

DEG_list<-intersections[["smo Down nic Up"]]
name<-"smoDown_nicUp"

option<-"smo"
nic_vs_smo_heatmaps(DEG_list, option, name)

option<-"nic"
nic_vs_smo_heatmaps(DEG_list, option, name)

option<-"nic_and_smo"
nic_vs_smo_heatmaps(DEG_list, option, name)



### Heatmap for DEG Up regulated in smoking and down in nicotine

DEG_list<-intersections[["smo Up nic Down"]]
name<-"smoUp_nicDown"

option<-"smo"
nic_vs_smo_heatmaps(DEG_list, option, name)

option<-"nic"
nic_vs_smo_heatmaps(DEG_list, option, name)

option<-"nic_and_smo"
nic_vs_smo_heatmaps(DEG_list, option, name)








## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
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
# date     2024-01-04
# rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
# pandoc   NA
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                * version   date (UTC) lib source
# AnnotationDbi          * 1.63.2    2023-07-03 [1] Bioconductor
# AnnotationHub            3.9.1     2023-06-14 [1] Bioconductor
# ape                      5.7-1     2023-03-13 [1] CRAN (R 4.3.0)
# aplot                    0.1.10    2023-03-08 [1] CRAN (R 4.3.0)
# Biobase                * 2.61.0    2023-06-02 [1] Bioconductor
# BiocFileCache            2.9.1     2023-07-14 [1] Bioconductor
# BiocGenerics           * 0.48.1    2023-11-02 [1] Bioconductor
# BiocManager              1.30.21.1 2023-07-18 [1] CRAN (R 4.3.0)
# BiocParallel             1.35.3    2023-07-07 [1] Bioconductor
# BiocVersion              3.18.0    2023-05-11 [1] Bioconductor
# biomaRt                  2.57.1    2023-06-14 [1] Bioconductor
# biomartr               * 1.0.7     2023-12-02 [1] CRAN (R 4.3.1)
# Biostrings               2.69.2    2023-07-05 [1] Bioconductor
# bit                      4.0.5     2022-11-15 [1] CRAN (R 4.3.0)
# bit64                    4.0.5     2020-08-30 [1] CRAN (R 4.3.0)
# bitops                   1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
# blob                     1.2.4     2023-03-17 [1] CRAN (R 4.3.0)
# cachem                   1.0.8     2023-05-01 [1] CRAN (R 4.3.0)
# cli                      3.6.1     2023-03-23 [1] CRAN (R 4.3.0)
# clusterProfiler        * 4.9.2     2023-07-14 [1] Bioconductor
# codetools                0.2-19    2023-02-01 [1] CRAN (R 4.3.0)
# colorspace               2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# cowplot                * 1.1.1     2020-12-30 [1] CRAN (R 4.3.0)
# crayon                   1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# curl                     5.0.1     2023-06-07 [1] CRAN (R 4.3.0)
# data.table               1.14.8    2023-02-17 [1] CRAN (R 4.3.0)
# DBI                      1.1.3     2022-06-18 [1] CRAN (R 4.3.0)
# dbplyr                   2.3.3     2023-07-07 [1] CRAN (R 4.3.0)
# DelayedArray             0.26.6    2023-07-02 [1] Bioconductor
# digest                   0.6.33    2023-07-07 [1] CRAN (R 4.3.0)
# DOSE                     3.27.2    2023-07-07 [1] Bioconductor
# downloader               0.4       2015-07-09 [1] CRAN (R 4.3.0)
# dplyr                    1.1.2     2023-04-20 [1] CRAN (R 4.3.0)
# ellipsis                 0.3.2     2021-04-29 [1] CRAN (R 4.3.0)
# enrichplot               1.21.1    2023-07-03 [1] Bioconductor
# fansi                    1.0.5     2023-10-08 [1] CRAN (R 4.3.1)
# farver                   2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
# fastmap                  1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
# fastmatch                1.1-3     2021-07-23 [1] CRAN (R 4.3.0)
# fgsea                    1.27.0    2023-05-20 [1] Bioconductor
# filelock                 1.0.2     2018-10-05 [1] CRAN (R 4.3.0)
# fs                       1.6.3     2023-07-20 [1] CRAN (R 4.3.0)
# gargle                   1.5.2     2023-07-20 [1] CRAN (R 4.3.0)
# generics                 0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb           * 1.37.2    2023-06-21 [1] Bioconductor
# GenomeInfoDbData         1.2.10    2023-05-28 [1] Bioconductor
# GenomicRanges          * 1.54.1    2023-10-30 [1] Bioconductor
# ggforce                  0.4.1     2022-10-04 [1] CRAN (R 4.3.0)
# ggfun                    0.1.1     2023-06-24 [1] CRAN (R 4.3.0)
# ggplot2                * 3.4.4     2023-10-12 [1] CRAN (R 4.3.1)
# ggplotify                0.1.1     2023-06-27 [1] CRAN (R 4.3.0)
# ggraph                   2.1.0     2022-10-09 [1] CRAN (R 4.3.0)
# ggrepel                  0.9.3     2023-02-03 [1] CRAN (R 4.3.0)
# ggtree                   3.9.0     2023-05-20 [1] Bioconductor
# glue                     1.6.2     2022-02-24 [1] CRAN (R 4.3.0)
# GO.db                    3.17.0    2023-05-28 [1] Bioconductor
# googledrive              2.1.1     2023-06-11 [1] CRAN (R 4.3.0)
# GOSemSim                 2.27.2    2023-07-14 [1] Bioconductor
# graphlayouts             1.0.0     2023-05-01 [1] CRAN (R 4.3.0)
# gridExtra                2.3       2017-09-09 [1] CRAN (R 4.3.0)
# gridGraphics             0.5-1     2020-12-13 [1] CRAN (R 4.3.0)
# gson                     0.1.0     2023-03-07 [1] CRAN (R 4.3.0)
# gtable                   0.3.4     2023-08-21 [1] CRAN (R 4.3.0)
# HDO.db                   0.99.1    2023-05-28 [1] Bioconductor
# here                   * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
# hms                      1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
# HPO.db                   0.99.2    2023-06-28 [1] Bioconductor
# htmltools                0.5.5     2023-03-23 [1] CRAN (R 4.3.0)
# httpuv                   1.6.11    2023-05-11 [1] CRAN (R 4.3.0)
# httr                     1.4.6     2023-05-08 [1] CRAN (R 4.3.0)
# igraph                   1.5.0     2023-06-16 [1] CRAN (R 4.3.0)
# interactiveDisplayBase   1.39.0    2023-06-02 [1] Bioconductor
# IRanges                * 2.36.0    2023-10-26 [1] Bioconductor
# jaffelab               * 0.99.32   2023-05-28 [1] Github (LieberInstitute/jaffelab@21e6574)
# jsonlite                 1.8.8     2023-12-04 [1] CRAN (R 4.3.1)
# KEGGREST                 1.41.0    2023-07-07 [1] Bioconductor
# labeling                 0.4.3     2023-08-29 [1] CRAN (R 4.3.0)
# later                    1.3.1     2023-05-02 [1] CRAN (R 4.3.0)
# lattice                  0.21-8    2023-04-05 [1] CRAN (R 4.3.0)
# lazyeval                 0.2.2     2019-03-15 [1] CRAN (R 4.3.0)
# lifecycle                1.0.3     2022-10-07 [1] CRAN (R 4.3.0)
# limma                    3.57.6    2023-06-21 [1] Bioconductor
# magrittr                 2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# MASS                     7.3-60    2023-05-04 [1] CRAN (R 4.3.0)
# Matrix                   1.6-4     2023-11-30 [1] CRAN (R 4.3.1)
# MatrixGenerics         * 1.13.0    2023-05-20 [1] Bioconductor
# matrixStats            * 1.0.0     2023-06-02 [1] CRAN (R 4.3.0)
# memoise                  2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
# mime                     0.12      2021-09-28 [1] CRAN (R 4.3.0)
# MPO.db                   0.99.7    2023-05-31 [1] Bioconductor
# munsell                  0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# nlme                     3.1-162   2023-01-31 [1] CRAN (R 4.3.0)
# org.Mm.eg.db           * 3.18.0    2024-01-01 [1] Bioconductor
# patchwork                1.1.2     2022-08-19 [1] CRAN (R 4.3.0)
# pheatmap               * 1.0.12    2019-01-04 [1] CRAN (R 4.3.0)
# pillar                   1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig                2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# plyr                     1.8.8     2022-11-11 [1] CRAN (R 4.3.0)
# png                      0.1-8     2022-11-29 [1] CRAN (R 4.3.0)
# polyclip                 1.10-4    2022-10-20 [1] CRAN (R 4.3.0)
# prettyunits              1.1.1     2020-01-24 [1] CRAN (R 4.3.0)
# progress                 1.2.2     2019-05-16 [1] CRAN (R 4.3.0)
# promises                 1.2.0.1   2021-02-11 [1] CRAN (R 4.3.0)
# purrr                    1.0.1     2023-01-10 [1] CRAN (R 4.3.0)
# qvalue                   2.33.0    2023-05-11 [1] Bioconductor
# R6                       2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# rafalib                * 1.0.0     2015-08-09 [1] CRAN (R 4.3.0)
# ragg                     1.2.5     2023-01-12 [1] CRAN (R 4.3.0)
# rappdirs                 0.3.3     2021-01-31 [1] CRAN (R 4.3.0)
# RColorBrewer             1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                     1.0.11    2023-07-06 [1] CRAN (R 4.3.0)
# RCurl                    1.98-1.12 2023-03-27 [1] CRAN (R 4.3.0)
# reshape2                 1.4.4     2020-04-09 [1] CRAN (R 4.3.0)
# rlang                  * 1.1.1     2023-04-28 [1] CRAN (R 4.3.0)
# rprojroot                2.0.3     2022-04-02 [1] CRAN (R 4.3.0)
# RSQLite                  2.3.1     2023-04-03 [1] CRAN (R 4.3.0)
# rstudioapi               0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
# S4Arrays                 1.1.4     2023-06-02 [1] Bioconductor
# S4Vectors              * 0.40.2    2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
# scales                   1.2.1     2022-08-20 [1] CRAN (R 4.3.0)
# scatterpie               0.2.1     2023-06-07 [1] CRAN (R 4.3.0)
# segmented                1.6-4     2023-04-13 [1] CRAN (R 4.3.0)
# sessioninfo            * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
# shadowtext               0.1.2     2022-04-22 [1] CRAN (R 4.3.0)
# shiny                    1.7.4.1   2023-07-06 [1] CRAN (R 4.3.0)
# stringi                  1.7.12    2023-01-11 [1] CRAN (R 4.3.0)
# stringr                  1.5.0     2022-12-02 [1] CRAN (R 4.3.0)
# SummarizedExperiment   * 1.30.2    2023-06-06 [1] Bioconductor
# systemfonts              1.0.4     2022-02-11 [1] CRAN (R 4.3.0)
# textshaping              0.3.6     2021-10-13 [1] CRAN (R 4.3.0)
# tibble                   3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidygraph                1.2.3     2023-02-01 [1] CRAN (R 4.3.0)
# tidyr                    1.3.0     2023-01-24 [1] CRAN (R 4.3.0)
# tidyselect               1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# tidytree                 0.4.4     2023-07-15 [1] CRAN (R 4.3.0)
# treeio                   1.25.1    2023-07-07 [1] Bioconductor
# tweenr                   2.0.2     2022-09-06 [1] CRAN (R 4.3.0)
# utf8                     1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                    0.6.4     2023-10-12 [1] CRAN (R 4.3.1)
# viridis                  0.6.3     2023-05-03 [1] CRAN (R 4.3.0)
# viridisLite              0.4.2     2023-05-02 [1] CRAN (R 4.3.0)
# withr                    2.5.2     2023-10-30 [1] CRAN (R 4.3.1)
# XML                      3.99-0.14 2023-03-19 [1] CRAN (R 4.3.0)
# xml2                     1.3.5     2023-07-06 [1] CRAN (R 4.3.0)
# xtable                   1.8-4     2019-04-21 [1] CRAN (R 4.3.0)
# XVector                  0.41.1    2023-06-02 [1] Bioconductor
# yaml                     2.3.8     2023-12-11 [1] CRAN (R 4.3.1)
# yulab.utils              0.0.6     2022-12-20 [1] CRAN (R 4.3.0)
# zlibbioc                 1.47.0    2023-05-20 [1] Bioconductor
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────