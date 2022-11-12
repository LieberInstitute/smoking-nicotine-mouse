
# 1. Visualize differentially expressed genes

library(here)
library(SummarizedExperiment)
library(stats)
library(jaffelab)
library(pheatmap)
library(rlang)
library(sessioninfo)


load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_nicotine.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_smoking.Rdata"))
load(here("processed-data/04_DEA/results_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/results_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/de_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/05_GO_KEGG/all_DEG.Rdata"))
load(here("processed-data/05_GO_KEGG/intersections.Rdata"))


## (For fitted models only)


## Manually determine coloring for plot annotation of all heatmaps
ann_colors = list()
ann_colors[["Sex"]]=c("F"="#FF99CC", "M"="#99CCFF")
ann_colors[["Group"]]=c("Control"="#FFFF66", "Experimental"="#FF9933")
ann_colors[["FDR nic"]]=c("0.01-0.05"="#9900CC", "<0.01"="#CC99FF")
ann_colors[["FDR smo"]]=c("0.01-0.05"="#FF33FF", "<0.01"="#FF99FF")
ann_colors[["FC smo"]]=c(">=1"="#00CC66", "<1"="#CCFF99")
ann_colors[["FC nic"]]=c(">=1"="#CC6600", "<1"="#FFCC99")
ann_colors[["Significance"]]=c("Nic only"="#9900CC", "Smo only"="#FF33FF", "Both"="#33FF33")
ann_colors[["Expt & Group"]]=c("Nicotine Control"="#FF9999", "Nicotine Experimental"="#FF0066", 
                               "Smoking Control"="#99CC99", "Smoking Experimental"="#669966")




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
  
  ## Remove technical variables' contributions 
  formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr +
    overallMapRate + mitoRate
  model<- model.matrix(formula, data=colData(rse))
  vGene_DEG<-cleaningY(vGene_DEG, model, P=2)
  
  ## Center the data to make differences more evident
  vGene_DEG<-(vGene_DEG-rowMeans(vGene_DEG))/rowSds(vGene_DEG)
  
  
  
  ## Samples' info
  ann_cols <- as.data.frame(vGene$targets[, c("Group", "Sex")])
  rownames(ann_cols) <- vGene$targets$SAMPLE_ID
  colnames(ann_cols) <- c("Group", "Sex")
  
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
  ann_rows<-cbind(FDRs, FCs)
  
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
    width=9
    height=8
  }
  else {
    width=13
    height=12
  }
  
  
  
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
    fontsize=8.5, 
    width = width,
    height = height,
    filename=paste("plots/06_Heatmap_DEG/Heatmap_DEG_", filename, ".pdf", sep="")
    
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
    
    ## Remove technical variables' contributions 
    formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr +
      overallMapRate + mitoRate
    model_nic<- model.matrix(formula, data=colData(rse_gene_brain_pups_nicotine))
    vGene_DEG_nic<-cleaningY(vGene_DEG_nic, model_nic, P=2)
    model_smo<- model.matrix(formula, data=colData(rse_gene_brain_pups_smoking))
    vGene_DEG_smo<-cleaningY(vGene_DEG_smo, model_smo, P=2)
    
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
    ann_rows<-cbind("Significance"=FDRs$Significance, FCs)
    
    ## Join annotation info
    anns<-list("Expt & Group"=ann_cols$`Expt & Group`, "Sex"=ann_cols$Sex, "Significance"=ann_rows$Significance,
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
      fontsize=8, 
      width = 16,
      height = 14,
      filename=paste("plots/06_Heatmap_DEG/Heatmap_", filename, "_DEG.pdf", sep="")
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
   
    ## Remove technical variables' contributions 
    formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr +
      overallMapRate + mitoRate
    model<- model.matrix(formula, data=colData(rse))
    vGene_DEG<-cleaningY(vGene_DEG, model, P=2)
    
    ## Center the data to make differences more evident
    vGene_DEG<-(vGene_DEG-rowMeans(vGene_DEG))/rowSds(vGene_DEG)
    
    
    
    ## Samples' info
    ann_cols <- as.data.frame(vGene$targets[, c("Group", "Sex")])
    rownames(ann_cols) <- vGene$targets$SAMPLE_ID
    
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
      width = 12,
      height = 11,
      filename=paste("plots/06_Heatmap_DEG/Heatmap_", filename, "_DEG_", option, ".pdf", sep="")
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
# setting  value
# version  R version 4.2.0 (2022-04-22 ucrt)
# os       Windows 10 x64 (build 19044)
# system   x86_64, mingw32
# ui       RStudio
# language (EN)
# collate  Spanish_Mexico.utf8
# ctype    Spanish_Mexico.utf8
# tz       America/Mexico_City
# date     2022-10-23
# rstudio  2022.07.2+576 Spotted Wakerobin (desktop)
# pandoc   2.19.2 @ C:/Program Files/RStudio/bin/quarto/bin/tools/ (via rmarkdown)