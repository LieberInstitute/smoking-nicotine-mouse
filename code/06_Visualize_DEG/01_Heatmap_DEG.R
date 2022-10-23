
# 1. Visualize differentially expressed genes

library(here)
library(SummarizedExperiment)
library(stats)
library(jaffelab)
library(pheatmap)


load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_nicotine.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_smoking.Rdata"))
load(here("processed-data/04_DEA/results_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/results_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/de_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/05_GO_KEGG/all_DEG.Rdata"))



## 1.1 Heatmaps of Nicotine and Smoking DEG 

DEG_heatmaps<- function(rse, results, de_genes, name){
  
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
  for (i in 1:nrow(vGene_DEG)){
    vGene_DEG[i,]<-vGene_DEG[i,]-mean(vGene_DEG[i,])
  }
  
  
  
  ## Samples' info
  df <- as.data.frame(vGene$targets[, c("Group", "Sex")])
  rownames(df) <- vGene$targets$SAMPLE_ID
  colnames(df) <- c("Group", "Sex")
  ## Genes' info
  FDRs<-data.frame("FDR"=signif(de_genes$adj.P.Val, digits = 3))
  rownames(FDRs)<-rownames(de_genes)
  
  
  
  ## Manually determine coloring for plot annotation
  ## Column annotation
  palette_names = c('Dark2', 'Accent')
  ann_colors = list()
  for (i in 1:ncol(df)) {
    col_name = colnames(df)[i]
    n_uniq_colors = length(unique(df[,col_name]))
    
    ## Use a unique palette with the correct number of levels, named with those levels
    ann_colors[[col_name]] = RColorBrewer::brewer.pal(n=6, palette_names[i])[4:3+n_uniq_colors]
    names(ann_colors[[col_name]]) = unique(df[,col_name])
    
  }
  ## Row annotation 
  ann_colors[["FDR"]]=colorRampPalette(colors = c("#FFB6C1", "#CD6090"))(n = nrow(vGene_DEG))
  
  
  
  ## Heatmap colors
  break1<-seq(min(vGene_DEG),0.001,by=0.009)
  break2<-seq(0.01,max(vGene_DEG),by=0.001)
  my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(break1)-1),
                  c(colorRampPalette(colors = c("darkred", "tomato1"))(n = length(break2)-1))
  )
  
  if (name=="nicotine"){
    width=9
    height=8
  }
  else {
    width=12
    height=11
  }
  
  
  ## Display heatmap
  pheatmap(
    vGene_DEG,
    breaks = c(break1, break2),
    color=my_palette,
    cluster_rows = TRUE,
    show_rownames = FALSE,
    cluster_cols = TRUE,
    annotation_col = df,
    annotation_row = FDRs,
    annotation_colors = ann_colors, 
    fontsize=8, 
    width = width,
    height = height,
    filename=paste("plots/06_Heatmap_DEG/Heatmap_DEG", name, ".pdf", sep="")
    
  )
}


## For fitted models only
## Nicotine DEG
DEG_heatmaps(rse_gene_brain_pups_nicotine, results_pups_nicotine_fitted, de_genes_pups_nicotine_fitted, "nicotine")
## Smoking DEG
DEG_heatmaps(rse_gene_brain_pups_smoking, results_pups_smoking_fitted, de_genes_pups_smoking_fitted, "smoking")







## 1.2 Heatmaps for Nicotine VS Smoking DEG


## All DEG in either smo or nic

## Extract lognorm counts of all DEG
vGene_DEG_all_nic <-vGene_nic$E[which(vGene_nic$genes$ensemblID %in% all$ensemblID), ]
vGene_DEG_all_smo <-vGene_smo$E[which(vGene_smo$genes$ensemblID %in% all$ensemblID), ]

## Remove technical variables' contributions 
formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr +
  overallMapRate + mitoRate
model_nic<- model.matrix(formula, data=colData(rse_gene_brain_pups_nicotine))
vGene_DEG_all_nic<-cleaningY(vGene_DEG_all_nic, model_nic, P=2)
model_smo<- model.matrix(formula, data=colData(rse_gene_brain_pups_smoking))
vGene_DEG_all_smo<-cleaningY(vGene_DEG_all_smo, model_smo, P=2)

## Join data of nic and smo samples 
vGene_DEG_all<-cbind(vGene_DEG_all_nic, vGene_DEG_all_smo)

## Center the data to make differences more evident
for (i in 1:nrow(vGene_DEG_all)){
  vGene_DEG_all[i,]<-vGene_DEG_all[i,]-mean(vGene_DEG_all[i,])
}




## Samples' info
df_nic <- as.data.frame(vGene_nic$targets[, c("Group", "Sex", "Expt")])
df_smo <- as.data.frame(vGene_smo$targets[, c("Group", "Sex", "Expt")])
df<-rbind(df_nic, df_smo)

## Manually determine coloring for plot annotation
palette_names = c('Dark2', 'Paired', 'RdPu')
ann_colors = list()
for (i in 1:ncol(df)) {
  col_name = colnames(df)[i]
  n_uniq_colors = length(unique(df[,col_name]))
  
  # # Use a unique palette with the correct number of levels, named with those levels
  ann_colors[[col_name]] = RColorBrewer::brewer.pal(n_uniq_colors, palette_names[i])[1:n_uniq_colors]
  names(ann_colors[[col_name]]) = unique(df[,col_name])
  
}

## Display heatmap
pheatmap(
  vGene_DEG_all,
  cluster_rows = TRUE,
  show_rownames = FALSE,
  cluster_cols = TRUE,
  annotation_col = df,
  annotation_colors = ann_colors
)


