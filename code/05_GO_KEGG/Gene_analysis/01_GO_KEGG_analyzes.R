
# 1. Gene Ontology and KEGG analyzes of DEG 


library(here)
library(SummarizedExperiment)
library(clusterProfiler)
library(org.Mm.eg.db)
library(jaffelab)
library(ggplot2)
library(cowplot)
library(rlang)
library(biomartr)
library(sessioninfo)


load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/results_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/results_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/DEG_fitted_smo_vs_nic_up.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/DEG_fitted_smo_vs_nic_down.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/DEG_fitted_smoDown_nicUp.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/DEG_fitted_smoUp_nicDown.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_smoking.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_nicotine.Rdata"))


## GO and KEGG analyzes for DEG from fitted model only

## Groups of DEG
up_nic<-de_genes_pups_nicotine_fitted[de_genes_pups_nicotine_fitted$logFC>0, c("EntrezID", "Symbol", "ensemblID")]
up_smo<-de_genes_pups_smoking_fitted[de_genes_pups_smoking_fitted$logFC>0, c("EntrezID", "Symbol", "ensemblID")]
all_up<-unique(rbind(up_nic, up_smo))

down_nic<-de_genes_pups_nicotine_fitted[de_genes_pups_nicotine_fitted$logFC<0, c("EntrezID", "Symbol", "ensemblID")]
down_smo<-de_genes_pups_smoking_fitted[de_genes_pups_smoking_fitted$logFC<0, c("EntrezID", "Symbol", "ensemblID")]
all_down<-unique(rbind(down_nic, down_smo))

all_nic<-unique(rbind(up_nic, down_nic))
all_smo<-unique(rbind(up_smo, down_smo))
all<-unique(rbind(all_nic, all_smo))
save(all, file="processed-data/05_GO_KEGG/Gene_analysis/all_DEG.Rdata")

## Intersections between groups
smoUp_nicDown<-merge(up_smo, down_nic)
smoDown_nicUp<-merge(down_smo, up_nic)
smoUp_nicUp<-merge(up_smo, up_nic)
smoDown_nicDown<-merge(down_smo, down_nic)
only_up_nic<-up_nic[which(! (up_nic$Symbol %in% smoDown_nicUp$Symbol | 
                               up_nic$Symbol %in% smoUp_nicUp$Symbol)),]
only_down_nic<-down_nic[which(! (down_nic$Symbol %in% smoDown_nicDown$Symbol | 
                                   down_nic$Symbol %in% smoUp_nicDown$Symbol)),]
only_up_smo<-up_smo[which(! (up_smo$Symbol %in% smoUp_nicUp$Symbol | 
                               up_smo$Symbol %in% smoUp_nicDown$Symbol)),]
only_down_smo<-down_smo[which(! (down_smo$Symbol %in% smoDown_nicDown$Symbol | 
                                   down_smo$Symbol %in% smoDown_nicUp$Symbol)),]
intersections<-list("only up nic"=only_up_nic, "only up smo"=only_up_smo, 
                    "only down nic"=only_down_nic, "only down smo"=only_down_smo, 
                    "smo Up nic Up"=smoUp_nicUp, "smo Down nic Down"=smoDown_nicDown, 
                    "smo Up nic Down"=smoUp_nicDown, "smo Down nic Up"=smoDown_nicUp)

save(intersections, file="processed-data/05_GO_KEGG/Gene_analysis/intersections.Rdata")


## Function to do GO and KEGG analyzes

GO_KEGG<- function(sigGeneList, geneUniverse, name){
  
  if (name=="intersections"){
    height=17
    width=15
  }
  else {
    height=10
    width=9
  }
  
  ## Do GO 
  ## Obtain biological processes 
  goBP_Adj <- compareCluster(
    sigGeneList,
    fun = "enrichGO",
    universe = geneUniverse,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  ## Save
  pdf(paste("plots/05_GO_KEGG/Gene_analysis/GO_BP_", name, ".pdf", sep=""), height = height, width = width)
  print(dotplot(goBP_Adj, title="GO Enrichment Analysis: Biological processes"))
  dev.off()
  
  
  ## Obtain molecular functions
  goMF_Adj <- compareCluster(
    sigGeneList,
    fun = "enrichGO",
    universe = geneUniverse,
    OrgDb = org.Mm.eg.db,
    ont = "MF",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  ## Save
  pdf(paste("plots/05_GO_KEGG/Gene_analysis/GO_MF_", name, ".pdf", sep=""), height = height, width = width)
  print(dotplot(goMF_Adj, title="GO Enrichment Analysis: Molecular function"))
  dev.off()
  
  
  ## Obtain cellular components
  goCC_Adj <- compareCluster(
    sigGeneList,
    fun = "enrichGO",
    universe = geneUniverse,
    OrgDb = org.Mm.eg.db,
    ont = "CC",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  ## Save
  pdf(paste("plots/05_GO_KEGG/Gene_analysis/GO_CC_", name, ".pdf", sep=""), height = height, width = width)
  print(dotplot(goCC_Adj, title="GO Enrichment Analysis: Cellular components"))
  dev.off()
  
  
  ## Do KEGG
  kegg_Adj <- compareCluster(
    sigGeneList,
    fun = "enrichKEGG",
    organism = "mmu",
    universe = geneUniverse,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  
  ## Save
  pdf(paste("plots/05_GO_KEGG/Gene_analysis/KEGG_", name, ".pdf", sep=""), height = height, width = width)
  print(dotplot(kegg_Adj, title="KEGG Enrichment Analysis"))
  dev.off()
  
  
  goList <- list(
    BP = goBP_Adj,
    MF = goMF_Adj,
    CC = goCC_Adj,
    KEGG = kegg_Adj
  )
  
  return(goList)
}




######################
# Up/Down DEG
######################

## List of DEG sets
sigGeneList <- list("up"=all_up$EntrezID, "down"=all_down$EntrezID) 
sigGeneList <-lapply(sigGeneList, function(x) {
  x[!is.na(x)]
})
## Background genes
geneUniverse <- as.character(union(top_genes_pups_nicotine_fitted$EntrezID,
                                   top_genes_pups_smoking_fitted$EntrezID))
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

goList_global<-GO_KEGG(sigGeneList, geneUniverse, "global")
save(goList_global, file="processed-data/05_GO_KEGG/Gene_analysis/goList_global.Rdata")



###################################
# Nicotine pups Up/Down DEG
###################################

## List of DEG sets
sigGeneList <- list("up"=up_nic$EntrezID, "down"=down_nic$EntrezID) 
sigGeneList <-lapply(sigGeneList, function(x) {
  x[!is.na(x)]
})
## Background genes
geneUniverse <- as.character(top_genes_pups_nicotine_fitted$EntrezID)
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

goList_nic<-GO_KEGG(sigGeneList, geneUniverse, "nicotine")
save(goList_nic, file="processed-data/05_GO_KEGG/Gene_analysis/goList_nic.Rdata")



###################################
# Smoking pups Up/Down DEG
###################################

## List of DEG sets
sigGeneList <- list("up"=up_smo$EntrezID, "down"=down_smo$EntrezID) 
sigGeneList <-lapply(sigGeneList, function(x) {
  x[!is.na(x)]
})
## Background genes
geneUniverse <- as.character(top_genes_pups_smoking_fitted$EntrezID)
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

goList_smo<-GO_KEGG(sigGeneList, geneUniverse, "smoking")
save(goList_smo, file="processed-data/05_GO_KEGG/Gene_analysis/goList_smo.Rdata")



##############################################
# Up/Down Nicotine VS Up/Down Smoking DEG
##############################################

sigGeneList <- list("Only up nic"=only_up_nic$EntrezID, "Only up smo"=only_up_smo$EntrezID,
                    "Only down nic"=only_down_nic$EntrezID, "Only down smo"=only_down_smo$EntrezID,
                    "Smo up, nic down"=smoUp_nicDown$EntrezID, "Smo down, nic up"=smoDown_nicUp$EntrezID,
                    "Smo up, nic up"=smoUp_nicUp$EntrezID, "Smo down, nic down"=smoDown_nicDown$EntrezID) 
sigGeneList <-lapply(sigGeneList, function(x) {
  x[!is.na(x)]
})
## Background genes
geneUniverse <- as.character(union(top_genes_pups_nicotine_fitted$EntrezID,
                                   top_genes_pups_smoking_fitted$EntrezID))
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

goList_intersections<-GO_KEGG(sigGeneList, geneUniverse, "intersections")
save(goList_intersections, file="processed-data/05_GO_KEGG/Gene_analysis/goList_intersections.Rdata")








## 1.1 Boxplots of top genes 

### 1.1.1 Top genes in Nic vs Smo and Up vs Down groups

## Lognorm counts of genes in both nicotine and smoking fitted models
vGene_smo<-results_pups_smoking_fitted[[1]][[2]]
vGene_nic<-results_pups_nicotine_fitted[[1]][[2]]
## Regress out residuals 
formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + 
  overallMapRate + mitoRate
model<- model.matrix(formula, data=colData(rse_gene_brain_pups_smoking))
vGene_smo$E<-cleaningY(vGene_smo$E, model, P=2)
rownames(vGene_nic$E)<-vGene_nic$genes$Symbol
model<- model.matrix(formula, data=colData(rse_gene_brain_pups_nicotine))
vGene_nic$E<-cleaningY(vGene_nic$E, model, P=2)
rownames(vGene_smo$E)<-vGene_smo$genes$Symbol




## Data frame for a single gene with its nicotine and smoking logcounts 
get_df_DEG<- function(gene, vGene_nic, vGene_smo) {
  
  ## Logcounts for that gene 
  logcounts_nic<-vGene_nic$E[which(rownames(vGene_nic)==gene),]
  logcounts_smo<-vGene_smo$E[which(rownames(vGene_smo)==gene),]
  df_nic<-data.frame("Gene_counts"=logcounts_nic, "Expt"=rep("Nicotine", length(logcounts_nic)),
                     "Group"=vGene_nic$targets$Group, "SampleID"=vGene_nic$targets$SAMPLE_ID)
  df_smo<-data.frame("Gene_counts"=logcounts_smo, "Expt"=rep("Smoking", length(logcounts_smo)),
                     "Group"=vGene_smo$targets$Group, "SampleID"=vGene_smo$targets$SAMPLE_ID)
  df<-rbind(df_nic, df_smo)
  return(df)
  
}



## Extract top 6 genes from a list based on their max q-values 
## in nicotine and smoking 

extract_top_genes <- function(DEG_list){
  
  ## If the genes are only from smoking 
  if (length(which(DEG_list %in% de_genes_pups_nicotine_fitted$Symbol))==0){
    
    ## q-values for the gene  
    q_values_smo<-top_genes_pups_smoking_fitted[which(top_genes_pups_smoking_fitted$Symbol %in% DEG_list),
                                                       c("Symbol", "adj.P.Val")]
    q_values_smo$adj.P.Val<-signif(q_values_smo$adj.P.Val, digits = 3)
    ## Genes with the lowest q-values
    genes<-q_values_smo[order(q_values_smo$adj.P.Val),1]
    
  }

  ## If the genes are only from nicotine
  else if (length(which(DEG_list %in% de_genes_pups_smoking_fitted$Symbol))==0){
    
    ## q-values for the gene  
    q_values_nic<-top_genes_pups_nicotine_fitted[which(top_genes_pups_nicotine_fitted$Symbol %in% DEG_list),
                                                c("Symbol", "adj.P.Val")]
    q_values_nic$adj.P.Val<-signif(q_values_nic$adj.P.Val, digits = 3)
    ## Genes with the lowest q-values
    genes<-q_values_nic[order(q_values_nic$adj.P.Val),1]
    
  }
  
  else {
    
    nic_qvals<-vector()
    smo_qvals<-vector()
    max_qvals<-vector()
    
    for (DEgene in DEG_list){
      ## q-values for the gene in nicotine and smoking 
      q_value_nic<-signif(top_genes_pups_nicotine_fitted[which(top_genes_pups_nicotine_fitted$Symbol==DEgene),
                                                         "adj.P.Val"], digits = 3)
      q_value_smo<-signif(top_genes_pups_smoking_fitted[which(top_genes_pups_smoking_fitted$Symbol==DEgene),
                                                        "adj.P.Val"], digits = 3)
      nic_qvals<-append(nic_qvals, q_value_nic)
      smo_qvals<-append(smo_qvals, q_value_smo)
      ## Max q-value
      max_qvals<-append(max_qvals, max(q_value_nic, q_value_smo))
    }
    
    q_vals<-data.frame(Gene=DEG_list, Nicotine=nic_qvals, Smoking=smo_qvals, Max=max_qvals)
    
    ## Genes with the lowest max q-values
    genes<-q_vals[order(q_vals$Max),1]
  }

  
  if (length(genes)<6){
    return(genes)
  }
  else {
    return(genes[1:6])
  }
}



## Boxplot to compare the nicotine vs smoking lognorm counts for a single gene
DEG_GO_boxplot <- function(DEgene){
  
  ## Extract necessary data
  df<-get_df_DEG(gene = DEgene, vGene_nic, vGene_smo)
  ## Extract Ensembl ID of the gene
  ensemblID<-top_genes_pups_nicotine_fitted[which(top_genes_pups_nicotine_fitted$Symbol==DEgene), "ensemblID"]
  
  ## q-value for the gene in nicotine and smoking 
  q_value_nic<-signif(top_genes_pups_nicotine_fitted[which(top_genes_pups_nicotine_fitted$Symbol==DEgene),
                                                    "adj.P.Val"], digits = 3)
  q_value_smo<-signif(top_genes_pups_smoking_fitted[which(top_genes_pups_smoking_fitted$Symbol==DEgene),
                                                   "adj.P.Val"], digits = 3)
  ## Log FC for the gene in nicotine and smoking 
  FC_nic<-signif(2**(top_genes_pups_nicotine_fitted[which(top_genes_pups_nicotine_fitted$Symbol==DEgene),
                                                     "logFC"]), digits = 3)
  FC_smo<-signif(2**(top_genes_pups_smoking_fitted[which(top_genes_pups_smoking_fitted$Symbol==DEgene),
                                                    "logFC"]), digits = 3)

  ## Boxplot for each DE gene
  p <-ggplot(data=as.data.frame(df), aes(x=Group,y=Gene_counts)) + 
      geom_boxplot(outlier.color = "#FFFFFFFF") +
      geom_jitter(aes(color=Group), position=position_jitter(0.2)) +
      theme_classic() +
      labs(x = "Experiment", y = "logcounts - covariates",
           title = paste(DEgene, ensemblID, sep=" - "),
           subtitle=" ") +
      theme(plot.margin=unit (c (1,1.5,1,1), 'cm'), legend.position = "none",
            plot.title = element_text(hjust=0.5, size=10, face="bold"), 
            plot.subtitle = element_text(size=17)) +
      scale_color_manual(values = c("orangered", "Dodgerblue")) +
      facet_wrap(~ Expt, scales = "free") +
      scale_x_discrete(labels=c("Ctrl", "Expt"))
  
  p <-ggdraw(p) + 
      draw_label(paste("FDR:", q_value_nic), x = 0.35, y = 0.87, size=9, color = "darkslategray") +
      draw_label(paste("FC:", FC_nic), x = 0.35, y = 0.84, size=9, color = "darkslategray") +
      draw_label(paste("FDR:", q_value_smo), x = 0.72, y = 0.87, size=9, color = "darkslategray") +
      draw_label(paste("FC:", FC_smo), x = 0.71, y = 0.84, size=9, color = "darkslategray") 
  
  return(p)
  
}



## Multiple plots for the top 6 genes in each group
GO_boxplots<- function (DEG_list, groups){
  plots<-list()
  i=1
  for (DEG in DEG_list){
    p<-DEG_GO_boxplot(DEG)
    plots[[i]]<-p
    i=i+1
  }
  plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], ncol=3)
  ggsave(here(paste("plots/05_GO_KEGG/Gene_analysis/Top", length(DEG_list), "_DEG_boxplots_", groups, ".pdf", sep="")), 
         width = 40, height = 25, units = "cm") 
}




## Extract symbols of DEG Up/Down in Nic/Smo
nic_smo_Up<-intersect(DEG_fitted_smo_vs_nic_up[[1]], DEG_fitted_smo_vs_nic_up[[2]])
nic_smo_Down<-intersect(DEG_fitted_smo_vs_nic_down[[1]], DEG_fitted_smo_vs_nic_down[[2]])
nicUp_smoDown<-intersect(DEG_fitted_smoDown_nicUp[[1]], DEG_fitted_smoDown_nicUp[[2]])
nicDown_smoUp<-intersect(DEG_fitted_smoUp_nicDown[[1]], DEG_fitted_smoUp_nicDown[[2]])



## Boxplots of the top 6 genes in each group

## DEG Up in Nic and Up in Smo
nic_smo_Up<-extract_top_genes(nic_smo_Up)
GO_boxplots(nic_smo_Up, "nic_smo_Up")

## DEG Down in Nic and Down in Smo
nic_smo_Down<-extract_top_genes(nic_smo_Down)
GO_boxplots(nic_smo_Down, "nic_smo_Down")

## DEG Up in Nic and Down in Smo
nicUp_smoDown<-extract_top_genes(nicUp_smoDown)
GO_boxplots(nicUp_smoDown, "nicUp_smoDown")

## DEG Down in Nic and Up in Smo
nicDown_smoUp<-extract_top_genes(nicDown_smoUp)
GO_boxplots(nicDown_smoUp, "nicDown_smoUp")








### 1.1.2 Genes in GO and KEGG descriptions

## Extract genes from each BP, CC, MF and KEGG
GO_KEGG_genes<- function(golist, term, cluster, description){
  
  GOdata<-as.data.frame(eval(parse_expr(paste(golist, "$", term, sep=""))))
  genes<-unique(GOdata[which(GOdata$Description==description & GOdata$Cluster==cluster), "geneID"])
  genes<-strsplit(genes, "/")
  genes<-unique(unlist(genes))
  
  if (term=="KEGG"){
    ## KEGG ids (entrez ids) to gene symbols
    symbols<-biomart(genes  = genes,
                     mart       = "ENSEMBL_MART_ENSEMBL",
                     dataset    = "mmusculus_gene_ensembl",
                     attributes = c("external_gene_name"),
                     filters    = "entrezgene_id")
    genes<-symbols$external_gene_name
  }
  return(genes)
  
}


## Boxplots of the first 6 genes involved in a process/pathway
GO_KEGG_boxplots<-function(DEG_list, description, cluster){
  plots<-list()
  i=1
  for (DEgene in DEG_list){
    plots[[i]]<-DEG_GO_boxplot(DEgene)
    i=i+1
  }
  if (length(DEG_list)<6){
    for (i in (length(plots)+1):6){
      plots[[i]]<-NA
    }
  }

  
  options(warn = - 1)   
  plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], ncol=3)
  ggsave(here(paste("plots/05_GO_KEGG/Gene_analysis/Top", length(DEG_list), "_", description,"_boxplots_",cluster, 
                    ".pdf", sep="")), width = 40, height = 25, units = "cm") 
  
}



## Boxplots 


## 1. Cellular components

## Genes of the SNARE complex
GO_genes<-GO_KEGG_genes("goList_intersections", "CC", "Only up smo", "SNARE complex")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "SNARE_complex", "Only_up_smo")

GO_genes<-GO_KEGG_genes("goList_intersections", "CC", "Smo up, nic down", "SNARE complex")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "SNARE_complex", "smoUp_nicDown")

GO_genes<-GO_KEGG_genes("goList_smo", "CC", "up", "SNARE complex")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "SNARE_complex", "up_smo")


## Genes in asymmetric synapses
GO_genes<-GO_KEGG_genes("goList_nic", "CC", "up", "asymmetric synapse")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Asymmetric_synapse", "up_nic")



## 2. Molecular functions

## Genes with heat shock protein binding activity
GO_genes<-GO_KEGG_genes("goList_intersections", "MF", "Smo down, nic up", "heat shock protein binding")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "HeatShock_protein_binding", "smoDown_nicUp")



## 3. Pathways

## Genes involved in Parkinson disease
GO_genes<-GO_KEGG_genes("goList_global", "KEGG", "up", "Parkinson disease")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Parkinson_disease", "up")

## Genes involved in dopaminergic synapses
GO_genes<-GO_KEGG_genes("goList_intersections", "KEGG", "Only up nic", "Dopaminergic synapse")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Dopaminergic_synapse", "Only_up_nic")

GO_genes<-GO_KEGG_genes("goList_nic", "KEGG", "up", "Dopaminergic synapse")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Dopaminergic_synapse", "up_nic")

## Genes involved in long−term depression
GO_genes<-GO_KEGG_genes("goList_intersections", "KEGG", "Only up nic", "Long-term depression")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Long-term_depression", "Only_up_nic")

## Genes involved in SNARE interactions in vesicle transport
GO_genes<-GO_KEGG_genes("goList_intersections", "KEGG", "Smo up, nic down", "SNARE interactions in vesicular transport")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "SNARE_int_ves_transport", "smoUp_nicDown")








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
# date     2022-10-02
# rstudio  2022.07.2+576 Spotted Wakerobin (desktop)
# pandoc   2.19.2 @ C:/Program Files/RStudio/bin/quarto/bin/tools/ (via rmarkdown)


