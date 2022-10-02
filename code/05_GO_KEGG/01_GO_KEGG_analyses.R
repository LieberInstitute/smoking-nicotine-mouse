
# 1. Gene Ontology and KEGG analyses


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


load(here("processed-data/04_DEA/top_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/top_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/de_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/results_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/results_pups_smoking_fitted.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_nicotine.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_smoking.Rdata"))
load(here("processed-data/04_DEA/DEG_fitted_smo_vs_nic_up.Rdata"))
load(here("processed-data/04_DEA/DEG_fitted_smo_vs_nic_down.Rdata"))
load(here("processed-data/04_DEA/DEG_fitted_smoDown_nicUp.Rdata"))
load(here("processed-data/04_DEA/DEG_fitted_smoUp_nicDown.Rdata"))


## GO and KEGG analyses for DEG from fitted model only

## Groups of DEG
up_nic<-de_genes_pups_nicotine_fitted[de_genes_pups_nicotine_fitted$logFC>0, c("EntrezID", "Symbol")]
up_smo<-de_genes_pups_smoking_fitted[de_genes_pups_smoking_fitted$logFC>0, c("EntrezID", "Symbol")]
all_up<-unique(rbind(up_nic, up_smo))

down_nic<-de_genes_pups_nicotine_fitted[de_genes_pups_nicotine_fitted$logFC<0, c("EntrezID", "Symbol")]
down_smo<-de_genes_pups_smoking_fitted[de_genes_pups_smoking_fitted$logFC<0, c("EntrezID", "Symbol")]
all_down<-unique(rbind(down_nic, down_smo))

all_nic<-unique(rbind(up_nic, down_nic))
all_smo<-unique(rbind(up_smo, down_smo))
all<-unique(rbind(all_nic, all_smo))

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




## Function to do GO and KEGG analyses

GO_KEGG<- function(sigGeneList, geneUniverse, name){
  
  if (name=="intersections"){
    height=26
    width=18
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
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  
  ## Save
  pdf(paste("plots/05_GO_KEGG/GO_BP_", name, ".pdf", sep=""), height = height, width = width)
  print(dotplot(goBP_Adj, title="GO Enrichment Analysis: Biological processes"))
  dev.off()
  pdf(paste("plots/05_GO_KEGG/GO_BP_", name, "_includeAllFalse.pdf", sep=""), height = height, width = width)
  print(dotplot(goBP_Adj, title="GO Enrichment Analysis: Biological processes", includeAll=FALSE))
  dev.off()
  
  
  ## Obtain molecular functions
  goMF_Adj <- compareCluster(
    sigGeneList,
    fun = "enrichGO",
    universe = geneUniverse,
    OrgDb = org.Mm.eg.db,
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  
  ## Save
  pdf(paste("plots/05_GO_KEGG/GO_MF_", name, ".pdf", sep=""), height = height, width = width)
  print(dotplot(goMF_Adj, title="GO Enrichment Analysis: Molecular function"))
  dev.off()
  pdf(paste("plots/05_GO_KEGG/GO_MF_", name, "_includeAllFalse.pdf", sep=""), height = height, width = width)
  print(dotplot(goMF_Adj, title="GO Enrichment Analysis: Molecular function", includeAll=FALSE))
  dev.off()
  
  
  ## Obtain cellular components
  goCC_Adj <- compareCluster(
    sigGeneList,
    fun = "enrichGO",
    universe = geneUniverse,
    OrgDb = org.Mm.eg.db,
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  
  ## Save
  pdf(paste("plots/05_GO_KEGG/GO_CC_", name, ".pdf", sep=""), height = height, width = width)
  print(dotplot(goCC_Adj, title="GO Enrichment Analysis: Cellular components"))
  dev.off()
  pdf(paste("plots/05_GO_KEGG/GO_CC_", name, "_includeAllFalse.pdf", sep=""), height = height, width = width)
  print(dotplot(goCC_Adj, title="GO Enrichment Analysis: Cellular components", includeAll=FALSE))
  dev.off()
  
  
  ## Do KEGG
  kegg_Adj <- compareCluster(
    sigGeneList,
    fun = "enrichKEGG",
    organism = "mmu",
    universe = geneUniverse,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  
  ## Save
  pdf(paste("plots/05_GO_KEGG/KEGG_", name, ".pdf", sep=""), height = height, width = width)
  print(dotplot(kegg_Adj, title="KEGG Enrichment Analysis"))
  dev.off()
  pdf(paste("plots/05_GO_KEGG/KEGG_", name, "_includeAllFalse.pdf", sep=""), height =height, width = width)
  print(dotplot(kegg_Adj, title="KEGG Enrichment Analysis", includeAll=FALSE))
  dev.off()
  
  
  goList <- list(
    BP = goBP_Adj,
    MF = goMF_Adj,
    CC = goCC_Adj,
    KEGG = kegg_Adj
  )
  save(goList, file = paste("processed-data/05_GO_KEGG/goList_", name, ".Rdata", sep=""))
  
  return(goList)
}




######################
# All/Up/Down DEG
######################

## List of DEG sets
sigGeneList <- list("all"=all$EntrezID, "up"=all_up$EntrezID, "down"=all_down$EntrezID) 
sigGeneList <-lapply(sigGeneList, function(x) {
  x[!is.na(x)]
})
## Background genes
geneUniverse <- as.character(union(top_genes_pups_nicotine_fitted$EntrezID,
                                   top_genes_pups_smoking_fitted$EntrezID))
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

goList_global<-GO_KEGG(sigGeneList, geneUniverse, "global")
save(goList_global, file="processed-data/05_GO_KEGG/goList_global.Rdata")



###################################
# Nicotine pups All/Up/Down DEG
###################################

## List of DEG sets
sigGeneList <- list("all"=all_nic$EntrezID, "up"=up_nic$EntrezID, "down"=down_nic$EntrezID) 
sigGeneList <-lapply(sigGeneList, function(x) {
  x[!is.na(x)]
})
## Background genes
geneUniverse <- as.character(top_genes_pups_nicotine_fitted$EntrezID)
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

goList_nic<-GO_KEGG(sigGeneList, geneUniverse, "nicotine")
save(goList_nic, file="processed-data/05_GO_KEGG/goList_nic.Rdata")



###################################
# Smoking pups All/Up/Down DEG
###################################

## List of DEG sets
sigGeneList <- list("all"=all_smo$EntrezID, "up"=up_smo$EntrezID, "down"=down_smo$EntrezID) 
sigGeneList <-lapply(sigGeneList, function(x) {
  x[!is.na(x)]
})
## Background genes
geneUniverse <- as.character(top_genes_pups_smoking_fitted$EntrezID)
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

goList_smo<-GO_KEGG(sigGeneList, geneUniverse, "smoking")
save(goList_smo, file="processed-data/05_GO_KEGG/goList_smo.Rdata")



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
save(goList_intersections, file="processed-data/05_GO_KEGG/goList_intersections.Rdata")








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


## Compare the Nicotine vs Smoking lognorm counts for a single gene
DEG_GO_boxplot <- function(DEgene){
  
  ## Extract necessary data
  df<-get_df_DEG(gene = DEgene, vGene_nic, vGene_smo)
  ## q-value for the gene in nicotine and smoking 
  q_value_nic<-signif(top_genes_pups_nicotine_fitted[which(top_genes_pups_nicotine_fitted$Symbol==DEgene),
                                                    "adj.P.Val"], digits = 3)
  q_value_smo<-signif(top_genes_pups_smoking_fitted[which(top_genes_pups_smoking_fitted$Symbol==DEgene),
                                                   "adj.P.Val"], digits = 3)
  
  ## Boxplot for each DE gene
  p<-ggplot(data=as.data.frame(df), aes(x=Expt,y=Gene_counts)) + 
    ## Hide outliers
    geom_boxplot(aes(fill=Group),  outlier.color = "#FFFFFFFF") +
    theme_classic() +
    labs(x = "Experiment", y = "lognorm counts",
         title = DEgene, 
         subtitle = paste("FDR in nicotine:", q_value_nic, "               ", "FDR in smoking:",
                          q_value_smo)) +
    theme(plot.margin=unit (c (1,1.5,1,1), 'cm'), legend.position = "right",
          plot.title = element_text(hjust=0.5, size=10, face="bold"),
          plot.subtitle = element_text(size = 9, hjust = 0.5), axis.text.x=element_blank()) +
    scale_fill_manual(values = c("Dodgerblue", "orangered")) +
    facet_wrap(~ Expt, scales = "free") 
  
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
  ggsave(here(paste("plots/05_GO_KEGG/Top_DEG_boxplots_", groups, ".pdf", sep="")), 
         width = 50, height = 20, units = "cm") 
}




## Extract symbols of DEG Up/Down in Nic/Smo
nic_smo_Up<-intersect(DEG_fitted_smo_vs_nic_up[[1]], DEG_fitted_smo_vs_nic_up[[2]])
nic_smo_Down<-intersect(DEG_fitted_smo_vs_nic_down[[1]], DEG_fitted_smo_vs_nic_down[[2]])
nicUp_smoDown<-intersect(DEG_fitted_smoDown_nicUp[[1]], DEG_fitted_smoDown_nicUp[[2]])
nicDown_smoUp<-intersect(DEG_fitted_smoUp_nicDown[[1]], DEG_fitted_smoUp_nicDown[[2]])



## Top 6 genes in nicotine

## DEG Up in Nic and Up in Smo
## Order genes by their FDRs  
nic_smo_Up_DEG<-de_genes_pups_nicotine_fitted[which(de_genes_pups_nicotine_fitted$Symbol %in% nic_smo_Up),]
nic_smo_Up_sorted<-nic_smo_Up_DEG[order(nic_smo_Up_DEG$adj.P.Val),"Symbol"]
GO_boxplots(nic_smo_Up_sorted[1:6], "nic_smo_Up_inNic")

## DEG Down in Nic and Down in Smo
nic_smo_Down_DEG<-de_genes_pups_nicotine_fitted[which(de_genes_pups_nicotine_fitted$Symbol %in% nic_smo_Down),]
nic_smo_Down_sorted<-nic_smo_Down_DEG[order(nic_smo_Down_DEG$adj.P.Val),"Symbol"]
GO_boxplots(nic_smo_Down_sorted[1:6], "nic_smo_Down_inNic")

## DEG Up in Nic and Down in Smo
nicUp_smoDown_DEG<-de_genes_pups_nicotine_fitted[which(de_genes_pups_nicotine_fitted$Symbol %in% nicUp_smoDown),]
nicUp_smoDown_sorted<-nicUp_smoDown_DEG[order(nicUp_smoDown_DEG$adj.P.Val),"Symbol"]
GO_boxplots(nicUp_smoDown_sorted[1:6], "nicUp_smoDown_inNic")

## DEG Down in Nic and Up in Smo
nicDown_smoUp_DEG<-de_genes_pups_nicotine_fitted[which(de_genes_pups_nicotine_fitted$Symbol %in% nicDown_smoUp),]
nicDown_smoUp_sorted<-nicDown_smoUp_DEG[order(nicDown_smoUp_DEG$adj.P.Val),"Symbol"]
GO_boxplots(nicDown_smoUp_sorted[1:6], "nicDown_smoUp_inNic")





## Top 6 genes in smoking

## DEG Up in Nic and Up in Smo
nic_smo_Up_DEG<-de_genes_pups_smoking_fitted[which(de_genes_pups_smoking_fitted$Symbol %in% nic_smo_Up),]
nic_smo_Up_sorted<-nic_smo_Up_DEG[order(nic_smo_Up_DEG$adj.P.Val),"Symbol"]
GO_boxplots(nic_smo_Up_sorted[1:6], "nic_smo_Up_inSmo")

## DEG Down in Nic and Down in Smo
nic_smo_Down_DEG<-de_genes_pups_smoking_fitted[which(de_genes_pups_smoking_fitted$Symbol %in% nic_smo_Down),]
nic_smo_Down_sorted<-nic_smo_Down_DEG[order(nic_smo_Down_DEG$adj.P.Val),"Symbol"]
GO_boxplots(nic_smo_Down_sorted[1:6], "nic_smo_Down_inSmo")

## DEG Up in Nic and Down in Smo
nicUp_smoDown_DEG<-de_genes_pups_smoking_fitted[which(de_genes_pups_smoking_fitted$Symbol %in% nicUp_smoDown),]
nicUp_smoDown_sorted<-nicUp_smoDown_DEG[order(nicUp_smoDown_DEG$adj.P.Val),"Symbol"]
GO_boxplots(nicUp_smoDown_sorted[1:6], "nicUp_smoDown_inSmo")

## DEG Down in Nic and Up in Smo
nicDown_smoUp_DEG<-de_genes_pups_smoking_fitted[which(de_genes_pups_smoking_fitted$Symbol %in% nicDown_smoUp),]
nicDown_smoUp_sorted<-nicDown_smoUp_DEG[order(nicDown_smoUp_DEG$adj.P.Val),"Symbol"]
GO_boxplots(nicDown_smoUp_sorted[1:6], "nicDown_smoUp_inSmo")









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
  ggsave(here(paste("plots/05_GO_KEGG/", description,"_boxplots_",cluster, ".pdf", sep="")), 
         width = 50, height = 20, units = "cm") 
  
}



## Boxplots 

## Genes involved in postsynapse organization
genes<-GO_KEGG_genes("goList_intersections", "BP", "Only up nic", "postsynapse organization")
GO_KEGG_boxplots(genes[1:6], "Postsynapse_organization", "Only_up_nic" )

genes<-GO_KEGG_genes("goList_intersections", "BP", "Only up smo", "postsynapse organization")
GO_KEGG_boxplots(genes[1:6], "Postsynapse_organization", "Only_up_smo" )

genes<-GO_KEGG_genes("goList_intersections", "BP", "Only down smo", "postsynapse organization")
GO_KEGG_boxplots(genes[1:6], "Postsynapse_organization", "Only_down_smo" )

genes<-GO_KEGG_genes("goList_intersections", "BP", "Smo down, nic up", "postsynapse organization")
GO_KEGG_boxplots(genes[1:4], "Postsynapse_organization", "SmoDown_nicUp" )

genes<-GO_KEGG_genes("goList_intersections", "BP", "Smo up, nic up", "postsynapse organization")
GO_KEGG_boxplots(genes[1:4], "Postsynapse_organization", "SmoUp_nicUp" )

genes<-GO_KEGG_genes("goList_intersections", "BP", "Only up nic", "neuronal stem cell population maintenance")
GO_KEGG_boxplots(genes[1], "Neuronal_stem_cell_population_maintenance", "Only_up_nic" )

genes<-GO_KEGG_genes("goList_intersections", "BP", "Only down smo", "neuronal stem cell population maintenance")
GO_KEGG_boxplots(genes[1:4], "Neuronal_stem_cell_population_maintenance", "Only_down_smo" )

genes<-GO_KEGG_genes("goList_intersections", "BP", "Smo down, nic up", "neuronal stem cell population maintenance")
GO_KEGG_boxplots(genes[1:2], "Neuronal_stem_cell_population_maintenance", "SmoDown_nicUp" )



## Genes involved in long-term synaptic potentiation
genes<-GO_KEGG_genes("goList_intersections", "BP", "Only up nic", "long-term synaptic potentiation")
GO_KEGG_boxplots(genes[1:3], "Long-term_synaptic_potentiation", "Only_up_nic" )

genes<-GO_KEGG_genes("goList_intersections", "BP", "Only up smo", "long-term synaptic potentiation")
GO_KEGG_boxplots(genes[1:6], "Long-term_synaptic_potentiation", "Only_up_smo" )

genes<-GO_KEGG_genes("goList_intersections", "BP", "Only down nic", "long-term synaptic potentiation")
GO_KEGG_boxplots(genes[1:2], "Long-term_synaptic_potentiation", "Only_down_nic" )

genes<-GO_KEGG_genes("goList_intersections", "BP", "Only down smo", "long-term synaptic potentiation")
GO_KEGG_boxplots(genes[1:6], "Long-term_synaptic_potentiation", "Only_down_smo" )

genes<-GO_KEGG_genes("goList_intersections", "BP", "Smo up, nic down", "long-term synaptic potentiation")
GO_KEGG_boxplots(genes[1], "Long-term_synaptic_potentiation", "SmoUp_nicDown" )

genes<-GO_KEGG_genes("goList_intersections", "BP", "Smo down, nic up", "long-term synaptic potentiation")
GO_KEGG_boxplots(genes[1:2], "Long-term_synaptic_potentiation", "SmoDown_nicUp" )

genes<-GO_KEGG_genes("goList_intersections", "BP", "Smo up, nic up", "long-term synaptic potentiation")
GO_KEGG_boxplots(genes[1], "Long-term_synaptic_potentiation", "SmoUp_nicUp" )



## Genes involved in response to inorganic substance
genes<-GO_KEGG_genes("goList_intersections", "BP", "Only up nic", "response to inorganic substance")
GO_KEGG_boxplots(genes[1:6], "Response_to_inorganic_substance", "Only_up_nic" )

genes<-GO_KEGG_genes("goList_intersections", "BP", "Only up smo", "response to inorganic substance")
GO_KEGG_boxplots(genes[1:6], "Response_to_inorganic_substance", "Only_up_smo" )

genes<-GO_KEGG_genes("goList_intersections", "BP", "Only down nic", "response to inorganic substance")
GO_KEGG_boxplots(genes[1:5], "Response_to_inorganic_substance", "Only_down_nic" )

genes<-GO_KEGG_genes("goList_intersections", "BP", "Only down smo", "response to inorganic substance")
GO_KEGG_boxplots(genes[1:6], "Response_to_inorganic_substance", "Only_down_smo" )

genes<-GO_KEGG_genes("goList_intersections", "BP", "Smo up, nic up", "response to inorganic substance")
GO_KEGG_boxplots(genes[1:5], "Response_to_inorganic_substance", "SmoUp_nicUp" )