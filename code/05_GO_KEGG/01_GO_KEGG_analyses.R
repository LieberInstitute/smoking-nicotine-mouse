
# 1. Gene Ontology and KEGG analyses


library(here)
library(SummarizedExperiment)
library(clusterProfiler)
library(org.Mm.eg.db)
library(jaffelab)

load(here("processed-data/04_DEA/top_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/top_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/de_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/results_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/results_pups_smoking_fitted.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_nicotine.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_smoking.Rdata"))


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
save(goList_global, "processed-data/05_GO_KEGG/goList_global.Rdata")



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
save(goList_nic, "processed-data/05_GO_KEGG/goList_nic.Rdata")



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
save(goList_smo, "processed-data/05_GO_KEGG/goList_smo.Rdata")



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
save(goList_intersections, "processed-data/05_GO_KEGG/goList_intersections.Rdata")

