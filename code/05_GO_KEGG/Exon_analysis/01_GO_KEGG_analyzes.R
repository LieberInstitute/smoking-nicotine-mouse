
# 1. Gene Ontology and KEGG analyzes of DE exons' genes


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


load(here("processed-data/04_DEA/Exon_analysis/de_exons_nic.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/de_exons_smo.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/top_exons_nic.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/top_exons_smo.Rdata"))



## Groups of DE exons' genes

nic_up_genes<-unique(de_exons_nic[which(de_exons_nic$logFC>0),c("EntrezID", "Symbol", "ensemblID")])
smo_up_genes<-unique(de_exons_smo[which(de_exons_smo$logFC>0),c("EntrezID", "Symbol", "ensemblID")])
all_up_genes<-unique(rbind(nic_up_genes, smo_up_genes))

nic_down_genes<-unique(de_exons_nic[which(de_exons_nic$logFC<0),c("EntrezID", "Symbol", "ensemblID")])
smo_down_genes<-unique(de_exons_smo[which(de_exons_smo$logFC<0),c("EntrezID", "Symbol", "ensemblID")])
all_down_genes<-unique(rbind(nic_down_genes, smo_down_genes))

all_nic_genes<-unique(rbind(nic_up_genes, nic_down_genes))
all_smo_genes<-unique(rbind(smo_up_genes, smo_down_genes))
all<-unique(rbind(all_nic_genes, all_smo_genes))
save(all, file="processed-data/05_GO_KEGG/Exon_analysis/all_DE_exons_genes.Rdata")

## Intersections between groups
smoUp_nicUp_genes<-merge(nic_up_genes, smo_up_genes)
smoDown_nicDown_genes<-merge(nic_down_genes, smo_down_genes)
smoUp_nicDown_genes<-merge(nic_down_genes, smo_up_genes)
smoDown_nicUp_genes<-merge(nic_up_genes, smo_down_genes)
smoUp_smoDown<-merge(smo_up_genes, smo_down_genes)
nicUp_nicDown<-merge(nic_up_genes, nic_down_genes)
only_up_nic_genes<-merge(merge(nic_up_genes[which(! nic_up_genes$ensemblID %in% nic_down_genes$ensemblID),],
                               nic_up_genes[which(! nic_up_genes$ensemblID %in% smo_up_genes$ensemblID),]),
                               nic_up_genes[which(! nic_up_genes$ensemblID %in% smo_down_genes$ensemblID),])
only_up_smo_genes<-merge(merge(smo_up_genes[which(! smo_up_genes$ensemblID %in% smo_down_genes$ensemblID),],
                               smo_up_genes[which(! smo_up_genes$ensemblID %in% nic_up_genes$ensemblID),]),
                               smo_up_genes[which(! smo_up_genes$ensemblID %in% nic_down_genes$ensemblID),])
only_down_nic_genes<-merge(merge(nic_down_genes[which(! nic_down_genes$ensemblID %in% nic_up_genes$ensemblID),],
                                 nic_down_genes[which(! nic_down_genes$ensemblID %in% smo_up_genes$ensemblID),]),
                                 nic_down_genes[which(! nic_down_genes$ensemblID %in% smo_down_genes$ensemblID),])
only_down_smo_genes<-merge(merge(smo_down_genes[which(! smo_down_genes$ensemblID %in% smo_up_genes$ensemblID),],
                                 smo_down_genes[which(! smo_down_genes$ensemblID %in% nic_up_genes$ensemblID),]),
                                 smo_down_genes[which(! smo_down_genes$ensemblID %in% nic_down_genes$ensemblID),])

intersections<-list("only up nic"=only_up_nic_genes, "only up smo"=only_up_smo_genes, 
                    "only down nic"=only_down_nic_genes, "only down smo"=only_down_smo_genes, 
                    "smo Up nic Up"=smoUp_nicUp_genes, "smo Down nic Down"=smoDown_nicDown_genes, 
                    "smo Up nic Down"=smoUp_nicDown_genes, "smo Down nic Up"=smoDown_nicUp_genes,
                    "smo Up smo Down"=smoUp_smoDown, "nic Up nic Down"=nicUp_nicDown)
save(intersections, file="processed-data/05_GO_KEGG/Exon_analysis/intersections_exons_genes.Rdata")



## Function to do GO and KEGG analyzes

GO_KEGG<- function(sigGeneList, geneUniverse, name){
  
  if (name=="intersections"){
    height=17
    width=15
  }
  else {
    height=7
    width=7
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
  if(!is.null(goBP_Adj)){
    pdf(paste("plots/05_GO_KEGG/Exon_analysis/GO_BP_", name, ".pdf", sep=""), height = height, width = width)
    print(dotplot(goBP_Adj, title="GO Enrichment Analysis: Biological processes"))
    dev.off()
  }

  
  
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
  if (!is.null(goMF_Adj)){
    pdf(paste("plots/05_GO_KEGG/Exon_analysis/GO_MF_", name, ".pdf", sep=""), height = height, width = width)
    print(dotplot(goMF_Adj, title="GO Enrichment Analysis: Molecular function"))
    dev.off()
  }

  
  
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
  if(!is.null(goCC_Adj)){
    pdf(paste("plots/05_GO_KEGG/Exon_analysis/GO_CC_", name, ".pdf", sep=""), height = height, width = width)
    print(dotplot(goCC_Adj, title="GO Enrichment Analysis: Cellular components"))
    dev.off()
  }

  
  
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
  if(!is.null(kegg_Adj)){
    pdf(paste("plots/05_GO_KEGG/Exon_analysis/KEGG_", name, ".pdf", sep=""), height = height, width = width)
    print(dotplot(kegg_Adj, title="KEGG Enrichment Analysis"))
    dev.off()
  }

  
  
  goList <- list(
    BP = goBP_Adj,
    MF = goMF_Adj,
    CC = goCC_Adj,
    KEGG = kegg_Adj
  )
  
  return(goList)
}





###########################
# Up/Down DE exons' genes
###########################

## List of sets 
sigGeneList <- list("up"=all_up_genes$EntrezID, "down"=all_down_genes$EntrezID) 
sigGeneList <-lapply(sigGeneList, function(x) {
  x[!is.na(x)]
})
## Background genes
geneUniverse <- as.character(union(top_exons_nic$EntrezID,
                                   top_exons_smo$EntrezID))
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

goList_global<-GO_KEGG(sigGeneList, geneUniverse, "global")
save(goList_global, file="processed-data/05_GO_KEGG/Exon_analysis/goList_global.Rdata")



####################################
# Nicotine Up/Down DE exons' genes
####################################

## List of sets
sigGeneList <- list("up"=nic_up_genes$EntrezID, "down"=nic_down_genes$EntrezID) 
sigGeneList <-lapply(sigGeneList, function(x) {
  x[!is.na(x)]
})
## Background genes
geneUniverse <- as.character(unique(top_exons_nic$EntrezID))
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

goList_nic<-GO_KEGG(sigGeneList, geneUniverse, "nicotine")
save(goList_nic, file="processed-data/05_GO_KEGG/Gene_analysis/goList_nic.Rdata")



####################################
# Smoking Up/Down DE exons' genes
####################################

## List of sets
sigGeneList <- list("up"=smo_up_genes$EntrezID, "down"=smo_down_genes$EntrezID) 
sigGeneList <-lapply(sigGeneList, function(x) {
  x[!is.na(x)]
})
## Background genes
geneUniverse <- as.character(unique(top_exons_smo$EntrezID))
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

goList_smo<-GO_KEGG(sigGeneList, geneUniverse, "smoking")
save(goList_smo, file="processed-data/05_GO_KEGG/Gene_analysis/goList_smo.Rdata")






























## GO and KEGG of Up/Down nic VS Up/Down smo DE exons' genes


## GO and KEGG of exons' genes not considered at the gene level or from non-DE genes and compare



## Boxplots of top 6 genes in each group

## Boxplots of top genes in certain pathways 




## Compare BP, MF, etc of DEG and DE exons' genes
## Compare DEG and DE exons' genes in certain pathways or processes