
# 1. Gene Ontology and KEGG analyses of DE txs' genes


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


load(here("processed-data/04_DEA/Tx_analysis/de_tx_nic.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/de_tx_smo.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/top_tx_nic.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/top_tx_smo.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/results_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/results_pups_smoking_fitted.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_smoking.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_nicotine.Rdata"))


## Add ensembl ID (without period) of the txs' genes
de_tx_nic$ensemblID<-sapply(de_tx_nic$ensembl_id, function(x){strsplit(x, "[.]")[[1]][1]})
de_tx_smo$ensemblID<-sapply(de_tx_smo$ensembl_id, function(x){strsplit(x, "[.]")[[1]][1]})
top_tx_nic$ensemblID<-sapply(top_tx_nic$ensembl_id, function(x){strsplit(x, "[.]")[[1]][1]})
top_tx_smo$ensemblID<-sapply(top_tx_smo$ensembl_id, function(x){strsplit(x, "[.]")[[1]][1]})

## Extract entrez ID 
extract_entrez <- function(gene_list){
  entrez_ids <-biomart(genes  = gene_list,
                   mart       = "ENSEMBL_MART_ENSEMBL",
                   dataset    = "mmusculus_gene_ensembl",
                   attributes = c("entrezgene_id"),
                   filters    = "ensembl_gene_id")
  return(entrez_ids)
}

## Add entrez ID of all txs' genes
all_genes_EntrezID<-extract_entrez(top_tx_nic$ensemblID)$entrezgene_id


## Groups of DE txs' genes

nic_up_genes<-unique(de_tx_nic[which(de_tx_nic$logFC>0),c("Symbol", "ensemblID")])
smo_up_genes<-unique(de_tx_smo[which(de_tx_smo$logFC>0),c("Symbol", "ensemblID")])
all_up_genes<-unique(rbind(nic_up_genes, smo_up_genes))

nic_down_genes<-unique(de_tx_nic[which(de_tx_nic$logFC<0),c("Symbol", "ensemblID")])
smo_down_genes<-unique(de_tx_smo[which(de_tx_smo$logFC<0),c("Symbol", "ensemblID")])
all_down_genes<-unique(rbind(nic_down_genes, smo_down_genes))

all_nic_genes<-unique(rbind(nic_up_genes, nic_down_genes))
all_smo_genes<-unique(rbind(smo_up_genes, smo_down_genes))
all<-unique(rbind(all_nic_genes, all_smo_genes))
save(all, file="processed-data/05_GO_KEGG/Tx_analysis/all_DE_txs_genes.Rdata")

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
save(intersections, file="processed-data/05_GO_KEGG/Tx_analysis/intersections_txs_genes.Rdata")



## Function to do GO and KEGG analyses

GO_KEGG<- function(sigGeneList, geneUniverse, name){
  
  if (name=="intersections"){
    height=17
    width=15
  }
  else {
    height=8.5
    width=9.5
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
    pdf(paste("plots/05_GO_KEGG/Tx_analysis/GO_BP_", name, ".pdf", sep=""), height = height, width = width)
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
    pdf(paste("plots/05_GO_KEGG/Tx_analysis/GO_MF_", name, ".pdf", sep=""), height = height, width = width)
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
    pdf(paste("plots/05_GO_KEGG/Tx_analysis/GO_CC_", name, ".pdf", sep=""), height = height, width = width)
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
    pdf(paste("plots/05_GO_KEGG/Tx_analysis/KEGG_", name, ".pdf", sep=""), height = height, width = width)
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
# Up/Down DE txs' genes
###########################

## List of sets 
all_up_genes_Entrez<-extract_entrez(all_up_genes$ensemblID)$entrezgene_id
all_down_genes_Entrez <-extract_entrez(all_down_genes$ensemblID)$entrezgene_id
sigGeneList <- list("up"=all_up_genes_Entrez, "down"=all_down_genes_Entrez) 
sigGeneList <-lapply(sigGeneList, function(x) {
  x[!is.na(x)]
})
## Background genes
geneUniverse <- all_genes_EntrezID
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

goList_global<-GO_KEGG(sigGeneList, geneUniverse, "global")
save(goList_global, file="processed-data/05_GO_KEGG/Tx_analysis/goList_global.Rdata")



####################################
# Nicotine Up/Down DE txs' genes
####################################

## List of sets
nic_up_genes_Entrez<-extract_entrez(nic_up_genes$ensemblID)$entrezgene_id
nic_down_genes_Entrez <-extract_entrez(nic_down_genes$ensemblID)$entrezgene_id
sigGeneList <- list("up"=nic_up_genes_Entrez, "down"=nic_down_genes_Entrez) 
sigGeneList <-lapply(sigGeneList, function(x) {
  x[!is.na(x)]
})
## Background genes
geneUniverse <- all_genes_EntrezID
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

goList_nic<-GO_KEGG(sigGeneList, geneUniverse, "nicotine")
save(goList_nic, file="processed-data/05_GO_KEGG/Tx_analysis/goList_nic.Rdata")



####################################
# Smoking Up/Down DE txs' genes
####################################

## List of sets
smo_up_genes_Entrez<-extract_entrez(smo_up_genes$ensemblID)$entrezgene_id
smo_down_genes_Entrez <-extract_entrez(smo_down_genes$ensemblID)$entrezgene_id
sigGeneList <- list("up"=smo_up_genes_Entrez, "down"=smo_down_genes_Entrez) 
sigGeneList <-lapply(sigGeneList, function(x) {
  x[!is.na(x)]
})
## Background genes
geneUniverse <- all_genes_EntrezID
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

goList_smo<-GO_KEGG(sigGeneList, geneUniverse, "smoking")
save(goList_smo, file="processed-data/05_GO_KEGG/Tx_analysis/goList_smo.Rdata")



#######################################################
# Up/Down Nicotine VS Up/Down Smoking DE txs' genes
#######################################################

only_up_nic_genes_Entrez <-extract_entrez(only_up_nic_genes$ensemblID)$entrezgene_id
only_up_smo_genes_Entrez <-extract_entrez(only_up_smo_genes$ensemblID)$entrezgene_id
only_down_nic_genes_Entrez <-extract_entrez(only_down_nic_genes$ensemblID)$entrezgene_id
only_down_smo_genes_Entrez <-extract_entrez(only_down_smo_genes$ensemblID)$entrezgene_id
smoUp_nicDown_genes_Entrez <-extract_entrez(smoUp_nicDown_genes$ensemblID)$entrezgene_id
smoDown_nicUp_genes_Entrez <-extract_entrez(smoDown_nicUp_genes$ensemblID)$entrezgene_id
smoUp_nicUp_genes_Entrez <-extract_entrez(smoUp_nicUp_genes$ensemblID)$entrezgene_id
smoDown_nicDown_genes_Entrez <-extract_entrez(smoDown_nicDown_genes$ensemblID)$entrezgene_id
smoUp_smoDown_genes_Entrez <-extract_entrez(smoUp_smoDown_genes$ensemblID)$entrezgene_id
nicUp_nicDown_genes_Entrez <-extract_entrez(nicUp_nicDown_genes$ensemblID)$entrezgene_id

sigGeneList <- list("Only up nic"=only_up_nic_genes_Entrez, "Only up smo"=only_up_smo_genes_Entrez,
                    "Only down nic"=only_down_nic_genes_Entrez, "Only down smo"=only_down_smo_genes_Entrez,
                    "Smo up, nic down"=smoUp_nicDown_genes_Entrez, "Smo down, nic up"=smoDown_nicUp_genes_Entrez,
                    "Smo up, nic up"=smoUp_nicUp_genes_Entrez, "Smo down, nic down"=smoDown_nicDown_genes_Entrez,
                    "Smo Up, smo Down"=smoUp_smoDown_Entrez, "Nic Up, nic Down"=nicUp_nicDown_Entrez) 
sigGeneList <-lapply(sigGeneList, function(x) {
  x[!is.na(x)]
})
## Background genes
geneUniverse <- all_genes_EntrezID
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

goList_intersections<-GO_KEGG(sigGeneList, geneUniverse, "intersections")
save(goList_intersections, file="processed-data/05_GO_KEGG/Tx_analysis/goList_intersections.Rdata")



##################################################################
# Txs' genes not considered at the gene or exon levels or non-DE genes or DE  
##################################################################

GO_KEGG_no_exons_genes<- function(expt){
  top_exons<-eval(parse_expr(paste("top_exons_", substr(expt,1,3), sep="")))
  top_genes<-eval(parse_expr(paste("top_genes_pups_", expt, "_fitted", sep="")))
  de_genes<-eval(parse_expr(paste("de_genes_pups_", expt, "_fitted", sep="")))
  de_exons<-eval(parse_expr(paste("de_exons_", substr(expt,1,3), sep="")))
  
  ## Exons' genes
  exons_genes<-unique(top_exons$ensemblID)
  
  ## Common genes: considered at the gene and exon level
  common_genes<-exons_genes[which(exons_genes %in% top_genes$ensemblID)]
  common_genes<-top_genes[which(top_genes$ensemblID %in% common_genes), c("ensemblID","EntrezID")]
  ## non-DE genes from the common genes containing DE exons
  non_DEG<-common_genes[which(common_genes$ensemblID %in% de_exons$ensemblID & ! common_genes$ensemblID %in% de_genes$ensemblID),
                        "EntrezID"]
  ## DEG with DE exons
  DEG<-common_genes[which(common_genes$ensemblID %in% de_exons$ensemblID & common_genes$ensemblID %in% de_genes$ensemblID),
                    "EntrezID"]
  ## DE exons' genes not present at the gene level
  no_present_genes<-unique(de_exons[which(!de_exons$ensemblID %in% top_genes$ensemblID), "EntrezID"])
  
  sigGeneList <- list("non-DEG with DEE"=non_DEG, "DEG with DEE"=DEG, "No at gene level"=no_present_genes) 
  sigGeneList <-lapply(sigGeneList, function(x) {
    x[!is.na(x)]
  })
  ## Background genes
  geneUniverse <- unique(top_exons$EntrezID)
  geneUniverse <- geneUniverse[!is.na(geneUniverse)]
  
  goList_noDEG<-GO_KEGG(sigGeneList, geneUniverse, paste("noDEG_", substr(expt,1,3), sep=""))
  save(goList_noDEG, file=paste("processed-data/05_GO_KEGG/Exon_analysis/goList_noDEG_", substr(expt,1,3), ".Rdata", sep=""))
  
}

## Analyses 
GO_KEGG_no_exons_genes("nicotine")
GO_KEGG_no_exons_genes("smoking")

