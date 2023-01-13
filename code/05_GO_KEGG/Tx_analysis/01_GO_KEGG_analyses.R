
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
all_genes_EntrezID<-as.character(extract_entrez(unique(top_tx_nic$ensemblID))$entrezgene_id)


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
smoUp_smoDown_genes<-merge(smo_up_genes, smo_down_genes)
nicUp_nicDown_genes<-merge(nic_up_genes, nic_down_genes)
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
                    "Smo Up, smo Down"=smoUp_smoDown_genes_Entrez, "Nic Up, nic Down"=nicUp_nicDown_genes_Entrez) 
sigGeneList <-lapply(sigGeneList, function(x) {
  x[!is.na(x)]
})
## Background genes
geneUniverse <- all_genes_EntrezID
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

goList_intersections<-GO_KEGG(sigGeneList, geneUniverse, "intersections")
save(goList_intersections, file="processed-data/05_GO_KEGG/Tx_analysis/goList_intersections.Rdata")



##############################################################
# Txs' genes not considered or non-DE at gene level 
##############################################################

GO_KEGG_no_txs_genes<- function(expt){
  
  top_genes<-eval(parse_expr(paste("top_genes_pups_", expt, "_fitted", sep="")))
  top_tx<-eval(parse_expr(paste("top_tx_", substr(expt,1,3), sep="")))
  de_genes<-eval(parse_expr(paste("de_genes_pups_", expt, "_fitted", sep="")))
  de_tx<-eval(parse_expr(paste("de_tx_", substr(expt,1,3), sep="")))
  
  ## Txs' genes
  tx_genes<-unique(top_tx$ensemblID)
  
  ## Common genes: considered at the gene and tx level
  common_genes<-tx_genes[which(tx_genes %in% top_genes$ensemblID)]
  common_genes<-top_genes[which(top_genes$ensemblID %in% common_genes), c("ensemblID","EntrezID")]
  ## non-DE genes from the common genes containing DE txs
  non_DEG<-unique(common_genes[which(common_genes$ensemblID %in% de_tx$ensemblID & ! common_genes$ensemblID %in% de_genes$ensemblID),
                        "EntrezID"])
  ## DEG with DE txs
  DEG<-unique(common_genes[which(common_genes$ensemblID %in% de_tx$ensemblID & common_genes$ensemblID %in% de_genes$ensemblID),
                    "EntrezID"])
  ## DE txs' genes not present at the gene level
  no_present_genes<-unique(de_tx[which(!de_tx$ensemblID %in% top_genes$ensemblID), "ensemblID"])
  no_present_genes<-extract_entrez(no_present_genes)$entrezgene_id
  
  sigGeneList <- list("non-DEG with DEtxs"=non_DEG, "DEG with DEtxs"=DEG, "No at gene level"=no_present_genes) 
  sigGeneList <-lapply(sigGeneList, function(x) {
    x[!is.na(x)]
  })
  ## Background genes
  geneUniverse <- all_genes_EntrezID
  geneUniverse <- geneUniverse[!is.na(geneUniverse)]
  
  goList_noDEG<-GO_KEGG(sigGeneList, geneUniverse, paste("noDEG_", substr(expt,1,3), sep=""))
  save(goList_noDEG, file=paste("processed-data/05_GO_KEGG/Tx_analysis/goList_noDEG_", substr(expt,1,3), ".Rdata", sep=""))
  
}

## Analyses 
GO_KEGG_no_txs_genes("nicotine")
GO_KEGG_no_txs_genes("smoking")



###########################################################
# GO and KEGG of DEG vs DE exons' genes vs DE txs' genes
###########################################################

compare_DE <- function(expt){
  
  top_exons<-eval(parse_expr(paste("top_exons_", substr(expt,1,3), sep="")))
  top_genes<-eval(parse_expr(paste("top_genes_pups_", expt, "_fitted", sep="")))
  top_tx<-eval(parse_expr(paste("top_tx", substr(expt,1,3), sep="")))
  de_genes<-eval(parse_expr(paste("de_genes_pups_", expt, "_fitted", sep="")))
  de_exons<-eval(parse_expr(paste("de_exons_", substr(expt,1,3), sep="")))
  de_tx<-eval(parse_expr(paste("de_tx", substr(expt,1,3), sep="")))
  
  de_tx_Entrez <- extract_entrez(de_tx$ensemblID)$entrezgene_id
  
  ## Define groups of genes
  
  ## DEG with DE txs and exons
  DE_all <- intersect(de_genes$EntrezID, de_exons$EntrezID, de_tx_Entrez)
  ## DEG with DE exons but no DE txs
  DE_genes_exons <- intersect(de_genes$EntrezID, de_exons$EntrezID)
  DE_genes_exons <- DE_genes_exons[which(! DE_genes_exons %in% de_tx_Entrez)]
  ## DEG with DE txs but no DE exons
  DE_genes_txs <- intersect(de_genes$EntrezID, de_tx_Entrez)
  DE_genes_txs <- DE_genes_txs[which(! DE_genes_txs %in% de_exons$EntrezID)]
  ## Non-DE genes with DE exons and txs
  DE_exons_txs <- intersect(de_exons$EntrezID, de_tx_Entrez)
  DE_exons_txs <- DE_exons_txs[which(! DE_exons_txs %in% de_genes$EntrezID)]
  ## DEG without DE exons or txs
  DE_genes <- union(de_genes[which(! de_genes$EntrezID %in% de_exons$EntrezID), "EntrezID"],
                    de_genes[which(! de_genes$EntrezID %in% de_tx_Entrez)])
  ## DEE from non-DE genes without DE txs
  DE_exons <- union(de_exons[which(! de_exons$EntrezID %in% de_genes$EntrezID), "EntrezID"],
                    de_exons[which(! de_exons$EntrezID %in% de_tx_Entrez)])
  ## DEtxs from non-DE genes without DE exons
  DE_txs <- union(de_tx[which(! de_tx_Entrez %in% de_genes$EntrezID), "EntrezID"],
                  de_tx[which(! de_tx_Entrez %in% de_exons$EntrezID)])
  
  sigGeneList <- list("3 levels"= DE_all, "DEG & DEE"=DE_genes_exons, 
                      "DEG & DEtxs"= DE_genes_txs, "DEE & DEtxs"= DE_exons_txs, 
                      "DEG only"= DE_genes, "DEE only"= DE_exons, "DEtxs only"= DE_txs) 
  sigGeneList <-lapply(sigGeneList, function(x) {
    x[!is.na(x)]
  })
  ## Background genes
  geneUniverse <- all_genes_EntrezID
  geneUniverse <- geneUniverse[!is.na(geneUniverse)]
  
  goList_DE_comparisons<-GO_KEGG(sigGeneList, geneUniverse, paste("DE_comparisons_", substr(expt,1,3), sep=""))
  save(goList_DE_comparisons, file=paste("processed-data/05_GO_KEGG/Tx_analysis/goList_DE_comparisons", substr(expt,1,3), ".Rdata", sep=""))
    
}

## Analyses 
compare_DE("nicotine")
compare_DE("smoking")





## 1.1 Boxplots of txs' genes 

### 1.1.1 Genes in GO and KEGG descriptions

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





