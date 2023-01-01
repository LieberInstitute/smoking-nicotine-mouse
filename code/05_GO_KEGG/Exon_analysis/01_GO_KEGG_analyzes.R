
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
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata"))


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
geneUniverse <- union(top_exons_nic$EntrezID,
                      top_exons_smo$EntrezID)
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
geneUniverse <- unique(top_exons_nic$EntrezID)
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

goList_nic<-GO_KEGG(sigGeneList, geneUniverse, "nicotine")
save(goList_nic, file="processed-data/05_GO_KEGG/Exon_analysis/goList_nic.Rdata")



####################################
# Smoking Up/Down DE exons' genes
####################################

## List of sets
sigGeneList <- list("up"=smo_up_genes$EntrezID, "down"=smo_down_genes$EntrezID) 
sigGeneList <-lapply(sigGeneList, function(x) {
  x[!is.na(x)]
})
## Background genes
geneUniverse <- unique(top_exons_smo$EntrezID)
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

goList_smo<-GO_KEGG(sigGeneList, geneUniverse, "smoking")
save(goList_smo, file="processed-data/05_GO_KEGG/Exon_analysis/goList_smo.Rdata")



#######################################################
# Up/Down Nicotine VS Up/Down Smoking DE exons' genes
#######################################################

sigGeneList <- list("Only up nic"=only_up_nic_genes$EntrezID, "Only up smo"=only_up_smo_genes$EntrezID,
                    "Only down nic"=only_down_nic_genes$EntrezID, "Only down smo"=only_down_smo_genes$EntrezID,
                    "Smo up, nic down"=smoUp_nicDown_genes$EntrezID, "Smo down, nic up"=smoDown_nicUp_genes$EntrezID,
                    "Smo up, nic up"=smoUp_nicUp_genes$EntrezID, "Smo down, nic down"=smoDown_nicDown_genes$EntrezID,
                    "Smo Up, smo Down"=smoUp_smoDown$EntrezID, "Nic Up, nic Down"=nicUp_nicDown$EntrezID) 
sigGeneList <-lapply(sigGeneList, function(x) {
  x[!is.na(x)]
})
## Background genes
geneUniverse <- union(top_exons_nic$EntrezID,
                      top_exons_smo$EntrezID)
geneUniverse <- geneUniverse[!is.na(geneUniverse)]

goList_intersections<-GO_KEGG(sigGeneList, geneUniverse, "intersections")
save(goList_intersections, file="processed-data/05_GO_KEGG/Exon_analysis/goList_intersections.Rdata")



##################################################################
# Exons' genes not considered at the gene level or non-DE genes 
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

## Analyzes 
GO_KEGG_no_exons_genes("nicotine")
GO_KEGG_no_exons_genes("smoking")





## 1.1 Boxplots of exons' genes 

### 1.1.1 Genes in GO and KEGG descriptions

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







## Reproducibility information

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
# date     2022-12-18
# rstudio  2022.07.2+576 Spotted Wakerobin (desktop)
# pandoc   NA








