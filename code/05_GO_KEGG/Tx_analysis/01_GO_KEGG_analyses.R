
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
load(here("processed-data/04_DEA/Exon_analysis/top_exons_nic.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/de_exons_nic.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/top_exons_smo.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/de_exons_smo.Rdata"))
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
                    "smo Up smo Down"=smoUp_smoDown_genes, "nic Up nic Down"=nicUp_nicDown_genes)
save(intersections, file="processed-data/05_GO_KEGG/Tx_analysis/intersections_txs_genes.Rdata")



## Function to do GO and KEGG analyses

GO_KEGG<- function(sigGeneList, geneUniverse, name){
  
  if (name=="intersections"){
    height=7
    width=7.5
  }
  else if (substr(name,1,14)=="DE_comparisons"){
    height=13
    width=11
  }
  else {
    height=8
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
    pdf(paste("plots/05_GO_KEGG/Tx_analysis/GO_BP_", name, ".pdf", sep=""), height = height, width = width)
    print(dotplot(goBP_Adj, title="GO Enrichment Analysis: Biological processes", font.size=9))
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
    print(dotplot(goMF_Adj, title="GO Enrichment Analysis: Molecular function", font.size=9))
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
    print(dotplot(goCC_Adj, title="GO Enrichment Analysis: Cellular components", font.size=9))
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
    print(dotplot(kegg_Adj, title="KEGG Enrichment Analysis", font.size=9))
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
  
  de_genes<-eval(parse_expr(paste("de_genes_pups_", expt, "_fitted", sep="")))
  de_exons<-eval(parse_expr(paste("de_exons_", substr(expt,1,3), sep="")))
  de_tx<-eval(parse_expr(paste("de_tx_", substr(expt,1,3), sep="")))
  
  ## Define groups of genes
  
  ## DEG with DE txs and exons
  DE_all <- intersect(intersect(de_genes$ensemblID, de_exons$ensemblID), de_tx$ensemblID)
  DE_all <- extract_entrez(unique(DE_all))$entrezgene_id
  ## DEG with DE exons but no DE txs
  DE_genes_exons <- intersect(de_genes$ensemblID, de_exons$ensemblID)
  DE_genes_exons <- DE_genes_exons[which(! DE_genes_exons %in% de_tx$ensemblID)]
  DE_genes_exons <- extract_entrez(unique(DE_genes_exons))$entrezgene_id
  ## DEG with DE txs but no DE exons
  DE_genes_txs <- intersect(de_genes$ensemblID, de_tx$ensemblID)
  DE_genes_txs <- DE_genes_txs[which(! DE_genes_txs %in% de_exons$ensemblID)]
  DE_genes_txs <- extract_entrez(unique(DE_genes_txs))$entrezgene_id
  ## Non-DE genes with DE exons and txs
  DE_exons_txs <- intersect(de_exons$ensemblID, de_tx$ensemblID)
  DE_exons_txs <- DE_exons_txs[which(! DE_exons_txs %in% de_genes$ensemblID)]
  DE_exons_txs <- extract_entrez(unique(DE_exons_txs))$entrezgene_id
  ## DEG without DE exons or txs
  DE_genes <- intersect(de_genes[which(! de_genes$ensemblID %in% de_exons$ensemblID), "ensemblID"],
                    de_genes[which(! de_genes$ensemblID %in% de_tx$ensemblID), "ensemblID"])
  DE_genes <- extract_entrez(unique(DE_genes))$entrezgene_id
  ## DEE from non-DE genes without DE txs
  DE_exons <- intersect(de_exons[which(! de_exons$ensemblID %in% de_genes$ensemblID), "ensemblID"],
                    de_exons[which(! de_exons$ensemblID %in% de_tx$ensemblID), "ensemblID"])
  DE_exons <- extract_entrez(unique(DE_exons))$entrezgene_id
  ## DEtxs from non-DE genes without DE exons
  DE_txs <- intersect(de_tx[which(! de_tx$ensemblID %in% de_genes$ensemblID), "ensemblID"],
                  de_tx[which(! de_tx$ensemblID %in% de_exons$ensemblID), "ensemblID"])
  DE_txs <- extract_entrez(unique(DE_txs))$entrezgene_id
  
  sigGeneList <- list("3 levels"= DE_all, "DEG & DEE"=DE_genes_exons, 
                      "DEG & DEtxs"= DE_genes_txs, "DEE & DEtxs"= DE_exons_txs, 
                      "DEG only"= DE_genes, "DEE only"= DE_exons, "DEtxs only"= DE_txs) 
  sigGeneList <-lapply(sigGeneList, function(x) {
    x[!is.na(x)]
  })
  ## Background genes
  geneUniverse <- union(union(top_exons_nic$EntrezID, top_genes_pups_nicotine_fitted$EntrezID),
                        all_genes_EntrezID)
  geneUniverse <- geneUniverse[!is.na(geneUniverse)]
  
  goList_DE_comparisons<-GO_KEGG(sigGeneList, geneUniverse, paste("DE_comparisons_", substr(expt,1,3), sep=""))
  save(goList_DE_comparisons, file=paste("processed-data/05_GO_KEGG/Tx_analysis/goList_DE_comparisons_", substr(expt,1,3), ".Rdata", sep=""))
    
}

## Analyses 
compare_DE("nicotine")
compare_DE("smoking")





## 1.1 Boxplots of txs' genes 

### 1.1.1 Genes in GO and KEGG descriptions

vGene_smo<-results_pups_smoking_fitted[[1]][[2]]
vGene_nic<-results_pups_nicotine_fitted[[1]][[2]]
rownames(vGene_nic$E)<-vGene_nic$genes$Symbol
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
    geom_boxplot(outlier.color = "#FFFFFFFF", width=0.35) +
    geom_jitter(aes(color=Group), shape=16, position=position_jitter(0.2), size=2.1) +
    theme_bw() +
    labs(x = "Experiment", y = "lognorm counts",
         title = paste(DEgene, ensemblID, sep=" - "),
         subtitle=" ") +
    scale_color_manual(values=c("Control" = "seashell3", "Experimental" = "orange3")) +
    scale_x_discrete(labels=c("Control"="Ctrl","Experimental"="Expt")) +
    facet_wrap(~ Expt, scales = "free") +
    scale_x_discrete(labels=c("Ctrl", "Expt")) +
    theme(plot.margin=unit (c (1,1.5,1,1), 'cm'), 
          legend.position = "none",
          plot.title = element_text(hjust=0.5, size=12, face="bold"), 
          plot.subtitle = element_text(size=17), 
          axis.title = element_text(size = (12)),
          axis.text = element_text(size = 10.5)) 
  
  p <-ggdraw(p) + 
    draw_label(paste("FDR:", q_value_nic), x = 0.35, y = 0.83, size=9, color = "darkslategray") +
    draw_label(paste("FC:", FC_nic), x = 0.35, y = 0.80, size=9, color = "darkslategray") +
    draw_label(paste("FDR:", q_value_smo), x = 0.72, y = 0.83, size=9, color = "darkslategray") +
    draw_label(paste("FC:", FC_smo), x = 0.71, y = 0.80, size=9, color = "darkslategray") 
  
  return(p)
  
}


## Extract genes from each BP, CC, MF and KEGG descriptions
GO_KEGG_genes<- function(golist, term, cluster, description){
  
  GOdata<-as.data.frame(eval(parse_expr(paste(golist, "$", term, sep=""))))
  genes<-unique(GOdata[which(GOdata$Description==description & GOdata$Cluster==cluster), "geneID"])
  genes<-strsplit(genes, "/")
  genes<-unique(unlist(genes))
  ## Retain only genes considered at the gene level
  if (! term=="KEGG"){
    genes<-genes[which(genes %in% top_genes_pups_nicotine_fitted$Symbol)]
  }
  
  else if (term=="KEGG"){
    genes<-genes[which(genes %in% top_genes_pups_nicotine_fitted$EntrezID)]
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
  ggsave(here(paste("plots/05_GO_KEGG/Tx_analysis/Top", length(DEG_list), "_", description,"_boxplots_",cluster, 
                    ".pdf", sep="")), width = 40, height = 20, units = "cm") 
  
}



## Boxplots 

## 1. Cellular components

## Genes in SNARE complex
GO_genes<-GO_KEGG_genes("goList_global", "CC", "up", "SNARE complex")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "SNARE_complex", "up")

GO_genes<-GO_KEGG_genes("goList_smo", "CC", "up", "SNARE complex")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "SNARE_complex", "smo_up")

GO_genes<-GO_KEGG_genes("goList_intersections", "CC", "Only up smo", "SNARE complex")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "SNARE_complex", "Only_up_smo")


## Genes in transport vesicle
GO_genes<-GO_KEGG_genes("goList_global", "CC", "up", "transport vesicle")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "transport_vesicle", "up")

GO_genes<-GO_KEGG_genes("goList_smo", "CC", "up", "transport vesicle")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "transport_vesicle", "smo_up")

GO_genes<-GO_KEGG_genes("goList_intersections", "CC", "Only up smo", "transport vesicle")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "transport_vesicle", "Only_up_smo")



## 2. Pathways

## Genes involved in Pathways of neurodegeneration − multiple diseases
GO_genes<-GO_KEGG_genes("goList_global", "KEGG", "up", "Pathways of neurodegeneration - multiple diseases - Mus musculus (house mouse)")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Neurodegeneration", "up")

GO_genes<-GO_KEGG_genes("goList_smo", "KEGG", "up", "Pathways of neurodegeneration - multiple diseases - Mus musculus (house mouse)")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Neurodegeneration", "smo_up")


## Genes involved in Parkinson disease
GO_genes<-GO_KEGG_genes("goList_global", "KEGG", "up", "Parkinson disease - Mus musculus (house mouse)")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Parkinson_disease", "up")

GO_genes<-GO_KEGG_genes("goList_smo", "KEGG", "up", "Parkinson disease - Mus musculus (house mouse)")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Parkinson_disease", "smo_up")
 
GO_genes<-GO_KEGG_genes("goList_intersections", "KEGG", "Only up smo", "Parkinson disease - Mus musculus (house mouse)")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Parkinson_disease", "Only_up_smo")


## Genes involved in Prion disease
GO_genes<-GO_KEGG_genes("goList_global", "KEGG", "up", "Prion disease - Mus musculus (house mouse)")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Prion_disease", "up")

GO_genes<-GO_KEGG_genes("goList_smo", "KEGG", "up", "Prion disease - Mus musculus (house mouse)")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Prion_disease", "smo_up")

GO_genes<-GO_KEGG_genes("goList_intersections", "KEGG", "Only up smo", "Prion disease - Mus musculus (house mouse)")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Prion_disease", "Only_up_smo")


## Genes involved in Huntington disease
GO_genes<-GO_KEGG_genes("goList_intersections", "KEGG", "Only up smo", "Huntington disease - Mus musculus (house mouse)")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Huntington_disease", "Only_up_smo")


## Genes involved in SNARE interactions in vesicular transport
GO_genes<-GO_KEGG_genes("goList_intersections", "KEGG", "Smo up, nic down", "SNARE interactions in vesicular transport")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "SNARE_interactions", "smoUp_nicDown")







## Reproducibility information

options(width = 120)
session_info()

# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.0 (2023-04-21)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Mexico_City
# date     2024-01-02
# rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
# pandoc   NA
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
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
# ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

