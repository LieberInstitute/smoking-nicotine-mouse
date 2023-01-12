
# 1. Gene Ontology and KEGG analyses of DE exons' genes


library(here)
library(SummarizedExperiment)
library(clusterProfiler)
library(org.Mm.eg.db)
library(jaffelab)
library(ggplot2)
library(cowplot)
library(rlang)
library(sessioninfo)


load(here("processed-data/04_DEA/Exon_analysis/de_exons_nic.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/de_exons_smo.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/top_exons_nic.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/top_exons_smo.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/results_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/results_pups_smoking_fitted.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_smoking.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_nicotine.Rdata"))


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

## Analyses 
GO_KEGG_no_exons_genes("nicotine")
GO_KEGG_no_exons_genes("smoking")





## 1.1 Boxplots of exons' genes 

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


## Extract genes from each BP, CC, MF and KEGG descriptions
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
  ggsave(here(paste("plots/05_GO_KEGG/Exon_analysis/Top", length(DEG_list), "_", description,"_boxplots_",cluster, 
                    ".pdf", sep="")), width = 40, height = 25, units = "cm") 
  
}



## Boxplots 

## 1. Biological processes

## Genes involved in regulation of synaptic plasticity
GO_genes<-GO_KEGG_genes("goList_global", "BP", "up", "regulation of synaptic plasticity")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "regulation_of_synaptic_plasticity", "up")

GO_genes<-GO_KEGG_genes("goList_smo", "BP", "up", "regulation of synaptic plasticity")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "regulation_of_synaptic_plasticity", "smo_up")



## 2. Cellular components

## Genes in secretory vesicles
GO_genes<-GO_KEGG_genes("goList_global", "CC", "up", "secretory vesicle")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "secretory_vesicle", "up")

GO_genes<-GO_KEGG_genes("goList_smo", "CC", "up", "secretory vesicle")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "secretory_vesicle", "smo_up")

GO_genes<-GO_KEGG_genes("goList_intersections", "CC", "Only up smo", "secretory vesicle")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "secretory_vesicle", "intersections_smo_up")


## Genes in transport vesicles
GO_genes<-GO_KEGG_genes("goList_global", "CC", "up", "transport vesicle")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "transport_vesicle", "up")

GO_genes<-GO_KEGG_genes("goList_smo", "CC", "up", "transport vesicle")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "transport_vesicle", "smo_up")

GO_genes<-GO_KEGG_genes("goList_intersections", "CC", "Only up smo", "transport vesicle")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "transport_vesicle", "intersections_smo_up")


## Genes in postsynaptic specialization
GO_genes<-GO_KEGG_genes("goList_global", "CC", "up", "postsynaptic specialization")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "postsynaptic_specialization", "up")

GO_genes<-GO_KEGG_genes("goList_smo", "CC", "up", "postsynaptic specialization")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "postsynaptic_specialization", "smo_up")


## Genes in asymmetric synapses
GO_genes<-GO_KEGG_genes("goList_global", "CC", "up", "asymmetric synapse")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "asymmetric_synapse", "up")

GO_genes<-GO_KEGG_genes("goList_smo", "CC", "up", "asymmetric synapse")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "asymmetric_synapse", "smo_up")


## Genes in exocytic vesicles
GO_genes<-GO_KEGG_genes("goList_global", "CC", "up", "exocytic vesicle")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "exocytic_vesicle", "up")

GO_genes<-GO_KEGG_genes("goList_smo", "CC", "up", "exocytic vesicle")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "exocytic_vesicle", "smo_up")

GO_genes<-GO_KEGG_genes("goList_intersections", "CC", "Only up smo", "exocytic vesicle")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "exocytic_vesicle", "intersections_smo_up")


## Genes in secretory granules
GO_genes<-GO_KEGG_genes("goList_intersections", "CC", "Only up smo", "secretory granule")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "secretory_granule", "intersections_smo_up")


## Genes in the SNARE complex
GO_genes<-GO_KEGG_genes("goList_intersections", "CC", "Only up smo", "SNARE complex")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "SNARE_complex", "intersections_smo_up")


## Genes in synaptic membranes
GO_genes<-GO_KEGG_genes("goList_intersections", "CC", "Smo Up, smo Down", "synaptic membrane")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "synaptic_membrane", "Smo_up_Smo_down")


## Genes in postsynaptic membranes
GO_genes<-GO_KEGG_genes("goList_intersections", "CC", "Smo Up, smo Down", "postsynaptic membrane")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "postsynaptic_membrane", "Smo_up_Smo_down")


## Genes in postsynaptic density
GO_genes<-GO_KEGG_genes("goList_intersections", "CC", "Smo Up, smo Down", "postsynaptic density")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "postsynaptic_density", "Smo_up_Smo_down")


## Genes in asymmetric synapse
GO_genes<-GO_KEGG_genes("goList_intersections", "CC", "Smo Up, smo Down", "asymmetric synapse")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "asymmetric_synapse", "Smo_up_Smo_down")


## Genes in postsynaptic specialization
GO_genes<-GO_KEGG_genes("goList_intersections", "CC", "Smo Up, smo Down", "postsynaptic specialization")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "postsynaptic_specialization", "Smo_up_Smo_down")



## 3. Molecular function

## Genes with voltage−gated ion channel activity
GO_genes<-GO_KEGG_genes("goList_smo", "MF", "up", "voltage-gated ion channel activity")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "voltage−gated_ion_channel_activity", "smo_up")

GO_genes<-GO_KEGG_genes("goList_intersections", "MF", "Only up smo", "voltage-gated ion channel activity")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "voltage−gated_ion_channel_activity", "only_smo_up")


## Genes with voltage−gated channel activity
GO_genes<-GO_KEGG_genes("goList_smo", "MF", "up", "voltage-gated channel activity")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "voltage−gated_channel_activity", "smo_up")

GO_genes<-GO_KEGG_genes("goList_intersections", "MF", "Only up smo", "voltage-gated channel activity")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "voltage−gated_channel_activity", "only_smo_up")


## Genes with gated channel activity
GO_genes<-GO_KEGG_genes("goList_smo", "MF", "up", "gated channel activity")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "gated_channel_activity", "smo_up")

GO_genes<-GO_KEGG_genes("goList_intersections", "MF", "Only up smo", "gated channel activity")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "gated_channel_activity", "only_smo_up")


## Genes with voltage−gated cation channel activity
GO_genes<-GO_KEGG_genes("goList_intersections", "MF", "Only up smo", "voltage-gated cation channel activity")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "voltage−gated_cation_channel_activity", "only_smo_up")







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








