
# 1. Gene Ontology and KEGG analyses of DEG 


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


## GO and KEGG analyses for DEG from fitted model only

## Groups of DEG
up_nic<-de_genes_pups_nicotine_fitted[de_genes_pups_nicotine_fitted$logFC>0, c("EntrezID", "Symbol", "ensemblID", "gencodeID")]
up_smo<-de_genes_pups_smoking_fitted[de_genes_pups_smoking_fitted$logFC>0, c("EntrezID", "Symbol", "ensemblID", "gencodeID")]
all_up<-unique(rbind(up_nic, up_smo))

down_nic<-de_genes_pups_nicotine_fitted[de_genes_pups_nicotine_fitted$logFC<0, c("EntrezID", "Symbol", "ensemblID", "gencodeID")]
down_smo<-de_genes_pups_smoking_fitted[de_genes_pups_smoking_fitted$logFC<0, c("EntrezID", "Symbol", "ensemblID", "gencodeID")]
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


## Function to do GO and KEGG analyses

GO_KEGG<- function(sigGeneList, geneUniverse, name){
  
  if (name=="intersections"){
    height=12.5
    width=10.5
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
  pdf(paste("plots/05_GO_KEGG/Gene_analysis/GO_BP_", name, ".pdf", sep=""), height = height, width = width)
  print(dotplot(goBP_Adj, title="GO Enrichment Analysis: Biological processes", font.size=9))
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
  print(dotplot(goMF_Adj, title="GO Enrichment Analysis: Molecular function", font.size=9))
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
  print(dotplot(goCC_Adj, title="GO Enrichment Analysis: Cellular components", font.size=9))
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
  print(dotplot(kegg_Adj, title="KEGG Enrichment Analysis", font.size=9))
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
         width = 40, height = 20, units = "cm") 
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
                    ".pdf", sep="")), width = 40, height = 20, units = "cm") 
  
}



## Boxplots for genes up/down in nic/smo clusters with enriched terms

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

## Genes of the SMN−Sm protein complex
GO_genes<-GO_KEGG_genes("goList_intersections", "CC", "Smo up, nic up", "SMN-Sm protein complex")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "SMN_Sm_protein_complex", "smoUp_nicUp")

## Genes of the postsynaptic endosome
GO_genes<-GO_KEGG_genes("goList_intersections", "CC", "Smo up, nic up", "postsynaptic endosome")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "postsynaptic_endosome", "smoUp_nicUp")


## Genes in asymmetric synapses
GO_genes<-GO_KEGG_genes("goList_nic", "CC", "up", "asymmetric synapse")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Asymmetric_synapse", "up_nic")



## 2. Molecular functions

## Genes with heat shock protein binding activity
GO_genes<-GO_KEGG_genes("goList_intersections", "MF", "Smo down, nic up", "heat shock protein binding")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "HeatShock_protein_binding", "smoDown_nicUp")

## Genes with SNAP receptor activity
GO_genes<-GO_KEGG_genes("goList_intersections", "MF", "Smo up, nic down", "SNAP receptor activity")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "SNAP_receptor_activity", "smoUp_nicDown")


## 3. Pathways

## Genes involved in Parkinson disease
GO_genes<-GO_KEGG_genes("goList_global", "KEGG", "up", "Parkinson disease - Mus musculus (house mouse)")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Parkinson_disease", "up")

## Genes involved in dopaminergic synapses
GO_genes<-GO_KEGG_genes("goList_intersections", "KEGG", "Only up nic", "Dopaminergic synapse - Mus musculus (house mouse)")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Dopaminergic_synapse", "Only_up_nic")

GO_genes<-GO_KEGG_genes("goList_nic", "KEGG", "up", "Dopaminergic synapse - Mus musculus (house mouse)")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Dopaminergic_synapse", "up_nic")

## Genes involved in long−term depression
GO_genes<-GO_KEGG_genes("goList_intersections", "KEGG", "Only up nic", "Long-term depression - Mus musculus (house mouse)")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "Long-term_depression", "Only_up_nic")

## Genes involved in SNARE interactions in vesicle transport
GO_genes<-GO_KEGG_genes("goList_intersections", "KEGG", "Smo up, nic down", "SNARE interactions in vesicular transport - Mus musculus (house mouse)")
top_DEG<-extract_top_genes(GO_genes)
GO_KEGG_boxplots(top_DEG, "SNARE_int_ves_transport", "smoUp_nicDown")








## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.0 (2023-04-21)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Mexico_City
# date     2024-01-01
# rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
# pandoc   NA
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                * version   date (UTC) lib source
# AcidBase                 0.7.3     2023-12-15 [1] Bioconductor
# AcidCLI                  0.3.0     2023-10-03 [1] Bioconductor
# AcidGenerics             0.7.6     2023-12-15 [1] Bioconductor
# AcidGenomes            * 0.7.2     2023-12-06 [1] Bioconductor
# AcidPlyr                 0.5.3     2023-12-13 [1] local
# AnnotationDbi          * 1.63.2    2023-07-03 [1] Bioconductor
# AnnotationFilter       * 1.26.0    2023-10-26 [1] Bioconductor
# AnnotationHub            3.9.1     2023-06-14 [1] Bioconductor
# ape                      5.7-1     2023-03-13 [1] CRAN (R 4.3.0)
# aplot                    0.1.10    2023-03-08 [1] CRAN (R 4.3.0)
# Biobase                * 2.61.0    2023-06-02 [1] Bioconductor
# BiocFileCache            2.9.1     2023-07-14 [1] Bioconductor
# BiocGenerics           * 0.48.1    2023-11-02 [1] Bioconductor
# BiocIO                   1.11.0    2023-06-02 [1] Bioconductor
# BiocManager            * 1.30.21.1 2023-07-18 [1] CRAN (R 4.3.0)
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
# edgeR                  * 3.43.7    2023-06-21 [1] Bioconductor
# ellipsis                 0.3.2     2021-04-29 [1] CRAN (R 4.3.0)
# enrichplot               1.21.1    2023-07-03 [1] Bioconductor
# ensembldb              * 2.26.0    2023-10-26 [1] Bioconductor
# fansi                    1.0.5     2023-10-08 [1] CRAN (R 4.3.1)
# farver                   2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
# fastmap                  1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
# fastmatch                1.1-3     2021-07-23 [1] CRAN (R 4.3.0)
# fgsea                    1.27.0    2023-05-20 [1] Bioconductor
# filelock                 1.0.2     2018-10-05 [1] CRAN (R 4.3.0)
# formatR                  1.14      2023-01-17 [1] CRAN (R 4.3.0)
# fs                       1.6.3     2023-07-20 [1] CRAN (R 4.3.0)
# futile.logger          * 1.4.3     2016-07-10 [1] CRAN (R 4.3.0)
# futile.options           1.0.1     2018-04-20 [1] CRAN (R 4.3.0)
# gargle                   1.5.2     2023-07-20 [1] CRAN (R 4.3.0)
# generics                 0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb           * 1.37.2    2023-06-21 [1] Bioconductor
# GenomeInfoDbData         1.2.10    2023-05-28 [1] Bioconductor
# GenomicAlignments        1.37.0    2023-07-07 [1] Bioconductor
# GenomicFeatures        * 1.53.1    2023-06-22 [1] Bioconductor
# GenomicRanges          * 1.54.1    2023-10-30 [1] Bioconductor
# ggforce                  0.4.1     2022-10-04 [1] CRAN (R 4.3.0)
# ggfun                    0.1.1     2023-06-24 [1] CRAN (R 4.3.0)
# ggplot2                * 3.4.4     2023-10-12 [1] CRAN (R 4.3.1)
# ggplotify                0.1.1     2023-06-27 [1] CRAN (R 4.3.0)
# ggraph                   2.1.0     2022-10-09 [1] CRAN (R 4.3.0)
# ggrepel                * 0.9.3     2023-02-03 [1] CRAN (R 4.3.0)
# ggtree                   3.9.0     2023-05-20 [1] Bioconductor
# glue                     1.6.2     2022-02-24 [1] CRAN (R 4.3.0)
# GO.db                    3.17.0    2023-05-28 [1] Bioconductor
# goalie                   0.7.7     2023-12-04 [1] Bioconductor
# googledrive              2.1.1     2023-06-11 [1] CRAN (R 4.3.0)
# GOSemSim                 2.27.2    2023-07-14 [1] Bioconductor
# graphlayouts             1.0.0     2023-05-01 [1] CRAN (R 4.3.0)
# gridExtra              * 2.3       2017-09-09 [1] CRAN (R 4.3.0)
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
# lambda.r                 1.2.4     2019-09-18 [1] CRAN (R 4.3.0)
# later                    1.3.1     2023-05-02 [1] CRAN (R 4.3.0)
# lattice                  0.21-8    2023-04-05 [1] CRAN (R 4.3.0)
# lazyeval                 0.2.2     2019-03-15 [1] CRAN (R 4.3.0)
# lifecycle                1.0.3     2022-10-07 [1] CRAN (R 4.3.0)
# limma                  * 3.57.6    2023-06-21 [1] Bioconductor
# locfit                   1.5-9.8   2023-06-11 [1] CRAN (R 4.3.0)
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
# pipette                  0.15.2    2023-12-15 [1] Bioconductor
# pkgconfig                2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# plyr                     1.8.8     2022-11-11 [1] CRAN (R 4.3.0)
# png                      0.1-8     2022-11-29 [1] CRAN (R 4.3.0)
# polyclip                 1.10-4    2022-10-20 [1] CRAN (R 4.3.0)
# prettyunits              1.1.1     2020-01-24 [1] CRAN (R 4.3.0)
# progress                 1.2.2     2019-05-16 [1] CRAN (R 4.3.0)
# promises                 1.2.0.1   2021-02-11 [1] CRAN (R 4.3.0)
# ProtGenerics             1.34.0    2023-10-26 [1] Bioconductor
# purrr                    1.0.1     2023-01-10 [1] CRAN (R 4.3.0)
# qvalue                   2.33.0    2023-05-11 [1] Bioconductor
# R.methodsS3            * 1.8.2     2022-06-13 [1] CRAN (R 4.3.0)
# R.oo                   * 1.25.0    2022-06-12 [1] CRAN (R 4.3.0)
# R.utils                * 2.12.2    2022-11-11 [1] CRAN (R 4.3.0)
# R6                       2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# rafalib                * 1.0.0     2015-08-09 [1] CRAN (R 4.3.0)
# ragg                     1.2.5     2023-01-12 [1] CRAN (R 4.3.0)
# rappdirs                 0.3.3     2021-01-31 [1] CRAN (R 4.3.0)
# RColorBrewer             1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                     1.0.11    2023-07-06 [1] CRAN (R 4.3.0)
# RCurl                    1.98-1.12 2023-03-27 [1] CRAN (R 4.3.0)
# reshape2                 1.4.4     2020-04-09 [1] CRAN (R 4.3.0)
# restfulr                 0.0.15    2022-06-16 [1] CRAN (R 4.3.0)
# rjson                    0.2.21    2022-01-09 [1] CRAN (R 4.3.0)
# rlang                  * 1.1.1     2023-04-28 [1] CRAN (R 4.3.0)
# rprojroot                2.0.3     2022-04-02 [1] CRAN (R 4.3.0)
# Rsamtools                2.17.0    2023-07-07 [1] Bioconductor
# RSQLite                  2.3.1     2023-04-03 [1] CRAN (R 4.3.0)
# rstudioapi               0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
# rtracklayer            * 1.61.0    2023-07-07 [1] Bioconductor
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
# syntactic                0.7.1     2023-10-27 [1] Bioconductor
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
# VennDiagram            * 1.7.3     2022-04-12 [1] CRAN (R 4.3.0)
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
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────


