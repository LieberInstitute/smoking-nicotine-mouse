
## 1.2 Comparison of DE jnxs 

load(here("processed-data/04_DEA/Jx_analysis/de_jxns_nic.Rdata"))
load(here("processed-data/04_DEA/Jx_analysis/top_jxns_nic.Rdata"))
load(here("processed-data/04_DEA/Jx_analysis/de_jxns_smo.Rdata"))
load(here("processed-data/04_DEA/Jx_analysis/top_jxns_smo.Rdata"))
load(here("processed-data/07_Jxn_anno/novel_jxns_foundGenes.Rdata"))

load(here("processed-data/04_DEA/Tx_analysis/de_tx_nic.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/top_tx_nic.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/de_tx_smo.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/top_tx_smo.Rdata"))

load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/results_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/results_pups_smoking_fitted.Rdata"))

load(here("processed-data/04_DEA/Exon_analysis/de_exons_nic.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/top_exons_nic.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/de_exons_smo.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/top_exons_smo.Rdata"))


### 1.2.1 Venn diagrams

## Function to create multiple Venn diagrams
venn_plot<-function(DE_lists, name, titles){
  
  if(name!='DEG_VS_txs_VS_exons_VS_jxns'){
    cat=1.05
    dist=0.2
  }
  else{
    cat=1.2
    dist=0.2
  }

  height=6 
  width=17.5
  cex=1.2
  main_cex =1.2 
  main_pos=c(0.5, 0.25)
  
  colors<-c("pink", "lightgoldenrod3", "turquoise", "lightsteelblue3")
  
  plots<-list()
  pdf(file = paste("plots/04_DEA/02_Comparisons/Jx_analysis/Venn_", name, ".pdf", sep=""), height = height, width = width)
  for (i in 1:length(DE_lists)){
    v<-venn.diagram(DE_lists[[i]], fill=colors, alpha = rep(0.5, length(DE_lists[[i]])), margin= 0.4,
                    lwd =0, cat.cex=cat, cex=cex, cat.dist=rep(dist, length(DE_lists[[i]])), filename=NULL, main = titles[i], 
                    main.cex = main_cex, main.pos = main_pos, disable.logging=TRUE)
    plots[[i]]<-v
  }
  
  gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], ncol=3)
  dev.off()

}


## Define groups of DEG
DEG_nic<-unique(de_genes_pups_nicotine_fitted$gencodeID)
DEG_smo<-unique(de_genes_pups_smoking_fitted$gencodeID)

## Define groups of DE txs' genes
DEtxs_genes_nic<-unique(de_tx_nic$ensembl_id)
DEtxs_genes_smo<-unique(de_tx_smo$ensembl_id)

## Define groups of DE exons' genes
DEexons_genes_nic<-unique(de_exons_nic$gencodeID)
DEexons_genes_smo<-unique(de_exons_smo$gencodeID)

## Define groups of all DE jxns' genes (from any jxn class)
DEjxns_genes_nic <- unique(de_jxns_nic$newGeneID)[which(!is.na(unique(de_jxns_nic$newGeneID)))]
DEjxns_genes_smo <- unique(de_jxns_smo$newGeneID)[which(!is.na(unique(de_jxns_smo$newGeneID)))]

## Define groups of genes associated with Novel DE jxns
novel_de_jxns_nic <- unique(de_jxns_nic[which(de_jxns_nic$Class=="Novel"), "newGeneID"])
Novel_DEjxns_genes_nic <- novel_de_jxns_nic[which(!is.na(novel_de_jxns_nic))]
novel_de_jxns_smo <- unique(de_jxns_smo[which(de_jxns_smo$Class=="Novel"), "newGeneID"])
Novel_DEjxns_genes_smo <- novel_de_jxns_smo[which(!is.na(novel_de_jxns_smo))]

## Define groups of the nearest genes of Novel DE jxns without associated gene
## (No NAs)
nearest_genes_nic <- unique(novel_jxns_foundGenes[["nearest_genes_nic"]])
nearest_genes_smo <- unique(novel_jxns_foundGenes[["nearest_genes_smo"]])

## Define groups of the following genes of Novel DE jxns without associated gene
## (No NAs)
following_genes_nic <- unique(novel_jxns_foundGenes[["following_genes_nic"]])
following_genes_smo <- unique(novel_jxns_foundGenes[["following_genes_smo"]])

## Define groups of the preceding genes of Novel DE jxns without associated gene
## (No NAs)
preceding_genes_nic <- unique(novel_jxns_foundGenes[["preceding_genes_nic"]])
preceding_genes_smo <- unique(novel_jxns_foundGenes[["preceding_genes_smo"]])

## Define groups of AltStartEnd DE jxns' genes
## (No NAs)
alt_de_jxns_nic <- unique(de_jxns_nic[which(de_jxns_nic$Class=="AltStartEnd"), "newGeneID"])
alt_de_jxns_smo <- unique(de_jxns_smo[which(de_jxns_smo$Class=="AltStartEnd"), "newGeneID"])

## Define groups of ExonSkip DE jxns' genes
## (No NAs)
exonSkip_de_jxns_nic <- unique(de_jxns_nic[which(de_jxns_nic$Class=="ExonSkip"), "newGeneID"])
exonSkip_de_jxns_smo <- unique(de_jxns_smo[which(de_jxns_smo$Class=="ExonSkip"), "newGeneID"])

## Define groups of InGen DE jxns' genes
## (No NAs)
inGen_de_jxns_nic <- unique(de_jxns_nic[which(de_jxns_nic$Class=="InGen"), "newGeneID"])
inGen_de_jxns_smo <- unique(de_jxns_smo[which(de_jxns_smo$Class=="InGen"), "newGeneID"])



## Venn diagrams

#################################################################################################
## Compare DEG vs DE txs' genes vs DE exons' genes vs All DE jxns' genes (from all jxn classes)
#################################################################################################

## All genes
DEG_vs_Txs_vs_Exons_vs_Jxns_all <-list(
  "DEG"= union(DEG_nic, DEG_smo),
  "DE txs' genes"=union(DEtxs_genes_nic, DEtxs_genes_smo),
  "DE exons' genes"=union(DEexons_genes_nic, DEexons_genes_smo),
  "DE jxns' genes"=union(DEjxns_genes_nic, DEjxns_genes_smo)                       
)

## Nic genes
DEG_vs_Txs_vs_Exons_vs_Jxns_nic <-list(
  "DEG"= DEG_nic,
  "DE txs' genes"=DEtxs_genes_nic,
  "DE exons' genes"=DEexons_genes_nic,
  "DE jxns' genes"=DEjxns_genes_nic
)

## Smo genes
DEG_vs_Txs_vs_Exons_vs_Jxns_smo <-list(
  "DEG"= DEG_smo,
  "DE txs' genes"=DEtxs_genes_smo,
  "DE exons' genes"=DEexons_genes_smo,
  "DE jxns' genes"=DEjxns_genes_smo
)

DE_lists<-list(DEG_vs_Txs_vs_Exons_vs_Jxns_all, DEG_vs_Txs_vs_Exons_vs_Jxns_nic, DEG_vs_Txs_vs_Exons_vs_Jxns_smo)
venn_plot(DE_lists, "DEG_VS_txs_VS_exons_VS_jxns", c("All", "Nicotine", "Smoking"))



##########################################################################################
## Compare DEG vs DE txs' genes vs DE exons' genes vs genes associated with Novel DE jxns
##########################################################################################

## All genes
DEG_vs_Txs_vs_Exons_vs_NovelAssignedJxns_all <-list(
  "DEG"= union(DEG_nic, DEG_smo),
  "DE txs' genes"=union(DEtxs_genes_nic, DEtxs_genes_smo),
  "DE exons' genes"=union(DEexons_genes_nic, DEexons_genes_smo),
  "Genes assigned to Novel DE jxns"=union(Novel_DEjxns_genes_nic, Novel_DEjxns_genes_smo)                       
)

## Nic genes
DEG_vs_Txs_vs_Exons_vs_NovelAssignedJxns_nic <-list(
  "DEG"= DEG_nic,
  "DE txs' genes"=DEtxs_genes_nic,
  "DE exons' genes"=DEexons_genes_nic,
  "Genes assigned to Novel DE jxns"=Novel_DEjxns_genes_nic
)

## Smo genes
DEG_vs_Txs_vs_Exons_vs_NovelAssignedJxns_smo <-list(
  "DEG"= DEG_smo,
  "DE txs' genes"=DEtxs_genes_smo,
  "DE exons' genes"=DEexons_genes_smo,
  "Genes assigned to Novel DE jxns"=Novel_DEjxns_genes_smo
)

DE_lists<-list(DEG_vs_Txs_vs_Exons_vs_NovelAssignedJxns_all, DEG_vs_Txs_vs_Exons_vs_NovelAssignedJxns_nic, DEG_vs_Txs_vs_Exons_vs_NovelAssignedJxns_smo)
venn_plot(DE_lists, "DEG_VS_txs_VS_exons_VS_jxns_NovelAssigned", c("All", "Nicotine", "Smoking"))



#########################################################################################################################
## Compare DEG vs DE txs' genes vs DE exons' genes vs nearest genes to Novel DE jxns (that were not assigned to a gene)
#########################################################################################################################

## All genes
DEG_vs_Txs_vs_Exons_vs_NovelNearestJxns_all <-list(
  "DEG"= union(DEG_nic, DEG_smo),
  "DE txs' genes"=union(DEtxs_genes_nic, DEtxs_genes_smo),
  "DE exons' genes"=union(DEexons_genes_nic, DEexons_genes_smo),
  "Nearest genes to Novel DE jxns"=union(nearest_genes_nic, nearest_genes_smo)                       
)

## Nic genes
DEG_vs_Txs_vs_Exons_vs_NovelNearestJxns_nic <-list(
  "DEG"= DEG_nic,
  "DE txs' genes"=DEtxs_genes_nic,
  "DE exons' genes"=DEexons_genes_nic,
  "Nearest genes to Novel DE jxns"= nearest_genes_nic
)

## Smo genes
DEG_vs_Txs_vs_Exons_vs_NovelNearestJxns_smo <-list(
  "DEG"= DEG_smo,
  "DE txs' genes"=DEtxs_genes_smo,
  "DE exons' genes"=DEexons_genes_smo,
  "Nearest genes to Novel DE jxns"= nearest_genes_smo
)

DE_lists<-list(DEG_vs_Txs_vs_Exons_vs_NovelNearestJxns_all, DEG_vs_Txs_vs_Exons_vs_NovelNearestJxns_nic, DEG_vs_Txs_vs_Exons_vs_NovelNearestJxns_smo)
venn_plot(DE_lists, "DEG_VS_txs_VS_exons_VS_jxns_NovelNearest", c("All", "Nicotine", "Smoking"))



###########################################################################################################################
## Compare DEG vs DE txs' genes vs DE exons' genes vs following genes to Novel DE jxns (that were not assigned to a gene)
###########################################################################################################################

## All genes
DEG_vs_Txs_vs_Exons_vs_NovelFollowingJxns_all <-list(
  "DEG"= union(DEG_nic, DEG_smo),
  "DE txs' genes"=union(DEtxs_genes_nic, DEtxs_genes_smo),
  "DE exons' genes"=union(DEexons_genes_nic, DEexons_genes_smo),
  "Following genes to Novel DE jxns"=union(following_genes_nic, following_genes_smo)                       
)

## Nic genes
DEG_vs_Txs_vs_Exons_vs_NovelFollowingJxns_nic <-list(
  "DEG"= DEG_nic,
  "DE txs' genes"=DEtxs_genes_nic,
  "DE exons' genes"=DEexons_genes_nic,
  "Following genes to Novel DE jxns"= following_genes_nic
)

## Smo genes
DEG_vs_Txs_vs_Exons_vs_NovelFollowingJxns_smo <-list(
  "DEG"= DEG_smo,
  "DE txs' genes"=DEtxs_genes_smo,
  "DE exons' genes"=DEexons_genes_smo,
  "Following genes to Novel DE jxns"= following_genes_smo
)

DE_lists<-list(DEG_vs_Txs_vs_Exons_vs_NovelFollowingJxns_all, DEG_vs_Txs_vs_Exons_vs_NovelFollowingJxns_nic, DEG_vs_Txs_vs_Exons_vs_NovelFollowingJxns_smo)
venn_plot(DE_lists, "DEG_VS_txs_VS_exons_VS_jxns_NovelFollowing", c("All", "Nicotine", "Smoking"))



###########################################################################################################################
## Compare DEG vs DE txs' genes vs DE exons' genes vs preceding genes to Novel DE jxns (that were not assigned to a gene)
###########################################################################################################################

## All genes
DEG_vs_Txs_vs_Exons_vs_NovelPrecedingJxns_all <-list(
  "DEG"= union(DEG_nic, DEG_smo),
  "DE txs' genes"=union(DEtxs_genes_nic, DEtxs_genes_smo),
  "DE exons' genes"=union(DEexons_genes_nic, DEexons_genes_smo),
  "Preceding genes to Novel DE jxns"=union(preceding_genes_nic, preceding_genes_smo)                       
)

## Nic genes
DEG_vs_Txs_vs_Exons_vs_NovelPrecedingJxns_nic <-list(
  "DEG"= DEG_nic,
  "DE txs' genes"=DEtxs_genes_nic,
  "DE exons' genes"=DEexons_genes_nic,
  "Preceding genes to Novel DE jxns"= preceding_genes_nic
)

## Smo genes
DEG_vs_Txs_vs_Exons_vs_NovelPrecedingJxns_smo <-list(
  "DEG"= DEG_smo,
  "DE txs' genes"=DEtxs_genes_smo,
  "DE exons' genes"=DEexons_genes_smo,
  "Preceding genes to Novel DE jxns"= preceding_genes_smo
)

DE_lists<-list(DEG_vs_Txs_vs_Exons_vs_NovelPrecedingJxns_all, DEG_vs_Txs_vs_Exons_vs_NovelPrecedingJxns_nic, DEG_vs_Txs_vs_Exons_vs_NovelPrecedingJxns_smo)
venn_plot(DE_lists, "DEG_VS_txs_VS_exons_VS_jxns_NovelPreceding", c("All", "Nicotine", "Smoking"))



###################################################################################
## Compare DEG vs DE txs' genes vs DE exons' genes vs AltStartEnd DE jxns' genes 
###################################################################################

## All genes
DEG_vs_Txs_vs_Exons_vs_AltJxns_all <-list(
  "DEG"= union(DEG_nic, DEG_smo),
  "DE txs' genes"=union(DEtxs_genes_nic, DEtxs_genes_smo),
  "DE exons' genes"=union(DEexons_genes_nic, DEexons_genes_smo),
  "AltStartEnd DE jxns' genes"=union(alt_de_jxns_nic, alt_de_jxns_smo)                       
)

## Nic genes
DEG_vs_Txs_vs_Exons_vs_AltJxns_nic <-list(
  "DEG"= DEG_nic,
  "DE txs' genes"=DEtxs_genes_nic,
  "DE exons' genes"=DEexons_genes_nic,
  "AltStartEnd DE jxns' genes"=alt_de_jxns_nic
)

## Smo genes
DEG_vs_Txs_vs_Exons_vs_AltJxns_smo <-list(
  "DEG"= DEG_smo,
  "DE txs' genes"=DEtxs_genes_smo,
  "DE exons' genes"=DEexons_genes_smo,
  "AltStartEnd DE jxns' genes"=alt_de_jxns_smo
)

DE_lists<-list(DEG_vs_Txs_vs_Exons_vs_AltJxns_all, DEG_vs_Txs_vs_Exons_vs_AltJxns_nic, DEG_vs_Txs_vs_Exons_vs_AltJxns_smo)
venn_plot(DE_lists, "DEG_VS_txs_VS_exons_VS_jxns_AltStartEnd", c("All", "Nicotine", "Smoking"))



###################################################################################
## Compare DEG vs DE txs' genes vs DE exons' genes vs InGen DE jxns' genes 
###################################################################################

## All genes
DEG_vs_Txs_vs_Exons_vs_InGenJxns_all <-list(
  "DEG"= union(DEG_nic, DEG_smo),
  "DE txs' genes"=union(DEtxs_genes_nic, DEtxs_genes_smo),
  "DE exons' genes"=union(DEexons_genes_nic, DEexons_genes_smo),
  "InGen DE jxns' genes"=union(inGen_de_jxns_nic, inGen_de_jxns_smo)                       
)

## Nic genes
DEG_vs_Txs_vs_Exons_vs_InGenJxns_nic <-list(
  "DEG"= DEG_nic,
  "DE txs' genes"=DEtxs_genes_nic,
  "DE exons' genes"=DEexons_genes_nic,
  "InGen DE jxns' genes"=inGen_de_jxns_nic
)

## Smo genes
DEG_vs_Txs_vs_Exons_vs_InGenJxns_smo <-list(
  "DEG"= DEG_smo,
  "DE txs' genes"=DEtxs_genes_smo,
  "DE exons' genes"=DEexons_genes_smo,
  "InGen DE jxns' genes"=inGen_de_jxns_smo
)

DE_lists<-list(DEG_vs_Txs_vs_Exons_vs_InGenJxns_all, DEG_vs_Txs_vs_Exons_vs_InGenJxns_nic, DEG_vs_Txs_vs_Exons_vs_InGenJxns_smo)
venn_plot(DE_lists, "DEG_VS_txs_VS_exons_VS_jxns_InGen", c("All", "Nicotine", "Smoking"))



###################################################################################
## Compare DEG vs DE txs' genes vs DE exons' genes vs ExonSkip DE jxns' genes 
###################################################################################

## All genes
DEG_vs_Txs_vs_Exons_vs_ExonSkipJxns_all <-list(
  "DEG"= union(DEG_nic, DEG_smo),
  "DE txs' genes"=union(DEtxs_genes_nic, DEtxs_genes_smo),
  "DE exons' genes"=union(DEexons_genes_nic, DEexons_genes_smo),
  "ExonSkip DE jxns' genes"=union(exonSkip_de_jxns_nic, exonSkip_de_jxns_smo)                       
)

## Nic genes
DEG_vs_Txs_vs_Exons_vs_ExonSkipJxns_nic <-list(
  "DEG"= DEG_nic,
  "DE txs' genes"=DEtxs_genes_nic,
  "DE exons' genes"=DEexons_genes_nic,
  "ExonSkip DE jxns' genes"=exonSkip_de_jxns_nic
)

## Smo genes
DEG_vs_Txs_vs_Exons_vs_ExonSkipJxns_smo <-list(
  "DEG"= DEG_smo,
  "DE txs' genes"=DEtxs_genes_smo,
  "DE exons' genes"=DEexons_genes_smo,
  "ExonSkip DE jxns' genes"=exonSkip_de_jxns_smo
)

DE_lists<-list(DEG_vs_Txs_vs_Exons_vs_ExonSkipJxns_all, DEG_vs_Txs_vs_Exons_vs_ExonSkipJxns_nic, DEG_vs_Txs_vs_Exons_vs_ExonSkipJxns_smo)
venn_plot(DE_lists, "DEG_VS_txs_VS_exons_VS_jxns_ExonSkip", c("All", "Nicotine", "Smoking"))





### 1.2.2 Explore expression levels of genes with DE features in nic and smo

## Create MA plots
MAplot <- function(top_genes, results, expt){
  
  DEG <- eval(parse_expr(paste("DEG", expt, sep="_")))
  DEtxs_genes <- eval(parse_expr(paste("DEtxs_genes", expt, sep="_")))
  DEjxns_genes <- eval(parse_expr(paste("DEjxns_genes", expt, sep="_")))
  DEexons_genes <- eval(parse_expr(paste("DEexons_genes", expt, sep="_")))
  de_genes <- results[[2]]
  
  ## Colors for plot
  cols <- c("DEG with DE txs, exons and jxns" = "orangered3", "DEG without DE txs, jxns and exons" = "cornflowerblue", 
            "Non-DE Gene with DE txs only"= "darkturquoise", "Non-DE Gene with DE exons only"= "chartreuse3", 
            "Non-DE Gene with DE jxns only"="deeppink1", "Rest of DEG"="wheat", "ns" = "grey") 
  sizes <- c("DEG with DE txs, exons and jxns" = 2, "DEG without DE txs, jxns and exons" = 2, 
             "Non-DE Gene with DE txs only"= 2, "Non-DE Gene with DE exons only"= 2, 
             "Non-DE Gene with DE jxns only"=2, "Rest of DEG"=1, "ns" = 0.5) 
  alphas <- c("DEG with DE txs, exons and jxns" = 1, "DEG without DE txs, jxns and exons" = 1, 
              "Non-DE Gene with DE txs only"= 1, "Non-DE Gene with DE exons only"= 1, 
              "Non-DE Gene with DE jxns only"=1, "Rest of DEG"=0.9, "ns" = 0.3)
  
  ## Obtain genes with DE features 
  
  ###### DEG without DE txs, jxns and exons ######
  DEG_only <- DEG[which(! (DEG %in% DEtxs_genes | 
                           DEG %in% DEexons_genes | 
                           DEG %in% DEjxns_genes ))]
  
  ###### DE txs without DE gene, jxns and exons ######
  DEtx_only <- DEtxs_genes[which(! (DEtxs_genes %in% DEG | 
                                    DEtxs_genes %in% DEexons_genes | 
                                    DEtxs_genes %in% DEjxns_genes ))]
  
  ###### DE exons without DE gene, txs and jxns ######
  DEexon_only <- DEexons_genes[which(! (DEexons_genes %in% DEG | 
                                        DEexons_genes %in% DEtxs_genes | 
                                        DEexons_genes %in% DEjxns_genes ))]
  
  ###### DE jxns without DE gene, tx and exon ######
  DEjxn_only <- DEjxns_genes[which(! (DEjxns_genes %in% DEG | 
                                      DEjxns_genes %in% DEtxs_genes | 
                                      DEjxns_genes %in% DEexons_genes))]
  
  ###### DEG with DE txs, exons and jxns ######
  DE_all_levels <- intersect(DEG, intersect(DEtxs_genes , intersect(DEjxns_genes, DEexons_genes)))
  
  ## Find those genes in the whole gene dataset
  DEfeatures_onlyOnelevel <- vector()
  for (gene in top_genes$gencodeID){
    if (gene %in% DEG_only){
      DEfeatures_onlyOnelevel <- append(DEfeatures_onlyOnelevel, "DEG without DE txs, jxns and exons")
    }
    else if (gene %in% DEtx_only){
      DEfeatures_onlyOnelevel <- append(DEfeatures_onlyOnelevel, "Non-DE Gene with DE txs only")
    }
    else if (gene %in% DEjxn_only){
      DEfeatures_onlyOnelevel <- append(DEfeatures_onlyOnelevel, "Non-DE Gene with DE jxns only")
    }
    else if (gene %in% DEexon_only){
      DEfeatures_onlyOnelevel <- append(DEfeatures_onlyOnelevel, "Non-DE Gene with DE exons only")
    }
    else if (gene %in% DE_all_levels){
      DEfeatures_onlyOnelevel <- append(DEfeatures_onlyOnelevel, "DEG with DE txs, exons and jxns")
    }
    else if(gene %in% de_genes$gencodeID){
      DEfeatures_onlyOnelevel <- append(DEfeatures_onlyOnelevel, "Rest of DEG")
    }
    else {
      DEfeatures_onlyOnelevel <- append(DEfeatures_onlyOnelevel, "ns")
    }
  }
  
  top_genes$DEfeatures_onlyOnelevel <- DEfeatures_onlyOnelevel
  top_genes$DEfeatures_onlyOnelevel <- factor(top_genes$DEfeatures_onlyOnelevel, levels=names(cols))
  
  ## Plots
  
  ## Mean expression values of cpm
  vGene <- results[[1]][[2]]
  top_genes$mean_log_expr<-apply(vGene$E, 1, mean)

  
  ## MA plots
  p <-ggplot(data = top_genes, 
            aes(x = mean_log_expr,y = logFC,
                 fill = DEfeatures_onlyOnelevel,    
                 size = DEfeatures_onlyOnelevel,
                 alpha = DEfeatures_onlyOnelevel)) + 
    geom_point(data=subset(top_genes, DEfeatures_onlyOnelevel!='ns'), shape = 21, colour="black") +
    geom_point(data=subset(top_genes, DEfeatures_onlyOnelevel=='ns'), shape = 21, colour="gray49") +
    theme_bw() +
    scale_fill_manual(values = cols) + 
    scale_size_manual(values = sizes) + 
    scale_alpha_manual(values = alphas) +
    labs(x="Mean of lognorm counts", fill="Genes and their DE features", size="Genes and their DE features", 
         alpha="Genes and their DE features", y='Log2FC (Exposed vs Ctrl)') +
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          axis.title = element_text(size = (11)),
          axis.text = element_text(size = (10)),
          legend.text = element_text(size=10),
          legend.title = element_text(size=11))

  ggsave(paste("plots/04_DEA/02_Comparisons/Jx_analysis/MAplots_DEfeatures_", expt, ".pdf", sep=""),
        width = 25, height = 15, units = "cm")
 
 
   ## Create facets of MA plots
   p <- p +  facet_wrap(~DEfeatures_onlyOnelevel, scales = "fixed", ncol=4) +
     theme(legend.position = 'none',
           strip.text = element_text(size = 9, face = "bold"),
           strip.background = element_blank())
   
   ggsave(paste("plots/04_DEA/02_Comparisons/Jx_analysis/facets_MAplots_DEfeatures_", expt, ".pdf", sep=""),
          width = 26, height = 11, units = "cm")
}

## MA plots for nicotine genes
MAplot(top_genes_pups_nicotine_fitted, results_pups_nicotine_fitted, "nic")

## MA plots for smoking genes
MAplot(top_genes_pups_smoking_fitted, results_pups_smoking_fitted, "smo")





### 1.2.3 Final results: obtain DE features of DEG

## The next data frame contains the following info for all the DEG:
##        - Gene: ID of the DEG
##        - DEG_in_nic: 1 if the gene is DE in nicotine, 0 if not
##        - DEG_in_smo: 1 if the gene is DE in smoking, 0 if not
##        - with_DE_Txs_Nic: 1 if the DEG has DE txs in nicotine, 0 if not
##        - with_DE_Exons_Nic: 1 if the DEG has DE exons in nicotine, 0 if not
##        - with_DE_Jxns_Nic: 1 if the DEG has DE jxns in nicotine, 0 if not
##        - with_DE_Txs_Smo: 1 if the DEG has DE txs in smoking, 0 if not
##        - with_DE_Exons_Smo: 1 if the DEG has DE exons in smoking, 0 if not
##        - with_DE_Jxns_Smo: 1 if the DEG has DE jxns in smoking, 0 if not

## All DEG 
DEG_all <- union(DEG_nic, DEG_smo)
final_DE_results <- data.frame(
  Gene=union(DEG_nic, DEG_smo),
  DEG_in_nic= unlist(sapply(DEG_all, function(x){if (x %in% DEG_nic) {"Yes"} else {"No"}})),
  DEG_in_smo= unlist(sapply(DEG_all, function(x){if (x %in% DEG_smo) {"Yes"} else {"No"}})),
  with_DE_Txs_Nic= unlist(sapply(DEG_all, function(x){if (x %in% DEtxs_genes_nic) {"Yes"} else {"No"}})),
  with_DE_Exons_Nic= unlist(sapply(DEG_all, function(x){if (x %in% DEexons_genes_nic) {"Yes"} else {"No"}})),
  with_DE_Jxns_Nic= unlist(sapply(DEG_all, function(x){if (x %in% DEjxns_genes_nic) {"Yes"} else {"No"}})),
  with_DE_Txs_Smo= unlist(sapply(DEG_all, function(x){if (x %in% DEtxs_genes_smo) {"Yes"} else {"No"}})),
  with_DE_Exons_Smo= unlist(sapply(DEG_all, function(x){if (x %in% DEexons_genes_smo) {"Yes"} else {"No"}})),
  with_DE_Jxns_Smo= unlist(sapply(DEG_all, function(x){if (x %in% DEjxns_genes_smo) {"Yes"} else {"No"}}))
)

save(final_DE_results, file="processed-data/04_DEA/Jx_analysis/final_DE_results.Rdata")







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
# date     2023-12-29
# rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
# pandoc   NA
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# AcidBase               0.7.3     2023-12-15 [1] Bioconductor
# AcidCLI                0.3.0     2023-10-03 [1] Bioconductor
# AcidGenerics           0.7.6     2023-12-15 [1] Bioconductor
# AcidGenomes          * 0.7.2     2023-12-06 [1] Bioconductor
# AcidPlyr               0.5.3     2023-12-13 [1] local
# AnnotationDbi          1.63.2    2023-07-03 [1] Bioconductor
# Biobase              * 2.61.0    2023-06-02 [1] Bioconductor
# BiocFileCache          2.9.1     2023-07-14 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-02 [1] Bioconductor
# BiocIO                 1.11.0    2023-06-02 [1] Bioconductor
# BiocParallel           1.35.3    2023-07-07 [1] Bioconductor
# biomaRt                2.57.1    2023-06-14 [1] Bioconductor
# biomartr             * 1.0.7     2023-12-02 [1] CRAN (R 4.3.1)
# Biostrings             2.69.2    2023-07-05 [1] Bioconductor
# bit                    4.0.5     2022-11-15 [1] CRAN (R 4.3.0)
# bit64                  4.0.5     2020-08-30 [1] CRAN (R 4.3.0)
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
# blob                   1.2.4     2023-03-17 [1] CRAN (R 4.3.0)
# cachem                 1.0.8     2023-05-01 [1] CRAN (R 4.3.0)
# cli                    3.6.1     2023-03-23 [1] CRAN (R 4.3.0)
# codetools              0.2-19    2023-02-01 [1] CRAN (R 4.3.0)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# cowplot              * 1.1.1     2020-12-30 [1] CRAN (R 4.3.0)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# curl                   5.0.1     2023-06-07 [1] CRAN (R 4.3.0)
# data.table             1.14.8    2023-02-17 [1] CRAN (R 4.3.0)
# DBI                    1.1.3     2022-06-18 [1] CRAN (R 4.3.0)
# dbplyr                 2.3.3     2023-07-07 [1] CRAN (R 4.3.0)
# DelayedArray           0.26.6    2023-07-02 [1] Bioconductor
# digest                 0.6.33    2023-07-07 [1] CRAN (R 4.3.0)
# dplyr                  1.1.2     2023-04-20 [1] CRAN (R 4.3.0)
# edgeR                * 3.43.7    2023-06-21 [1] Bioconductor
# fansi                  1.0.5     2023-10-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
# fastmap                1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
# filelock               1.0.2     2018-10-05 [1] CRAN (R 4.3.0)
# formatR                1.14      2023-01-17 [1] CRAN (R 4.3.0)
# fs                     1.6.3     2023-07-20 [1] CRAN (R 4.3.0)
# futile.logger        * 1.4.3     2016-07-10 [1] CRAN (R 4.3.0)
# futile.options         1.0.1     2018-04-20 [1] CRAN (R 4.3.0)
# gargle                 1.5.2     2023-07-20 [1] CRAN (R 4.3.0)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.37.2    2023-06-21 [1] Bioconductor
# GenomeInfoDbData       1.2.10    2023-05-28 [1] Bioconductor
# GenomicAlignments      1.37.0    2023-07-07 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-30 [1] Bioconductor
# ggplot2              * 3.4.4     2023-10-12 [1] CRAN (R 4.3.1)
# ggrepel              * 0.9.3     2023-02-03 [1] CRAN (R 4.3.0)
# glue                   1.6.2     2022-02-24 [1] CRAN (R 4.3.0)
# goalie                 0.7.7     2023-12-04 [1] Bioconductor
# googledrive            2.1.1     2023-06-11 [1] CRAN (R 4.3.0)
# gridExtra            * 2.3       2017-09-09 [1] CRAN (R 4.3.0)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.0)
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
# hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
# httr                   1.4.6     2023-05-08 [1] CRAN (R 4.3.0)
# IRanges              * 2.36.0    2023-10-26 [1] Bioconductor
# jaffelab             * 0.99.32   2023-05-28 [1] Github (LieberInstitute/jaffelab@21e6574)
# jsonlite               1.8.8     2023-12-04 [1] CRAN (R 4.3.1)
# KEGGREST               1.41.0    2023-07-07 [1] Bioconductor
# labeling               0.4.3     2023-08-29 [1] CRAN (R 4.3.0)
# lambda.r               1.2.4     2019-09-18 [1] CRAN (R 4.3.0)
# lattice                0.21-8    2023-04-05 [1] CRAN (R 4.3.0)
# lifecycle              1.0.3     2022-10-07 [1] CRAN (R 4.3.0)
# limma                * 3.57.6    2023-06-21 [1] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [1] CRAN (R 4.3.0)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# MASS                   7.3-60    2023-05-04 [1] CRAN (R 4.3.0)
# Matrix                 1.6-4     2023-11-30 [1] CRAN (R 4.3.1)
# MatrixGenerics       * 1.13.0    2023-05-20 [1] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.3.0)
# memoise                2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# nlme                   3.1-162   2023-01-31 [1] CRAN (R 4.3.0)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pipette                0.15.2    2023-12-15 [1] Bioconductor
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# png                    0.1-8     2022-11-29 [1] CRAN (R 4.3.0)
# prettyunits            1.1.1     2020-01-24 [1] CRAN (R 4.3.0)
# progress               1.2.2     2019-05-16 [1] CRAN (R 4.3.0)
# purrr                  1.0.1     2023-01-10 [1] CRAN (R 4.3.0)
# R.methodsS3          * 1.8.2     2022-06-13 [1] CRAN (R 4.3.0)
# R.oo                 * 1.25.0    2022-06-12 [1] CRAN (R 4.3.0)
# R.utils              * 2.12.2    2022-11-11 [1] CRAN (R 4.3.0)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.3.0)
# ragg                   1.2.5     2023-01-12 [1] CRAN (R 4.3.0)
# rappdirs               0.3.3     2021-01-31 [1] CRAN (R 4.3.0)
# RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.11    2023-07-06 [1] CRAN (R 4.3.0)
# RCurl                  1.98-1.12 2023-03-27 [1] CRAN (R 4.3.0)
# restfulr               0.0.15    2022-06-16 [1] CRAN (R 4.3.0)
# rjson                  0.2.21    2022-01-09 [1] CRAN (R 4.3.0)
# rlang                * 1.1.1     2023-04-28 [1] CRAN (R 4.3.0)
# rprojroot              2.0.3     2022-04-02 [1] CRAN (R 4.3.0)
# Rsamtools              2.17.0    2023-07-07 [1] Bioconductor
# RSQLite                2.3.1     2023-04-03 [1] CRAN (R 4.3.0)
# rstudioapi             0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
# rtracklayer          * 1.61.0    2023-07-07 [1] Bioconductor
# S4Arrays               1.1.4     2023-06-02 [1] Bioconductor
# S4Vectors            * 0.40.2    2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
# scales                 1.2.1     2022-08-20 [1] CRAN (R 4.3.0)
# segmented              1.6-4     2023-04-13 [1] CRAN (R 4.3.0)
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
# stringi                1.7.12    2023-01-11 [1] CRAN (R 4.3.0)
# stringr                1.5.0     2022-12-02 [1] CRAN (R 4.3.0)
# SummarizedExperiment * 1.30.2    2023-06-06 [1] Bioconductor
# syntactic              0.7.1     2023-10-27 [1] Bioconductor
# systemfonts            1.0.4     2022-02-11 [1] CRAN (R 4.3.0)
# textshaping            0.3.6     2021-10-13 [1] CRAN (R 4.3.0)
# tibble                 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.4     2023-10-12 [1] CRAN (R 4.3.1)
# VennDiagram          * 1.7.3     2022-04-12 [1] CRAN (R 4.3.0)
# withr                  2.5.2     2023-10-30 [1] CRAN (R 4.3.1)
# XML                    3.99-0.14 2023-03-19 [1] CRAN (R 4.3.0)
# xml2                   1.3.5     2023-07-06 [1] CRAN (R 4.3.0)
# XVector                0.41.1    2023-06-02 [1] Bioconductor
# yaml                   2.3.8     2023-12-11 [1] CRAN (R 4.3.1)
# zlibbioc               1.47.0    2023-05-20 [1] Bioconductor
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────



