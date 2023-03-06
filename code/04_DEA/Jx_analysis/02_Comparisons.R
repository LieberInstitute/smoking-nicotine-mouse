
## 1.2 Comparison of DE jnxs 

load(here("processed-data/04_DEA/Jx_analysis/de_jxns_nic.Rdata"))
load(here("processed-data/04_DEA/Jx_analysis/de_jxns_smo.Rdata"))
load(here("processed-data/07_Jxn_anno/novel_jxns_foundGenes.Rdata"))

load(here("processed-data/04_DEA/Tx_analysis/de_tx_nic.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/de_tx_smo.Rdata"))

load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata"))

load(here("processed-data/04_DEA/Exon_analysis/de_exons_nic.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/de_exons_smo.Rdata"))



### 1.2.1 Venn diagrams

## Function to create multiple Venn diagrams
venn_plot<-function(DE_lists, colors, name, titles){

  height=5
  width=14
  dist=0.1
  cat=0.75
  cex=0.8
  main_cex =1 
  main_pos=c(0.5, 0.09)
  
  plots<-list()
  pdf(file = paste("plots/04_DEA/02_Comparisons/Jx_analysis/Venn_", name, ".pdf", sep=""), height = height, width = width)
  for (i in 1:length(DE_lists)){
    v<-venn.diagram(DE_lists[[i]], fill=colors[[i]], alpha = rep(0.5, length(DE_lists[[i]])), 
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

## Define groups of genes associated to Novel DE jxns
novel_de_jxns_nic <- unique(de_jxns_nic[which(de_jxns_nic$Class=="Novel"), "newGeneID"])
Novel_DEjxns_genes_nic <- novel_de_jxns_nic[which(!is.na(novel_de_jxns_nic))]
novel_de_jxns_smo <- unique(de_jxns_smo[which(de_jxns_smo$Class=="Novel"), "newGeneID"])
Novel_DEjxns_genes_smo <- novel_de_jxns_smo[which(!is.na(novel_de_jxns_smo))]

## Define groups of the nearest genes of Novel DE jxns without associated gene
## (No NAs)
nearest_genes_nic <- unique(novel_jxns_foundGenes[["nearest_genes_nic"]]$geneIdVersion)
nearest_genes_smo <- unique(novel_jxns_foundGenes[["nearest_genes_smo"]]$geneIdVersion)

## Define groups of the following genes of Novel DE jxns without associated gene
## (No NAs)
following_genes_nic <- unique(novel_jxns_foundGenes[["following_genes_nic"]]$geneIdVersion)
following_genes_smo <- unique(novel_jxns_foundGenes[["following_genes_smo"]]$geneIdVersion)

## Define groups of the preceding genes of Novel DE jxns without associated gene
## (No NAs)
preceding_genes_nic <- unique(novel_jxns_foundGenes[["preceding_genes_nic"]]$geneIdVersion)
preceding_genes_smo <- unique(novel_jxns_foundGenes[["preceding_genes_smo"]]$geneIdVersion)

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
colors<-list(c("rosybrown2", "navajowhite2", "lightblue3", "darkolivegreen2"), c("plum", "lightgoldenrod2", "lightcyan3", "darkolivegreen3"), 
             c("pink1", "khaki2", "lightsteelblue3", "darkolivegreen4"))
venn_plot(DE_lists, colors, "DEG_VS_txs_VS_exons_VS_jxns", c("All", "Nicotine", "Smoking"))



##########################################################################################
## Compare DEG vs DE txs' genes vs DE exons' genes vs genes associated to Novel DE jxns
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
colors<-list(c("rosybrown2", "navajowhite2", "lightblue3", "darkolivegreen2"), c("plum", "lightgoldenrod2", "lightcyan3", "darkolivegreen3"), 
             c("pink1", "khaki2", "lightsteelblue3", "darkolivegreen4"))
venn_plot(DE_lists, colors, "DEG_VS_txs_VS_exons_VS_jxns_NovelAssigned", c("All", "Nicotine", "Smoking"))



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
colors<-list(c("rosybrown2", "navajowhite2", "lightblue3", "darkolivegreen2"), c("plum", "lightgoldenrod2", "lightcyan3", "darkolivegreen3"), 
             c("pink1", "khaki2", "lightsteelblue3", "darkolivegreen4"))
venn_plot(DE_lists, colors, "DEG_VS_txs_VS_exons_VS_jxns_NovelNearest", c("All", "Nicotine", "Smoking"))



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
colors<-list(c("rosybrown2", "navajowhite2", "lightblue3", "darkolivegreen2"), c("plum", "lightgoldenrod2", "lightcyan3", "darkolivegreen3"), 
             c("pink1", "khaki2", "lightsteelblue3", "darkolivegreen4"))
venn_plot(DE_lists, colors, "DEG_VS_txs_VS_exons_VS_jxns_NovelFollowing", c("All", "Nicotine", "Smoking"))



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
colors<-list(c("rosybrown2", "navajowhite2", "lightblue3", "darkolivegreen2"), c("plum", "lightgoldenrod2", "lightcyan3", "darkolivegreen3"), 
             c("pink1", "khaki2", "lightsteelblue3", "darkolivegreen4"))
venn_plot(DE_lists, colors, "DEG_VS_txs_VS_exons_VS_jxns_NovelPreceding", c("All", "Nicotine", "Smoking"))



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
colors<-list(c("rosybrown2", "navajowhite2", "lightblue3", "darkolivegreen2"), c("plum", "lightgoldenrod2", "lightcyan3", "darkolivegreen3"), 
             c("pink1", "khaki2", "lightsteelblue3", "darkolivegreen4"))
venn_plot(DE_lists, colors, "DEG_VS_txs_VS_exons_VS_jxns_AltStartEnd", c("All", "Nicotine", "Smoking"))



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
colors<-list(c("rosybrown2", "navajowhite2", "lightblue3", "darkolivegreen2"), c("plum", "lightgoldenrod2", "lightcyan3", "darkolivegreen3"), 
             c("pink1", "khaki2", "lightsteelblue3", "darkolivegreen4"))
venn_plot(DE_lists, colors, "DEG_VS_txs_VS_exons_VS_jxns_InGen", c("All", "Nicotine", "Smoking"))



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
colors<-list(c("rosybrown2", "navajowhite2", "lightblue3", "darkolivegreen2"), c("plum", "lightgoldenrod2", "lightcyan3", "darkolivegreen3"), 
             c("pink1", "khaki2", "lightsteelblue3", "darkolivegreen4"))
venn_plot(DE_lists, colors, "DEG_VS_txs_VS_exons_VS_jxns_ExonSkip", c("All", "Nicotine", "Smoking"))





### 1.2.1.1 Explore expression levels of DE features at only one level in nic and smo
## MA plots

## Obtain DE features at only one level
## DEG
DEG_only <- 



cols <- c("Up" = "#ffad73", "Down" = "#26b3ff", "ns" = "grey") 
sizes <- c("Up" = 2, "Down" = 2, "ns" = 1) 
alphas <- c("Up" = 1, "Down" = 1, "ns" = 0.5)
top_genes$mean_log_expr<-apply(vGene$E, 1, mean)
p1<-ggplot(data = top_genes, 
           aes(x = mean_log_expr,y = logFC,
               fill = DE,    
               size = DE,
               alpha = DE)) + 
  geom_point(shape = 21,    
             colour = "black") +
  scale_fill_manual(values = cols) + 
  scale_size_manual(values = sizes) + 
  scale_alpha_manual(values = alphas) +
  labs(x="Mean of normalized counts")







### 1.2.2 Final results: obtain DE features of DEG

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

# setting  value
# version  R version 4.2.2 (2022-10-31)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Mexico_City
# date     2023-03-02
# rstudio  2022.12.0+353 Elsbeth Geranium (desktop)
# pandoc   NA



