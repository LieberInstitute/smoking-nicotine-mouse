
## 1.2 Comparison of DE jnxs 

load(here("processed-data/04_DEA/Jx_analysis/de_jxns_nic.Rdata"))
load(here("processed-data/04_DEA/Jx_analysis/de_jxns_smo.Rdata"))

load(here("processed-data/04_DEA/Tx_analysis/de_tx_nic.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/de_tx_smo.Rdata"))

load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata"))

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

## Define groups of all DE jxn's genes
DEjxns_genes_nic <- unique(de_jxns_nic$newGeneID)[which(!is.na(unique(de_jxns_nic$newGeneID)))]
DEjxns_genes_smo <- unique(de_jxns_smo$newGeneID)[which(!is.na(unique(de_jxns_smo$newGeneID)))]



## Compare DE txs' genes vs DE exons' genes vs DE jxns' genes vs DEG

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


## Venn diagrams
DE_lists<-list(DEG_vs_Txs_vs_Exons_vs_Jxns_all, DEG_vs_Txs_vs_Exons_vs_Jxns_nic, DEG_vs_Txs_vs_Exons_vs_Jxns_smo)
colors<-list(c("rosybrown2", "navajowhite2", "lightblue3", "darkolivegreen2"), c("plum", "lightgoldenrod2", "lightcyan3", "darkolivegreen3"), 
             c("pink1", "khaki2", "lightsteelblue3", "darkolivegreen4"))
venn_plot(DE_lists, colors, "DEG_VS_txs_VS_exons_VS_jxns", c("All", "Nicotine", "Smoking"))

