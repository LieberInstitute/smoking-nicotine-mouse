
## 1.2 Comparison of DE transcripts


load(here("processed-data/04_DEA/Tx_analysis/top_tx_nic.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/de_tx_nic.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/top_tx_smo.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/de_tx_smo.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/rse_tx_brain_pups_nicotine.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/rse_tx_brain_pups_smoking.Rdata"))

load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/results_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/results_pups_smoking_fitted.Rdata"))
load(here("processed-data/05_GO_KEGG/Gene_analysis/intersections.Rdata"))     

load(here("processed-data/04_DEA/Exon_analysis/de_exons_nic.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/de_exons_smo.Rdata"))


### 1.2.1 T-stats plots

### 1.2.1.1 T-stats of tx in nic vs smo

## Function to add DE info of tx in both groups of samples

add_DE_info <-function(top_tx) {
  
  DE<-vector()
  for (i in 1:dim(top_tx)[1]) {
    
    if(top_tx$transcript_id[i] %in% de_tx_nic$transcript_id) {
      ## DE tx in both groups
      if(top_tx$transcript_id[i] %in% de_tx_smo$transcript_id){
        DE<-append(DE, "sig Both")
      }
      ## DE tx in nic only
      else{
        DE<-append(DE, "sig nic")
      }
    }
    
    else if (top_tx$transcript_id[i] %in% de_tx_smo$transcript_id){
      if(!top_tx$transcript_id[i] %in% de_tx_nic$transcript_id){
        ## DE tx in smo only
        DE<-append(DE, "sig smo")
      }
    }
    else {
      ## No DE tx in any group
      DE<-append(DE, "None")
    }
  }
  return(DE)
}



## Compare t-stats of tx from different groups of samples
t_stat_plot <- function(top_tx1, top_tx2, name_1, name_2, title){
  
  ## Correlation coeff
  rho <- cor(top_tx1$t, top_tx2$t, method = "spearman")
  rho_anno = paste0("rho = ", format(round(rho, 2), nsmall = 2))
  
  cols <- cols <- c("deeppink3", "thistle3","navajowhite2", "darkgrey") 
  names(cols)<-c("sig Both", "sig nic", "sig smo", "None")
  alphas <- c( 1, 1, 1,0.5)  
  names(alphas)<-c("sig Both", "sig nic", "sig smo", "None")
  
  ## Merge data
  t_stats<-data.frame(t1=top_tx1$t, t2=top_tx2$t)
  ## Add DE info for both groups 
  t_stats$DE<-add_DE_info(top_tx1)
  t_stats$DE <- factor(t_stats$DE, levels=names(cols))
  
  plot <- ggplot(t_stats, aes(x = t1, y = t2, color=DE, alpha=DE)) +
    geom_point(size = 1.5) +
    labs(x = paste("t-stats", name_1), 
         y = paste("t-stats", name_2),
         title = title, 
         subtitle = rho_anno, 
         parse = T) +
    scale_color_manual(values = cols, labels=names(cols), drop = FALSE) + 
    scale_alpha_manual(values = alphas, labels=names(alphas), drop=FALSE) +
    guides(alpha = 'none') + 
    theme_bw() +
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.text = element_text(size=11),
          legend.title = element_text(size=12))
  
  plot
  ggsave(filename=paste("plots/04_DEA/02_Comparisons/Tx_analysis/t_stats_", gsub(" ", "_", title), 
                        ".pdf", sep=""), height = 20, width = 25, units = "cm")
}


#####################################
# Smoking vs nicotine (transcripts)
#####################################

t_stat_plot(top_tx_nic, top_tx_smo, "Nicotine pups", "Smoking pups", "Nic vs Smo DE tx")





### 1.2.1.2 T-stats of genes vs transcripts

## Function to add DE info of tx and genes

add_DE_info_tx_vs_genes <-function(t_stats) {
  
  DE<-vector()
  for (i in 1:dim(t_stats)[1]) {
    
    ## DE tx from DEG 
    if (t_stats$adj.P.Val_tx[i]<0.05 && t_stats$adj.P.Val_genes[i]<0.05){
      DE<-append(DE, "sig Both")
    }
    
    ## DE tx only
    else if (t_stats$adj.P.Val_tx[i]<0.05 && t_stats$adj.P.Val_genes[i]>=0.05){
      DE<-append(DE, "sig Tx")
    }
    
    ## DE genes only
    else if(t_stats$adj.P.Val_tx[i]>=0.05 && t_stats$adj.P.Val_genes[i]<0.05){
      DE<-append(DE, "sig Gene")
    }
    
    ## Neither DE tx nor DE genes
    else {
      DE<-append(DE, "None")      
    }
  }
  return(DE)
}



## Create plots of t-stats of tx vs genes

t_stat_tx_vs_genes<- function(expt){
  
  top_tx<-eval(parse_expr(paste("top_tx_", substr(expt,1,3), sep="")))
  top_genes<-eval(parse_expr(paste("top_genes_pups_", expt, "_fitted", sep="")))
  de_tx<-eval(parse_expr(paste("de_tx_", substr(expt,1,3), sep="")))
  
  if (expt=="nicotine"){
    abs_t_tx=6
  }
  else{
    abs_t_tx=4
  }
  
  ## Transcripts' genes
  tx_genes<-unique(top_tx$ensembl_id)
  
  ## Common genes
  tx_genes<-tx_genes[which(tx_genes %in% top_genes$gencodeID)]
  
  ## Extract transcripts' info
  t_stats<-top_tx[which(top_tx$ensembl_id %in% tx_genes), c("transcript_id", "adj.P.Val", "Symbol", 
                                                                    "t", "ensembl_id", "logFC")]
  colnames(t_stats)[2]<-"adj.P.Val_tx"
  colnames(t_stats)[4]<-"t_tx"
  colnames(t_stats)[6]<-"logFC_tx"
  
  ## Add t-stats and FDRs of transcripts' genes 
  t_genes<-vector()
  FDRs<-vector()
  logFCs <- vector()
  
  for (i in 1:dim(t_stats)[1]){
    t<-top_genes[which(top_genes$gencodeID==t_stats$ensembl_id[i]), "t"]
    FDR<-top_genes[which(top_genes$gencodeID==t_stats$ensembl_id[i]), "adj.P.Val"]
    logFC <- top_genes[which(top_genes$gencodeID==t_stats$ensembl_id[i]), "logFC"]
    
    t_genes<-append(t_genes, t)
    FDRs<-append(FDRs, FDR)
    logFCs <- append(logFCs, logFC)
  }
  t_stats$t_genes<-t_genes
  t_stats$adj.P.Val_genes<-FDRs
  t_stats$logFC_genes <- logFCs
  
  ## Correlation coeff between t-stats of genes and transcripts
  rho <- cor(t_stats$t_tx, t_stats$t_genes, method = "spearman")
  rho_anno = paste0("rho = ", format(round(rho, 2), nsmall = 2))
  
  cols <- c("coral2","pink", "lightgoldenrod3", "grey") 
  names(cols) <- c("sig Both","sig Gene", "sig Tx", "None")
  alphas <- c(1, 1, 1, 0.2)  
  names(alphas) <- names(cols)
  
  ## Add DE info for both groups 
  t_stats$DE<-add_DE_info_tx_vs_genes(t_stats)
  t_stats$DE <- factor(t_stats$DE, levels=names(cols))
  
  
  ## Gene-tx symbols of:
  
  # ------- Top 3 DE txs ------- 
  paste0(de_tx[order(de_tx$adj.P.Val, decreasing = FALSE), c('Symbol')][1:3], '-', 
                       de_tx[order(de_tx$adj.P.Val, decreasing = FALSE), c('transcript_id')][1:3])
  top_de_txs <- de_tx[order(de_tx$adj.P.Val, decreasing = FALSE), c('transcript_id')][1:3]
  # [1] "Phf3-ENSMUST00000185521.1"    "Pnisr-ENSMUST00000029911.11"  "Ankrd11-ENSMUST00000212937.1"     <- for nic
  # [1] "Ccnb2-ENSMUST00000034742.7" "Top2a-ENSMUST00000068031.7" "Mt2-ENSMUST00000034214.7"    <- for smo
  
  
  #  ------- DE tx from non-DEG with |t-stat|>abs_t_tx ------- 
  if (expt == 'nicotine'){
    paste0(t_stats[which(t_stats$DE=='sig Tx' & (abs(t_stats$t_tx)>abs_t_tx)), 'Symbol'], '-', 
           t_stats[which(t_stats$DE=='sig Tx' & (abs(t_stats$t_tx)>abs_t_tx)), 'transcript_id'])
    de_txs_nonDEGs <- t_stats[which(t_stats$DE=='sig Tx' & (abs(t_stats$t_tx)>abs_t_tx)), 'transcript_id']
    # [1] "Phf3-ENSMUST00000185521.1"    "Trpc4-ENSMUST00000199359.1"   "Ankrd11-ENSMUST00000212937.1" "Bcl11a-ENSMUST00000124148.1" 
    # [5] "Galnt10-ENSMUST00000108846.1" "Scaf11-ENSMUST00000227133.1"    <- for nic
    
    ## Ignore Galnt10-ENSMUST00000108846.1 in nic
    de_txs_nonDEGs <- de_txs_nonDEGs[-which(de_txs_nonDEGs=='ENSMUST00000108846.1')]
  }
  
  ## Take only Btf3, Srsf6, Cyhr1 and H13 txs for smo
  else {
   nonDEGs <- c("Btf3", "Srsf6", "Cyhr1", "H13")
   paste0(t_stats[which(t_stats$DE=='sig Tx' & (abs(t_stats$t_tx)>abs_t_tx) & t_stats$Symbol %in% nonDEGs), 'Symbol'], '-', 
          t_stats[which(t_stats$DE=='sig Tx' & (abs(t_stats$t_tx)>abs_t_tx) & t_stats$Symbol %in% nonDEGs), 'transcript_id'])
   # [1] "H13-ENSMUST00000125366.7"    "H13-ENSMUST00000089059.8"    "Srsf6-ENSMUST00000017065.14" "Btf3-ENSMUST00000022163.14"  "Cyhr1-ENSMUST00000081291.12"
   # [6] "Cyhr1-ENSMUST00000176274.1"    <- for smo
   de_txs_nonDEGs <- t_stats[which(t_stats$DE=='sig Tx' & (abs(t_stats$t_tx)>abs_t_tx) & t_stats$Symbol %in% nonDEGs), 'transcript_id']
  }
  
  
  #  ------- DE tx whose DEG have an opposite sign in logFC ------- 
  if (expt == 'nicotine'){
    paste0(t_stats[which(t_stats$DE=='sig Both' & (sign(t_stats$logFC_tx)!=sign(t_stats$logFC_genes))), 'Symbol'], '-',
           t_stats[which(t_stats$DE=='sig Both' & (sign(t_stats$logFC_tx)!=sign(t_stats$logFC_genes))), 'transcript_id'])
    de_txs_DEGs <- t_stats[which(t_stats$DE=='sig Both' & (sign(t_stats$logFC_tx)!=sign(t_stats$logFC_genes))), 'transcript_id']
    # [1] "Pnisr-ENSMUST00000148561.1"   "Dcun1d5-ENSMUST00000216770.1" "Dgcr8-ENSMUST00000115633.2"   <- for nic
  }
  else{
    ## Ignore for smo
    de_txs_DEGs <- NULL
  }

  
  #  ------- DE tx from DEGs with up and down txs ------- 
  ## Up and down transcripts' genes
  tx_up_genes<-unique(de_tx[which(de_tx$logFC>0),"ensembl_id"])
  tx_down_genes<-unique(de_tx[which(de_tx$logFC<0),"ensembl_id"])
  ## Genes with up and down tx
  interest_genes<-intersect(tx_up_genes, tx_down_genes)
  ## Retain only the transcripts' genes that were considered at the gene level
  interest_genes<-intersect(interest_genes, tx_genes)
  
  if (expt == 'nicotine'){
    ## Only DEGs and DE txs
    paste0(t_stats[which(t_stats$ensembl_id %in% interest_genes & t_stats$adj.P.Val_genes< 0.05 & t_stats$adj.P.Val_tx<0.05), 'Symbol'], '-',
           t_stats[which(t_stats$ensembl_id %in% interest_genes & t_stats$adj.P.Val_genes< 0.05 & t_stats$adj.P.Val_tx<0.05), 'transcript_id'])
    de_txs_up_down_DEGs <- t_stats[which(t_stats$ensembl_id %in% interest_genes & t_stats$adj.P.Val_genes< 0.05 & t_stats$adj.P.Val_tx<0.05), 'transcript_id']
    # [1] "Pnisr-ENSMUST00000029911.11"  "Pnisr-ENSMUST00000148561.1"   "Dcun1d5-ENSMUST00000215683.1" "Dcun1d5-ENSMUST00000216770.1"    <- for nic
  }
  else{
    ## DE txs from the DEGs Meaf6, Ivns1abp, Morf4l2, Sin3b and Ppp2r5c for smo
    DEGs <- c("Meaf6", "Ivns1abp", "Morf4l2", "Sin3b", "Ppp2r5c")
    paste0(t_stats[which(t_stats$ensembl_id %in% interest_genes & t_stats$adj.P.Val_genes< 0.05 & t_stats$adj.P.Val_tx<0.05 & t_stats$Symbol %in% DEGs), 'Symbol'], '-',
           t_stats[which(t_stats$ensembl_id %in% interest_genes & t_stats$adj.P.Val_genes< 0.05 & t_stats$adj.P.Val_tx<0.05 & t_stats$Symbol %in% DEGs), 'transcript_id'])
    # [1] "Ivns1abp-ENSMUST00000023918.12" "Ivns1abp-ENSMUST00000111887.9"  "Meaf6-ENSMUST00000154689.7"     "Meaf6-ENSMUST00000184205.7"    
    # [5] "Sin3b-ENSMUST00000109950.4"     "Sin3b-ENSMUST00000004494.15"    "Sin3b-ENSMUST00000212095.1"     "Ppp2r5c-ENSMUST00000221715.1"  
    # [9] "Ppp2r5c-ENSMUST00000109832.2"   "Morf4l2-ENSMUST00000080411.12"  "Morf4l2-ENSMUST00000169418.7"   "Morf4l2-ENSMUST00000113095.7"  
    de_txs_up_down_DEGs <- t_stats[which(t_stats$ensembl_id %in% interest_genes & t_stats$adj.P.Val_genes< 0.05 & t_stats$adj.P.Val_tx<0.05 & t_stats$Symbol %in% DEGs), 'transcript_id']
  }
  
  tx_symbols<-vector()
  
  for (i in 1:dim(t_stats)[1]) {
    if (t_stats$transcript_id[i] %in% c(top_de_txs, de_txs_nonDEGs, de_txs_DEGs, de_txs_up_down_DEGs)){
      tx_symbols <- append(tx_symbols, paste0(t_stats$Symbol[i], '-', t_stats$transcript_id[i]))
    }
    else {
      tx_symbols <- append(tx_symbols, NA)
    }
  }
  t_stats$tx_symbols<-tx_symbols
  
  ## Plot
  
  if(expt == 'nicotine'){
    plot <- ggplot(t_stats, aes(x = t_genes, y = t_tx, color=DE, alpha=DE, label= tx_symbols)) +
      geom_point(size = 1.2) +
      labs(x = "t-stats genes", 
           y = "t-stats tx",
           title = paste(capitalize(expt),"genes vs tx", sep=" "), 
           subtitle = rho_anno, 
           parse = T) +
      geom_label_repel(data=subset(t_stats, !tx_symbols %in% c('Trpc4-ENSMUST00000199359.1', 'Dgcr8-ENSMUST00000115633.2')), 
                       aes(fontface = 'bold'), fill='white',
                       size=2.8,
                       max.overlaps = Inf,
                       min.segment.length = unit(0, "cm"),
                       point.padding = unit(0.1, "cm"),
                       box.padding = 0.2,
                       label.padding = 0.2,
                       label.size = 0.2,
                       nudge_y=0.8,
                       nudge_x = 1,
                       show.legend=FALSE) +
      geom_label_repel(data=subset(t_stats, tx_symbols %in% c('Trpc4-ENSMUST00000199359.1', 'Dgcr8-ENSMUST00000115633.2')), 
                       aes(fontface = 'bold'), fill='white',
                       size=2.8,
                       max.overlaps = Inf,
                       min.segment.length = unit(0, "cm"),
                       point.padding = unit(0.1, "cm"),
                       box.padding = 0.2,
                       label.padding = 0.2,
                       label.size = 0.2,
                       nudge_y=-1.9,
                       nudge_x = -1,
                       show.legend=FALSE) +
      theme_bw() +
      guides(alpha = 'none') + 
      scale_color_manual(values = cols, labels=names(cols), drop = FALSE) + 
      scale_alpha_manual(values = alphas, labels=names(alphas), drop=FALSE) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10),
            legend.text = element_text(size=11),
            legend.title = element_text(size=12))
  }
  
  else{
    plot <- ggplot(t_stats, aes(x = t_genes, y = t_tx, color=DE, alpha=DE, label= tx_symbols)) +
      geom_point(size = 1.2) +
      labs(x = "t-stats genes", 
           y = "t-stats tx",
           title = paste(capitalize(expt),"genes vs tx", sep=" "), 
           subtitle = rho_anno, 
           parse = T) +
      ## Label only the most significant DE txs per gene
      geom_label_repel(data=subset(t_stats, tx_symbols %in% c('Btf3-ENSMUST00000022163.14', 'Cyhr1-ENSMUST00000176274.1')),
                       aes(fontface = 'bold'), fill='white',
                       size=2.8,
                       max.overlaps = Inf,
                       min.segment.length = unit(0, "cm"),
                       point.padding = unit(0.1, "cm"),
                       box.padding = 0.2,
                       label.padding = 0.2,
                       label.size = 0.2,
                       nudge_y= 3.8,
                       nudge_x = -3,
                       show.legend=FALSE) +
      geom_label_repel(data=subset(t_stats, tx_symbols %in% c('Srsf6-ENSMUST00000017065.14')),
                       aes(fontface = 'bold'), fill='white',
                       size=2.8,
                       max.overlaps = Inf,
                       min.segment.length = unit(0, "cm"),
                       point.padding = unit(0.1, "cm"),
                       box.padding = 0.2,
                       label.padding = 0.2,
                       label.size = 0.2,
                       nudge_y= -3,
                       nudge_x = 3.5,
                       show.legend=FALSE) +
    geom_label_repel(data=subset(t_stats, tx_symbols %in% c('Ppp2r5c-ENSMUST00000221715.1')),
                     aes(fontface = 'bold'), fill='white',
                     size=2.8,
                     max.overlaps = Inf,
                     min.segment.length = unit(0, "cm"),
                     point.padding = unit(0.1, "cm"),
                     box.padding = 0.2,
                     label.padding = 0.2,
                     label.size = 0.2,
                     nudge_y= 2,
                     nudge_x = -2,
                     show.legend=FALSE) +
      geom_label_repel(data=subset(t_stats, tx_symbols %in% c('Ivns1abp-ENSMUST00000111887.9')),
                       aes(fontface = 'bold'), fill='white',
                       size=2.8,
                       max.overlaps = Inf,
                       min.segment.length = unit(0, "cm"),
                       point.padding = unit(0.1, "cm"),
                       box.padding = 0.2,
                       label.padding = 0.2,
                       label.size = 0.2,
                       nudge_y= 0,
                       nudge_x = -2,
                       show.legend=FALSE) +
      geom_label_repel(data=subset(t_stats, tx_symbols %in% c('Meaf6-ENSMUST00000184205.7')),
                       aes(fontface = 'bold'), fill='white',
                       size=2.8,
                       max.overlaps = Inf,
                       min.segment.length = unit(0, "cm"),
                       point.padding = unit(0.1, "cm"),
                       box.padding = 0.2,
                       label.padding = 0.2,
                       label.size = 0.2,
                       nudge_y= 2,
                       nudge_x = 1,
                       show.legend=FALSE) +
      geom_label_repel(data=subset(t_stats, tx_symbols %in% c('Mt2-ENSMUST00000034214.7', 'Sin3b-ENSMUST00000109950.4',
                                                              'Ppp2r5c-ENSMUST00000109832.2')),
                       aes(fontface = 'bold'), fill='white',
                       size=2.8,
                       max.overlaps = Inf,
                       min.segment.length = unit(0, "cm"),
                       point.padding = unit(0.1, "cm"),
                       box.padding = 0.2,
                       label.padding = 0.2,
                       label.size = 0.2,
                       nudge_y= 0.2,
                       nudge_x = 0.1,
                       show.legend=FALSE) +
    geom_label_repel(data=subset(t_stats, tx_symbols %in% c('Meaf6-ENSMUST00000154689.7')),
                     aes(fontface = 'bold'), fill='white',
                     size=2.8,
                     max.overlaps = Inf,
                     min.segment.length = unit(0, "cm"),
                     point.padding = unit(0.1, "cm"),
                     box.padding = 0.2,
                     label.padding = 0.2,
                     label.size = 0.2,
                     nudge_y= -0.5,
                     nudge_x = 0.1,
                     show.legend=FALSE) +
      geom_label_repel(data=subset(t_stats, tx_symbols %in% c('Morf4l2-ENSMUST00000169418.7')),
                       aes(fontface = 'bold'), fill='white',
                       size=2.8,
                       max.overlaps = Inf,
                       min.segment.length = unit(0, "cm"),
                       point.padding = unit(0.1, "cm"),
                       box.padding = 0.2,
                       label.padding = 0.2,
                       label.size = 0.2,
                       nudge_y= -0.01,
                       nudge_x = 0.5,
                       show.legend=FALSE) +
      geom_label_repel(data=subset(t_stats, tx_symbols %in% c('Morf4l2-ENSMUST00000080411.12', 
                                                              'Sin3b-ENSMUST00000004494.15')),
                       aes(fontface = 'bold'), fill='white',
                       size=2.8,
                       max.overlaps = Inf,
                       min.segment.length = unit(0, "cm"),
                       point.padding = unit(0.1, "cm"),
                       box.padding = 0.2,
                       label.padding = 0.2,
                       label.size = 0.2,
                       nudge_y= -0.9,
                       nudge_x = 0.5,
                       show.legend=FALSE) +
      geom_label_repel(data=subset(t_stats, tx_symbols %in% c('Top2a-ENSMUST00000068031.7', 'Ccnb2-ENSMUST00000034742.7',
                                                              'Ivns1abp-ENSMUST00000023918.12')),
                       aes(fontface = 'bold'), fill='white',
                       size=2.8,
                       max.overlaps = Inf,
                       min.segment.length = unit(0, "cm"),
                       point.padding = unit(0.1, "cm"),
                       box.padding = 0.2,
                       label.padding = 0.2,
                       label.size = 0.2,
                       nudge_y= -1.8,
                       nudge_x = -3,
                       show.legend=FALSE) +
      geom_label_repel(data=subset(t_stats, tx_symbols %in% c('H13-ENSMUST00000089059.8')),
                       aes(fontface = 'bold'), fill='white',
                       size=2.8,
                       max.overlaps = Inf,
                       min.segment.length = unit(0, "cm"),
                       point.padding = unit(0.1, "cm"),
                       box.padding = 0.2,
                       label.padding = 0.2,
                       label.size = 0.2,
                       nudge_y= -0.5,
                       nudge_x = -6,
                       show.legend=FALSE) +
      theme_bw() +
      guides(alpha = 'none') + 
      scale_color_manual(values = cols, labels=names(cols), drop = FALSE) + 
      scale_alpha_manual(values = alphas, labels=names(alphas), drop=FALSE) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10),
            legend.text = element_text(size=11),
            legend.title = element_text(size=12))
  }

  
  plot
  ggsave(filename=paste("plots/04_DEA/02_Comparisons/Tx_analysis/t_stats_tx_vs_genes_", substr(expt,1,3), 
                        ".pdf", sep=""), height = 14.5, width = 18, units = "cm")
  
}


##################################### 
# Nicotine genes vs nicotine tx
#####################################
expt<-"nicotine"
t_stat_tx_vs_genes(expt)


##################################### 
# Smoking genes vs smoking tx
#####################################
expt<-"smoking"
t_stat_tx_vs_genes(expt)







## 1.2.2 Boxplots of relevant genes and their tx

## Each boxplot
create_boxplot<- function(counts, y, title, q_value, FC, tx_tpm){
  
  ## Boxplot
  p<-ggplot(data=counts, 
            aes(x=Group,y=eval(parse_expr(y)))) + 
    geom_boxplot(outlier.color = "#FFFFFFFF", width=0.65) +
    geom_jitter(aes(colour=Group),shape=16, 
                position=position_jitter(0.2), size=2.7) +
    theme_bw() +
    labs(x = "Group", y = "lognorm counts",
         title = title,
         ## Add FDR, FC and % of tpm
         subtitle=paste(" FDR:", q_value, "       ", tx_tpm,  
                        "\n", "FC:", FC)) +
    scale_color_manual(values=c("Control" = "seashell3", "Experimental" = "orange3")) +
    scale_x_discrete(labels=c("Control"="Ctrl","Experimental"="Expt")) +
    theme(plot.margin=unit (c(0.4,0.4,0.4,0.4), 'cm'), 
          legend.position = "none",
          plot.title = element_text(hjust=0.5, size=12, face="bold"), 
          plot.subtitle = element_text(size = 10),
          axis.title = element_text(size = (12)),
          axis.text = element_text(size = 10.5)) 
  
  print(p)
}



## Boxplots
gene_tx_boxplots<- function(expt, gene, tx1, tx2){
  
  top_tx<-eval(parse_expr(paste("top_tx_", substr(expt,1,3), sep="")))
  top_genes<-eval(parse_expr(paste("top_genes_pups_", expt, "_fitted", sep="")))
  results_genes<-eval(parse_expr(paste("results_pups_", expt, "_fitted", sep="")))
  RSE<-eval(parse_expr((paste("rse_tx_brain_pups_", expt, sep=""))))
  
  ## Expression values of all genes
  vGene<-results_genes[[1]][[2]]
  
  ## Get ensembl ID and symbol of the gene
  if(! gene %in% rownames(vGene)){
    gene_symbol<-gene
    gene<-vGene$genes[which(vGene$genes$Symbol==gene), "gencodeID"]
  }
  else {
    gene_symbol<-vGene$genes[which(vGene$genes$ensemblID==gene), "Symbol"]
  }
  
  ## Gene ID for plot: Symbol-ensemblID
  gene_ID<-paste(gene_symbol, "-", gene, sep="")
  
  ## Extract expression values of the gene
  counts_gene<-vGene$E[rownames(vGene)==gene,]
  
  ## Gene's transcripts 
  gene_txs<-top_tx[which(top_tx$ensembl_id==gene),]
  ## Order tx by FDR and extract the top 3
  top_gene_txs<-gene_txs[order(gene_txs$adj.P.Val),"transcript_id"][1:3]
  
  ## Obtain total tpm of the gene's transcripts
  gene_txs_tpm<-assays(RSE)$tpm[which(rowData(RSE)$gene_id==gene),]
  if(is.null(dim(gene_txs_tpm))){
    gene_txs_tpm<-sum(gene_txs_tpm)
  }
  else{
    gene_txs_tpm<-sum(apply(gene_txs_tpm, 1, sum)) 
  }
 
  
  ## Add specific transcripts
  if (!is.null(tx1) | !is.null(tx2)){
    
    if (!is.null(tx1) &  is.null(tx2)){
      tx<-tx1
      ## Extract tx ID 
      tx_ID<-strsplit(tx, "−")[[1]][2]
      if (! tx_ID %in% top_gene_txs[1:3]){
        top_gene_txs<-append(tx_ID, top_gene_txs) 
      }
    }
    else if(!is.null(tx2) & is.null(tx1)){
      tx<-tx2
      tx_ID<-strsplit(tx, "−")[[1]][2]
      if (! tx_ID %in% top_gene_txs[1:3]){
        top_gene_txs<-append(tx_ID, top_gene_txs) 
      }
    }
    
    ## Add two txs
    else {
      ## Tx IDs
      for (i in 1:2){
        tx<-eval(parse_expr(paste("tx", i, sep="")))
        tx_ID<-strsplit(tx, "−")[[1]][2]
        if (! tx_ID %in% top_gene_txs[1:3]){
          top_gene_txs<-append(tx_ID, top_gene_txs) 
        }
      }
    }
  }
  
  ## Log-expression values of the txs
  logcounts<-assays(RSE)$logcounts
  counts_txs<-logcounts[which(rownames(logcounts) %in% top_gene_txs[1:3]),]
  if(!is.null(dim(counts_txs))){
    counts_txs<-t(counts_txs)
  }
  
  ## Data frame with all expression values and samples' group
  counts<-cbind(counts_gene, counts_txs, "Group"=vGene$targets$Group)
  counts<-as.data.frame(counts)
  counts$counts_gene<-as.numeric(counts$counts_gene)
  colnames(counts)[1]<-gene
  if(!is.null(dim(counts_txs))){
    counts[,2]<-as.numeric(counts[,2])
    counts[,3]<-as.numeric(counts[,3])
    counts[,4]<-as.numeric(counts[,4])
  }
  else{
    counts$counts_txs <- as.numeric(counts$counts_txs)
    colnames(counts)[2] <- top_gene_txs[1]
  }

  
  plots<-list()
  for (i in 1:(dim(counts)[2]-1)){
    
    ## Boxplot of the gene
    if (i==1){
      ## q-value for the gene
      q_value=signif(top_genes[which(top_genes$gencodeID==gene), "adj.P.Val"], digits = 3)
      FC=signif(2**(top_genes[which(top_genes$gencodeID==gene), "logFC"]), digits = 3)
      y<-gene
      title<-gene_ID
      tx_tpm<-NULL
    }
    
    ## Boxplot of the txs
    else {
      q_value=signif(top_tx[which(top_tx$transcript_id==colnames(counts)[i]), "adj.P.Val"], digits = 3)
      FC=signif(2**(top_tx[which(top_tx$transcript_id==colnames(counts)[i]), "logFC"]), digits = 3)
      y<-colnames(counts)[i]
      ## Tx name for plot
      title<-paste(gene_symbol, "-", y, sep="")
      tx_tpm<-sum(assays(RSE)$tpm[which(rowData(RSE)$transcript_id==colnames(counts)[i]),])
      tx_tpm<-paste("tpm prop: ", signif(100*tx_tpm/gene_txs_tpm, 4), "%", sep="")
      
    }
    
    ## Plots
    p<-create_boxplot(counts, y, title, q_value, FC, tx_tpm)
    plots[[i]]<-p
  }
  
  if(!is.null(dim(counts_txs))){
    plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol = 4)
    ggsave(here(paste("plots/04_DEA/02_Comparisons/Tx_analysis/Boxplots_", gene_symbol, ".pdf", sep="")), 
           width = 34, height = 10, units = "cm")
  }
  else{
    plot_grid(plots[[1]], plots[[2]], ncol = 2)
    ggsave(here(paste("plots/04_DEA/02_Comparisons/Tx_analysis/Boxplots_", gene_symbol, ".pdf", sep="")), 
           width = 17, height = 10, units = "cm")
  }
 
  
}



#########################################
## Boxplots of nicotine genes and txs
#########################################

#### Non-DE genes but with DE txs (including top DE txs already) ####

## Phf3: predicted to be involved in transcription, DNA-templated;
## expressed in CNS
gene_tx_boxplots("nicotine", "Phf3", "Phf3−ENSMUST00000185521.1", NULL)

## Ankrd11: acts in head morphogenesis; expressed in cerebral cortex
gene_tx_boxplots("nicotine", "Ankrd11", "Ankrd11−ENSMUST00000212937.1", NULL)

## Trpc4: acts upstream of or within gamma-aminobutyric acid secretion and 
## oligodendrocyte differentiation; expressed in brain
gene_tx_boxplots("nicotine", "Trpc4", "Trpc4−ENSMUST00000199359.1", NULL)

## Scaf11: predicted to be involved in spliceosomal complex assembly;
## expressed in diencephalon lateral wall ventricular layer; heart ventricle; 
## midbrain ventricular layer; and telencephalon ventricular layer
gene_tx_boxplots("nicotine", "Scaf11", "Scaf11−ENSMUST00000227133.1", NULL)

## Bcl11a: expressed in NS and sensoty organ; human ortholog(s) of this gene 
## are implicated in autism spectrum disorder and schizophrenia
gene_tx_boxplots("nicotine", "Bcl11a", "Bcl11a−ENSMUST00000124148.1", NULL)



#### DE tx whose DEG have opposite direction of regulation ####

## Dgcr8: enables primary miRNA binding activity; located in postsynaptic density;
## expressed in CNS, brain, branchial arch and facial prominence
gene_tx_boxplots("nicotine", "Dgcr8", "Dgcr8−ENSMUST00000115633.2", NULL)



#### DEG with Up and Down DE tx ####

## Pnisr: predicted to be active in presynaptic active zone; expressed in NS
gene_tx_boxplots("nicotine", "Pnisr", "Pnisr−ENSMUST00000029911.11", 
                 "Pnisr−ENSMUST00000148561.1")

## Dcun1d5: predicted to be involved in protein modification by small protein 
## conjugation or removal, protein neddylation, and regulation of cell growth;
## expressed in NS
gene_tx_boxplots("nicotine", "Dcun1d5", "Dcun1d5−ENSMUST00000215683.1", 
                 "Dcun1d5−ENSMUST00000216770.1")



#########################################
## Boxplots of smoking genes and txs
#########################################

#### Top DE txs ####
gene_tx_boxplots('smoking', 'Top2a', 'Top2a−ENSMUST00000068031.7', NULL)
gene_tx_boxplots('smoking', 'Ccnb2', 'Ccnb2−ENSMUST00000034742.7', NULL)
gene_tx_boxplots('smoking', 'Mt2', 'Mt2−ENSMUST00000034214.7', NULL)



#### Non-DE genes but with DE txs ####

## Btf3: acts upstream of or within in utero embryonic development; expressed in brain
gene_tx_boxplots("smoking", "Btf3", "Btf3−ENSMUST00000022163.14", NULL)

## Srsf6: critical for mRNA splicing, involved in mRNA export from the nucleus; 
## expressed in CNS
gene_tx_boxplots("smoking", "Srsf6", "Srsf6−ENSMUST00000017065.14", NULL)

## Cyhr1: predicted to enable zinc ion binding activity; located in cytoplasm and nuclear envelope
gene_tx_boxplots("smoking", "Cyhr1", "Cyhr1−ENSMUST00000176274.1", 
                 "Cyhr1−ENSMUST00000081291.12")

## H13: acts upstream of or within in utero embryonic development; expressed in  brain, 
## embryo ectoderm, extraembryonic component, hemolymphoid system and intervertebral disc
gene_tx_boxplots("smoking", "H13", "H13−ENSMUST00000089059.8", 
                 "H13−ENSMUST00000125366.7")


#### DEG with Up and Down DE tx ####

## Ivns1abp: acts upstream of or within negative regulation of intrinsic apoptotic 
## signaling pathway; expressed in NS
gene_tx_boxplots("smoking", "Ivns1abp", "Ivns1abp−ENSMUST00000111887.9", 
                 "Ivns1abp−ENSMUST00000023918.12")

## Morf4l2: predicted to act upstream of or within DNA repair and regulation of growth
## biased expression in placenta adult and CNS
gene_tx_boxplots("smoking", "Morf4l2", "Morf4l2−ENSMUST00000169418.7", 
                 "Morf4l2−ENSMUST00000080411.12")

## Sin3b: involved in skeletal muscle tissue development; expressed in CNS
gene_tx_boxplots("smoking", "Sin3b", "Sin3b−ENSMUST00000109950.4", 
                 "Sin3b−ENSMUST00000004494.15")

## Ppp2r5c: predicted to be involved in negative regulation of cell population 
## proliferation, protein dephosphorylation, and signal transduction by p53 class mediator;
## expressed in CNS
gene_tx_boxplots("smoking", "Ppp2r5c", "Ppp2r5c−ENSMUST00000221715.1", 
                 "Ppp2r5c−ENSMUST00000109832.2")

## Meaf6: predicted to act upstream of or within chromatin organization;
## expressed in CNS and whole brain
gene_tx_boxplots("smoking", "Meaf6", "Meaf6−ENSMUST00000184205.7", "Meaf6−ENSMUST00000154689.7")








### 1.2.2 Venn diagrams

## Function to create multiple Venn diagrams
venn_plot<-function(DE_lists, colors, name, titles){
  
  if (name=="smo_VS_nic_DE_txs"){
    height=10
    width=15
    margin=0.2
    dist=0.1
    cat=0.75
    cex=0.8
    main_pos = c(0.5,0.476)
    main_cex =1
  }
  
  else if (name=="DEG_VS_txs_genes"){
    height=12
    width=16
    margin=5
    dist=0.07
    cat=0.75
    cex=0.8
    main_pos = c(0.5,0.476)
    main_cex =1
  }
  
  else if (name=="smo_VS_nic_DE_txs_genes"){
    height=12
    width=16
    margin=0.2
    dist=0.09
    cat=0.75
    cex=0.8
    main_pos = c(0.5,0.476)
    main_cex =1
  }
  
  else if (name=="intersections_DEG_VS_txs_genes"){
    height=10
    width=20
    margin=7
    dist=0.06
    cat=0.9
    cex=1
    main_pos = c(0.5,0.476)
    main_cex =1
  }
  
  else if (name=="DEG_VS_txs_VS_exons"){
    height=15
    width=20
    margin=7
    dist=0.08
    cat=0.9
    cex=1
    main_pos = c(0.5,0.476)
    main_cex =1
  }
  
  else if (name=="intersections_DEG_VS_txs_VS_exons"){
    height=15
    width=30
    margin=1.15
    dist=0.09
    cat=1.3
    cex=1.3
    main_pos = c(0.5,0.375)
    main_cex =1.5
  }
  
  plots<-list()
  pdf(file = paste("plots/04_DEA/02_Comparisons/Tx_analysis/Venn_", name, ".pdf", sep=""), 
      height = height, width = width)
  for (i in 1:length(DE_lists)){
    v<-venn.diagram(DE_lists[[i]], fill=colors[[i]], alpha = rep(0.5, length(DE_lists[[i]])), 
                    lwd =0, margin=margin, cat.cex=cat, cex=cex, height = 10, width = 12.5, units = "cm", 
                    cat.dist=rep(dist, length(DE_lists[[i]])), filename=NULL, main = titles[i], 
                    main.pos = main_pos, main.cex = main_cex, disable.logging=TRUE)
    plots[[i]]<-v
  }
  
  if (i==6){
    gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], 
                            plots[[6]], ncol=3)
    dev.off()
  }
  else if (i==9) {
    gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
                            plots[[6]], plots[[7]], plots[[8]], plots[[9]], ncol=4)
    dev.off()
  }
  else {
    gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
                            plots[[6]], plots[[7]], plots[[8]], ncol=4)
    dev.off()    
  }
}



## Define groups of DE txs
nic_up<-de_tx_nic[which(de_tx_nic$logFC>0),"transcript_id"]
smo_up<-de_tx_smo[which(de_tx_smo$logFC>0),"transcript_id"]
nic_down<-de_tx_nic[which(de_tx_nic$logFC<0),"transcript_id"]
smo_down<-de_tx_smo[which(de_tx_smo$logFC<0),"transcript_id"]


## Define groups of DE txs' genes
nic_up_genes<-unique(de_tx_nic[which(de_tx_nic$logFC>0),"ensembl_id"])
smo_up_genes<-unique(de_tx_smo[which(de_tx_smo$logFC>0),"ensembl_id"])
nic_down_genes<-unique(de_tx_nic[which(de_tx_nic$logFC<0),"ensembl_id"])
smo_down_genes<-unique(de_tx_smo[which(de_tx_smo$logFC<0),"ensembl_id"])

smoUp_nicUp_genes<-intersect(nic_up_genes, smo_up_genes)
smoDown_nicDown_genes<-intersect(nic_down_genes, smo_down_genes)
smoUp_nicDown_genes<-intersect(nic_down_genes, smo_up_genes)
smoDown_nicUp_genes<-intersect(nic_up_genes, smo_down_genes)
only_up_nic_genes<-nic_up_genes[which(! (nic_up_genes %in% smo_up_genes | 
                                           nic_up_genes %in% smo_down_genes | 
                                           nic_up_genes %in% nic_down_genes))]
only_up_smo_genes<-smo_up_genes[which(! (smo_up_genes %in% smo_down_genes | 
                                           smo_up_genes %in% nic_down_genes | 
                                           smo_up_genes %in% nic_up_genes))]
only_down_nic_genes<-nic_down_genes[which(! (nic_down_genes %in% smo_up_genes | 
                                               nic_down_genes %in% smo_down_genes | 
                                               nic_down_genes %in% nic_up_genes))]
only_down_smo_genes<-smo_down_genes[which(! (smo_down_genes %in% smo_up_genes | 
                                               smo_down_genes %in% nic_down_genes | 
                                               smo_down_genes %in% nic_up_genes))]


## Define groups of DEG
nic_DEG_up<-de_genes_pups_nicotine_fitted[which(de_genes_pups_nicotine_fitted$logFC>0),"gencodeID"]
nic_DEG_down<-de_genes_pups_nicotine_fitted[which(de_genes_pups_nicotine_fitted$logFC<0),"gencodeID"]
smo_DEG_up<-de_genes_pups_smoking_fitted[which(de_genes_pups_smoking_fitted$logFC>0),"gencodeID"]
smo_DEG_down<-de_genes_pups_smoking_fitted[which(de_genes_pups_smoking_fitted$logFC<0),"gencodeID"]

only_up_nic_DEG<-intersections[[1]]$gencodeID
only_up_smo_DEG<-intersections[[2]]$gencodeID
only_down_nic_DEG<-intersections[[3]]$gencodeID
only_down_smo_DEG<-intersections[[4]]$gencodeID
smoUp_nicUp_DEG<-intersections[[5]]$gencodeID
smoDown_nicDown_DEG<-intersections[[6]]$gencodeID
smoUp_nicDown_DEG<-intersections[[7]]$gencodeID
smoDown_nicUp_DEG<-intersections[[8]]$gencodeID


## Define groups of DE exons' genes
nic_up_exons_genes<-unique(de_exons_nic[which(de_exons_nic$logFC>0),"gencodeID"])
smo_up_exons_genes<-unique(de_exons_smo[which(de_exons_smo$logFC>0),"gencodeID"])
nic_down_exons_genes<-unique(de_exons_nic[which(de_exons_nic$logFC<0),"gencodeID"])
smo_down_exons_genes<-unique(de_exons_smo[which(de_exons_smo$logFC<0),"gencodeID"])

smoUp_nicUp_exons_genes<-intersect(nic_up_exons_genes, smo_up_exons_genes)
smoDown_nicDown_exons_genes<-intersect(nic_down_exons_genes, smo_down_exons_genes)
smoUp_nicDown_exons_genes<-intersect(nic_down_exons_genes, smo_up_exons_genes)
smoDown_nicUp_exons_genes<-intersect(nic_up_exons_genes, smo_down_exons_genes)
only_up_nic_exons_genes<-nic_up_exons_genes[which(! (nic_up_exons_genes %in% smo_up_exons_genes | 
                                           nic_up_exons_genes %in% smo_down_exons_genes | 
                                           nic_up_exons_genes %in% nic_down_exons_genes))]
only_up_smo_exons_genes<-smo_up_exons_genes[which(! (smo_up_exons_genes %in% smo_down_exons_genes | 
                                           smo_up_exons_genes %in% nic_down_exons_genes | 
                                           smo_up_exons_genes %in% nic_up_exons_genes))]
only_down_nic_exons_genes<-nic_down_exons_genes[which(! (nic_down_exons_genes %in% smo_up_exons_genes | 
                                               nic_down_exons_genes %in% smo_down_exons_genes | 
                                               nic_down_exons_genes %in% nic_up_exons_genes))]
only_down_smo_exons_genes<-smo_down_exons_genes[which(! (smo_down_exons_genes %in% smo_up_exons_genes | 
                                               smo_down_exons_genes %in% nic_down_exons_genes | 
                                               smo_down_exons_genes %in% nic_up_exons_genes))]





## Compare smoking VS nicotine DE txs

## Smo vs nic DE txs
DE_txs_smo_vs_nic <-list(
  "Nicotine"=de_tx_nic$transcript_id,
  "Smoking"=de_tx_smo$transcript_id)

## Smo vs nic Up DE txs
DE_txs_smo_vs_nic_Up <-list(
  "Nicotine up"=nic_up,
  "Smoking up"=smo_up)

## Smo vs nic Down DE txs
DE_txs_smo_vs_nic_Down <-list(
  "Nicotine down"=nic_down,
  "Smoking down"=smo_down)

## Smo Up vs nic Down DE txs
DE_txs_smoUp_vs_nicDown <-list(
  "Nicotine down"=nic_down,
  "Smoking up"=smo_up)

## Smo Down vs nic Up DE txs
DE_txs_smoDown_vs_nicUp <-list(
  "Nicotine up"=nic_up,
  "Smoking down"=smo_down)

## All DE txs
DE_txs_all<-list(
  "Nicotine up"=nic_up,
  "Nicotine down"=nic_down,
  "Smoking up"=smo_up, 
  "Smoking down"=smo_down)



## Venn diagrams
DE_lists<-list(DE_txs_smo_vs_nic, DE_txs_smo_vs_nic_Up, DE_txs_smo_vs_nic_Down, 
               DE_txs_smoUp_vs_nicDown, DE_txs_smoDown_vs_nicUp, DE_txs_all)
colors<-list(c("olivedrab1", "rosybrown2"), c("brown3", "coral"), c("cyan3", "cyan4"), 
             c("cyan3", "coral"), c("brown3", "cyan4"), c("brown3", "cyan3", "coral", "cyan4"))
venn_plot(DE_lists, colors, "smo_VS_nic_DE_txs", NULL)





## Compare DE txs' genes 

## Smo vs nic genes
DE_txs_genes_smo_vs_nic <-list(
  "Nicotine"=unique(de_tx_nic$ensembl_id),
  "Smoking"=unique(de_tx_smo$ensembl_id))

## Smo vs nic Up genes
DE_txs_genes_smo_vs_nic_Up <-list(
  "Nicotine up"=nic_up_genes,
  "Smoking up"=smo_up_genes)

## Smo vs nic Down genes
DE_txs_genes_smo_vs_nic_Down <-list(
  "Nicotine down"=nic_down_genes,
  "Smoking down"=smo_down_genes)

## Smo Up vs nic Down genes
DE_txs_genes_smoUp_vs_nicDown <-list(
  "Nicotine down"=nic_down_genes,
  "Smoking up"=smo_up_genes)

## Smo Down vs nic Up genes
DE_txs_genes_smoDown_vs_nicUp <-list(
  "Nicotine up"=nic_up_genes,
  "Smoking down"=smo_down_genes)

## Up vs Down genes
DE_txs_genes_Up_vs_Down <-list(
  "Up"=union(nic_up_genes, smo_up_genes),
  "Down"=union(nic_down_genes, smo_down_genes))

## Nic Up vs nic Down genes
DE_txs_genes_nicUp_vs_nicDown <-list(
  "Nicotine up"=nic_up_genes,
  "Nicotine down"=nic_down_genes)

## Smo Up vs smo Down genes
DE_txs_genes_smoUp_vs_smoDown <-list(
  "Smoking up"=smo_up_genes,
  "Smoking down"=smo_down_genes)

## All DE exons
DE_txs_genes_all<-list(
  "Nicotine up"=nic_up_genes,
  "Nicotine down"=nic_down_genes,
  "Smoking up"=smo_up_genes, 
  "Smoking down"=smo_down_genes)



## Venn diagrams
DE_lists<-list(DE_txs_genes_smo_vs_nic, DE_txs_genes_smo_vs_nic_Up, 
               DE_txs_genes_smo_vs_nic_Down, DE_txs_genes_smoUp_vs_nicDown, 
               DE_txs_genes_smoDown_vs_nicUp, DE_txs_genes_Up_vs_Down,
               DE_txs_genes_nicUp_vs_nicDown, DE_txs_genes_smoUp_vs_smoDown, 
               DE_txs_genes_all)
colors<-list(c("olivedrab1", "rosybrown2"), c("brown3", "coral"), 
             c("cyan3", "cyan4"), c("cyan3", "coral"), c("brown3", "cyan4"), 
             c("firebrick", "dodgerblue"),c("brown3", "cyan3"), c("coral", "cyan4"),
             c("brown3", "cyan3", "coral", "cyan4"))
venn_plot(DE_lists, colors, "smo_VS_nic_DE_txs_genes", NULL)





## Compare DE txs' genes with DEG

## All DEG vs all txs' genes
DEG_vs_Txs_all <-list(
  "DEG"= union(union(nic_DEG_down, nic_DEG_up), union(smo_DEG_down, smo_DEG_up)),
  "DE txs' genes"=union(union(nic_down_genes, nic_up_genes), union(smo_down_genes, smo_up_genes))
)

## Nic DEGs vs nic txs' genes 
DEG_nic_vs_Txs_nic <-list(
  "DEG"= union(nic_DEG_down, nic_DEG_up),
  "DE txs' genes"=union(nic_down_genes, nic_up_genes)
)

## Smo DEGs vs smo txs' genes
DEG_smo_vs_Txs_smo <-list(
  "DEG"= union(smo_DEG_down, smo_DEG_up),
  "DE txs' genes"=union(smo_down_genes, smo_up_genes)
)

## Up DEGs vs up txs' genes
DEG_up_vs_Txs_up <-list(
  "DEG"= union(nic_DEG_up, smo_DEG_up),
  "DE txs' genes"=union(nic_up_genes, smo_up_genes)
)

## Down DEGs vs down txs' genes
DEG_down_vs_Txs_down <-list(
  "DEG"= union(nic_DEG_down, smo_DEG_down),
  "DE txs' genes"=union(nic_down_genes, smo_down_genes)
)

## Nic up DEGs vs nic up txs' genes
DEG_nicUp_vs_Txs_nicUp <-list(
  "DEG"= nic_DEG_up,
  "DE txs' genes"= nic_up_genes
)

## Nic down DEGs vs nic down txs' genes
DEG_nicDown_vs_Txs_nicDown <-list(
  "DEG"= nic_DEG_down,
  "DE txs' genes"= nic_down_genes
)

## Smo up DEGs vs smo up txs' genes
DEG_smoUp_vs_Txs_smoUp <-list(
  "DEG"= smo_DEG_up,
  "DE txs' genes"= smo_up_genes
)

## Smo down DEGs vs smo down txs' genes
DEG_smoDown_vs_Txs_smoDown <-list(
  "DEG"= smo_DEG_down,
  "DE txs' genes"= smo_down_genes
)



## Venn diagrams
DE_lists<-list(DEG_vs_Txs_all, DEG_nic_vs_Txs_nic, DEG_smo_vs_Txs_smo,
               DEG_up_vs_Txs_up, DEG_down_vs_Txs_down, DEG_nicUp_vs_Txs_nicUp,
               DEG_nicDown_vs_Txs_nicDown, DEG_smoUp_vs_Txs_smoUp, DEG_smoDown_vs_Txs_smoDown)
colors<-list(c("rosybrown2", "navajowhite2"), c("plum", "lightgoldenrod2"), 
             c("pink1", "khaki2"), c("rosybrown4", "navajowhite4"), c("rosybrown1", "navajowhite1"), 
             c("plum3", "lightgoldenrod4"),c("plum1", "lightgoldenrod1"), c("pink3", "khaki3"),
             c("pink", "khaki1"))
venn_plot(DE_lists, colors, "DEG_VS_txs_genes", c("All", "Nicotine", "Smoking", "Up", "Down", 
                                                    "Nicotine Up", "Nicotine Down", "Smoking Up", "Smoking Down"))





## Compare DEG intersections with those of txs' genes

## Only up nic DEGs vs txs' genes 
DEG_vs_Txs_only_up_nic <-list(
  "DEG"= only_up_nic_DEG,
  "DE txs' genes"= only_up_nic_genes
)

## Only down nic DEGs vs txs' genes 
DEG_vs_Txs_only_down_nic <-list(
  "DEG"= only_down_nic_DEG,
  "DE txs' genes"= only_down_nic_genes
)

## Only up smo DEGs vs txs' genes 
DEG_vs_Txs_only_up_smo <-list(
  "DEG"= only_up_smo_DEG,
  "DE txs' genes"= only_up_smo_genes
)

## Only down smo DEGs vs txs' genes 
DEG_vs_Txs_only_down_smo <-list(
  "DEG"= only_down_smo_DEG,
  "DE txs' genes"= only_down_smo_genes
)

## Smo and nic up DEGs vs txs' genes
DEG_vs_Txs_smoUp_nicUp <-list(
  "DEG"=smoUp_nicUp_DEG,
  "DE txs' genes"=smoUp_nicUp_genes
)

## Smo and nic down DEGs vs txs' genes
DEG_vs_Txs_smoDown_nicDown <-list(
  "DEG"=smoDown_nicDown_DEG,
  "DE txs' genes"=smoDown_nicDown_genes
)

## Smo up and nic down DEGs vs txs' genes
DEG_vs_Txs_smoUp_nicDown <-list(
  "DEG"=smoUp_nicDown_DEG,
  "DE txs' genes"=smoUp_nicDown_genes
)

## Smo down and nic up DEGs vs txs' genes
DEG_vs_Txs_smoDown_nicUp <-list(
  "DEG"=smoDown_nicUp_DEG,
  "DE txs' genes"=smoDown_nicUp_genes
)



## Venn diagrams
DE_lists<-list(DEG_vs_Txs_only_up_nic, DEG_vs_Txs_only_down_nic, DEG_vs_Txs_only_up_smo,
               DEG_vs_Txs_only_down_smo, DEG_vs_Txs_smoUp_nicUp, DEG_vs_Txs_smoDown_nicDown,
               DEG_vs_Txs_smoUp_nicDown, DEG_vs_Txs_smoDown_nicUp)
colors<-list(c("darkslategray3", "lightgoldenrod3"), c("darkslategray1", "lightgoldenrod1"), 
             c("lightpink3", "olivedrab3"), c("pink", "olivedrab1"), c("palegreen3", "navajowhite3"),
             c("palegreen", "navajowhite"), c("palegreen","navajowhite3"), c("palegreen3", "navajowhite"))
venn_plot(DE_lists, colors, "intersections_DEG_VS_txs_genes", c("Only up in nic", "Only down in nic",
                                                                  "Only up in smo", "Only down in smo", 
                                                                  "Smo up, nic up", "Smo down, nic down",
                                                                  "Smo up, nic down", "Smo down, nic up"))





## Compare DE txs' genes vs DE exons' genes vs DEG

## All genes
DEG_vs_Txs_vs_Exons_all <-list(
  "DEG"= union(union(nic_DEG_down, nic_DEG_up), union(smo_DEG_down, smo_DEG_up)),
  "DE txs' genes"=union(union(nic_down_genes, nic_up_genes), union(smo_down_genes, smo_up_genes)),
  "DE exons' genes"=union(union(nic_down_exons_genes, nic_up_exons_genes), union(smo_down_exons_genes, smo_up_exons_genes))
)

## Nic genes
DEG_vs_Txs_vs_Exons_nic <-list(
  "DEG"= union(nic_DEG_down, nic_DEG_up),
  "DE txs' genes"=union(nic_down_genes, nic_up_genes),
  "DE exons' genes"=union(nic_down_exons_genes, nic_up_exons_genes)
)

## Smo genes
DEG_vs_Txs_vs_Exons_smo <-list(
  "DEG"= union(smo_DEG_down, smo_DEG_up),
  "DE txs' genes"=union(smo_down_genes, smo_up_genes),
  "DE exons' genes"=union(smo_down_exons_genes, smo_up_exons_genes)
)

## Up genes
DEG_vs_Txs_vs_Exons_up <-list(
  "DEG"= union(nic_DEG_up, smo_DEG_up),
  "DE txs' genes"=union(nic_up_genes, smo_up_genes),
  "DE exons' genes"=union(nic_up_exons_genes, smo_up_exons_genes)
)

## Down genes
DEG_vs_Txs_vs_Exons_down <-list(
  "DEG"= union(nic_DEG_down, smo_DEG_down),
  "DE txs' genes"=union(nic_down_genes, smo_down_genes),
  "DE exons' genes"=union(nic_down_exons_genes, smo_down_exons_genes)
)

## Nic up genes
DEG_vs_Txs_vs_Exons_nicUp <-list(
  "DEG"= nic_DEG_up,
  "DE txs' genes"= nic_up_genes,
  "DE exons' genes"= nic_up_exons_genes
)

## Nic down genes
DEG_vs_Txs_vs_Exons_nicDown <-list(
  "DEG"= nic_DEG_down,
  "DE txs' genes"= nic_down_genes,
  "DE exons' genes"= nic_down_exons_genes
)

## Smo up genes
DEG_vs_Txs_vs_Exons_smoUp <-list(
  "DEG"= smo_DEG_up,
  "DE txs' genes"= smo_up_genes,
  "DE exons' genes"= smo_up_exons_genes
)

## Smo down genes
DEG_vs_Txs_vs_Exons_smoDown <-list(
  "DEG"= smo_DEG_down,
  "DE txs' genes"= smo_down_genes,
  "DE exons' genes"=smo_down_exons_genes
)



## Venn diagrams
DE_lists<-list(DEG_vs_Txs_vs_Exons_all, DEG_vs_Txs_vs_Exons_nic, DEG_vs_Txs_vs_Exons_smo,
               DEG_vs_Txs_vs_Exons_up, DEG_vs_Txs_vs_Exons_down, DEG_vs_Txs_vs_Exons_nicUp,
               DEG_vs_Txs_vs_Exons_nicDown, DEG_vs_Txs_vs_Exons_smoUp, DEG_vs_Txs_vs_Exons_smoDown)
colors<-list(c("rosybrown2", "navajowhite2", "lightblue3"), c("plum", "lightgoldenrod2", "lightcyan3"), 
             c("pink1", "khaki2", "lightsteelblue3"), c("rosybrown4", "navajowhite4", "lightblue4"), 
             c("rosybrown1", "navajowhite1", "lightblue1"), c("plum3", "lightgoldenrod4", "lightcyan4"),
             c("plum1", "lightgoldenrod1", "lightcyan1"), c("pink3", "khaki3", "lightsteelblue4"),
             c("pink", "khaki1", "lightsteelblue1"))
venn_plot(DE_lists, colors, "DEG_VS_txs_VS_exons", c("All", "Nicotine", "Smoking", "Up", "Down", 
                                                  "Nicotine Up", "Nicotine Down", "Smoking Up", "Smoking Down"))





## Compare intersections of txs' genes vs DE exons' genes vs DEG

## Only up nic genes 
DEG_vs_Txs_vs_Exons_only_up_nic <-list(
  "DEG"= only_up_nic_DEG,
  "DE txs' genes"= only_up_nic_genes,
  "DE exons' genes"= only_up_nic_exons_genes
)

## Only down nic genes 
DEG_vs_Txs_vs_Exons_only_down_nic <-list(
  "DEG"= only_down_nic_DEG,
  "DE txs' genes"= only_down_nic_genes,
  "DE exons' genes"= only_down_nic_exons_genes
)

## Only up smo DEGs vs txs' genes 
DEG_vs_Txs_vs_Exons_only_up_smo <-list(
  "DEG"= only_up_smo_DEG,
  "DE txs' genes"= only_up_smo_genes,
  "DE exons' genes"= only_up_smo_exons_genes
)

## Only down smo genes 
DEG_vs_Txs_vs_Exons_only_down_smo <-list(
  "DEG"= only_down_smo_DEG,
  "DE txs' genes"= only_down_smo_genes,
  "DE exons' genes"= only_down_smo_exons_genes
)

## Smo and nic up genes
DEG_vs_Txs_vs_Exons_smoUp_nicUp <-list(
  "DEG"=smoUp_nicUp_DEG,
  "DE txs' genes"=smoUp_nicUp_genes,
  "DE exons' genes"= smoUp_nicUp_exons_genes
)

## Smo and nic down genes
DEG_vs_Txs_vs_Exons_smoDown_nicDown <-list(
  "DEG"=smoDown_nicDown_DEG,
  "DE txs' genes"=smoDown_nicDown_genes,
  "DE exons' genes"= smoDown_nicDown_exons_genes
)

## Smo up and nic down genes
DEG_vs_Txs_vs_Exons_smoUp_nicDown <-list(
  "DEG"=smoUp_nicDown_DEG,
  "DE txs' genes"=smoUp_nicDown_genes,
  "DE exons' genes"= smoUp_nicDown_exons_genes
)

## Smo down and nic up genes
DEG_vs_Txs_vs_Exons_smoDown_nicUp <-list(
  "DEG"=smoDown_nicUp_DEG,
  "DE txs' genes"=smoDown_nicUp_genes,
  "DE exons' genes"= smoDown_nicUp_exons_genes
)



## Venn diagrams
DE_lists<-list(DEG_vs_Txs_vs_Exons_only_up_nic, DEG_vs_Txs_vs_Exons_only_down_nic, DEG_vs_Txs_vs_Exons_only_up_smo,
               DEG_vs_Txs_vs_Exons_only_down_smo, DEG_vs_Txs_vs_Exons_smoUp_nicUp, DEG_vs_Txs_vs_Exons_smoDown_nicDown,
               DEG_vs_Txs_vs_Exons_smoUp_nicDown, DEG_vs_Txs_vs_Exons_smoDown_nicUp)
colors<-list(c("darkslategray3", "lightgoldenrod3", "orangered"), c("darkslategray1", "lightgoldenrod1", "sienna1"), 
             c("lightpink3", "olivedrab3", "indianred"), c("pink", "olivedrab1", "indianred1"), 
             c("palegreen3", "navajowhite3", "orange3"), c("palegreen", "navajowhite", "orange1"), 
             c("palegreen2","navajowhite2", "orange2"), c("palegreen2", "navajowhite2", "orange2"))
venn_plot(DE_lists, colors, "intersections_DEG_VS_txs_VS_exons", c("Only up in nic", "Only down in nic",
                                                                "Only up in smo", "Only down in smo", 
                                                                "Smo up, nic up", "Smo down, nic down",
                                                                "Smo up, nic down", "Smo down, nic up"))







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
# date     2023-01-11
# rstudio  2022.07.2+576 Spotted Wakerobin (desktop)
# pandoc   NA
