
## 1.2 Comparison of DE transcripts


load(here("processed-data/04_DEA/Tx_analysis/top_tx_nic.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/de_tx_nic.Rdata"))
#load(here("processed-data/04_DEA/Tx_analysis/results_nic.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/top_tx_smo.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/de_tx_smo.Rdata"))
#load(here("processed-data/04_DEA/Tx_analysis/results_smo.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/rse_tx_brain_pups_nicotine.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/rse_tx_brain_pups_smoking.Rdata"))

load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/results_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/results_pups_smoking_fitted.Rdata"))
#load(here("processed-data/05_GO_KEGG/Gene_analysis/intersections.Rdata"))     

# load(here("processed-data/04_DEA/Exon_analysis/top_exons_nic.Rdata"))
# load(here("processed-data/04_DEA/Exon_analysis/de_exons_nic.Rdata"))
# load(here("processed-data/04_DEA/Exon_analysis/top_exons_smo.Rdata"))
# load(here("processed-data/04_DEA/Exon_analysis/de_exons_smo.Rdata"))


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
  
  ## Merge data
  t_stats<-data.frame(t1=top_tx1$t, t2=top_tx2$t)
  ## Add DE info for both groups 
  t_stats$DE<-add_DE_info(top_tx1)
  
  cols <- c("red", "#ffad73","#26b3ff", "dark grey") 
  names(cols)<-c("sig Both", "sig nic", "sig smo", "None")
  alphas <- c( 1, 1, 1,0.5)  
  names(alphas)<-c("sig Both", "sig nic", "sig smo", "None")
  
  plot <- ggplot(t_stats, aes(x = t1, y = t2, color=DE, alpha=DE)) +
    geom_point(size = 1) +
    labs(x = paste("t-stats", name_1), 
         y = paste("t-stats", name_2),
         title = title, 
         subtitle = rho_anno, 
         parse = T) +
    theme_bw() +
    scale_color_manual(values = cols) + 
    scale_alpha_manual(values = alphas)
  
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
      DE<-append(DE, "sig tx")
    }
    
    ## DE genes only
    else if(t_stats$adj.P.Val_tx[i]>=0.05 && t_stats$adj.P.Val_genes[i]<0.05){
      DE<-append(DE, "sig gene")
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
    both_t=c(3,-4)
    FDR=1
  }
  else{
    abs_t_tx=7
    both_t=c(4,-5)
    FDR=0.001
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
  for (i in 1:dim(t_stats)[1]){
    t<-top_genes[which(top_genes$gencodeID==t_stats$ensembl_id[i]), "t"]
    FDR<-top_genes[which(top_genes$gencodeID==t_stats$ensembl_id[i]), "adj.P.Val"]
    t_genes<-append(t_genes, t)
    FDRs<-append(FDRs, FDR)
  }
  t_stats$t_genes<-t_genes
  t_stats$adj.P.Val_genes<-FDRs
  
  ## Correlation coeff between t-stats of genes and transcripts
  rho <- cor(t_stats$t_tx, t_stats$t_genes, method = "spearman")
  rho_anno = paste0("rho = ", format(round(rho, 2), nsmall = 2))
  
  ## Add DE info for both groups 
  t_stats$DE<-add_DE_info_tx_vs_genes(t_stats)
  
  
  ## Gene-tx symbols of DE tx with no DEG, 
  ## DE tx whose DEG have an opposite sign in logFC
  ## and also label DE tx from genes with up and down tx
  
  ## Up and down transcripts' genes
  tx_up_genes<-unique(de_tx[which(de_tx$logFC>0),"ensembl_id"])
  tx_down_genes<-unique(de_tx[which(de_tx$logFC<0),"ensembl_id"])
  ## Genes with up and down tx
  interest_genes<-intersect(tx_up_genes, tx_down_genes)
  ## Retain only the transcripts' genes that were considered at the gene level
  interest_genes<-intersect(interest_genes, tx_genes)
  
  tx_symbols<-vector()
  for (i in 1:dim(t_stats)[1]) {
    
    if (t_stats$DE[i]=="sig tx" & (abs(t_stats$t_tx[i])>abs_t_tx | 
       ((t_stats$ensembl_id[i] %in% interest_genes) & (t_stats$adj.P.Val_tx[i]< FDR)))) {
      tx_symbols<-append(tx_symbols, paste(t_stats$Symbol[i], "-", t_stats$transcript_id[i], sep=""))
    }
    else if(t_stats$DE[i]=="sig Both"){
      if ((t_stats$t_genes[i]> both_t[1] & t_stats$t_tx[i]< both_t[2]) | 
          (t_stats$t_genes[i]< -both_t[1] & t_stats$t_tx[i]> -both_t[2])){
        tx_symbols<-append(tx_symbols, paste(t_stats$Symbol[i], "-", t_stats$transcript_id[i], sep=""))
      }
      else if((t_stats$ensembl_id[i] %in% interest_genes) & (t_stats$adj.P.Val_tx[i]< FDR)){
        tx_symbols<-append(tx_symbols, paste(t_stats$Symbol[i], "-", t_stats$transcript_id[i], sep=""))
      }      
      else{
        tx_symbols<-append(tx_symbols, NA)
      }
    }
    else {
      tx_symbols<-append(tx_symbols, NA)
    }
  }
  t_stats$tx_symbols<-tx_symbols
  
  ## Plot
  cols <- c("red", "#ffad73","#26b3ff", "dark grey") 
  names(cols)<-c("sig Both","sig tx", "sig gene", "None")
  alphas <- c(1, 1, 1,0.5)  
  names(alphas)<-c("sig Both", "sig tx", "sig gene", "None")
  
  plot <- ggplot(t_stats, aes(x = t_genes, y = t_tx, color=DE, alpha=DE, label= tx_symbols)) +
    geom_point(size = 1) +
    labs(x = "t-stats genes", 
         y = "t-stats tx",
         title = paste(capitalize(expt),"genes vs tx", sep=" "), 
         subtitle = rho_anno, 
         parse = T) +
    geom_label_repel(fill="white", size=2, max.overlaps = Inf,  
                     box.padding = 0.2, 
                     show.legend=FALSE) +
    theme_bw() +
    scale_color_manual(values = cols) + 
    scale_alpha_manual(values = alphas)
  
  plot
  ggsave(filename=paste("plots/04_DEA/02_Comparisons/Tx_analysis/t_stats_tx_vs_genes_", substr(expt,1,3), 
                        ".pdf", sep=""), height = 20, width = 25, units = "cm")
  
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





### 1.2.1.3 T-stats of exons vs transcripts

## Function to add DE info of tx and exons

add_DE_info_tx_vs_exons <-function(t_stats) {
  
  DE<-vector()
  for (i in 1:dim(t_stats)[1]) {
    
    ## DE tx with DE exons
    if (t_stats$adj.P.Val_tx[i]<0.05 && t_stats$adj.P.Val_genes[i]<0.05){
      DE<-append(DE, "sig Both")
    }
    
    ## DE tx only
    else if (t_stats$adj.P.Val_tx[i]<0.05 && t_stats$adj.P.Val_genes[i]>=0.05){
      DE<-append(DE, "sig tx")
    }
    
    ## DE genes only
    else if(t_stats$adj.P.Val_tx[i]>=0.05 && t_stats$adj.P.Val_genes[i]<0.05){
      DE<-append(DE, "sig gene")
    }
    
    ## Neither DE tx nor DE genes
    else {
      DE<-append(DE, "None")      
    }
  }
  return(DE)
}



## Create plots of t-stats of tx vs exons

t_stat_tx_vs_exons<- function(expt){
  
  top_tx<-eval(parse_expr(paste("top_tx_", substr(expt,1,3), sep="")))
  top_exons<-eval(parse_expr(paste("top_exons_", substr(expt,1,3), sep="")))
  de_tx<-eval(parse_expr(paste("de_tx_", substr(expt,1,3), sep="")))
  de_exons<-eval(parse_expr(paste("de_exons_", substr(expt,1,3), sep="")))
  
  if (expt=="nicotine"){
    abs_t_tx=6
    both_t=c(3,-4)
    FDR=1
  }
  else{
    abs_t_tx=7
    both_t=c(4,-5)
    FDR=0.001
  }
  
  ## Transcripts' 
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
  for (i in 1:dim(t_stats)[1]){
    t<-top_genes[which(top_genes$gencodeID==t_stats$ensembl_id[i]), "t"]
    FDR<-top_genes[which(top_genes$gencodeID==t_stats$ensembl_id[i]), "adj.P.Val"]
    t_genes<-append(t_genes, t)
    FDRs<-append(FDRs, FDR)
  }
  t_stats$t_genes<-t_genes
  t_stats$adj.P.Val_genes<-FDRs
  
  ## Correlation coeff between t-stats of genes and transcripts
  rho <- cor(t_stats$t_tx, t_stats$t_genes, method = "spearman")
  rho_anno = paste0("rho = ", format(round(rho, 2), nsmall = 2))
  
  ## Add DE info for both groups 
  t_stats$DE<-add_DE_info_tx_vs_genes(t_stats)
  
  
  ## Gene-tx symbols of DE tx with no DEG, 
  ## DE tx whose DEG have an opposite sign in logFC
  ## and also label DE tx from genes with up and down tx
  
  ## Up and down transcripts' genes
  tx_up_genes<-unique(de_tx[which(de_tx$logFC>0),"ensembl_id"])
  tx_down_genes<-unique(de_tx[which(de_tx$logFC<0),"ensembl_id"])
  ## Genes with up and down tx
  interest_genes<-intersect(tx_up_genes, tx_down_genes)
  ## Retain only the transcripts' genes that were considered at the gene level
  interest_genes<-intersect(interest_genes, tx_genes)
  
  tx_symbols<-vector()
  for (i in 1:dim(t_stats)[1]) {
    
    if (t_stats$DE[i]=="sig tx" & (abs(t_stats$t_tx[i])>abs_t_tx | 
                                   ((t_stats$ensembl_id[i] %in% interest_genes) & (t_stats$adj.P.Val_tx[i]< FDR)))) {
      tx_symbols<-append(tx_symbols, paste(t_stats$Symbol[i], "-", t_stats$transcript_id[i], sep=""))
    }
    else if(t_stats$DE[i]=="sig Both"){
      if ((t_stats$t_genes[i]> both_t[1] & t_stats$t_tx[i]< both_t[2]) | 
          (t_stats$t_genes[i]< -both_t[1] & t_stats$t_tx[i]> -both_t[2])){
        tx_symbols<-append(tx_symbols, paste(t_stats$Symbol[i], "-", t_stats$transcript_id[i], sep=""))
      }
      else if((t_stats$ensembl_id[i] %in% interest_genes) & (t_stats$adj.P.Val_tx[i]< FDR)){
        tx_symbols<-append(tx_symbols, paste(t_stats$Symbol[i], "-", t_stats$transcript_id[i], sep=""))
      }      
      else{
        tx_symbols<-append(tx_symbols, NA)
      }
    }
    else {
      tx_symbols<-append(tx_symbols, NA)
    }
  }
  t_stats$tx_symbols<-tx_symbols
  
  ## Plot
  cols <- c("red", "#ffad73","#26b3ff", "dark grey") 
  names(cols)<-c("sig Both","sig tx", "sig gene", "None")
  alphas <- c(1, 1, 1,0.5)  
  names(alphas)<-c("sig Both", "sig tx", "sig gene", "None")
  
  plot <- ggplot(t_stats, aes(x = t_genes, y = t_tx, color=DE, alpha=DE, label= tx_symbols)) +
    geom_point(size = 1) +
    labs(x = "t-stats genes", 
         y = "t-stats tx",
         title = paste(capitalize(expt),"genes vs tx", sep=" "), 
         subtitle = rho_anno, 
         parse = T) +
    geom_label_repel(fill="white", size=2, max.overlaps = Inf,  
                     box.padding = 0.2, 
                     show.legend=FALSE) +
    theme_bw() +
    scale_color_manual(values = cols) + 
    scale_alpha_manual(values = alphas)
  
  plot
  ggsave(filename=paste("plots/04_DEA/02_Comparisons/Tx_analysis/t_stats_tx_vs_genes_", substr(expt,1,3), 
                        ".pdf", sep=""), height = 20, width = 25, units = "cm")
  
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
    geom_boxplot(outlier.color = "#FFFFFFFF") +
    geom_jitter(aes(colour=Group),shape=16, 
                position=position_jitter(0.2)) +
    theme_classic() +
    labs(x = "Group", y = "norm counts",
         title = title,
         ## Add FDR and FC of genes and txs
         subtitle=paste(" FDR:", q_value, "       ", tx_tpm,  
                        "\n", "FC:", FC)) +
    theme(plot.margin=unit (c (1,1.5,1,1), 'cm'), legend.position = "none",
          plot.title = element_text(hjust=0.5, size=10, face="bold"), 
          plot.subtitle = element_text(size = 9)) 
  
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
  
  ## Regress out residuals
  formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
  model=model.matrix(formula, data=vGene$targets)
  vGene$E<-cleaningY(vGene$E, model, P=2)
  
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
  gene_txs_tpm<-sum(apply(gene_txs_tpm, 1, sum))
  
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
  
  ## Log-expression values of that txs
  ## Regress out residuals
  formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
  model=model.matrix(formula, data=colData(RSE))
  logcounts<-cleaningY(assays(RSE)$logcounts, model, P=2)
  
  counts_txs<-logcounts[which(rownames(logcounts) %in% top_gene_txs[1:3]),]
  counts_txs<-t(counts_txs)
  
  ## Data frame with all expression values and samples' group
  counts<-cbind(counts_gene, counts_txs, "Group"=vGene$targets$Group)
  counts<-as.data.frame(counts)
  counts$counts_gene<-as.numeric(counts$counts_gene)
  colnames(counts)[1]<-gene
  counts[,2]<-as.numeric(counts[,2])
  counts[,3]<-as.numeric(counts[,3])
  counts[,4]<-as.numeric(counts[,4])
  
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
  
  plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol = 4)
  ggsave(here(paste("plots/04_DEA/02_Comparisons/Tx_analysis/Boxplots_", gene_symbol, ".pdf", sep="")), 
         width = 35, height = 10, units = "cm")
  
}



#########################################
## Boxplots of nicotine genes and txs
#########################################

#### Non-DE genes but with DE txs ####

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



#### Non-DE genes with Up and Down DE tx ####

## Phf3: predicted to be involved in transcription, DNA-templated;
## expressed in CNS
gene_tx_boxplots("nicotine", "Phf3", "Phf3−ENSMUST00000185521.1", 
                 "Phf3−ENSMUST00000191329.1")



#########################################
## Boxplots of smoking genes and txs
#########################################

#### Non-DE genes but with DE txs ####

## Btf3: acts upstream of or within in utero embryonic development; expressed in brain
gene_tx_boxplots("smoking", "Btf3", "Btf3−ENSMUST00000022163.14", NULL)

## Srsf6: critical for mRNA splicing, involved in mRNA export from the nucleus; 
## expressed in CNS
gene_tx_boxplots("smoking", "Srsf6", "Srsf6−ENSMUST00000017065.14", NULL)



#### DE tx whose DEG have opposite direction of regulation ####

## Meaf6: predicted to act upstream of or within chromatin organization;
## expressed in CNS and whole brain
gene_tx_boxplots("smoking", "Meaf6", "Meaf6−ENSMUST00000184205.7", NULL)



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



#### Non-DE genes with Up and Down DE tx ####

## Cyhr1: predicted to enable zinc ion binding activity; located in cytoplasm and nuclear envelope
gene_tx_boxplots("smoking", "Cyhr1", "Cyhr1−ENSMUST00000176274.1", 
                 "Cyhr1−ENSMUST00000081291.12")

## H13: acts upstream of or within in utero embryonic development; expressed in  brain, 
## embryo ectoderm, extraembryonic component, hemolymphoid system and intervertebral disc
gene_tx_boxplots("smoking", "H13", "H13−ENSMUST00000089059.8", 
                 "H13−ENSMUST00000125366.7")







### 1.2.2 Venn diagrams

## Function to create multiple Venn diagrams
venn_plot<-function(DE_lists, colors, name, titles){
  
  if (name=="smo_VS_nic_DE_exons"){
    height=10
    width=15
    margin=0.2
    dist=0.1
    cat=0.75
    cex=0.8
  }
  
  else if (name=="DEG_VS_exons_genes"){
    height=12
    width=16
    margin=5
    dist=0.07
    cat=0.75
    cex=0.8
  }
  
  else if (name=="smo_VS_nic_DE_exons_genes"){
    height=12
    width=16
    margin=0.2
    dist=0.09
    cat=0.75
    cex=0.8
  }
  
  else if (name=="intersections_DEG_VS_exons_genes"){
    height=10
    width=20
    margin=5
    dist=0.06
    cat=0.9
    cex=1
  }
  
  plots<-list()
  pdf(file = paste("plots/04_DEA/02_Comparisons/Exon_analysis/Venn_", name, ".pdf", sep=""), 
      height = height, width = width)
  for (i in 1:length(DE_lists)){
    v<-venn.diagram(DE_lists[[i]], fill=colors[[i]], alpha = rep(0.5, length(DE_lists[[i]])), 
                    lwd =0, margin=margin, cat.cex=cat, cex=cex, height = 10, width = 12.5, units = "cm", 
                    cat.dist=rep(dist, length(DE_lists[[i]])), filename=NULL, main = titles[i], 
                    main.pos = c(0.5,0.476), disable.logging=TRUE)
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

