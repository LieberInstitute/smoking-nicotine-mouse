
## 1.2 Comparison of DE transcripts


load(here("processed-data/04_DEA/Tx_analysis/top_tx_nic.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/de_tx_nic.Rdata"))
#load(here("processed-data/04_DEA/Tx_analysis/results_nic.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/top_tx_smo.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/de_tx_smo.Rdata"))
#load(here("processed-data/04_DEA/Tx_analysis/results_smo.Rdata"))

load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata"))
#load(here("processed-data/04_DEA/Gene_analysis/results_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata"))
#load(here("processed-data/04_DEA/Gene_analysis/results_pups_smoking_fitted.Rdata"))
#load(here("processed-data/05_GO_KEGG/Gene_analysis/intersections.Rdata"))     



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







