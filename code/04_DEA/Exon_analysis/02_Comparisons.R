
## 1.2 Comparison of DE exons


load(here("processed-data/04_DEA/Exon_analysis/top_exons_nic.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/de_exons_nic.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/results_nic.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/top_exons_smo.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/de_exons_smo.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/results_smo.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/results_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/results_pups_smoking_fitted.Rdata"))
load(here("processed-data/05_GO_KEGG/Gene_analysis/intersections.Rdata"))     



### 1.2.1 T-stats plots

### 1.2.1.1 T-stats of exons in nic vs smo

## Function to add DE info of exons in both groups of samples

add_DE_info <-function(top_exons1, top_exons2) {
  
  DE<-vector()
  for (i in 1:dim(top_exons1)[1]) {
    
    if(top_exons1$exon_libdID[i] %in% de_exons_nic$exon_libdID) {
      ## DE exons in both groups
      if(top_exons2$exon_libdID[i] %in% de_exons_smo$exon_libdID){
        DE<-append(DE, "sig Both")
      }
      ## DE exons in nic only
      else{
        DE<-append(DE, "sig nic")
      }
    }
    
      else if (top_exons2$exon_libdID[i] %in% de_exons_smo$exon_libdID){
      if(!top_exons1$exon_libdID[i] %in% de_exons_nic$exon_libdID){
        ## DE exons in smo only
        DE<-append(DE, "sig smo")
      }
    }
    else {
      ## No DE exons in any group
      DE<-append(DE, "None")
    }
  }
  return(DE)
}



## Compare t-stats of exons from different groups of samples
t_stat_plot <- function(top_exons1, top_exons2, name_1, name_2, title){
  
  ## Correlation coeff
  rho <- cor(top_exons1$t, top_exons2$t, method = "spearman")
  rho_anno = paste0("rho = ", format(round(rho, 2), nsmall = 2))
  
  ## Merge data
  t_stats<-data.frame(t1=top_exons1$t, t2=top_exons2$t)
  ## Add DE info for both groups
  t_stats$DE<-add_DE_info(top_exons1, top_exons2)
  
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
  ggsave(filename=paste("plots/04_DEA/02_Comparisons/Exon_analysis/t_stats_", gsub(" ", "_", title), 
                        ".pdf", sep=""), height = 20, width = 25, units = "cm")
}


################################
# Smoking vs nicotine (exons)
################################

t_stat_plot(top_exons_nic, top_exons_smo, "Nicotine pups", "Smoking pups", "Nic vs Smo DE exons")





### 1.2.1.2 T-stats of genes vs exons

## Function to add DE info of exons and genes

add_DE_info_exons_vs_genes <-function(t_stats) {
  
  DE<-vector()
  for (i in 1:dim(t_stats)[1]) {

    ## DE exons in both genes and exons
    if (t_stats$adj.P.Val_exons[i]<0.05 && t_stats$adj.P.Val_genes[i]<0.05){
      
      if(abs(t_stats$logFC_exons[i])>0.25) {
        DE<-append(DE, "sig Both")
      }
      ## DE genes only
      else {
        DE<-append(DE, "sig gene")
      } 
    }
    
    else if(t_stats$adj.P.Val_exons[i]<0.05 && t_stats$adj.P.Val_genes[i]>=0.05){
      
      ## DE exons only
      if(abs(t_stats$logFC_exons[i])>0.25) {
        DE<-append(DE, "sig exon")
      }
      else{
        DE<-append(DE, "None")
      }
    }
    
    else if(t_stats$adj.P.Val_exons[i]>=0.05 && t_stats$adj.P.Val_genes[i]<0.05){
      ## DE genes only
      DE<-append(DE, "sig gene")
    }
      
    else {
      ## Neither DE exons nor DE genes
      DE<-append(DE, "None")      
    }
  }
  return(DE)
}



## Create plots of t-stats of exons vs genes

t_stat_exons_vs_genes<- function(expt){
  
  top_exons<-eval(parse_expr(paste("top_exons_", substr(expt,1,3), sep="")))
  top_genes<-eval(parse_expr(paste("top_genes_pups_", expt, "_fitted", sep="")))
  
  ## Exons' genes
  exons_genes<-unique(top_exons$ensemblID)
  
  ## Common genes
  exons_genes<-exons_genes[which(exons_genes %in% top_genes$ensemblID)]
  
  ## Extract exons' info
  t_stats<-top_exons[which(top_exons$ensemblID %in% exons_genes), c("exon_libdID", "adj.P.Val", "Symbol", 
                                                                    "t", "ensemblID", "logFC", "seqnames", 
                                                                    "start", "end")]
  colnames(t_stats)[2]<-"adj.P.Val_exons"
  colnames(t_stats)[4]<-"t_exons"
  colnames(t_stats)[6]<-"logFC_exons"
  
  ## Add t-stats and FDRs of exons' genes
  t_genes<-vector()
  FDRs<-vector()
  for (i in 1:dim(t_stats)[1]){
    t<-top_genes[which(top_genes$ensemblID==t_stats$ensemblID[i]), "t"]
    FDR<-top_genes[which(top_genes$ensemblID==t_stats$ensemblID[i]), "adj.P.Val"]
    t_genes<-append(t_genes, t)
    FDRs<-append(FDRs, FDR)
  }
  t_stats$t_genes<-t_genes
  t_stats$adj.P.Val_genes<-FDRs
  
  ## Correlation coeff between t-stats of genes and exons
  rho <- cor(t_stats$t_exons, t_stats$t_genes, method = "spearman")
  rho_anno = paste0("rho = ", format(round(rho, 2), nsmall = 2))
  
  ## Add DE info for both groups
  t_stats$DE<-add_DE_info_exons_vs_genes(t_stats)
  
  ## Gene-exon symbols of DE exons with no DEG 
  ## and for DE exons whose DEG have an opposite sign in logFC
  exon_symbols<-vector()
  for (i in 1:dim(t_stats)[1]) {
    if (t_stats$DE[i]=="sig exon" & abs(t_stats$t_exons[i])>6) {
      exon_symbols<-append(exon_symbols, paste(t_stats$Symbol[i], "-", t_stats$seqname[i], ":", t_stats$start[i], "-",
                                               t_stats$end[i], sep=""))
    }
    else if(t_stats$DE[i]=="sig Both"){
      if ((t_stats$t_genes[i]< -4 & t_stats$t_exons[i]> 3) | (t_stats$t_genes[i]> 3 & t_stats$t_exons[i]< -4)){
        exon_symbols<-append(exon_symbols, paste(t_stats$Symbol[i], "-", t_stats$seqname[i], ":", t_stats$start[i], "-",
                                                 t_stats$end[i], sep=""))
      }
      else{
        exon_symbols<-append(exon_symbols, NA)
      }

    }
    else {
      exon_symbols<-append(exon_symbols, NA)
    }
  }
  t_stats$exon_symbols<- exon_symbols
  
  ## Plot
  cols <- c("red", "#ffad73","#26b3ff", "dark grey") 
  names(cols)<-c("sig Both","sig exon", "sig gene", "None")
  alphas <- c( 1, 1, 1,0.5)  
  names(alphas)<-c("sig Both", "sig exon", "sig gene", "None")
  
  plot <- ggplot(t_stats, aes(x = t_genes, y = t_exons, color=DE, alpha=DE, label= exon_symbols)) +
    geom_point(size = 1) +
    labs(x = "t-stats genes", 
         y = "t-stats exons",
         title = paste(capitalize(expt),"genes vs exons", sep=" "), 
         subtitle = rho_anno, 
         parse = T) +
    geom_label_repel(fill="white", size=2, max.overlaps = Inf,  
                     box.padding = 0.2, 
                     show.legend=FALSE) +
    theme_bw() +
    scale_color_manual(values = cols) + 
    scale_alpha_manual(values = alphas)
  
  plot
  ggsave(filename=paste("plots/04_DEA/02_Comparisons/Exon_analysis/t_stats_global_", substr(expt,1,3), 
                        ".pdf", sep=""), height = 20, width = 25, units = "cm")

}


#################################### 
# Nicotine genes vs nicotine exons
#####################################
expt<-"nicotine"
t_stat_exons_vs_genes(expt)


#################################### 
# Smoking genes vs smoking exons
#####################################
expt<-"smoking"
t_stat_exons_vs_genes(expt)





### 1.2.1.3 T-stats of genes vs mean of |t-stat of gene - t-stats of gene's exons|

t_stat_tg<-function(expt){

  top_exons<-eval(parse_expr(paste("top_exons_", substr(expt,1,3), sep="")))
  top_genes<-eval(parse_expr(paste("top_genes_pups_", expt, "_fitted", sep="")))
  
  ## Exons' genes
  exons_genes<-unique(top_exons$ensemblID)
  ## Common genes (genes with exons)
  exons_genes<-exons_genes[which(exons_genes %in% top_genes$ensemblID)]
  
  ## Extract t-stats of each gene and its exons
  tg_te<-data.frame("ensemblID", "tg", "t")
  for (i in 1:length(exons_genes)){
    ## t-stats of the gene
    tg<-top_genes[which(top_genes$ensemblID==exons_genes[i]), "t"]
    ## t-stats of the gene's exons
    t_exons<-top_exons[which(top_exons$ensemblID==exons_genes[i]),  "t"]
    ## mean(|tg-te|)
    t<-mean(abs(tg-t_exons))
    tg_te[i,]<-c(exons_genes[i], tg, t)
  }
  colnames(tg_te)<-c("ensemblID", "tg", "t")
  tg_te$tg<-as.numeric(tg_te$tg)
  tg_te$t<-as.numeric(tg_te$t)
  
  ## Add DE info of genes
  DE<-vector()
  for (gene in tg_te$ensemblID){
    FDR=top_genes[which(top_genes$ensemblID==gene),"adj.P.Val"]
    if(FDR<0.05){
      DE<-append(DE, "sig gene")
    }
    else{
      DE<-append(DE, "ns")
    }
  }
  tg_te$DE<-DE
  
  ## Add gene symbols of DEG with mean(|tg-te|)>4
  gene_symbols<-vector()
  for (i in 1:length(tg_te$ensemblID)){
    if(tg_te$DE[i]=="sig gene" && tg_te$t[i]>4){
      symbol<-top_genes[which(top_genes$ensemblID==tg_te$ensemblID[i]), "Symbol"]
      gene_symbols<-append(gene_symbols, paste(symbol, "-", tg_te$ensemblID[i]))
    }
    else{
      gene_symbols<-append(gene_symbols, NA)
    }
  }
  tg_te$gene_symbols<-gene_symbols
  
  ## Plot
  cols <- c("red", "dark grey") 
  names(cols)<-c("sig gene", "ns")
  alphas <- c( 1,0.5)  
  names(alphas)<-c("sig gene", "ns")
  
  plot <- ggplot(tg_te, aes(x =tg, y = t, color=DE, alpha=DE, label=gene_symbols)) +
    geom_point(size = 1) +
    labs(x = "t-stats of gene", 
         y = "mean of |tg-te|",
         title = paste(capitalize(expt),"genes", sep=" ")) +
    geom_label_repel(fill="white", size=2, max.overlaps = Inf,
                     box.padding = 0.2,
                     show.legend=FALSE) +
    theme_bw() +
    scale_color_manual(values = cols) + 
    scale_alpha_manual(values = alphas)
  
  plot
  ggsave(filename=paste("plots/04_DEA/02_Comparisons/Exon_analysis/t_stats_tg_", substr(expt,1,3), 
                        ".pdf", sep=""), height = 20, width = 25, units = "cm")
}


##################
# Nicotine genes 
##################
expt<-"nicotine"
t_stat_tg(expt)


###################
# Smoking genes 
###################
expt<-"smoking"
t_stat_tg(expt)





### 1.2.1.4 T-stats of exons vs |t-stat of exon's gene - t-stat of exon|

t_stat_te<- function(expt){
  top_exons<-eval(parse_expr(paste("top_exons_", substr(expt,1,3), sep="")))
  top_genes<-eval(parse_expr(paste("top_genes_pups_", expt, "_fitted", sep="")))
  
  ## Exons' genes
  exons_genes<-unique(top_exons$ensemblID)
  exons_genes<-exons_genes[which(exons_genes %in% top_genes$ensemblID)]
  
  ## Extract exons' info
  t_stats<-top_exons[which(top_exons$ensemblID %in% exons_genes), c("exon_libdID", "adj.P.Val", "Symbol", 
                                                                    "t", "ensemblID", "logFC", "seqnames", 
                                                                    "start", "end")]
  colnames(t_stats)[2]<-"adj.P.Val_exons"
  colnames(t_stats)[4]<-"te"
  colnames(t_stats)[6]<-"logFC_exons"
  
  ## Add t-stats of exons' genes
  t_genes<-vector()
  for (i in 1:dim(t_stats)[1]){
    tg<-top_genes[which(top_genes$ensemblID==t_stats$ensemblID[i]), "t"]
    t_genes<-append(t_genes, tg)
  }
  t_stats$tg<-t_genes
  
  ## |te-tg| of each exon
  t_stats$t<-abs(t_stats$te - t_stats$tg)
  
  ## Add DE info of exons
  DE<-vector()
  for (i in 1:dim(t_stats)[1]){
    if(t_stats$adj.P.Val_exons[i]<0.05 && abs(t_stats$logFC_exons[i])>0.25){
      DE<-append(DE, "sig exon")
    }
    else{
      DE<-append(DE, "ns")
    }
  }
  t_stats$DE<-DE
  
  ## Add symbols of DE exons with |te-tg|>8.5
  exon_symbols<-vector()
  for (i in 1:length(t_stats$ensemblID)){
    if(t_stats$DE[i]=="sig exon" && t_stats$t[i]>8.5){
      symbol<-top_genes[which(top_genes$ensemblID==t_stats$ensemblID[i]), "Symbol"]
      exon_symbols<-append(exon_symbols, paste(symbol, "-", t_stats$seqnames[i], ":", 
                                               t_stats$start[i], "-", t_stats$end[i], sep=""))
    }
    else{
      exon_symbols<-append(exon_symbols, NA)
    }
  }
  t_stats$exon_symbols<-exon_symbols
  
  ## Plot
  cols <- c("red", "dark grey") 
  names(cols)<-c("sig exon", "ns")
  alphas <- c( 1,0.5)  
  names(alphas)<-c("sig exon", "ns")
  
  plot <- ggplot(t_stats, aes(x =te, y = t, color=DE, alpha=DE, label=exon_symbols)) +
    geom_point(size = 1) +
    labs(x = "t-stats of exon", 
         y = "|te-tg|",
         title = paste(capitalize(expt),"exons", sep=" ")) +
    geom_label_repel(fill="white", size=2, max.overlaps = Inf,
                     box.padding = 0.2,
                     show.legend=FALSE) +
    theme_bw() +
    scale_color_manual(values = cols) + 
    scale_alpha_manual(values = alphas)
  
  plot
  ggsave(filename=paste("plots/04_DEA/02_Comparisons/Exon_analysis/t_stats_te_", substr(expt,1,3), 
                        ".pdf", sep=""), height = 20, width = 25, units = "cm")
}


##################
# Nicotine exons
##################
expt<-"nicotine"
t_stat_te(expt)


###################
# Smoking exons 
###################
expt<-"smoking"
t_stat_te(expt)





## 1.2.2 Boxplots of relevant DEG and their exons

## Each boxplot
create_boxplot<- function(counts, y, title, q_value){
  
  ## Boxplot
  p<-ggplot(data=counts, 
            aes(x=Group,y=eval(parse_expr(y)))) + 
    geom_boxplot(outlier.color = "#FFFFFFFF") +
    geom_jitter(aes(colour=Group),shape=16, 
                position=position_jitter(0.2)) +
    theme_classic() +
    labs(x = "Group", y = "norm counts",
         title = title,
         subtitle = paste("FDR:", q_value)) +
    theme(plot.margin=unit (c (1,1.5,1,1), 'cm'), legend.position = "none",
          plot.title = element_text(hjust=0.5, size=10, face="bold"), 
          plot.subtitle = element_text(size = 9)) 
  print(p)
}



## Boxplots
gene_exons_boxplots<- function(expt, gene, exon){

  top_exons<-eval(parse_expr(paste("top_exons_", substr(expt,1,3), sep="")))
  top_genes<-eval(parse_expr(paste("top_genes_pups_", expt, "_fitted", sep="")))
  results_genes<-eval(parse_expr(paste("results_pups_", expt, "_fitted", sep="")))
  
  ## Expression values of all genes
  vGene<-results_genes[[1]][[2]]
  rownames(vGene)<- vGene$genes$ensemblID
  
  ## Regress out residuals
  formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
  model=model.matrix(formula, data=vGene$targets)
  vGene$E<-cleaningY(vGene$E, model, P=2)
  
  ## Get ensembl ID and symbol of the gene
  if(! gene %in% rownames(vGene)){
    gene_symbol<-gene
    gene<-vGene$genes[which(vGene$genes$Symbol==gene), "ensemblID"]
  }
  else {
    gene_symbol<-vGene$genes[which(vGene$genes$ensemblID==gene), "Symbol"]
  }
  
  ## Gene ID for plot: Symbol-ensemblID
  gene_ID<-paste(gene_symbol, "-", gene, sep="")
  
  ## Extract expression values of the gene
  counts_gene<-vGene$E[rownames(vGene)==gene,]
  
  ## Gene's exons
  gene_exons<-top_exons[which(top_exons$ensemblID==gene),]
  ## Order them by FDR and extract the top 3
  top_gene_exons<-gene_exons[order(gene_exons$adj.P.Val),"exon_libdID"][1:3]
  ## Add specific exons
  if (!is.null(exon)){
    ## Exon ID to exon unique ID
    exon_gene<-strsplit(exon, "−")[[1]][1]
    exon_chr<-strsplit(strsplit(exon, "−")[[1]][2], ":")[[1]][1]
    exon_start<-strsplit(strsplit(exon, "−")[[1]][2], ":")[[1]][2]
    exon_end<-strsplit(exon, "−")[[1]][3]
    exon_ID<-top_exons[which(top_exons$Symbol==exon_gene & top_exons$seqnames==exon_chr
                             & top_exons$start==exon_start & top_exons$end==exon_end), "exon_libdID"]
    if (! exon_ID %in% top_gene_exons){
      top_gene_exons<-append(exon_ID, top_gene_exons) 
    }
  }
  
  ## Expression values of that exons
  vExon<-results_nic[[1]][[2]]
  ## Regress out residuals
  formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
  model=model.matrix(formula, data=vExon$targets)
  vExon$E<-cleaningY(vExon$E, model, P=2)
  counts_exons<-vExon$E[which(rownames(vExon) %in% top_gene_exons[1:3]),]
  counts_exons<-t(counts_exons)
  
  ## Data frame with all expression values and samples' group
  counts<-cbind(counts_gene, counts_exons, "Group"=vGene$targets$Group)
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
      q_value=signif(top_genes[which(top_genes$ensemblID==gene), "adj.P.Val"], digits = 3)
      y<-gene
      title<-gene_ID
    }
    
    ## Boxplot of the exons
    else {
      q_value<-signif(top_exons[which(top_exons$exon_libdID==colnames(counts)[i]), "adj.P.Val"], digits = 3)
      y<-colnames(counts)[i]
      ## Exon ID for plot
      title<-paste(top_exons[which(top_exons$exon_libdID==colnames(counts)[i]), "Symbol"], "-",
                   top_exons[which(top_exons$exon_libdID==colnames(counts)[i]), "seqnames"], ":", 
                   top_exons[which(top_exons$exon_libdID==colnames(counts)[i]), "start"], "-", 
                   top_exons[which(top_exons$exon_libdID==colnames(counts)[i]), "end"], sep="")
    }
    
    ## Plots
    p<-create_boxplot(counts, y, title, q_value)
    plots[[i]]<-p
  }
  
  plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol = 4)
  ggsave(here(paste("plots/04_DEA/02_Comparisons/Exon_analysis/Boxplots_", gene_symbol, ".pdf", sep="")), 
         width = 35, height = 10, units = "cm")

}

#########################################
## Nicotine boxplots for genes and exons
#########################################

#### Genes no DE but with DE exons ####

## Ankrd11: acts in head morphogenesis; expressed in cerebral cortex
gene_exons_boxplots("nicotine", "Ankrd11", "Ankrd11−chr8:122931509−122932058")

## Scaf11: predicted to be involved in spliceosomal complex assembly
gene_exons_boxplots("nicotine", "Scaf11", "Scaf11−chr15:96449740−96452139")

## Pcdh7: predicted to be involved in cell adhesion; expressed in CNS
gene_exons_boxplots("nicotine", "Pcdh7", "Pcdh7−chr5:57728114−57729465")

## Rb1cc1: involved in autophagosome assembly; expressed in NS
gene_exons_boxplots("nicotine", "Rb1cc1", "Rb1cc1−chr1:6250931−6251258")



#### DE genes with high mean(|tg-te|) ####

## Foxn3: acts upstream of or within craniofacial suture morphogenesis; expressed
## in CNS, PNS and sensory organs
gene_exons_boxplots("nicotine", "Foxn3", NULL)

## Kat7: enables histone acetyltransferase activity; expressed in brain
gene_exons_boxplots("nicotine", "Kat7", NULL)

## Dis3l2: enables 3'-5'-exoribonuclease activity; magnesium ion binding activity; 
## and poly(U) RNA binding activity; expressed in CNS
gene_exons_boxplots("nicotine", "Dis3l2", NULL)



#### DE exons with high |tg-te| values ####

## Ttc17: predicted to be involved in actin filament polymerization and cilium organization;
## expressed in cerebellum adult 
gene_exons_boxplots("nicotine", "Ttc17", "Ttc17−chr2:94327937−94328485")





#########################################
## Smoking boxplots for genes and exons
#########################################

#### Genes no DE but with DE exons ####

## Cd200: involved in negative regulation of neuroinflammatory response; expressed in NS
gene_exons_boxplots("smoking", "Cd200", "Cd200−chr16:45407194−45407311")

## Rps6kb1: acts upstream of or within behavioral fear response; expressed in CNS
gene_exons_boxplots("smoking", "Rps6kb1", "Rps6kb1−chr11:86543931−86544045")

## N4bp2l1: expressed in frontal lobe adult and cortex adult
gene_exons_boxplots("smoking", "N4bp2l1", "N4bp2l1−chr5:150591755−150592559")

## Nsg2: predicted to be involved in clathrin coat assembly and endosomal transport;
## expressed in NS
gene_exons_boxplots("smoking", "Nsg2", "Nsg2−chr11:32031770−32033414")

## Robo2: involved in neuron development; expressed in NS
gene_exons_boxplots("smoking", "Robo2", "Robo2−chr16:73898912−73899184")



#### DE genes with high mean(|tg-te|) ####

## Rpl13a: involved in macrophage chemotaxis; expressed in CNS
gene_exons_boxplots("smoking", "Rpl13a", NULL)

## Kntc1: involved in mitotic sister chromatid segregation; biased expression in CNS
gene_exons_boxplots("smoking", "Kntc1", NULL)

## Ndufb2: involved in mitochondrial respiratory chain complex I assembly; expressed in cranium
gene_exons_boxplots("smoking", "Ndufb2", NULL)

## 5031439G07Rik: is expressed in nervous system and neural retina
gene_exons_boxplots("smoking", "5031439G07Rik", NULL)

## Ttf2: involved in RNA splicing and mRNA processing; expressed in CNS
gene_exons_boxplots("smoking", "Ttf2", NULL)



#### DE exons with high |tg-te| values ####

## Sin3b: involved in skeletal muscle tissue development; expressed in CNS
gene_exons_boxplots("smoking", "Sin3b", "Sin3b−chr8:72744458−72744665")

## Tcaf1: involved in negative regulation of cell migration; broad expression in CNS
gene_exons_boxplots("smoking", "Tcaf1", "Tcaf1−chr6:42690299−42690523")

## Exoc4: involved in cilium assembly and protein localization to ciliary membrane;
## expressed in CNS
gene_exons_boxplots("smoking", "Exoc4", "Exoc4−chr6:33475997−33477149")

## Ganab: involved in N-glycan processing; expressed in CNS
gene_exons_boxplots("smoking", "Ganab", "Ganab−chr19:8899196−8899394")





### 1.2.3 Venn diagrams

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



## Define groups of DE exons
nic_up<-de_exons_nic[which(de_exons_nic$logFC>0),"exon_libdID"]
smo_up<-de_exons_smo[which(de_exons_smo$logFC>0),"exon_libdID"]
nic_down<-de_exons_nic[which(de_exons_nic$logFC<0),"exon_libdID"]
smo_down<-de_exons_smo[which(de_exons_smo$logFC<0),"exon_libdID"]


## Define groups of DE exons' genes
nic_up_genes<-unique(de_exons_nic[which(de_exons_nic$logFC>0),"ensemblID"])
smo_up_genes<-unique(de_exons_smo[which(de_exons_smo$logFC>0),"ensemblID"])
nic_down_genes<-unique(de_exons_nic[which(de_exons_nic$logFC<0),"ensemblID"])
smo_down_genes<-unique(de_exons_smo[which(de_exons_smo$logFC<0),"ensemblID"])

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
nic_DEG_up<-de_genes_pups_nicotine_fitted[which(de_genes_pups_nicotine_fitted$logFC>0),"ensemblID"]
nic_DEG_down<-de_genes_pups_nicotine_fitted[which(de_genes_pups_nicotine_fitted$logFC<0),"ensemblID"]
smo_DEG_up<-de_genes_pups_smoking_fitted[which(de_genes_pups_smoking_fitted$logFC>0),"ensemblID"]
smo_DEG_down<-de_genes_pups_smoking_fitted[which(de_genes_pups_smoking_fitted$logFC<0),"ensemblID"]

only_up_nic_DEG<-intersections[[1]]$ensemblID
only_up_smo_DEG<-intersections[[2]]$ensemblID
only_down_nic_DEG<-intersections[[3]]$ensemblID
only_down_smo_DEG<-intersections[[4]]$ensemblID
smoUp_nicUp_DEG<-intersections[[5]]$ensemblID
smoDown_nicDown_DEG<-intersections[[6]]$ensemblID
smoUp_nicDown_DEG<-intersections[[7]]$ensemblID
smoDown_nicUp_DEG<-intersections[[8]]$ensemblID





## Compare smoking VS nicotine DE exons

## Smo vs nic DE exons
DE_exons_smo_vs_nic <-list(
  "Nicotine"=de_exons_nic$exon_libdID,
  "Smoking"=de_exons_smo$exon_libdID)

## Smo vs nic Up DE exons
DE_exons_smo_vs_nic_Up <-list(
  "Nicotine up"=nic_up,
  "Smoking up"=smo_up)

## Smo vs nic Down DE exons
DE_exons_smo_vs_nic_Down <-list(
  "Nicotine down"=nic_down,
  "Smoking down"=smo_down)

## Smo Up vs nic Down DE exons
DE_exons_smoUp_vs_nicDown <-list(
  "Nicotine down"=nic_down,
  "Smoking up"=smo_up)

## Smo Down vs nic Up DE exons
DE_exons_smoDown_vs_nicUp <-list(
  "Nicotine up"=nic_up,
  "Smoking down"=smo_down)

## All DE exons
DE_exons_all<-list(
  "Nicotine up"=nic_up,
  "Nicotine down"=nic_down,
  "Smoking up"=smo_up, 
  "Smoking down"=smo_down)



## Venn diagrams
DE_lists<-list(DE_exons_smo_vs_nic, DE_exons_smo_vs_nic_Up, DE_exons_smo_vs_nic_Down, 
              DE_exons_smoUp_vs_nicDown, DE_exons_smoDown_vs_nicUp, DE_exons_all)
colors<-list(c("olivedrab1", "rosybrown2"), c("brown3", "coral"), c("cyan3", "cyan4"), 
             c("cyan3", "coral"), c("brown3", "cyan4"), c("brown3", "cyan3", "coral", "cyan4"))
venn_plot(DE_lists, colors, "smo_VS_nic_DE_exons", NULL)





## Compare DE exons' genes 

## Smo vs nic genes
DE_exons_genes_smo_vs_nic <-list(
  "Nicotine"=unique(de_exons_nic$ensemblID),
  "Smoking"=unique(de_exons_smo$ensemblID))

## Smo vs nic Up genes
DE_exons_genes_smo_vs_nic_Up <-list(
  "Nicotine up"=nic_up_genes,
  "Smoking up"=smo_up_genes)

## Smo vs nic Down genes
DE_exons_genes_smo_vs_nic_Down <-list(
  "Nicotine down"=nic_down_genes,
  "Smoking down"=smo_down_genes)

## Smo Up vs nic Down genes
DE_exons_genes_smoUp_vs_nicDown <-list(
  "Nicotine down"=nic_down_genes,
  "Smoking up"=smo_up_genes)

## Smo Down vs nic Up genes
DE_exons_genes_smoDown_vs_nicUp <-list(
  "Nicotine up"=nic_up_genes,
  "Smoking down"=smo_down_genes)

## Up vs Down genes
DE_exons_genes_Up_vs_Down <-list(
  "Up"=union(nic_up_genes, smo_up_genes),
  "Down"=union(nic_down_genes, smo_down_genes))

## Nic Up vs nic Down genes
DE_exons_genes_nicUp_vs_nicDown <-list(
  "Nicotine up"=nic_up_genes,
  "Nicotine down"=nic_down_genes)

## Smo Up vs smo Down genes
DE_exons_genes_smoUp_vs_smoDown <-list(
  "Smoking up"=smo_up_genes,
  "Smoking down"=smo_down_genes)

## All DE exons
DE_exons_genes_all<-list(
  "Nicotine up"=nic_up_genes,
  "Nicotine down"=nic_down_genes,
  "Smoking up"=smo_up_genes, 
  "Smoking down"=smo_down_genes)



## Venn diagrams
DE_lists<-list(DE_exons_genes_smo_vs_nic, DE_exons_genes_smo_vs_nic_Up, 
               DE_exons_genes_smo_vs_nic_Down, DE_exons_genes_smoUp_vs_nicDown, 
               DE_exons_genes_smoDown_vs_nicUp, DE_exons_genes_Up_vs_Down,
               DE_exons_genes_nicUp_vs_nicDown, DE_exons_genes_smoUp_vs_smoDown, 
               DE_exons_genes_all)
colors<-list(c("olivedrab1", "rosybrown2"), c("brown3", "coral"), 
             c("cyan3", "cyan4"), c("cyan3", "coral"), c("brown3", "cyan4"), 
             c("firebrick", "dodgerblue"),c("brown3", "cyan3"), c("coral", "cyan4"),
             c("brown3", "cyan3", "coral", "cyan4"))
venn_plot(DE_lists, colors, "smo_VS_nic_DE_exons_genes", NULL)





## Compare DE exons' genes with DEG

## All DEG vs all exons' genes
DEG_vs_Exons_all <-list(
  "DEG"= union(union(nic_DEG_down, nic_DEG_up), union(smo_DEG_down, smo_DEG_up)),
  "DE exons' genes"=union(union(nic_down_genes, nic_up_genes), union(smo_down_genes, smo_up_genes))
)

## Nic DEGs vs nic exons' genes 
DEG_nic_vs_Exons_nic <-list(
  "DEG"= union(nic_DEG_down, nic_DEG_up),
  "DE exons' genes"=union(nic_down_genes, nic_up_genes)
)

## Smo DEGs vs mo exons' genes
DEG_smo_vs_Exons_smo <-list(
  "DEG"= union(smo_DEG_down, smo_DEG_up),
  "DE exons' genes"=union(smo_down_genes, smo_up_genes)
)

## Up DEGs vs up exons' genes
DEG_up_vs_Exons_up <-list(
  "DEG"= union(nic_DEG_up, smo_DEG_up),
  "DE exons' genes"=union(nic_up_genes, smo_up_genes)
)

## Down DEGs vs down exons' genes
DEG_down_vs_Exons_down <-list(
  "DEG"= union(nic_DEG_down, smo_DEG_down),
  "DE exons' genes"=union(nic_down_genes, smo_down_genes)
)

## Nic up DEGs vs nic up exons' genes
DEG_nicUp_vs_Exons_nicUp <-list(
  "DEG"=nic_up_genes,
  "DE exons' genes"=nic_DEG_up
)

## Nic down DEGs vs nic down exons' genes
DEG_nicDown_vs_Exons_nicDown <-list(
  "DEG"=nic_down_genes,
  "DE exons' genes"=nic_DEG_down
)

## Smo up DEGs vs smo up exons' genes
DEG_smoUp_vs_Exons_smoUp <-list(
  "DEG"=smo_up_genes,
  "DE exons' genes"=smo_DEG_up
)

## Smo down DEGs vs smo down exons' genes
DEG_smoDown_vs_Exons_smoDown <-list(
  "DEG"=smo_down_genes,
  "DE exons' genes"=smo_DEG_down
)



## Venn diagrams
DE_lists<-list(DEG_vs_Exons_all, DEG_nic_vs_Exons_nic, DEG_smo_vs_Exons_smo,
               DEG_up_vs_Exons_up, DEG_down_vs_Exons_down, DEG_nicUp_vs_Exons_nicUp,
               DEG_nicDown_vs_Exons_nicDown, DEG_smoUp_vs_Exons_smoUp, DEG_smoDown_vs_Exons_smoDown)
colors<-list(c("rosybrown2", "navajowhite2"), c("plum", "lightgoldenrod2"), 
             c("pink1", "khaki2"), c("rosybrown4", "navajowhite4"), c("rosybrown1", "navajowhite1"), 
             c("plum3", "lightgoldenrod4"),c("plum1", "lightgoldenrod1"), c("pink3", "khaki3"),
             c("pink", "khaki1"))
venn_plot(DE_lists, colors, "DEG_VS_exons_genes", c("All", "Nicotine", "Smoking", "Up", "Down", 
                                                    "Nicotine Up", "Nicotine Down", "Smoking Up", "Smoking Down"))





## Compare DEG intersections with those of exons' genes

## Only up nic DEGs vs exons' genes 
DEG_vs_Exons_only_up_nic <-list(
  "DEG"= only_up_nic_DEG,
  "DE exons' genes"= only_up_nic_genes
)

## Only down nic DEGs vs exons' genes 
DEG_vs_Exons_only_down_nic <-list(
  "DEG"= only_down_nic_DEG,
  "DE exons' genes"= only_down_nic_genes
)

## Only up smo DEGs vs exons' genes 
DEG_vs_Exons_only_up_smo <-list(
  "DEG"= only_up_smo_DEG,
  "DE exons' genes"= only_up_smo_genes
)

## Only down smo DEGs vs exons' genes 
DEG_vs_Exons_only_down_smo <-list(
  "DEG"= only_down_smo_DEG,
  "DE exons' genes"= only_down_smo_genes
)

## Smo and nic up DEGs vs exons' genes
DEG_vs_Exons_smoUp_nicUp <-list(
  "DEG"=smoUp_nicUp_DEG,
  "DE exon's genes"=smoUp_nicUp_genes
)

## Smo and nic down DEGs vs exons' genes
DEG_vs_Exons_smoDown_nicDown <-list(
  "DEG"=smoDown_nicDown_DEG,
  "DE exon's genes"=smoDown_nicDown_genes
)

## Smo up and nic down DEGs vs exons' genes
DEG_vs_Exons_smoUp_nicDown <-list(
  "DEG"=smoUp_nicDown_DEG,
  "DE exon's genes"=smoUp_nicDown_genes
)

## Smo down and nic up DEGs vs exons' genes
DEG_vs_Exons_smoDown_nicUp <-list(
  "DEG"=smoDown_nicUp_DEG,
  "DE exon's genes"=smoDown_nicUp_genes
)



## Venn diagrams
DE_lists<-list(DEG_vs_Exons_only_up_nic, DEG_vs_Exons_only_down_nic, DEG_vs_Exons_only_up_smo,
               DEG_vs_Exons_only_down_smo, DEG_vs_Exons_smoUp_nicUp, DEG_vs_Exons_smoDown_nicDown,
               DEG_vs_Exons_smoUp_nicDown, DEG_vs_Exons_smoDown_nicUp)
colors<-list(c("darkslategray3", "lightgoldenrod3"), c("darkslategray1", "lightgoldenrod1"), 
             c("lightpink3", "olivedrab3"), c("pink", "olivedrab1"), c("palegreen3", "navajowhite3"),
             c("palegreen", "navajowhite"), c("palegreen","navajowhite3"), c("palegreen3", "navajowhite"))
venn_plot(DE_lists, colors, "intersections_DEG_VS_exons_genes", c("Only up in nic", "Only down in nic",
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
# date     2022-11-19
# rstudio  2022.07.2+576 Spotted Wakerobin (desktop)
# pandoc   2.19.2 @ C:/Program Files/RStudio/bin/quarto/bin/tools/ (via rmarkdown)


