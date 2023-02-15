
## 1.2 Comparison of DEG 

load(here("processed-data/04_DEA/Gene_analysis/top_genes_blood_smoking_naive.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_blood_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_blood_smoking_interaction.Rdata"))

load(here("processed-data/04_DEA/Gene_analysis/top_genes_adults_nicotine_naive.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_adults_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_adults_nicotine_interaction.Rdata"))

load(here("processed-data/04_DEA/Gene_analysis/top_genes_adults_smoking_naive.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_adults_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_adults_smoking_interaction.Rdata"))

load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_naive.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_naive.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_interaction.Rdata"))

load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_naive.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_naive.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_interaction.Rdata"))

## Pup brain Tx data 
load(here("processed-data/04_DEA/Tx_analysis/top_tx_nic.Rdata"))
load(here("processed-data/04_DEA/Tx_analysis/top_tx_smo.Rdata"))

## Data from prenatal and adult human postmortem prefrontal cortices 
## (from Steve's paper https://github.com/LieberInstitute/Smoking_DLPFC_Devel/blob/master/README.md)
load(here("raw-data/Genes_DE_sva.rda"))



### 1.2.1 T-stats plots

## Function to add DE info of genes in both groups
add_DE_info <-function(top_genes1, top_genes2, name_1, name_2) {
  DE<-vector()
  for (i in 1:dim(top_genes1)[1]) {
    ## DE genes in both groups
    if (top_genes1$adj.P.Val[i]<0.05 && top_genes2$adj.P.Val[i]<0.05) {
      DE<-append(DE, "sig Both")
    }
    ## DE genes in only one of the groups
    else if (top_genes1$adj.P.Val[i]<0.05 && !top_genes2$adj.P.Val[i]<0.05) {
      DE<-append(DE, paste("sig", name_1))
    }
    
    else if (top_genes2$adj.P.Val[i]<0.05 && !top_genes1$adj.P.Val[i]<0.05) {
      DE<-append(DE, paste("sig",name_2))
    }
    ## No DE genes in neither group
    else {
      DE<-append(DE, "None")
    }
  }
  return(DE)
}


## Compare t-stats of genes from different groups of samples
t_stat_plot <- function(top_genes1, top_genes2, name_1, name_2, model_name){
  
  ## Correlation coeff
  rho <- cor(top_genes1$t, top_genes2$t, method = "spearman")
  rho_anno = paste0("rho = ", format(round(rho, 2), nsmall = 2))
  
  ## Merge data
  t_stats<-data.frame(t1=top_genes1$t, t2=top_genes2$t)
  ## Add DE info for both groups
  t_stats$DEG<-add_DE_info(top_genes1, top_genes2, name_1, name_2)
  
  cols <- c("red", "#ffad73","#26b3ff", "dark grey") 
  names(cols)<-c("sig Both", paste0("sig ", name_1), paste0("sig ", name_2), "None")
  alphas <- c( 1, 1, 1,0.5)  
  names(alphas)<-c("sig Both", paste0("sig ", name_1), paste0("sig ", name_2), "None")
   
  plot <- ggplot(t_stats, aes(x = t1, y = t2, color=DEG, alpha=DEG)) +
     geom_point(size = 1) +
     labs(x = paste("t-stats", name_1), 
          y = paste("t-stats", name_2),
          title = model_name, 
          subtitle = rho_anno, 
          parse = T) +
     theme_bw() +
     scale_color_manual(values = cols) + 
     scale_alpha_manual(values = alphas) 
 return(plot)
}


## Function to create multiple t-stats plots
tstats_plots<-function(top_genes_pairs, name_1, name_2, models){
  
  plots<-list()
  for (i in 1:length(top_genes_pairs)){
    p<-t_stat_plot(top_genes_pairs[[i]][[1]], top_genes_pairs[[i]][[2]], name_1, name_2, models[i])
    plots[[i]]<-p
  }
  plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol=2)
  ggsave(filename=paste("plots/04_DEA/02_Comparisons/Gene_analysis/t_stats_",gsub(" ", "_", name_1), "_VS_",
                        gsub(" ", "_", name_2), ".pdf", sep=""), 
                        height = 20, width = 25, units = "cm")
}



## Analyses

################################################################################
##       1. Compare tissues - Analysis of blood vs brain biomarkers
################################################################################

#################### 1.1 Compare t-stats in the 3 models #######################

####### Smoking adult blood vs Smoking adult brain #######
top_genes_pairs<-list(list(top_genes_adults_smoking_naive, top_genes_blood_smoking_naive), 
                      list(top_genes_adults_smoking_fitted, top_genes_blood_smoking_fitted),
                      list(top_genes_adults_smoking_interaction, top_genes_blood_smoking_interaction))
models<-c("Naive model", "Fitted model", "Interaction model")
tstats_plots(top_genes_pairs,  "Smoking adult brain", "Smoking adult blood", models)


####### Smoking adult blood vs Nicotine adult brain #######
top_genes_pairs<-list(list(top_genes_adults_nicotine_naive, top_genes_blood_smoking_naive), 
                      list(top_genes_adults_nicotine_fitted, top_genes_blood_smoking_fitted),
                      list(top_genes_adults_nicotine_interaction, top_genes_blood_smoking_interaction))
models<-c("Naive model", "Fitted model", "Interaction model")
tstats_plots(top_genes_pairs,  "Nicotine adult brain", "Smoking adult blood", models)


####### Smoking adult blood vs Smoking pup brain #######
top_genes_pairs<-list(list(top_genes_pups_smoking_naive, top_genes_blood_smoking_naive), 
                      list(top_genes_pups_smoking_fitted, top_genes_blood_smoking_fitted),
                      list(top_genes_pups_smoking_interaction, top_genes_blood_smoking_interaction))
models<-c("Naive model", "Fitted model", "Interaction model")
tstats_plots(top_genes_pairs,  "Smoking pup brain", "Smoking adult blood", models)


####### Smoking adult blood vs Nicotine pup brain #######
top_genes_pairs<-list(list(top_genes_pups_nicotine_naive, top_genes_blood_smoking_naive), 
                      list(top_genes_pups_nicotine_fitted, top_genes_blood_smoking_fitted),
                      list(top_genes_pups_nicotine_interaction, top_genes_blood_smoking_interaction))
models<-c("Naive model", "Fitted model", "Interaction model")
tstats_plots(top_genes_pairs,  "Nicotine pup brain", "Smoking adult blood", models)



#######################  1.2 Search mouse brain genes/txs that replicate in mouse blood ####################### 
## (With p<0.05 in blood, FDR<0.05 in pup brain/p<0.05 in adult brain, and the same logFC sign in both tissues)
## Use mouse genes and txs from fitted models only

t_stat_plot_brain_blood_replication <- function(age_mouse, expt_mouse, feature){
  
  ## Define blood dataset
  top_genes_blood <- top_genes_blood_smoking_fitted
  
  ## Compare blood genes vs brain genes
  if (feature=="genes"){
    
    ## Define brain dataset
    top_genes_brain <-eval(parse_expr(paste("top_genes", age_mouse, expt_mouse, "fitted", sep="_")))
    
    ## Merge brain and blood data
    ## (Brain and blood datasets contain the same genes in the same order)
    t_stats<-data.frame(ensemblID=rownames(top_genes_brain), gene_symbol=top_genes_brain$Symbol, t_blood=top_genes_blood$t, t_brain=top_genes_brain$t,
                        FDR_blood=top_genes_blood$adj.P.Val, FDR_brain=top_genes_brain$adj.P.Val, P.Val_blood=top_genes_blood$P.Value, P.Val_brain=top_genes_brain$P.Value,
                        logFC_blood=top_genes_blood$logFC, logFC_brain=top_genes_brain$logFC)
   
    ## Add DE and replication info for each gene
    DE <-vector()
    for (i in 1:dim(t_stats)[1]) {
      
      ## Pup brain genes
      if (age_mouse=="pups"){
        
        ## Pup brain genes (FDR<0.05) that replicate in blood (p<0.05) 
        if ((t_stats$P.Val_blood[i]<0.05 && t_stats$FDR_brain[i]<0.05) &&
            (sign(t_stats$logFC_blood[i])==sign(t_stats$logFC_brain[i]))) {
          DE<-append(DE, "Replicating genes (p<0.05 in blood, FDR<0.05 in brain)")
        }
        
        ## DEG in brain
        ## (Note that all replicating genes of pup brain are also DEG (FDR<0.05))
        else if (t_stats$FDR_brain[i]<0.05){
          DE<-append(DE, "Signif in brain (FDR<0.05)")
        }
        
        ## Non-significant genes 
        else {
          DE<-append(DE, "n.s. genes")      
        }
      }
      
      ## Adult brain genes 
      else {
        
        ## Adult brain genes (p<0.05) that replicate in blood (p<0.05) 
        if ((t_stats$P.Val_blood[i]<0.05 && t_stats$P.Val_brain[i]<0.05) &&
            (sign(t_stats$logFC_blood[i])==sign(t_stats$logFC_brain[i]))) {
          DE<-append(DE, "Replicating genes (p<0.05 in blood and brain)")
        }
        else {
          DE<-append(DE, "n.s. genes")      
        }
      }
    }
    
    t_stats$DE<- DE
    
    ## Correlation coeff between t-stats of blood vs brain genes 
    rho <- cor(t_stats$t_blood, t_stats$t_brain, method = "spearman")
    rho_anno = paste0("rho = ", format(round(rho, 2), nsmall = 2))
    
    ## Colors and alphas for plot
    if (age_mouse=="pups"){
      cols <- c("red",   "#26b3ff", "dark grey") 
      names(cols)<-c("Replicating genes (p<0.05 in blood, FDR<0.05 in brain)", "Signif in brain (FDR<0.05)", "n.s. genes")
      alphas <- c(1, 1, 0.5)  
      names(alphas)<-c("Replicating genes (p<0.05 in blood, FDR<0.05 in brain)", "Signif in brain (FDR<0.05)", "n.s. genes")
    }
    else {
      cols <- c("red","dark grey") 
      names(cols)<-c("Replicating genes (p<0.05 in blood and brain)", "n.s. genes")
      alphas <- c(1, 0.5)  
      names(alphas)<-c("Replicating genes (p<0.05 in blood and brain)", "n.s. genes")
    }
    
    ## Add labels of interesting genes 
    
    ## 3 most significant replicating genes in blood
    rep_genes <- t_stats[which(t_stats$DE==names(alphas)[1]),]
    top_rep_genes <- rep_genes[order(rep_genes$FDR_blood),"gene_symbol"][1:3]
    
    label <- vector()
    for (i in 1:dim(t_stats)[1]){
      ## Labels of the top 3 replicating genes
      if (t_stats$gene_symbol[i] %in% top_rep_genes){
        label <- append(label, paste(t_stats$ensemblID[i], t_stats$gene_symbol[i], sep="-"))
      }
      else{
        label <- append(label, NA)
      }
    }
    
    t_stats$label <- label
    
    ## Plot
    plot <- ggplot(t_stats, aes(x = t_brain, y = t_blood, color=DE, alpha=DE, label=label)) +
      geom_point(size = 1) +
      labs(x = paste("t-stats in", capitalize(expt_mouse), substr(age_mouse, 1, nchar(age_mouse)-1), "brain"), 
           y = "t-stats in Smoking adult blood",
           title = paste(capitalize(expt_mouse), "mouse brain vs Smoking mouse blood", sep=" "), 
           subtitle = rho_anno, 
           parse = T) +
      geom_label_repel(fill="white", size=2, max.overlaps = Inf,  
                       box.padding = 0.2, 
                       show.legend=FALSE) +
      theme_bw() +
      scale_color_manual(values = cols) + 
      scale_alpha_manual(values = alphas)
    
    plot + theme(legend.text = element_text(size=8))
    plot
    ggsave(filename=paste("plots/04_DEA/02_Comparisons/Gene_analysis/t_stats_replication_", capitalize(expt_mouse), "_", substr(age_mouse, 1, nchar(age_mouse)-1),
                          "_Brain_vs_Smoking_adult_Blood.pdf", sep=""), height = 12, width = 20, units = "cm")
    
    
    ## Quantify the number of brain genes that replicate in blood
    
    ## Total unique pup brain DEG (FDR<0.05)
    total_pups_DEG=length(unique(t_stats[which(t_stats$FDR_brain<0.05), "ensemblID"]))
    ## Total unique adult brain genes with p<0.05
    total_adults_P_val_genes=length(unique(t_stats[which(t_stats$P.Val_brain<0.05), "ensemblID"]))
    ## Unique replicating genes 
    rep_genes=length(unique(t_stats[which(t_stats$DE==names(alphas)[1]),"ensemblID"]))
  
    if (age_mouse=="pups"){
      ## Percentage 
      percentage=signif(rep_genes / total_pups_DEG *100, 3)
      print(paste(rep_genes, "out of", total_pups_DEG, "DEG in", expt_mouse, "pup brain (FDR<0.05) replicate in smoking adult blood (with p<0.05 and same logFC direction) -", 
                  paste(percentage, "%", sep="")))
    }
    else {
      ## Percentage 
      percentage=signif(rep_genes / total_adults_P_val_genes *100, 3)
      print(paste(rep_genes, "out of", total_adults_P_val_genes, "genes in", expt_mouse, "adult brain (p<0.05) replicate in smoking adult blood (also p<0.05 and same logFC direction) -", 
                  paste(percentage, "%", sep="")))
    }

  }
  
  ## Compare blood genes vs brain txs
  else {
    
    top_tx_brain <-eval(parse_expr(paste("top_tx", substr(expt_mouse, 1, 3), sep="_")))
    
    ## Common genes in gene and tx datasets
    tx_genes<-unique(top_tx_brain$ensembl_id)
    common_tx_genes<-tx_genes[which(tx_genes %in% top_genes_blood$gencodeID)]
    ## Txs info
    t_stats<-top_tx_brain[which(top_tx_brain$ensembl_id %in% common_tx_genes), c("transcript_id", "P.Value", "adj.P.Val", "Symbol", 
                                                              "t", "ensembl_id", "logFC")]
    colnames(t_stats)[c(2,3,5,7)] <- paste(colnames(t_stats[c(2,3,5,7)]), "tx", sep="_")
    
    ## Add t-stats and FDRs of the transcripts' genes 
    t_tx_genes<-vector()
    P_val_tx_genes<-vector()
    logFC_tx_genes<-vector()
    for (i in 1:dim(t_stats)[1]){
      t<-top_genes_blood[which(top_genes_blood$gencodeID==t_stats$ensembl_id[i]), "t"]
      P_val<-top_genes_blood[which(top_genes_blood$gencodeID==t_stats$ensembl_id[i]), "P.Value"]
      logFC<-top_genes_blood[which(top_genes_blood$gencodeID==t_stats$ensembl_id[i]), "logFC"]
      t_tx_genes<-append(t_tx_genes, t)
      P_val_tx_genes<-append(P_val_tx_genes, P_val)
      logFC_tx_genes<-append(logFC_tx_genes, logFC)
    }
    t_stats$t_gene<-t_tx_genes
    t_stats$P.Val_gene<-P_val_tx_genes
    t_stats$logFC_gene<-logFC_tx_genes
  
  
    ## Add DE and replication info for each tx
    DE <-vector()
    for (i in 1:dim(t_stats)[1]) {
      
      ## Pup brain txs (FDR<0.05) that replicate in mouse blood genes (p<0.05) 
      if ((t_stats$P.Val_gene[i]<0.05 && t_stats$adj.P.Val_tx[i]<0.05) &&
               (sign(t_stats$logFC_tx[i])==sign(t_stats$logFC_gene[i]))) {
        DE<-append(DE, "Replicating txs (blood gene with p<0.05, brain tx with FDR<0.05)")
      }
      
      ## DE txs in brain
      ## (Note that all replicating txs of pup brain are also significant (FDR<0.05))
      else if (t_stats$adj.P.Val_tx[i]<0.05){
        DE<-append(DE, "Signif txs in brain (FDR<0.05)")
      }
      
      ## Non-significant txs
      else {
        DE<-append(DE, "n.s. txs")      
      }
    }
    
    t_stats$DE<- DE
    
    ## Correlation coeff between t-stats of blood genes vs brain txs 
    rho <- cor(t_stats$t_gene, t_stats$t_tx, method = "spearman")
    rho_anno = paste0("rho = ", format(round(rho, 2), nsmall = 2))
    
    ## Colors and alphas for plot
    cols <- c("red", "#26b3ff", "dark grey") 
    names(cols)<-c("Replicating txs (blood gene with p<0.05, brain tx with FDR<0.05)", "Signif txs in brain (FDR<0.05)", "n.s. txs")
    alphas <- c(1, 1, 0.5)  
    names(alphas)<-c("Replicating txs (blood gene with p<0.05, brain tx with FDR<0.05)", "Signif txs in brain (FDR<0.05)", "n.s. txs")
    
    
    ## Add labels of interesting txs and their genes
    
    ## 3 most significant replicating txs 
    rep_txs <- t_stats[which(t_stats$DE==names(alphas)[1]),]
    top_rep_txs <- rep_txs[order(rep_txs$adj.P.Val_tx), "transcript_id"][1:3]
    
    label <- vector()
    for (i in 1:dim(t_stats)[1]){
      ## Labels of the top 3 replicating txs
      if (t_stats$transcript_id[i] %in% top_rep_txs) {
        label <- append(label, paste(t_stats$Symbol[i], t_stats$transcript_id[i], sep="-"))
      }
      else{
        label <- append(label, NA)
      }
    }
    
    t_stats$label <- label
    
    ## Plot
    plot <- ggplot(t_stats, aes(x = t_tx, y = t_gene, color=DE, alpha=DE, label=label)) +
      geom_point(size = 1) +
      labs(x = paste("t-stats of txs from", capitalize(expt_mouse), "pup brain"), 
           y = "t-stats of genes from Smoking adult blood",
           title = paste(capitalize(expt_mouse), "mouse brain vs Smoking mouse blood", sep=" "), 
           subtitle = rho_anno, 
           parse = T) +
      geom_label_repel(fill="white", size=2, max.overlaps = Inf,  
                       box.padding = 0.2, 
                       show.legend=FALSE) +
      theme_bw() +
      scale_color_manual(values = cols) + 
      scale_alpha_manual(values = alphas)
    
    plot + theme(legend.text = element_text(size=8))
    plot
    ggsave(filename=paste("plots/04_DEA/02_Comparisons/Gene_analysis/t_stats_replication_", capitalize(substr(expt_mouse, 1, 3)), 
                          "PupBrain_Tx_vs_SmoAdultBlood_Genes.pdf", sep=""), height = 12, width = 20, units = "cm")
    
    
    ## Quantify the number of brain txs that replicate in blood

    ## Total unique brain DE txs
    total_pups_DEtxs=length(top_tx_brain[which(top_tx_brain$adj.P.Val<0.05), "transcript_id"])
    ## Unique replicating txs
    rep_txs=length(unique(t_stats[which(t_stats$DE==names(alphas)[1]),"transcript_id"]))
    ## Percentage 
    percentage=signif(rep_txs / total_pups_DEtxs *100, 3)
    print(paste(rep_txs, "out of", total_pups_DEtxs, "DE txs in", expt_mouse, 
                "pup brain (FDR<0.05) replicate in smoking adult blood genes (with p<0.05 and same logFC direction) -", paste(percentage, "%", sep="")))

  }
 
  return(t_stats)
}


################### 1.2.1 Brain genes vs blood genes ###################

####### Smoking adult blood vs Smoking adult brain #######
Blood_vs_SmoAdultBrain_data <- t_stat_plot_brain_blood_replication(age_mouse = "adults", expt_mouse = "smoking", feature = "genes")
## "37 out of 772 genes in smoking adult brain (p<0.05) replicate in smoking adult blood (also p<0.05 and same logFC direction) - 4.79%"
save(Blood_vs_SmoAdultBrain_data, file="processed-data/04_DEA/Gene_analysis/Blood_vs_SmoAdultBrain_data.Rdata")
## Add replication info to top_table data
top_genes_adults_smoking_fitted$replication_in_blood <- Blood_vs_SmoAdultBrain_data$DE
save(top_genes_adults_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_adults_smoking_fitted.Rdata")


####### Smoking adult blood vs Nicotine adult brain #######
Blood_vs_NicAdultBrain_data <- t_stat_plot_brain_blood_replication(age_mouse = "adults", expt_mouse = "nicotine", feature = "genes")
## "33 out of 679 genes in nicotine adult brain (p<0.05) replicate in smoking adult blood (also p<0.05 and same logFC direction) - 4.86%"
save(Blood_vs_NicAdultBrain_data, file="processed-data/04_DEA/Gene_analysis/Blood_vs_NicAdultBrain_data.Rdata")
top_genes_adults_nicotine_fitted$replication_in_blood <- Blood_vs_NicAdultBrain_data$DE
save(top_genes_adults_nicotine_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_adults_nicotine_fitted.Rdata")


####### Smoking adult blood vs Smoking pup brain #######
Blood_vs_SmoPupBrain_data <- t_stat_plot_brain_blood_replication(age_mouse = "pups", expt_mouse = "smoking", feature = "genes")
## "126 out of 4165 DEG in smoking pup brain (FDR<0.05) replicate in smoking adult blood (with p<0.05 and same logFC direction) - 3.03%"
save(Blood_vs_SmoPupBrain_data, file="processed-data/04_DEA/Gene_analysis/Blood_vs_SmoPupBrain_data.Rdata")
top_genes_pups_smoking_fitted$replication_in_blood <- Blood_vs_SmoPupBrain_data$DE
save(top_genes_pups_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata")


####### Smoking adult blood vs Nicotine pup brain #######
Blood_vs_NicPupBrain_data <- t_stat_plot_brain_blood_replication(age_mouse = "pups", expt_mouse = "nicotine", feature = "genes")
## "15 out of 1010 DEG in nicotine pup brain (FDR<0.05) replicate in smoking adult blood (with p<0.05 and same logFC direction) - 1.49%"
save(Blood_vs_NicPupBrain_data, file="processed-data/04_DEA/Gene_analysis/Blood_vs_NicPupBrain_data.Rdata")
top_genes_pups_nicotine_fitted$replication_in_blood <- Blood_vs_NicPupBrain_data$DE
save(top_genes_pups_nicotine_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata")




################### 1.2.2 Brain txs vs blood genes ###################

####### Smoking adult blood (genes) vs Smoking pup brain Txs #######
Blood_vs_SmoPupBrainTx_data <- t_stat_plot_brain_blood_replication(age_mouse = NULL, expt_mouse = "smoking", feature = "txs")
## "112 out of 4059 DE txs in smoking pup brain (FDR<0.05) replicate in smoking adult blood genes (with p<0.05 and same logFC direction) - 2.76%"
save(Blood_vs_SmoPupBrainTx_data, file="processed-data/04_DEA/Gene_analysis/Blood_vs_SmoPupBrainTx_data.Rdata")
                                  
####### Smoking adult blood (genes) vs Nicotine pup brain Txs #######
Blood_vs_NicPupBrainTx_data <- t_stat_plot_brain_blood_replication(age_mouse = NULL, expt_mouse = "nicotine", feature = "txs")
## "9 out of 232 DE txs in nicotine pup brain (FDR<0.05) replicate in smoking adult blood genes (with p<0.05 and same logFC direction) - 3.88%"
save(Blood_vs_NicPupBrainTx_data, file="processed-data/04_DEA/Gene_analysis/Blood_vs_NicPupBrainTx_data.Rdata")



## (See code below to search for mouse blood genes that replicate in human brain and human brain genes that replicate in mouse blood)




######################################################
# 2. Compare experiments - Smoking vs Nicotine genes
######################################################

####### Smoking adults vs nicotine adults #######
top_genes_pairs<-list(list(top_genes_adults_smoking_naive, top_genes_adults_nicotine_naive), 
                      list(top_genes_adults_smoking_fitted, top_genes_adults_nicotine_fitted),
                      list(top_genes_adults_smoking_interaction, top_genes_adults_nicotine_interaction))
models<-c("Naive model", "Fitted model", "Interaction model")
tstats_plots(top_genes_pairs,  "Smoking adults", "Nicotine adults", models)


####### Smoking pups vs nicotine pups #######
top_genes_pairs<-list(list(top_genes_pups_smoking_naive, top_genes_pups_nicotine_naive), 
                      list(top_genes_pups_smoking_fitted, top_genes_pups_nicotine_fitted),
                      list(top_genes_pups_smoking_interaction, top_genes_pups_nicotine_interaction))
models<-c("Naive model", "Fitted model", "Interaction model")
tstats_plots(top_genes_pairs,  "Smoking pups", "Nicotine pups", models)


########################################
# 3. Compare ages - Pup vs adult genes
########################################

####### Smoking pups vs smoking adults #######
top_genes_pairs<-list(list(top_genes_pups_smoking_naive, top_genes_adults_smoking_naive), 
                      list(top_genes_pups_smoking_fitted, top_genes_adults_smoking_fitted),
                      list(top_genes_pups_smoking_interaction, top_genes_adults_smoking_interaction))
models<-c("Naive model", "Fitted model", "Interaction model")
tstats_plots(top_genes_pairs,  "Smoking pups", "Smoking adults", models)


####### Nicotine pups vs nicotine adults #######
top_genes_pairs<-list(list(top_genes_pups_nicotine_naive, top_genes_adults_nicotine_naive), 
                      list(top_genes_pups_nicotine_fitted, top_genes_adults_nicotine_fitted),
                      list(top_genes_pups_nicotine_interaction, top_genes_adults_nicotine_interaction))
models<-c("Naive model", "Fitted model", "Interaction model")
tstats_plots(top_genes_pairs,  "Nicotine pups", "Nicotine adults", models)


####################################################
# 4. Compare models - Naive vs fitted model genes
####################################################

####### Naive vs fitted models for nicotine pups #######
t<-t_stat_plot(top_genes_pups_nicotine_naive, top_genes_pups_nicotine_fitted, 
               "Naive model", "Fitted model", "Nicotine pups")
ggsave("plots/04_DEA/02_Comparisons/Gene_analysis/t_stats_Naive_VS_Fitted_Nicotine.pdf", t, 
       height = 10, width = 12, units = "cm")


####### Naive vs fitted models for smoking pups #######
t<-t_stat_plot(top_genes_pups_smoking_naive, top_genes_pups_smoking_fitted, 
               "Naive model", "Fitted model", "Smoking pups")
ggsave("plots/04_DEA/02_Comparisons/Gene_analysis/t_stats_Naive_VS_Fitted_Smoking.pdf", t, 
       height = 10, width = 12, units = "cm")



################################################################################
##           5. Compare human brain vs mouse brain/blood genes
################################################################################
## Samples from prenatal and adult human brain were exposed to smoking 
## Compare mouse genes from fitted models only


################### 5.1 Mouse genes that replicate in human #################### 

## Genes in prenatal and adult human brain are the same
setdiff(rownames(fetalGene), rownames(adultGene))
# character(0)

## " " to NA 
fetalGene[fetalGene == ""] <- NA
adultGene[adultGene == ""] <- NA

## Ensembl IDs of human genes 
fetalGene$ensemblID <- rownames(fetalGene)
adultGene$ensemblID <- rownames(adultGene)

## Find the homologous genes for human in mouse
human_mouse_ids<-biomart(genes  = fetalGene$ensemblID,
                 mart       = "ENSEMBL_MART_ENSEMBL",
                 dataset    = "hsapiens_gene_ensembl",
                 attributes = c("mmusculus_homolog_ensembl_gene", "mmusculus_homolog_associated_gene_name"),
                 filters    = "ensembl_gene_id")

## Common genes in human homologs and mouse datasets
## (All mouse datasets contain the same genes in the same order)
common_genes <- human_mouse_ids[which(human_mouse_ids$mmusculus_homolog_ensembl_gene %in% top_genes_pups_nicotine_fitted$ensemblID),]
common_genes$human_ensembl_gene_id <- common_genes$ensembl_gene_id
common_genes$ensembl_gene_id <- NULL

## Create plots to check if mouse genes replicate (FDR<5% for pups and p-value<5% for adults) in human brain (with p-value<5% and same logFC sign)
t_stat_plot_mouse_in_human <- function(age_mouse, expt_mouse, tissue_mouse, age_human){
  
  ## Define mouse dataset
  if (tissue_mouse=="blood"){
    top_genes <-eval(parse_expr(paste("top_genes", tissue_mouse, expt_mouse, "fitted", sep="_")))
    signif_measure_mouse="p-value"
  }
  else {
    top_genes <-eval(parse_expr(paste("top_genes", age_mouse, expt_mouse, "fitted", sep="_")))
    if (age_mouse=="pups"){
      signif_measure_mouse="FDR"
    }
    else {
      signif_measure_mouse="p-value"
    }
  }
  
  ## Define human dataset
  if (age_human=="prenatal"){
    humanGene <-fetalGene
  }
  else {
    humanGene <-adultGene
  }
  
  ## Extract mouse and human data of those common genes
  human_mouse_data <- data.frame(matrix(nrow = nrow(common_genes), ncol = 12))
  colnames(human_mouse_data) <- c("mmusculus_homolog_ensembl_gene", "mmusculus_homolog_associated_gene_name", "human_ensembl_gene_id",
                                  "t_mouse", "adj.P.Val_mouse", "P.Value_mouse", "logFC_mouse", "gene_symbol_human", "t_human", "adj.P.Val_human", 
                                  "P.Value_human", "logFC_human")
  for (i in 1:nrow(common_genes)){
    ## Find and extract info of mouse gene in mouse dataset
    mouse_data <-top_genes[which(top_genes$ensemblID==common_genes[i,1]), c("t", "adj.P.Val", "P.Value", "logFC")]
    ## Find and extract info of human gene in human dataset
    human_data <-humanGene[which(humanGene$ensemblID==common_genes[i, 3]), c("Symbol", "t", "adj.P.Val", "P.Value", "logFC")]
    human_mouse_data[i,] <- cbind(common_genes[i,], mouse_data, human_data)
  }
  
  ## Add DE and replication info of each gene
  DE<-vector()
  for (i in 1:dim(human_mouse_data)[1]) {
    
    ## Pup mouse genes with FDR<5% and human genes with p-value<5%
    if (age_mouse=="pups"){
      
      ## DEG in human (FDR<0.1) and mouse (FDR<0.05)
      if(human_mouse_data$adj.P.Val_human[i]<0.1 && human_mouse_data$adj.P.Val_mouse[i]<0.05) {
        DE<-append(DE, "Signif in both")
      }
      
      ## DEG in human 
      else if (human_mouse_data$adj.P.Val_human[i]<0.1){
        DE<-append(DE, "Signif in human (FDR<0.1)")
      }

      ## Replicating mouse genes
      else if ((human_mouse_data$adj.P.Val_mouse[i]<0.05 && human_mouse_data$P.Value_human[i]<0.05) &&
          (sign(human_mouse_data$logFC_mouse[i])==sign(human_mouse_data$logFC_human[i]))) {
        DE<-append(DE, "Replicating genes (p<0.05 in human, FDR<0.05 in mouse)")
      }
      
      ## DEG in mouse
      else if (human_mouse_data$adj.P.Val_mouse[i]<0.05){
        DE<-append(DE, "Signif in mouse (FDR<0.05)")
      }
  
      ## Non-significant genes 
      else {
        DE<-append(DE, "n.s. genes")      
      }
    }
    
    ## Adult mouse genes with p-value<5% and human genes with p-value<5%
    else {
      
      ## DEG in human (FDR<0.1) and mouse (FDR<0.05) <- there weren't DEG in adult mice
      if(human_mouse_data$adj.P.Val_human[i]<0.1 && human_mouse_data$adj.P.Val_mouse[i]<0.05) {
        DE<-append(DE, "Signif in both")
      }
      else if (human_mouse_data$adj.P.Val_human[i]<0.1){
        DE<-append(DE, "Signif in human (FDR<0.1)")
      }
      else if ((human_mouse_data$P.Value_mouse[i]<0.05 && human_mouse_data$P.Value_human[i]<0.05) &&
          (sign(human_mouse_data$logFC_mouse[i])==sign(human_mouse_data$logFC_human[i]))) {
        DE<-append(DE, "Replicating genes (p<0.05 in human and mouse)")
      }
      else {
        DE<-append(DE, "n.s. genes")      
      }
    }
  }
  
  human_mouse_data$DE<- DE
  
  ## Correlation coeff between t-stats of genes in human and mouse
  rho <- cor(human_mouse_data$t_human, human_mouse_data$t_mouse, method = "spearman")
  rho_anno = paste0("rho = ", format(round(rho, 2), nsmall = 2))
  
  ## Colors and alphas for plot
  if (age_mouse=="pups"){
    cols <- c("yellow3", "#ffad73", "red", "#26b3ff", "dark grey") 
    names(cols)<-c("Signif in both", "Signif in human (FDR<0.1)", "Replicating genes (p<0.05 in human, FDR<0.05 in mouse)","Signif in mouse (FDR<0.05)", "n.s. genes")
    alphas <- c(1, 1, 1, 1, 0.5)  
    names(alphas)<-c("Signif in both", "Signif in human (FDR<0.1)", "Replicating genes (p<0.05 in human, FDR<0.05 in mouse)","Signif in mouse (FDR<0.05)", "n.s. genes")
  }
  else {
    cols <- c("yellow3", "#ffad73", "red","dark grey") 
    names(cols)<-c("Signif in both", "Signif in human (FDR<0.1)", "Replicating genes (p<0.05 in human and mouse)","n.s. genes")
    alphas <- c(1, 1, 1, 0.5)  
    names(alphas)<-c("Signif in both", "Signif in human (FDR<0.1)", "Replicating genes (p<0.05 in human and mouse)","n.s. genes")
  }
  
  
  ## Add labels of interesting genes 

  ## 3 most significant replicating genes in human
  rep_genes <- human_mouse_data[which(human_mouse_data$DE==names(alphas)[3]),]
  top_rep_genes <- rep_genes[order(rep_genes$adj.P.Val_human),"gene_symbol_human"][1:3]
  
  label <- vector()
  for (i in 1:dim(human_mouse_data)[1]){
    ## DEG in both human and mouse
    if (human_mouse_data$DE[i]=="Signif in both"){
      ## Label of the form: [gene name in mouse]-[gene name in human]
      label <- append(label, paste(human_mouse_data$mmusculus_homolog_associated_gene_name[i], "-", human_mouse_data$gene_symbol_human[i], sep=""))
    }
    ## Labels of the top 3 replicating genes
    else if (human_mouse_data$gene_symbol_human[i] %in% top_rep_genes){
      label <- append(label, paste(human_mouse_data$mmusculus_homolog_associated_gene_name[i], "-", human_mouse_data$gene_symbol_human[i], sep=""))
    }
    else{
      label <- append(label, NA)
    }
  }
  
  human_mouse_data$label <- label
  
  ## Plot
  plot <- ggplot(human_mouse_data, aes(x = t_mouse, y = t_human, color=DE, alpha=DE, label=label)) +
    geom_point(size = 1) +
    labs(x = paste("t-stats in", substr(age_mouse, 1, nchar(age_mouse)-1), "mouse", tissue_mouse), 
         y = paste("t-stats in", age_human, "human brain"),
         title = paste(capitalize(expt_mouse),"mouse vs Smoking human", sep=" "), 
         subtitle = rho_anno, 
         parse = T) +
    geom_label_repel(fill="white", size=2, max.overlaps = Inf,  
                     box.padding = 0.2, 
                     show.legend=FALSE) +
    theme_bw() +
    scale_color_manual(values = cols) + 
    scale_alpha_manual(values = alphas)
  
  plot + theme(legend.text = element_text(size=8))
  plot
  ggsave(filename=paste("plots/04_DEA/02_Comparisons/Gene_analysis/t_stats_", age_human, "_Human_vs_Mouse_", age_mouse, "_", substr(expt_mouse,1,3), "_", 
                        tissue_mouse, ".pdf", sep=""), height = 12, width = 20, units = "cm")
  
  
  ## Quantify mouse genes that replicate in human
  
  ## Total unique DEG in pup brain (FDR<0.05)
  total_pups_DEG=length(unique(top_genes[which(top_genes$adj.P.Val<0.05),"ensemblID"]))
  ## Total unique adult brain genes with p<0.05
  total_adults_P_val_genes=length(unique(top_genes[which(top_genes$P.Val<0.05),"ensemblID"]))
  ## Unique replicating genes 
  rep_genes=length(unique(human_mouse_data[which(human_mouse_data$DE==names(alphas)[3]),"mmusculus_homolog_ensembl_gene"]))
  
  if (age_mouse=="pups"){
    ## Percentage 
    percentage=signif(rep_genes / total_pups_DEG *100, 3)
    print(paste(rep_genes, "out of", total_pups_DEG, "DEG in", expt_mouse, "mouse pup", tissue_mouse, "(FDR<0.05) replicate in smoking human", age_human, 
                "brain (with p<0.05 and same logFC direction) -", paste(percentage, "%", sep="")))
  }
  else {
    ## Percentage 
    percentage=signif(rep_genes / total_adults_P_val_genes *100, 3)
    print(paste(rep_genes, "out of", total_adults_P_val_genes, "genes in", expt_mouse, "adult mouse", tissue_mouse, "(p<0.05) replicate in smoking human",
                age_human, "brain (also p<0.05 and same logFC direction) -", paste(percentage, "%", sep="")))
  }
  
  return(human_mouse_data)
}



################################################################
##  Nicotine mouse pup brain vs Smoking human prenatal brain 
################################################################
prenatalHuman_pupNicMouse_data <- t_stat_plot_mouse_in_human(age_mouse = "pups", tissue_mouse = "brain", expt_mouse = "nicotine", age_human = "prenatal")
## "78 out of 1010 DEG in nicotine mouse pup brain (FDR<0.05) replicate in smoking human prenatal brain (with p<0.05 and same logFC direction) - 7.72%"
save(prenatalHuman_pupNicMouse_data, file="processed-data/04_DEA/Gene_analysis/prenatalHuman_pupNicMouse_data.Rdata")
## Add replication info to all mouse genes (NA if it has no human homolog or the homolog is not present in the human dataset)
## (Add the info of the first occurrence)
top_genes_pups_nicotine_fitted$replication_in_prenatalHumanBrain <- apply(top_genes_pups_nicotine_fitted, 1, function(x){prenatalHuman_pupNicMouse_data[match(x["ensemblID"], prenatalHuman_pupNicMouse_data$mmusculus_homolog_ensembl_gene), "DE"]})
save(top_genes_pups_nicotine_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata")


################################################################
##  Smoking mouse pup brain vs Smoking human prenatal brain 
################################################################
prenatalHuman_pupSmoMouse_data <- t_stat_plot_mouse_in_human(age_mouse = "pups", tissue_mouse = "brain", expt_mouse = "smoking", age_human = "prenatal")
## "267 out of 4165 DEG in smoking mouse pup brain (FDR<0.05) replicate in smoking human prenatal brain (with p<0.05 and same logFC direction) - 6.41%"
save(prenatalHuman_pupSmoMouse_data, file="processed-data/04_DEA/Gene_analysis/prenatalHuman_pupSmoMouse_data.Rdata")
top_genes_pups_smoking_fitted$replication_in_prenatalHumanBrain <- apply(top_genes_pups_smoking_fitted, 1, function(x){prenatalHuman_pupSmoMouse_data[match(x["ensemblID"], prenatalHuman_pupSmoMouse_data$mmusculus_homolog_ensembl_gene), "DE"]})
save(top_genes_pups_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata")


################################################################
##  Nicotine adult mouse brain vs Smoking human prenatal brain 
################################################################
prenatalHuman_adultNicMouse_data <- t_stat_plot_mouse_in_human(age_mouse = "adults", tissue_mouse = "brain", expt_mouse = "nicotine", age_human = "prenatal")
## "30 out of 679 genes in nicotine adult mouse brain (p<0.05) replicate in smoking human prenatal brain (also p<0.05 and same logFC direction) - 4.42%"
save(prenatalHuman_adultNicMouse_data, file="processed-data/04_DEA/Gene_analysis/prenatalHuman_adultNicMouse_data.Rdata")
top_genes_adults_nicotine_fitted$replication_in_prenatalHumanBrain <- apply(top_genes_adults_nicotine_fitted, 1, function(x){prenatalHuman_adultNicMouse_data[match(x["ensemblID"], prenatalHuman_adultNicMouse_data$mmusculus_homolog_ensembl_gene), "DE"]})
save(top_genes_adults_nicotine_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_adults_nicotine_fitted.Rdata")


################################################################
##  Smoking adult mouse brain vs Smoking human prenatal brain 
################################################################
prenatalHuman_adultSmoMouse_data <- t_stat_plot_mouse_in_human(age_mouse = "adults", tissue_mouse = "brain", expt_mouse = "smoking", age_human = "prenatal")
## "40 out of 772 genes in smoking adult mouse brain (p<0.05) replicate in smoking human prenatal brain (also p<0.05 and same logFC direction) - 5.18%"
save(prenatalHuman_adultSmoMouse_data, file="processed-data/04_DEA/Gene_analysis/prenatalHuman_adultSmoMouse_data.Rdata")
top_genes_adults_smoking_fitted$replication_in_prenatalHumanBrain <- apply(top_genes_adults_smoking_fitted, 1, function(x){prenatalHuman_adultSmoMouse_data[match(x["ensemblID"], prenatalHuman_adultSmoMouse_data$mmusculus_homolog_ensembl_gene), "DE"]})
save(top_genes_adults_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_adults_smoking_fitted.Rdata")


################################################################
##  Smoking adult mouse blood vs Smoking human prenatal brain 
################################################################
prenatalHuman_bloodMouse_data <- t_stat_plot_mouse_in_human(age_mouse = "adults", tissue_mouse = "blood", expt_mouse = "smoking", age_human = "prenatal")
## "101 out of 1499 genes in smoking adult mouse blood (p<0.05) replicate in smoking human prenatal brain (also p<0.05 and same logFC direction) - 6.74%"
save(prenatalHuman_bloodMouse_data, file="processed-data/04_DEA/Gene_analysis/prenatalHuman_bloodMouse_data.Rdata")
top_genes_blood_smoking_fitted$replication_in_prenatalHumanBrain <- apply(top_genes_blood_smoking_fitted, 1, function(x){prenatalHuman_bloodMouse_data[match(x["ensemblID"], prenatalHuman_bloodMouse_data$mmusculus_homolog_ensembl_gene), "DE"]})
save(top_genes_blood_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_blood_smoking_fitted.Rdata")


################################################################
##  Nicotine mouse pup brain vs Smoking human adult brain 
################################################################
adultHuman_pupNicMouse_data <- t_stat_plot_mouse_in_human(age_mouse = "pups", tissue_mouse = "brain", expt_mouse = "nicotine", age_human = "adult")
## "18 out of 1010 DEG in nicotine mouse pup brain (FDR<0.05) replicate in smoking human adult brain (with p<0.05 and same logFC direction) - 1.78%"
save(adultHuman_pupNicMouse_data, file="processed-data/04_DEA/Gene_analysis/adultHuman_pupNicMouse_data.Rdata")
top_genes_pups_nicotine_fitted$replication_in_adultHumanBrain <- apply(top_genes_pups_nicotine_fitted, 1, function(x){adultHuman_pupNicMouse_data[match(x["ensemblID"], adultHuman_pupNicMouse_data$mmusculus_homolog_ensembl_gene), "DE"]})
save(top_genes_pups_nicotine_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata")


################################################################
##  Smoking mouse pup brain vs Smoking human adult brain 
################################################################
adultHuman_pupSmoMouse_data<- t_stat_plot_mouse_in_human(age_mouse = "pups", tissue_mouse = "brain", expt_mouse = "smoking", age_human = "adult")
## "74 out of 4165 DEG in smoking mouse pup brain (FDR<0.05) replicate in smoking human adult brain (with p<0.05 and same logFC direction) - 1.78%"
save(adultHuman_pupSmoMouse_data, file="processed-data/04_DEA/Gene_analysis/adultHuman_pupSmoMouse_data.Rdata")
top_genes_pups_smoking_fitted$replication_in_adultHumanBrain <- apply(top_genes_pups_smoking_fitted, 1, function(x){adultHuman_pupSmoMouse_data[match(x["ensemblID"], adultHuman_pupSmoMouse_data$mmusculus_homolog_ensembl_gene), "DE"]})
save(top_genes_pups_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata")


################################################################
##  Nicotine adult mouse brain vs Smoking human adult brain 
################################################################
adultHuman_adultNicMouse_data <- t_stat_plot_mouse_in_human(age_mouse = "adults", tissue_mouse = "brain", expt_mouse = "nicotine", age_human = "adult")
## "13 out of 679 genes in nicotine adult mouse brain (p<0.05) replicate in smoking human adult brain (also p<0.05 and same logFC direction) - 1.91%"
save(adultHuman_adultNicMouse_data, file="processed-data/04_DEA/Gene_analysis/adultHuman_adultNicMouse_data.Rdata")
top_genes_adults_nicotine_fitted$replication_in_adultHumanBrain <- apply(top_genes_adults_nicotine_fitted, 1, function(x){adultHuman_adultNicMouse_data[match(x["ensemblID"], adultHuman_adultNicMouse_data$mmusculus_homolog_ensembl_gene), "DE"]})
save(top_genes_adults_nicotine_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_adults_nicotine_fitted.Rdata")


################################################################
##  Smoking adult mouse brain vs Smoking human adult brain 
################################################################
adultHuman_adultSmoMouse_data <- t_stat_plot_mouse_in_human(age_mouse = "adults", tissue_mouse = "brain", expt_mouse = "smoking", age_human = "adult")
## "9 out of 772 genes in smoking adult mouse brain (p<0.05) replicate in smoking human adult brain (also p<0.05 and same logFC direction) - 1.17%"
save(adultHuman_adultSmoMouse_data, file="processed-data/04_DEA/Gene_analysis/adultHuman_adultSmoMouse_data.Rdata")
top_genes_adults_smoking_fitted$replication_in_adultHumanBrain <- apply(top_genes_adults_smoking_fitted, 1, function(x){adultHuman_adultSmoMouse_data[match(x["ensemblID"], adultHuman_adultSmoMouse_data$mmusculus_homolog_ensembl_gene), "DE"]})
save(top_genes_adults_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_adults_smoking_fitted.Rdata")


################################################################
##  Smoking adult mouse blood vs Smoking human adult brain 
################################################################
adultHuman_bloodMouse_data <- t_stat_plot_mouse_in_human(age_mouse = "adults", tissue_mouse = "blood", expt_mouse = "smoking", age_human = "adult")
## "16 out of 1499 genes in smoking adult mouse blood (p<0.05) replicate in smoking human adult brain (also p<0.05 and same logFC direction) - 1.07%"
save(adultHuman_bloodMouse_data, file="processed-data/04_DEA/Gene_analysis/adultHuman_bloodMouse_data.Rdata")
top_genes_blood_smoking_fitted$replication_in_adultHumanBrain <- apply(top_genes_blood_smoking_fitted, 1, function(x){adultHuman_bloodMouse_data[match(x["ensemblID"], adultHuman_bloodMouse_data$mmusculus_homolog_ensembl_gene), "DE"]})
save(top_genes_blood_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_blood_smoking_fitted.Rdata")




################### 5.2 Human genes that replicate in mouse #################### 

## Obtain human brain genes that replicate (FDR<10%) in mouse blood or brain (with p-value<5% and same logFC sign)
replication_human_in_mouse<- function(age_mouse, expt_mouse, tissue_mouse, age_human){
  
  ## Define mouse dataset
  if (tissue_mouse=="blood"){
    top_genes <-eval(parse_expr(paste("top_genes", tissue_mouse, expt_mouse, "fitted", sep="_")))
    signif_measure_mouse="p-value"
  }
  else {
    top_genes <-eval(parse_expr(paste("top_genes", age_mouse, expt_mouse, "fitted", sep="_")))
  }
  
  ## Define human dataset
  if (age_human=="prenatal"){
    humanGene <-fetalGene
  }
  else {
    humanGene <-adultGene
  }
  
  ## Extract human DEG (with FDR<0.1)
  de_genes_human <- humanGene[which(humanGene$adj.P.Val<0.1),]
  colnames(de_genes_human) <- paste(colnames(de_genes_human), "human", sep="_")
  
  ## Add the IDs of the homologous genes in mouse
  mouse_ensemblIDs <- vector()
  mouse_gene_names <- vector()
  for (i in 1:dim(de_genes_human)[1]){
    mouse_ensemblID<- common_genes[which(common_genes$human_ensembl_gene_id==de_genes_human$ensemblID_human[i]),"mmusculus_homolog_ensembl_gene"]
    mouse_gene_name<- common_genes[which(common_genes$human_ensembl_gene_id==de_genes_human$ensemblID_human[i]),"mmusculus_homolog_associated_gene_name"]
    mouse_ensemblIDs <- append(mouse_ensemblIDs, mouse_ensemblID)
    mouse_gene_names <- append(mouse_gene_names, mouse_gene_name)
  }
  de_genes_human <- cbind(de_genes_human, "ensemblID_mouse"=mouse_ensemblIDs, "gene_name_mouse"=mouse_gene_names)
  
  ## Mouse data of human DEG homologs
  mouse_data <- data.frame(matrix(ncol = 4, nrow = nrow(de_genes_human)))
  colnames(mouse_data) <- c("t_mouse", "adj.P.Val_mouse", "P.Value_mouse", "logFC_mouse")
  for (i in 1:nrow(de_genes_human)){
    ## Extract info of mouse genes
    mouse_data[i,] <-top_genes[which(top_genes$ensemblID==de_genes_human[i,"ensemblID_mouse"]), c("t", "adj.P.Val", "P.Value", "logFC")]
  }
  
  human_mouse_data <- cbind(de_genes_human, mouse_data)
  
  ## Add replication info of each gene
  replication<-vector()
  for (i in 1:dim(human_mouse_data)[1]) {
    
    ## Human genes with FDR<10% and mouse genes with p-value<5%
    if (human_mouse_data$P.Value_mouse[i]<0.05 && (sign(human_mouse_data$logFC_mouse[i])==sign(human_mouse_data$logFC_human[i]))) {
      replication<-append(replication, "Replicating gene (FDR<0.1 in human, p<0.05 in mouse)")
    }
      
    ## Non-replicating genes
    else {
      replication<-append(replication, "Non-replicating gene")      
    }
  }
  
  human_mouse_data$replication_in_mouse <- replication
  
  
  ## Quantify human genes that replicate in mouse
  
  ## Total DEG in human brain (FDR<0.1)
  total_human_DEG=nrow(de_genes_human)
  rep_genes_ids <- human_mouse_data[which(human_mouse_data$replication_in_mouse=="Replicating gene (FDR<0.1 in human, p<0.05 in mouse)"),"Symbol_human"]
  ## Unique replicating genes 
  rep_genes=length(rep_genes_ids)
  ## Percentage 
  percentage=signif(rep_genes / total_human_DEG *100, 3)
  print(paste(rep_genes, "out of", total_human_DEG, "DEG in smoking human", age_human, "brain (FDR<0.1) replicate in", expt_mouse, substr(age_mouse, 1, nchar(age_mouse)-1), 
         "mouse", tissue_mouse, "(with p<0.05 and same logFC direction) -", paste(percentage, "%. Genes:", sep="")))
  print(rep_genes_ids)

  return(human_mouse_data)
}



## Plots
## Create plots to verify which human genes replicate in mouse 
t_stat_plot_human_in_mouse <- function(age_mouse, expt_mouse, tissue_mouse, age_human){
  
  ## Define mouse dataset
  if (tissue_mouse=="blood"){
    top_genes <-eval(parse_expr(paste("top_genes", tissue_mouse, expt_mouse, "fitted", sep="_")))
  }
  else {
    top_genes <-eval(parse_expr(paste("top_genes", age_mouse, expt_mouse, "fitted", sep="_")))
  }
  
  ## Define human dataset
  if (age_human=="prenatal"){
    humanGene <-fetalGene
  }
  else {
    humanGene <-adultGene
  }
  
  ## Extract mouse and human data of common genes
  human_mouse_data <- data.frame(matrix(nrow = nrow(common_genes), ncol = 12))
  colnames(human_mouse_data) <- c("mmusculus_homolog_ensembl_gene", "mmusculus_homolog_associated_gene_name", "human_ensembl_gene_id",
                                  "t_mouse", "adj.P.Val_mouse", "P.Value_mouse", "logFC_mouse", "gene_symbol_human", "t_human", "adj.P.Val_human", 
                                  "P.Value_human", "logFC_human")
  for (i in 1:nrow(common_genes)){
    ## Find and extract info of mouse gene in mouse dataset
    mouse_data <-top_genes[which(top_genes$ensemblID==common_genes[i,1]), c("t", "adj.P.Val", "P.Value", "logFC")]
    ## Find and extract info of human gene in human dataset
    human_data <-humanGene[which(humanGene$ensemblID==common_genes[i, 3]), c("Symbol", "t", "adj.P.Val", "P.Value", "logFC")]
    human_mouse_data[i,] <- cbind(common_genes[i,], mouse_data, human_data)
  }
  
  ## Add DE and replication info of each gene
  DE<-vector()
  for (i in 1:dim(human_mouse_data)[1]) {
    
    ## DEG in human (FDR<0.1) and mouse (FDR<0.05)
    if(human_mouse_data$adj.P.Val_human[i]<0.1 && human_mouse_data$adj.P.Val_mouse[i]<0.05) {
      DE<-append(DE, "Signif in both")
    }
    
    ## Replicating human genes (with FDR<10% and p-value<5% in mouse)
    else if ((human_mouse_data$adj.P.Val_human[i]<0.1 && human_mouse_data$P.Value_mouse[i]<0.05) &&
             (sign(human_mouse_data$logFC_mouse[i])==sign(human_mouse_data$logFC_human[i]))) {
      DE<-append(DE, "Replicating genes (FDR<0.1 in human, p<0.05 in mouse)")
    }
    
    ## DEG in human 
    else if (human_mouse_data$adj.P.Val_human[i]<0.1){
      DE<-append(DE, "Signif in human (FDR<0.1)")
    }
    
    ## DEG in mouse
    else if (human_mouse_data$adj.P.Val_mouse[i]<0.05){
      DE<-append(DE, "Signif in mouse (FDR<0.05)")
    }
    
    ## Non-significant genes 
    else {
      DE<-append(DE, "n.s. genes")      
    }
  }
  
  human_mouse_data$DE<- DE
  
  ## Correlation coeff between t-stats of genes in human and mouse
  rho <- cor(human_mouse_data$t_human, human_mouse_data$t_mouse, method = "spearman")
  rho_anno = paste0("rho = ", format(round(rho, 2), nsmall = 2))
  
  ## Colors and alphas for plot
  cols <- c("yellow3", "red", "#ffad73", "#26b3ff", "dark grey") 
  names(cols)<-c("Signif in both", "Replicating genes (FDR<0.1 in human, p<0.05 in mouse)", "Signif in human (FDR<0.1)", 
                 "Signif in mouse (FDR<0.05)", "n.s. genes")
  alphas <- c(1, 1, 1, 1, 0.5)  
  names(alphas)<-c("Signif in both", "Replicating genes (FDR<0.1 in human, p<0.05 in mouse)", "Signif in human (FDR<0.1)", 
                   "Signif in mouse (FDR<0.05)", "n.s. genes")
  
  
  ## Add labels of replicating human genes 
  
  rep_genes <- human_mouse_data[which(human_mouse_data$DE==names(alphas)[2]), "gene_symbol_human"]
  
  label <- vector()
  for (i in 1:dim(human_mouse_data)[1]){
    ## DEG in both human and mouse
    if (human_mouse_data$DE[i]=="Signif in both"){
      ## Label of the form: [gene name in mouse]-[gene name in human]
      label <- append(label, paste(human_mouse_data$mmusculus_homolog_associated_gene_name[i], "-", human_mouse_data$gene_symbol_human[i], sep=""))
    }
    ## Labels of replicating genes
    else if (human_mouse_data$gene_symbol_human[i] %in% rep_genes){
      label <- append(label, paste(human_mouse_data$mmusculus_homolog_associated_gene_name[i], "-", human_mouse_data$gene_symbol_human[i], sep=""))
    }
    else{
      label <- append(label, NA)
    }
  }
  
  human_mouse_data$label <- label
  
  ## Plot
  plot <- ggplot(human_mouse_data, aes(x = t_human, y = t_mouse, color=DE, alpha=DE, label=label)) +
    geom_point(size = 1) +
    labs(x = paste("t-stats in", age_human, "human brain"), 
         y = paste("t-stats in", substr(age_mouse, 1, nchar(age_mouse)-1), "mouse", tissue_mouse),
         title = paste("Smoking human vs", capitalize(expt_mouse), "mouse", sep=" "), 
         subtitle = rho_anno, 
         parse = T) +
    geom_label_repel(fill="white", size=2, max.overlaps = Inf,  
                     box.padding = 0.2, 
                     show.legend=FALSE) +
    theme_bw() +
    scale_color_manual(values = cols) + 
    scale_alpha_manual(values = alphas)
  
  plot + theme(legend.text = element_text(size=8))
  plot
  ggsave(filename=paste("plots/04_DEA/02_Comparisons/Gene_analysis/t_stats_replication_", age_human, "_Human_in_Mouse_", age_mouse, "_", substr(expt_mouse,1,3), "_", 
                        tissue_mouse, ".pdf", sep=""), height = 12, width = 20, units = "cm")
  
}


################################################################
##  Smoking human prenatal brain vs Nicotine mouse pup brain   
################################################################
t_stat_plot_human_in_mouse(age_mouse = "pups", tissue_mouse = "brain", expt_mouse = "nicotine", age_human = "prenatal")
prenatalHuman_in_NicPupMouse <- replication_human_in_mouse(age_mouse = "pups", tissue_mouse = "brain", expt_mouse = "nicotine", age_human = "prenatal")
## "2 out of 13 DEG in smoking human prenatal brain (FDR<0.1) replicate in nicotine pup mouse brain (with p<0.05 and same logFC direction) - 15.4%. Genes:"
## [1] "MPPED1" "SDC1" 

################################################################
##  Smoking human prenatal brain vs Smoking mouse pup brain   
################################################################
t_stat_plot_human_in_mouse(age_mouse = "pups", tissue_mouse = "brain", expt_mouse = "smoking", age_human = "prenatal")
prenatalHuman_in_SmoPupMouse <- replication_human_in_mouse(age_mouse = "pups", tissue_mouse = "brain", expt_mouse = "smoking", age_human = "prenatal")
## "1 out of 13 DEG in smoking human prenatal brain (FDR<0.1) replicate in smoking pup mouse brain (with p<0.05 and same logFC direction) - 7.69%. Genes:"
## [1] "NRCAM"

################################################################
##  Smoking human prenatal brain vs Nicotine adult mouse brain   
################################################################
t_stat_plot_human_in_mouse(age_mouse = "adults", tissue_mouse = "brain", expt_mouse = "nicotine", age_human = "prenatal")
prenatalHuman_in_NicAdultMouse <- replication_human_in_mouse(age_mouse = "adults", tissue_mouse = "brain", expt_mouse = "nicotine", age_human = "prenatal")
## "0 out of 13 DEG in smoking human prenatal brain (FDR<0.1) replicate in nicotine adult mouse brain (with p<0.05 and same logFC direction) - 0%. Genes:"
## character(0)

################################################################
##  Smoking human prenatal brain vs Smoking adult mouse brain   
################################################################
t_stat_plot_human_in_mouse(age_mouse = "adults", tissue_mouse = "brain", expt_mouse = "smoking", age_human = "prenatal")
prenatalHuman_in_SmoAdultMouse <- replication_human_in_mouse(age_mouse = "adults", tissue_mouse = "brain", expt_mouse = "smoking", age_human = "prenatal")
## "0 out of 13 DEG in smoking human prenatal brain (FDR<0.1) replicate in smoking adult mouse brain (with p<0.05 and same logFC direction) - 0%. Genes:"
## character(0)

################################################################
##  Smoking human prenatal brain vs Smoking adult mouse blood   
################################################################
t_stat_plot_human_in_mouse(age_mouse = "adults", tissue_mouse = "blood", expt_mouse = "smoking", age_human = "prenatal")
prenatalHuman_in_BloodMouse<- replication_human_in_mouse(age_mouse = "adults", tissue_mouse = "blood", expt_mouse = "smoking", age_human = "prenatal")
## "1 out of 13 DEG in smoking human prenatal brain (FDR<0.1) replicate in smoking adult mouse blood (with p<0.05 and same logFC direction) - 7.69%. Genes:"
## [1] "KCNN2"

################################################################
##  Smoking human adult brain vs Nicotine mouse pup brain   
################################################################
t_stat_plot_human_in_mouse(age_mouse = "pups", tissue_mouse = "brain", expt_mouse = "nicotine", age_human = "adult")
adultHuman_in_NicPupMouse <- replication_human_in_mouse(age_mouse = "pups", tissue_mouse = "brain", expt_mouse = "nicotine", age_human = "adult")
## "1 out of 1 DEG in smoking human adult brain (FDR<0.1) replicate in nicotine pup mouse brain (with p<0.05 and same logFC direction) - 100%. Genes:"
## [1] "MARCO"

################################################################
##  Smoking human adult brain vs Smoking mouse pup brain   
################################################################
t_stat_plot_human_in_mouse(age_mouse = "pups", tissue_mouse = "brain", expt_mouse = "smoking", age_human = "adult")
adultHuman_in_SmoPupMouse <- replication_human_in_mouse(age_mouse = "pups", tissue_mouse = "brain", expt_mouse = "smoking", age_human = "adult")
## "0 out of 1 DEG in smoking human adult brain (FDR<0.1) replicate in smoking pup mouse brain (with p<0.05 and same logFC direction) - 0%. Genes:"
## character(0)

################################################################
##  Smoking human adult brain vs Nicotine adult mouse brain   
################################################################
t_stat_plot_human_in_mouse(age_mouse = "adults", tissue_mouse = "brain", expt_mouse = "nicotine", age_human = "adult")
adultHuman_in_NicAdultMouse <- replication_human_in_mouse(age_mouse = "adults", tissue_mouse = "brain", expt_mouse = "nicotine", age_human = "adult")
## "0 out of 1 DEG in smoking human adult brain (FDR<0.1) replicate in nicotine adult mouse brain (with p<0.05 and same logFC direction) - 0%. Genes:"
## character(0)

################################################################
##  Smoking human adult brain vs Smoking adult mouse brain   
################################################################
t_stat_plot_human_in_mouse(age_mouse = "adults", tissue_mouse = "brain", expt_mouse = "smoking", age_human = "adult")
adultHuman_in_SmoAdultMouse <- replication_human_in_mouse(age_mouse = "adults", tissue_mouse = "brain", expt_mouse = "smoking", age_human = "adult")
## "0 out of 1 DEG in smoking human adult brain (FDR<0.1) replicate in smoking adult mouse brain (with p<0.05 and same logFC direction) - 0%. Genes:"
## character(0)

################################################################
##  Smoking human adult brain vs Smoking adult mouse blood   
################################################################
t_stat_plot_human_in_mouse(age_mouse = "adults", tissue_mouse = "blood", expt_mouse = "smoking", age_human = "adult")
adultHuman_in_BloodMouse<- replication_human_in_mouse(age_mouse = "adults", tissue_mouse = "blood", expt_mouse = "smoking", age_human = "adult")
## "0 out of 1 DEG in smoking human adult brain (FDR<0.1) replicate in smoking adult mouse blood (with p<0.05 and same logFC direction) - 0%. Genes:"
## character(0)


## Create table with the results 

human_genes <- rbind(cbind(fetalGene[which(fetalGene$adj.P.Val<0.1), c("Symbol", "ensemblID")], Age="Prenatal"), 
                     cbind(adultGene[which(adultGene$adj.P.Val<0.1), c("Symbol", "ensemblID")], Age="Adult"))
human_genes_in_mouse <- cbind(human_genes, Replicate_in_NicPupMouseBrain=c("No","No","No","No","Yes","No","No","Yes","No","No","No","No","No","Yes"))
human_genes_in_mouse <- cbind(human_genes_in_mouse, Replicate_in_SmoPupMouseBrain=c("No","No","No","No","No","No","No","No","No","No","No","No","Yes","No"))
human_genes_in_mouse <- cbind(human_genes_in_mouse, Replicate_in_NicAdultMouseBrain=c("No","No","No","No","No","No","No","No","No","No","No","No","No","No"))
human_genes_in_mouse <- cbind(human_genes_in_mouse, Replicate_in_SmoAdultMouseBrain=c("No","No","No","No","No","No","No","No","No","No","No","No","No","No"))
human_genes_in_mouse <- cbind(human_genes_in_mouse, Replicate_in_SmoAdultMouseBlood=c("Yes","No","No","No","No","No","No","No","No","No","No","No","No","No"))
save(human_genes_in_mouse, file="processed-data/04_DEA/Gene_analysis/human_genes_in_mouse.Rdata")







### 1.2.3 Venn diagrams

## Function to create multiple Venn diagrams
venn_plot<-function(DEG_lists, colors, filename){
  
  plots<-list()
  pdf(file = paste("plots/04_DEA/02_Comparisons/Gene_analysis/Venn_", filename, ".pdf", sep=""))
  for (i in 1:length(DEG_lists)){
     v<-venn.diagram(DEG_lists[[i]], fill=colors[[i]], alpha = rep(0.5, length(DEG_lists[[i]])), 
                     lwd =0, margin=0.2, cat.cex=0.6, cex=0.6, height = 35, width = 35, units = "cm", 
                     cat.dist=rep(0.09, length(DEG_lists[[i]])),filename=NULL)
     plots[[i]]<-v
  }
  
  if (i==4){
    gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol=2)
    dev.off()
  }
  else if (i==3){
    gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], ncol=2)
    dev.off()
  }
  else {
    gridExtra::grid.arrange(plots[[1]], plots[[2]], ncol=2)
    dev.off()
  }
}



################################################################################
##                   Compare naive VS fitted models DEG 
################################################################################

## Nicotine DEG
DEG_nic_naive_vs_fitted <-list(
  "Nicotine naive"=de_genes_pups_nicotine_naive$Symbol,
  "Nicotine fitted"=de_genes_pups_nicotine_fitted$Symbol)
## Smoking DEG
DEG_smo_naive_vs_fitted <-list(
  "Smoking naive"=de_genes_pups_smoking_naive$Symbol,
  "Smoking fitted"=de_genes_pups_smoking_fitted$Symbol)

## Compare smoking VS nicotine DEG
## Naive model DEG
DEG_naive_smo_vs_nic <-list(
  "Nicotine naive"=de_genes_pups_nicotine_naive$Symbol,
  "Smoking naive"=de_genes_pups_smoking_naive$Symbol)
## Fitted model DEG
DEG_fitted_smo_vs_nic <-list(
  "Nicotine fitted"=de_genes_pups_nicotine_fitted$Symbol,
  "Smoking fitted"=de_genes_pups_smoking_fitted$Symbol)


DEG_lists<-list(DEG_nic_naive_vs_fitted, DEG_smo_naive_vs_fitted, 
                DEG_naive_smo_vs_nic, DEG_fitted_smo_vs_nic  )
colors<-list(c("palegreen2", "yellow3"), c("lightsalmon", "slategray3"), 
             c("palegreen2", "lightsalmon"), c("yellow3", "slategray3"))
venn_plot(DEG_lists, colors, "model_smo_and_nic")



################################################################################
##                Venn diagrams for Up and Down regulated DEG
################################################################################

## Genes of each group
nic_naive_up<-de_genes_pups_nicotine_naive[which(de_genes_pups_nicotine_naive$logFC>0),"Symbol"]
nic_naive_down<-de_genes_pups_nicotine_naive[which(de_genes_pups_nicotine_naive$logFC<0),"Symbol"]
nic_fitted_up<-de_genes_pups_nicotine_fitted[which(de_genes_pups_nicotine_fitted$logFC>0),"Symbol"]
nic_fitted_down<-de_genes_pups_nicotine_fitted[which(de_genes_pups_nicotine_fitted$logFC<0),"Symbol"]
smo_naive_up<-de_genes_pups_smoking_naive[which(de_genes_pups_smoking_naive$logFC>0),"Symbol"]
smo_naive_down<-de_genes_pups_smoking_naive[which(de_genes_pups_smoking_naive$logFC<0),"Symbol"]
smo_fitted_up<-de_genes_pups_smoking_fitted[which(de_genes_pups_smoking_fitted$logFC>0),"Symbol"]
smo_fitted_down<-de_genes_pups_smoking_fitted[which(de_genes_pups_smoking_fitted$logFC<0),"Symbol"]


## Compare Up and down Regulated nicotine/smoking genes 
## Naive VS fitted models

## Nicotine DEG
DEG_nic_up<-list(
  "Nicotine naive up"=nic_naive_up,
  "Nicotine fitted up"=nic_fitted_up)
DEG_nic_down<-list(
  "Nicotine naive down"=nic_naive_down,
  "Nicotine fitted down"=nic_fitted_down)

## Smoking DEG
DEG_smo_up<-list(
  "Smoking naive up"=smo_naive_up,
  "Smoking fitted up"=smo_fitted_up)
DEG_smo_down<-list(
  "Smoking naive down"=smo_naive_down,
  "Smoking fitted down"=smo_fitted_down)

DEG_lists<-list(DEG_nic_up, DEG_nic_down, DEG_smo_up, DEG_smo_down)
colors<-list(c("firebrick3", "brown1"), c("dodgerblue2", "deepskyblue1"), 
             c("coral1", "darksalmon"), c("darkturquoise", "darkslategray2"))
venn_plot(DEG_lists, colors, "nic_and_smo_Up_and_Down")


## Compare smoking VS nicotine Up and Down regulated genes
## For genes from either naive or fitted model

DEG_smo_vs_nic_Up<-list(
  "Smoking up"=union(smo_fitted_up, smo_naive_up),
  "Nicotine up"=union(nic_fitted_up, nic_naive_up))

DEG_smo_vs_nic_Down<-list(
  "Smoking down"=union(smo_fitted_down, smo_naive_down),
  "Nicotine down"=union(nic_fitted_down, nic_naive_down))

DEG_smoUp_vs_nicDown<-list(
  "Smoking up"=union(smo_fitted_up, smo_naive_up),
  "Nicotine down"=union(nic_fitted_down, nic_naive_down))

DEG_smoDown_vs_nicUp<-list(
  "Smoking down"=union(smo_fitted_down, smo_naive_down),
  "Nicotine up"=union(nic_fitted_up, nic_naive_up))

DEG_lists<-list(DEG_smo_vs_nic_Up, DEG_smo_vs_nic_Down, DEG_smoUp_vs_nicDown, DEG_smoDown_vs_nicUp)
colors<-list(c("coral", "brown3"), c("cyan4", "cyan3"), c("coral", "cyan3"), c("cyan4", "brown3"))
venn_plot(DEG_lists, colors, "smo_VS_nic_Up_and_Down")


## Compare smoking VS nicotine Up and Down regulated genes 
## For naive/fitted genes only

## Naive model DEG
DEG_naive_smo_vs_nic_up<-list(
  "Smoking naive up"=smo_naive_up,
  "Nicotine naive up"=nic_naive_up)
DEG_naive_smo_vs_nic_down<-list(
  "Smoking naive down"=smo_naive_down,
  "Nicotine naive down"=nic_naive_down)
DEG_naive_smoUp_nicDown<-list(
  "Smoking naive up"=smo_naive_up,
  "Nicotine naive down"=nic_naive_down)
DEG_naive_smoDown_nicUp<-list(
  "Smoking naive down"=smo_naive_down,
  "Nicotine naive up"=nic_naive_up)

DEG_lists<-list(DEG_naive_smo_vs_nic_up, DEG_naive_smo_vs_nic_down, 
                DEG_naive_smoUp_nicDown, DEG_naive_smoDown_nicUp)
colors=list(c("coral1", "firebrick3"), c("darkturquoise", "dodgerblue2"), 
            c("coral1", "dodgerblue2"), c("darkturquoise", "firebrick3"))
venn_plot(DEG_lists, colors, "Naive_smo_VS_nic_Up_and_Down")


## Fitted model DEG
DEG_fitted_smo_vs_nic_up<-list(
  "Smoking fitted up"=smo_fitted_up,
  "Nicotine fitted up"=nic_fitted_up)
save(DEG_fitted_smo_vs_nic_up, file="processed-data/04_DEA/Gene_analysis/DEG_fitted_smo_vs_nic_up.Rdata")

DEG_fitted_smo_vs_nic_down<-list(
  "Smoking fitted down"=smo_fitted_down,
  "Nicotine fitted down"=nic_fitted_down)
save(DEG_fitted_smo_vs_nic_down, file="processed-data/04_DEA/Gene_analysis/DEG_fitted_smo_vs_nic_down.Rdata")

DEG_fitted_smoUp_nicDown<-list(
  "Smoking fitted up"=smo_fitted_up,
  "Nicotine fitted down"=nic_fitted_down)
save(DEG_fitted_smoUp_nicDown, file = "processed-data/04_DEA/Gene_analysis/DEG_fitted_smoUp_nicDown.Rdata")

DEG_fitted_smoDown_nicUp<-list(
  "Smoking fitted down"=smo_fitted_down,
  "Nicotine fitted up"=nic_fitted_up)
save(DEG_fitted_smoDown_nicUp, file="processed-data/04_DEA/Gene_analysis/DEG_fitted_smoDown_nicUp.Rdata")

DEG_lists<-list(DEG_fitted_smo_vs_nic_up, DEG_fitted_smo_vs_nic_down, 
                DEG_fitted_smoUp_nicDown, DEG_fitted_smoDown_nicUp)
colors=list(c("darksalmon", "brown1"), c("darkslategray2", "deepskyblue1"),
            c("darksalmon", "deepskyblue1"), c("darkslategray2", "brown1"))
venn_plot(DEG_lists, colors, "Fitted_smo_VS_nic_Up_and_Down")



################################################################################
##                Compare all 4 groups of DEG (by expt and model)
################################################################################

DEG_all<-list(
  "Nicotine naive"=de_genes_pups_nicotine_naive$Symbol,
  "Nicotine fitted"=de_genes_pups_nicotine_fitted$Symbol,
  "Smoking naive"=de_genes_pups_smoking_naive$Symbol,
  "Smoking fitted"=de_genes_pups_smoking_fitted$Symbol)

## Compare DEG from the fitted model
## Compare all Up/Down regulated genes
DEG_all_Up_and_Down<-list(
  "Nicotine up"=nic_fitted_up,
  "Nicotine down"=nic_fitted_down,
  "Smoking up"=smo_fitted_up, 
  "Smoking down"=smo_fitted_down)

## Compare all smoking VS all nicotine DEG
DEG_all_smo_vs_all_nic<-list(
  "Nicotine"=union(de_genes_pups_nicotine_naive$Symbol, de_genes_pups_nicotine_fitted$Symbol),
  "Smoking"=union(de_genes_pups_smoking_naive$Symbol, de_genes_pups_smoking_fitted$Symbol))


DEG_lists<-list(DEG_all, DEG_all_Up_and_Down, DEG_all_smo_vs_all_nic)
colors<-list(c("palegreen2", "yellow3", "lightsalmon", "slategray3"), 
             c("brown3", "cyan3", "coral", "cyan4"), c("olivedrab1", "rosybrown2"))
venn_plot(DEG_lists, colors, "all_DEG")



################################################################################
##          Venn diagrams of mouse replicating genes in human
################################################################################

## Note: each circle represents the number of unique gene pairs (mouse-human homologs), with p<0.05 in adult mouse / FDR<0.05 in mouse pups, 
## p<0.05 in prenatal/adult human brain and the same logFC sign in both species. Therefore, the intersection contains the replicating genes. 
## (Consider that some human genes could be associated with >1 mouse homolog gene and the other way around)

## Compare genes from different mouse and human datasets

venn_human_vs_mouse <- function(dataset, age_mouse, expt_mouse, tissue_mouse, age_human){
  
  data <- eval(parse_expr(dataset))
  
  if (age_mouse=="pup"){
    signif_measure <- "adj.P.Val_mouse"
    signif_measure_name <- "FDR"
  }
  else{
    signif_measure <- "P.Value_mouse"
    signif_measure_name <- "p"
  }

  
  ## 1. Up mouse and human genes
  
  ## Extract ensembl IDs of interest genes in mouse (FDR<0.05 in pups and p<0.05 in adults, with logFC>0) and their counterparts in human
  up_interest_genes_mouse <- data[which(eval(parse_expr(paste(dataset, signif_measure, sep="$")))<0.05 & data$logFC_mouse>0), 
                                 c("mmusculus_homolog_ensembl_gene", "human_ensembl_gene_id")]
  up_interest_genes_mouse_human <- apply(up_interest_genes_mouse, 1, function(x){paste(x[1], x[2], sep="-")})
  
  ## Extract ensembl IDs of interest genes in human (p<0.05 and logFC>0) and their counterparts in mouse
  up_interest_genes_human <- data[which(data$P.Value_human<0.05 & data$logFC_human>0), c("mmusculus_homolog_ensembl_gene", "human_ensembl_gene_id")]
  up_interest_genes_human_mouse <- apply(up_interest_genes_human, 1, function(x){paste(x[1], x[2], sep="-")})
  
  ## Define category names of the diagram
  mouse_cat_name <- paste("Up in", capitalize(expt_mouse), age_mouse, tissue_mouse, paste("(", signif_measure_name, "<0.05)", sep=""))
  human_cat_name <- paste("Up in Smoking", age_human, "human brain (p<0.05)")
  
  up_mouse_human_genes<-list(up_interest_genes_mouse_human, up_interest_genes_human_mouse)
  names(up_mouse_human_genes) <- c(mouse_cat_name, human_cat_name)
  
  
  ## 2. Down mouse and human genes
  down_interest_genes_mouse <- data[which(eval(parse_expr(paste(dataset, signif_measure, sep="$")))<0.05 & data$logFC_mouse<0), 
                                  c("mmusculus_homolog_ensembl_gene", "human_ensembl_gene_id")]
  down_interest_genes_mouse_human <- apply(down_interest_genes_mouse, 1, function(x){paste(x[1], x[2], sep="-")})
  
  down_interest_genes_human <- data[which(data$P.Value_human<0.05 & data$logFC_human<0), c("mmusculus_homolog_ensembl_gene", "human_ensembl_gene_id")]
  down_interest_genes_human_mouse <- apply(down_interest_genes_human, 1, function(x){paste(x[1], x[2], sep="-")})
  
  ## Define category names of the diagram
  mouse_cat_name <- paste("Down in", capitalize(expt_mouse), age_mouse, tissue_mouse, paste("(", signif_measure_name, "<0.05)", sep=""))
  human_cat_name <- paste("Down in Smoking", age_human, "human brain (p<0.05)")
  
  down_mouse_human_genes<-list(down_interest_genes_mouse_human, down_interest_genes_human_mouse)
  names(down_mouse_human_genes) <- c(mouse_cat_name, human_cat_name)
  
  gene_pairs_list<-list(up_mouse_human_genes, down_mouse_human_genes)
  colors<-list(c("coral2", "firebrick2"), c("dodgerblue2", "deepskyblue3" ))
  
  ## Diagrams
  plots<-list()
  pdf(file = paste("plots/04_DEA/02_Comparisons/Gene_analysis/Venn_", dataset, ".pdf", sep=""), height = 8, width = 15)
  for (i in 1:length(gene_pairs_list)){
    v<-venn.diagram(gene_pairs_list[[i]], fill=colors[[i]], alpha = rep(0.5, length(gene_pairs_list[[i]])), 
                    lwd =0, margin=0.3, cat.cex=0.9, cex=1, height = 20, width = 50, units = "cm", 
                    cat.dist=rep(0.15, length(gene_pairs_list[[i]])),filename=NULL)
    plots[[i]]<-v
  }
  gridExtra::grid.arrange(plots[[1]], plots[[2]], ncol=2)
  dev.off()
}


## Plots

## Nicotine mouse pup brain vs Smoking human prenatal brain 
venn_human_vs_mouse("prenatalHuman_pupNicMouse_data", "pup", "nicotine", "brain", "prenatal")

## Smoking mouse pup brain vs Smoking human prenatal brain 
venn_human_vs_mouse("prenatalHuman_pupSmoMouse_data", "pup", "smoking", "brain", "prenatal")

## Nicotine adult mouse brain vs Smoking human prenatal brain 
venn_human_vs_mouse("prenatalHuman_adultNicMouse_data", "adult", "nicotine", "brain", "prenatal")

## Smoking adult mouse brain vs Smoking human prenatal brain
venn_human_vs_mouse("prenatalHuman_adultSmoMouse_data", "adult", "smoking", "brain", "prenatal")

## Smoking adult mouse blood vs Smoking human prenatal brain
venn_human_vs_mouse("prenatalHuman_bloodMouse_data", "adult", "smoking", "blood", "prenatal")

## Nicotine mouse pup brain vs Smoking human adult brain 
venn_human_vs_mouse("adultHuman_pupNicMouse_data", "pup", "nicotine", "brain", "adult")

## Smoking mouse pup brain vs Smoking human adult brain 
venn_human_vs_mouse("adultHuman_pupSmoMouse_data", "pup", "smoking", "brain", "adult")

## Nicotine adult mouse brain vs Smoking human adult brain 
venn_human_vs_mouse("adultHuman_adultNicMouse_data", "adult", "nicotine", "brain", "adult")

## Smoking adult mouse brain vs Smoking human adult brain
venn_human_vs_mouse("adultHuman_adultSmoMouse_data", "adult", "smoking", "brain", "adult")

## Smoking adult mouse blood vs Smoking human adult brain
venn_human_vs_mouse("adultHuman_bloodMouse_data", "adult", "smoking", "blood", "adult")










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
# date     2022-09-22
# rstudio  2022.02.3+492 Prairie Trillium (desktop)
# pandoc   2.17.1.1 @ C:/Program Files/RStudio/bin/quarto/bin/ (via rmarkdown)






