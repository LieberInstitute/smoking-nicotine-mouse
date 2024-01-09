
## All libraries used for the whole DEA

library(here)
library(SummarizedExperiment)
library(stats)
library(edgeR)
library(limma)
library(ggplot2)
library(rlang)
library(cowplot)
library(ggrepel)
library(jaffelab)
library(VennDiagram) 
library(gridExtra)
library(biomartr)
library(Hmisc)
library(readxl)
library(sessioninfo)


# 1. Differential Expression Analysis at the gene level

load(here("raw-data/rse_gene_smoking_mouse_n208.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_gene_blood_qc.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_adults_nicotine.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_adults_smoking.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_nicotine.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_smoking.Rdata"))


# ------------------------------------------------------------------------------
## 1.1 Modeling
# ------------------------------------------------------------------------------

## Extract previous output from calcNormFactors for all samples
norm_factors<-calcNormFactors(rse_gene, method = "TMM")
samples_factors<-data.frame(SAMPLE_ID=norm_factors$samples$SAMPLE_ID,
                            norm.factors=norm_factors$samples$norm.factors,
                            lib.size=norm_factors$samples$lib.size)


## DEA for experimental (nicotine/smoking) vs ctrls
DEA_expt_vs_ctl<- function(RSE, formula, name, coef){
  
  ## Previous lib sizes of each sample
  match_samples <- match(RSE$SAMPLE_ID, samples_factors$SAMPLE_ID)
  stopifnot(all(!is.na(match_samples)))
  factors<-samples_factors[match_samples, ]
  
  pdf(file = paste("plots/04_DEA/01_Modeling/Gene_analysis/DEA_plots_", name, ".pdf", sep="" ))
  par(mfrow=c(2,2))
  
  ## Model matrix
  model=model.matrix(formula, data=colData(RSE))
  ## Use previous norm factors to scale the raw library sizes
  RSE_scaled = calcNormFactors(RSE)
  RSE_scaled$samples$lib.size<-factors$lib.size
  RSE_scaled$samples$norm.factors<-factors$norm.factors
  
  ## Transform counts to log2(CPM) 
  ## Estimate mean-variance relationship for each gene
  vGene = voom(RSE_scaled, design=model, plot=TRUE)
  
  ## Fit linear model for each gene
  fitGene = lmFit(vGene)
  
  ## Empirical Bayesian calculation to compute moderated t-statistics
  eBGene = eBayes(fitGene)

  ## Plot average log expression vs logFC
  limma::plotMA(eBGene, coef = coef, xlab = "Mean of normalized counts", 
         ylab="logFC")
  ## Plot -log(p-value) vs logFC
  volcanoplot(eBGene, coef = coef)
  
  ## Select top-ranked genes for Group (expt vs ctrl)
  top_genes = topTable(eBGene, coef=coef, p.value = 1, number=nrow(RSE), sort.by="none")
  
  ## Histogram of adjusted p values
  hist(top_genes$adj.P.Val, xlab="FDR", main="")
  
  dev.off()
  
  return(list(top_genes, vGene, eBGene))
  
}  



## Plots for DE genes
plots_DE<-function(top_genes, vGene, FDR=0.05, name) {
  
  ## NS/Down/Upregulated genes
  DE<-vector()
  for (i in 1:dim(top_genes)[1]) {
    if (top_genes$adj.P.Val[i]>FDR) {
      DE<-append(DE, "ns")
    }
    else {
      if (top_genes$logFC[i]>0) {
        DE<-append(DE, "Up")
      }
      else {
        DE<-append(DE, "Down")
      }
    }
  }
  top_genes$DE<- DE
  
  ## Gene symbols for DEG with |logFC|>1 and top 2 most significant up and down DEGs
  DEG_symbol<-vector()
  
  ## Top most significant up and down genes
  top2up <- subset(top_genes[order(top_genes$adj.P.Val, decreasing = FALSE), ], DE=='Up')[1:2, 'Symbol']
  top2down <- subset(top_genes[order(top_genes$adj.P.Val, decreasing = FALSE), ], DE=='Down')[1:2, 'Symbol']

  for (i in 1:dim(top_genes)[1]) {
    if ((top_genes$DE[i]!="ns" & abs(top_genes$logFC[i])>1) | top_genes$Symbol[i] %in% c(top2up, top2down)) {
      DEG_symbol<-append(DEG_symbol, top_genes$Symbol[i])
    }
    else {
      DEG_symbol<-append(DEG_symbol, NA)
    }
  }
  top_genes$DEG_symbol<- DEG_symbol
  
  
  ## MA plot for DE genes
  cols <- c("Up" = "red", "Down" = "#26b3ff", "ns" = "grey") 
  sizes <- c("Up" = 2, "Down" = 2, "ns" = 1) 
  alphas <- c("Up" = 1, "Down" = 1, "ns" = 0.3)
  top_genes$mean_log_expr<-apply(vGene$E, 1, mean)
  
  p1<-ggplot(data = top_genes,
             aes(x = mean_log_expr,y = logFC,
                 fill = DE,
                 size = DE,
                 alpha = DE)) +
    sm_hgrid(legends = TRUE) +
    geom_point(shape = 21) +
    theme_bw() +
    scale_fill_manual(values = cols, name=NULL) +
    scale_size_manual(values = sizes, name=NULL) +
    scale_alpha_manual(values = alphas, name=NULL) +
    labs(x="Mean of log-normalized counts", y="Log2FC (Exposed vs Ctrl)") +
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          legend.position = c(0.82, 0.15),
          legend.background = element_rect(fill=NA, color='black'),
          legend.key.height = unit(0.15,"cm"),
          axis.title = element_text(size = (10)),
          legend.text = element_text(size=11))
  
  
  ## Volcano plot
  p2<-ggplot(data = top_genes,
             aes(x = logFC,y = -log10(adj.P.Val),
                 fill = DE,
                 size = DE,
                 alpha = DE,
                 label= DEG_symbol)) +
    theme_bw() +
    geom_point(shape =21) +
    geom_hline(yintercept = -log10(FDR),
               linetype = "dashed", color = 'gray65', linewidth=0.5) +
    geom_vline(xintercept = c(-1,1),
               linetype = "dashed", color = 'gray65', linewidth=0.5) +
    geom_label_repel(aes(fontface = 'bold', fill=NULL, alpha=NULL),
                    size=3.3,
                    max.overlaps = 15,
                    min.segment.length = unit(0.2, "cm"),
                    point.padding = unit(0.1, "cm"),
                    box.padding = 0.2,
                    label.padding = 0.2,
                    label.size = 0.2,
                    show.legend=FALSE) +
    labs(y="-log10(FDR)", x="Log2FC (Exposed vs Ctrl)")+
    scale_fill_manual(values = cols, name=NULL) +
    scale_size_manual(values = sizes, name=NULL) +
    scale_alpha_manual(values = alphas, name=NULL) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = c(0.13, 0.85),
          legend.background = element_rect(fill='white', color='black'),
          legend.key.height = unit(0.15,"cm"),
          legend.text = element_text(size=11),
          legend.title = element_text(size = 12),
          axis.title = element_text(size = (12)),
          axis.text = element_text(size = 10)) 
    
  
  plot_grid(p1, p2, ncol=2)
  ggsave(paste("plots/04_DEA/01_Modeling/Gene_analysis/DEG_plots_", name, ".pdf", sep=""), 
         width = 27, height = 14, units = "cm")
}



## Boxplot of a single gene
DE_one_boxplot <- function (de_genes, lognorm_DE, DEgene){
  
  ## q-value for the gene
  q_value<-signif(de_genes[which(de_genes$Symbol==DEgene), "adj.P.Val"], digits = 2)
  
  ## FC
  FC<-signif(2**(de_genes[which(de_genes$Symbol==DEgene), "logFC"]), digits=2)
  
  ## Boxplot for each DE gene
  DEgene<-chartr(":", ".",DEgene)  
   
  ggplot(data=as.data.frame(lognorm_DE), 
    aes(x=Group,y=eval(parse_expr(DEgene)))) + 
    ## Hide outliers
    geom_boxplot(outlier.color = "#FFFFFFFF", width=0.35) +
    ## Samples colored by Group + noise
    geom_jitter(aes(colour=Group),shape=16, 
               position=position_jitter(0.2), size=2.1) +
    theme_bw() +
    scale_color_manual(values=c("Control" = "seashell3", "Experimental" = "orange3")) +
    scale_x_discrete(labels=c("Control"="Ctrl","Experimental"="Expt")) +
    labs(x = "Group", y = "lognorm counts",
        title = DEgene, 
        subtitle = paste("FDR:", q_value, '    ', 'FC:', FC)) +
    theme(plot.margin=unit (c(0.4,0.4,0.4,0.4), 'cm'), 
          legend.position = "none",
          plot.title = element_text(hjust=0.5, size=12, face="bold"), 
          plot.subtitle = element_text(size = 10),
          axis.title = element_text(size = (12)),
          axis.text = element_text(size = 10.5)) 
  
}


 
## Boxplot of lognorm counts for top 2 up and down DE genes
## Obtain lognorm counts of DE genes
DE_boxplots <- function(RSE, vGene, model, de_genes){
  
  ## Order genes by q-value
  de_genes<-de_genes[order(de_genes$adj.P.Val),]
  lognorm_DE<-vGene$E[rownames(de_genes),]
  ## Samples as rows and genes as columns
  lognorm_DE<-t(lognorm_DE)
  colnames(lognorm_DE)<-rowData(RSE)[rownames(de_genes),"Symbol"]
  ## Add samples' Group information
  lognorm_DE<-data.frame(lognorm_DE, "Group"=colData(RSE)$Group)
  
  ## Top most significant up and down genes
  top2up <- subset(de_genes[order(de_genes$adj.P.Val, decreasing = FALSE), ], logFC>0)[1:2, 'Symbol']
  top2down <- subset(de_genes[order(de_genes$adj.P.Val, decreasing = FALSE), ], logFC<0)[1:2, 'Symbol']
  
  plots<-list()
  for (i in 1:4){
    DEgene<-c(top2up, top2down)[i]
    p<-DE_one_boxplot(de_genes, lognorm_DE, DEgene)
    plots[[i]]<-p
  }
  plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol = 2, align = 'hv')
  ggsave(here(paste("plots/04_DEA/01_Modeling/Gene_analysis/DE_boxplots_",name, ".pdf", sep="")), 
         width = 13, height = 14, units = "cm")
}



## Perform DEA for each group of samples
apply_DEA<-function(RSE, formula, name, coef){
  ## DEA
  results<-DEA_expt_vs_ctl(RSE, formula, name, coef)
  top_genes<-results[[1]]
  
  ## If there are DEG
  if (length(which(top_genes$adj.P.Val<0.05))>0){
    model=model.matrix(formula, data=colData(RSE))
    vGene<-results[[2]]
    de_genes<-top_genes[top_genes$adj.P.Val < 0.05,]
    ## Plots for DE genes
    plots_DE(top_genes, vGene, 0.05, name)
    DE_boxplots(RSE, vGene, model, de_genes)
    print(paste(length(which(top_genes$adj.P.Val<0.05)), "differentially expressed genes", sep=" "))
    return(list(results, de_genes))
  }
  else {
    print("No differentially expressed genes")
    return(results)
  }
}



###############################
#         Blood DEA 
#   Smoking vs ctrls in blood 
###############################
 
RSE<-rse_gene_blood_qc

## Naive model
formula<- ~ Group + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate 
name<-"blood_smoking_naive"
coef<-"GroupExperimental"
results_blood_smoking_naive<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes"
top_genes_blood_smoking_naive<-results_blood_smoking_naive[[1]]
save(results_blood_smoking_naive, file="processed-data/04_DEA/Gene_analysis/results_blood_smoking_naive.Rdata")
save(top_genes_blood_smoking_naive, file="processed-data/04_DEA/Gene_analysis/top_genes_blood_smoking_naive.Rdata")



## Model fitted by Pregnancy
formula<- ~ Group + Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate
name<-"blood_smoking_fitted"
coef<-"GroupExperimental"
results_blood_smoking_fitted<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes"
top_genes_blood_smoking_fitted<-results_blood_smoking_fitted[[1]]
save(results_blood_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/results_blood_smoking_fitted.Rdata")
save(top_genes_blood_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_blood_smoking_fitted.Rdata")



## Model with interaction between Group and Pregnancy
formula<- ~ Group*Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate 
name<-"blood_smoking_interaction"
coef<-"GroupExperimental:PregnancyYes"
results_blood_smoking_interaction<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes" 
top_genes_blood_smoking_interaction<-results_blood_smoking_interaction[[1]]
save(results_blood_smoking_interaction, 
     file="processed-data/04_DEA/Gene_analysis/results_blood_smoking_interaction.Rdata")
save(top_genes_blood_smoking_interaction,
     file="processed-data/04_DEA/Gene_analysis/top_genes_blood_smoking_interaction.Rdata")






##################################
#     Nicotine adults DEA
#    Nicotine vs ctrls in adults
##################################

RSE<-rse_gene_brain_adults_nicotine

## Naive model
formula<- ~ Group + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr 
name<-"adults_nicotine_naive"
coef<-"GroupExperimental"
results_adults_nicotine_naive<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes" 
top_genes_adults_nicotine_naive<-results_adults_nicotine_naive[[1]]
save(results_adults_nicotine_naive, file="processed-data/04_DEA/Gene_analysis/results_adults_nicotine_naive.Rdata")
save(top_genes_adults_nicotine_naive, file="processed-data/04_DEA/Gene_analysis/top_genes_adults_nicotine_naive.Rdata")
 


## Model fitted by Pregnancy
formula<- ~ Group + Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate 
name<-"adults_nicotine_fitted"
coef<-"GroupExperimental"
results_adults_nicotine_fitted<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes" 
top_genes_adults_nicotine_fitted<-results_adults_nicotine_fitted[[1]]
save(results_adults_nicotine_fitted, file="processed-data/04_DEA/Gene_analysis/results_adults_nicotine_fitted.Rdata")
save(top_genes_adults_nicotine_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_adults_nicotine_fitted.Rdata")



## Model with interaction between Group and Pregnancy
formula<- ~ Group*Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate
name<-"adults_nicotine_interaction"
coef<-"GroupExperimental:PregnancyYes"
results_adults_nicotine_interaction<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes" 
top_genes_adults_nicotine_interaction<-results_adults_nicotine_interaction[[1]]
save(results_adults_nicotine_interaction, 
     file="processed-data/04_DEA/Gene_analysis/results_adults_nicotine_interaction.Rdata")
save(top_genes_adults_nicotine_interaction,
     file="processed-data/04_DEA/Gene_analysis/top_genes_adults_nicotine_interaction.Rdata")





####################################
#       Smoking adults DEA
#    Smoking vs ctrls in adults
####################################

RSE<-rse_gene_brain_adults_smoking

## Naive model
formula<- ~ Group + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr 
name<-"adults_smoking_naive"
coef<-"GroupExperimental"
results_adults_smoking_naive<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes"
top_genes_adults_smoking_naive<-results_adults_smoking_naive[[1]]
save(results_adults_smoking_naive, file="processed-data/04_DEA/Gene_analysis/results_adults_smoking_naive.Rdata")
save(top_genes_adults_smoking_naive, file="processed-data/04_DEA/Gene_analysis/top_genes_adults_smoking_naive.Rdata")



## Model fitted by Pregnancy
formula<- ~ Group + Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate 
name<-"adults_smoking_fitted"
coef<-"GroupExperimental"
results_adults_smoking_fitted<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes"
top_genes_adults_smoking_fitted<-results_adults_smoking_fitted[[1]]
save(results_adults_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/results_adults_smoking_fitted.Rdata")
save(top_genes_adults_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_adults_smoking_fitted.Rdata")



## Model with interaction between Group and Pregnancy
formula<- ~ Group*Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate 
name<-"adults_smoking_interaction"
coef<-"GroupExperimental:PregnancyYes"
results_adults_smoking_interaction<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes"
top_genes_adults_smoking_interaction<-results_adults_smoking_interaction[[1]]
save(results_adults_smoking_interaction,
     file="processed-data/04_DEA/Gene_analysis/results_adults_smoking_interaction.Rdata")
save(top_genes_adults_smoking_interaction,
     file="processed-data/04_DEA/Gene_analysis/top_genes_adults_smoking_interaction.Rdata")





##################################
#     Nicotine pups DEA
#    Nicotine vs ctrls in pups
##################################

RSE<-rse_gene_brain_pups_nicotine

## Naive model
formula<- ~ Group + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr + mitoRate
name<-"pups_nicotine_naive"
coef<-"GroupExperimental"
results_pups_nicotine_naive<-apply_DEA(RSE, formula, name, coef)
# "1038 differentially expressed genes"
top_genes_pups_nicotine_naive<-results_pups_nicotine_naive[[1]][[1]]
## DE genes
de_genes_pups_nicotine_naive<-results_pups_nicotine_naive[[2]]
save(results_pups_nicotine_naive, file="processed-data/04_DEA/Gene_analysis/results_pups_nicotine_naive.Rdata")
save(top_genes_pups_nicotine_naive, file="processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_naive.Rdata")
save(de_genes_pups_nicotine_naive, file="processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_naive.Rdata")



## Model fitted by Sex
formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
name<-"pups_nicotine_fitted"
coef<-"GroupExperimental"
results_pups_nicotine_fitted<-apply_DEA(RSE, formula, name, coef)
# "1010 differentially expressed genes"
top_genes_pups_nicotine_fitted<-results_pups_nicotine_fitted[[1]][[1]]
## DE genes
de_genes_pups_nicotine_fitted<-results_pups_nicotine_fitted[[2]]
save(results_pups_nicotine_fitted, file="processed-data/04_DEA/Gene_analysis/results_pups_nicotine_fitted.Rdata")
save(top_genes_pups_nicotine_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata")
save(de_genes_pups_nicotine_fitted, file="processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata")
de_genes_brain_pup_nicotine <- de_genes_pups_nicotine_fitted[order(de_genes_pups_nicotine_fitted$adj.P.Val),]
write.table(de_genes_pups_nicotine_fitted, file = "processed-data/04_DEA/Gene_analysis/de_genes_brain_pup_nicotine.csv", row.names = FALSE, col.names = TRUE, sep = '\t')




## Model with interaction between Group and Sex
formula<- ~ Group*Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
name<-"pups_nicotine_interaction"
coef<-"GroupExperimental:SexM"
results_pups_nicotine_interaction<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes"
top_genes_pups_nicotine_interaction<-results_pups_nicotine_interaction[[1]]
save(results_pups_nicotine_interaction, 
     file="processed-data/04_DEA/Gene_analysis/results_pups_nicotine_interaction.Rdata")
save(top_genes_pups_nicotine_interaction,
     file="processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_interaction.Rdata")





##################################
#     Smoking pups DEA
#    Smoking vs ctrls in pups
##################################

RSE<-rse_gene_brain_pups_smoking

## Naive model
formula<- ~ Group + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr + mitoRate
name<-"pups_smoking_naive"
coef<-"GroupExperimental"
results_pups_smoking_naive<-apply_DEA(RSE, formula, name, coef)
# "4108 differentially expressed genes"
top_genes_pups_smoking_naive<-results_pups_smoking_naive[[1]][[1]]
## DE genes
de_genes_pups_smoking_naive<-results_pups_smoking_naive[[2]]
save(results_pups_smoking_naive, file="processed-data/04_DEA/Gene_analysis/results_pups_smoking_naive.Rdata")
save(top_genes_pups_smoking_naive, file="processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_naive.Rdata")
save(de_genes_pups_smoking_naive, file="processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_naive.Rdata")



## Model fitted by Sex
formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
name<-"pups_smoking_fitted"
coef<-"GroupExperimental"
results_pups_smoking_fitted<-apply_DEA(RSE, formula, name, coef)
# "4165 differentially expressed genes"
top_genes_pups_smoking_fitted<-results_pups_smoking_fitted[[1]][[1]]
## DE genes
de_genes_pups_smoking_fitted<-results_pups_smoking_fitted[[2]]
save(results_pups_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/results_pups_smoking_fitted.Rdata")
save(top_genes_pups_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata")
save(de_genes_pups_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata")
de_genes_brain_pup_smoking <- de_genes_pups_smoking_fitted[order(de_genes_pups_smoking_fitted$adj.P.Val),]
write.table(de_genes_brain_pup_smoking, file = "processed-data/04_DEA/Gene_analysis/de_genes_brain_pup_smoking.csv", row.names = FALSE, col.names = TRUE, sep = '\t')



## Model with interaction between Group and Sex
formula<- ~ Group*Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
name<-"pups_smoking_interaction"
coef<-"GroupExperimental:SexM"
results_pups_smoking_interaction<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes"
top_genes_pups_smoking_interaction<-results_pups_smoking_interaction[[1]]
save(results_pups_smoking_interaction, 
     file="processed-data/04_DEA/Gene_analysis/results_pups_smoking_interaction.Rdata")
save(top_genes_pups_smoking_interaction, 
     file="processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_interaction.Rdata")






## Condense results for all genes (don't take replication in human brain)

top_genes_results <- cbind(top_genes_blood_smoking_fitted[, 1:13], 
                           top_genes_blood_smoking_fitted[,c("logFC", "t", "P.Value", "adj.P.Val")],
                           top_genes_adults_nicotine_fitted[,c("logFC", "t", "P.Value", "adj.P.Val", "replication_in_blood")],
                           top_genes_adults_smoking_fitted[,c("logFC", "t", "P.Value", "adj.P.Val", "replication_in_blood")],
                           top_genes_pups_nicotine_fitted[,c("logFC", "t", "P.Value", "adj.P.Val", "replication_in_blood")],
                           top_genes_pups_smoking_fitted[,c("logFC", "t", "P.Value", "adj.P.Val", "replication_in_blood")])

colnames(top_genes_results)[14:37] <- c(paste0(c("logFC", "t", "P.Value", "adj.P.Val"), '_blood_adult_smoking'),
                                        paste0(c("logFC", "t", "P.Value", "adj.P.Val", "replication_in_blood"), '_brain_adult_nicotine'),
                                        paste0(c("logFC", "t", "P.Value", "adj.P.Val", "replication_in_blood"),'_brain_adult_smoking'),
                                        paste0(c("logFC", "t", "P.Value", "adj.P.Val", "replication_in_blood"), '_brain_pup_nicotine'),
                                        paste0(c("logFC", "t", "P.Value", "adj.P.Val", "replication_in_blood"),'_brain_pup_smoking'))

save(top_genes_results, file="processed-data/04_DEA/Gene_analysis/top_genes_results.Rdata")
write.table(top_genes_results, file = "processed-data/04_DEA/Gene_analysis/top_genes_results.csv", row.names = FALSE, col.names = TRUE, sep = '\t')







## Reproducibility information
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
# date     2023-12-18
# rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
# pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version    date (UTC) lib source
# abind                  1.4-5      2016-07-21 [1] CRAN (R 4.3.0)
# aod                    1.3.2      2022-04-02 [1] CRAN (R 4.3.0)
# backports              1.4.1      2021-12-13 [1] CRAN (R 4.3.0)
# base64enc              0.1-3      2015-07-28 [1] CRAN (R 4.3.0)
# Biobase              * 2.61.0     2023-06-02 [1] Bioconductor
# BiocGenerics         * 0.47.0     2023-06-02 [1] Bioconductor
# BiocParallel         * 1.35.3     2023-07-07 [1] Bioconductor
# bitops                 1.0-7      2021-04-24 [1] CRAN (R 4.3.0)
# boot                   1.3-28.1   2022-11-22 [1] CRAN (R 4.3.0)
# broom                  1.0.5      2023-06-09 [1] CRAN (R 4.3.0)
# car                    3.1-2      2023-03-30 [1] CRAN (R 4.3.0)
# carData                3.0-5      2022-01-06 [1] CRAN (R 4.3.0)
# caTools                1.18.2     2021-03-28 [1] CRAN (R 4.3.0)
# cellranger             1.1.0      2016-07-27 [1] CRAN (R 4.3.0)
# checkmate              2.2.0      2023-04-27 [1] CRAN (R 4.3.0)
# cli                    3.6.1      2023-03-23 [1] CRAN (R 4.3.0)
# cluster                2.1.4      2022-08-22 [1] CRAN (R 4.3.0)
# codetools              0.2-19     2023-02-01 [1] CRAN (R 4.3.0)
# colorspace             2.1-0      2023-01-23 [1] CRAN (R 4.3.0)
# corpcor                1.6.10     2021-09-16 [1] CRAN (R 4.3.0)
# cowplot              * 1.1.1      2020-12-30 [1] CRAN (R 4.3.0)
# crayon                 1.5.2      2022-09-29 [1] CRAN (R 4.3.0)
# data.table             1.14.8     2023-02-17 [1] CRAN (R 4.3.0)
# DelayedArray           0.26.6     2023-07-02 [1] Bioconductor
# digest                 0.6.33     2023-07-07 [1] CRAN (R 4.3.0)
# dplyr                  1.1.2      2023-04-20 [1] CRAN (R 4.3.0)
# edgeR                * 3.43.7     2023-06-21 [1] Bioconductor
# EnvStats               2.8.0      2023-07-08 [1] CRAN (R 4.3.0)
# evaluate               0.21       2023-05-05 [1] CRAN (R 4.3.0)
# fANCOVA                0.6-1      2020-11-13 [1] CRAN (R 4.3.0)
# fansi                  1.0.5      2023-10-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1      2022-07-06 [1] CRAN (R 4.3.0)
# fastmap                1.1.1      2023-02-24 [1] CRAN (R 4.3.0)
# foreign                0.8-84     2022-12-06 [1] CRAN (R 4.3.0)
# Formula                1.2-5      2023-02-24 [1] CRAN (R 4.3.0)
# fs                     1.6.3      2023-07-20 [1] CRAN (R 4.3.0)
# gargle                 1.5.2      2023-07-20 [1] CRAN (R 4.3.0)
# generics               0.1.3      2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.37.2     2023-06-21 [1] Bioconductor
# GenomeInfoDbData       1.2.10     2023-05-28 [1] Bioconductor
# GenomicRanges        * 1.53.1     2023-06-02 [1] Bioconductor
# gghalves               0.1.4      2022-11-20 [1] CRAN (R 4.3.0)
# ggplot2              * 3.4.4      2023-10-12 [1] CRAN (R 4.3.1)
# ggpubr                 0.6.0      2023-02-10 [1] CRAN (R 4.3.0)
# ggrepel              * 0.9.3      2023-02-03 [1] CRAN (R 4.3.0)
# ggsignif               0.6.4      2022-10-13 [1] CRAN (R 4.3.0)
# glue                   1.6.2      2022-02-24 [1] CRAN (R 4.3.0)
# googledrive            2.1.1      2023-06-11 [1] CRAN (R 4.3.0)
# gplots                 3.1.3      2022-04-25 [1] CRAN (R 4.3.0)
# gridExtra            * 2.3        2017-09-09 [1] CRAN (R 4.3.0)
# gtable                 0.3.4      2023-08-21 [1] CRAN (R 4.3.0)
# gtools                 3.9.4      2022-11-27 [1] CRAN (R 4.3.0)
# here                 * 1.0.1      2020-12-13 [1] CRAN (R 4.3.0)
# Hmisc                * 5.1-0      2023-05-08 [1] CRAN (R 4.3.0)
# htmlTable              2.4.1      2022-07-07 [1] CRAN (R 4.3.0)
# htmltools              0.5.5      2023-03-23 [1] CRAN (R 4.3.0)
# htmlwidgets            1.6.2      2023-03-17 [1] CRAN (R 4.3.0)
# IRanges              * 2.35.2     2023-06-23 [1] Bioconductor
# iterators              1.0.14     2022-02-05 [1] CRAN (R 4.3.0)
# jaffelab             * 0.99.32    2023-05-28 [1] Github (LieberInstitute/jaffelab@21e6574)
# KernSmooth             2.23-22    2023-07-10 [1] CRAN (R 4.3.0)
# knitr                  1.43       2023-05-25 [1] CRAN (R 4.3.0)
# labeling               0.4.3      2023-08-29 [1] CRAN (R 4.3.0)
# lattice                0.21-8     2023-04-05 [1] CRAN (R 4.3.0)
# lifecycle              1.0.3      2022-10-07 [1] CRAN (R 4.3.0)
# limma                * 3.57.6     2023-06-21 [1] Bioconductor
# lme4                   1.1-34     2023-07-04 [1] CRAN (R 4.3.0)
# lmerTest               3.1-3      2020-10-23 [1] CRAN (R 4.3.0)
# locfit                 1.5-9.8    2023-06-11 [1] CRAN (R 4.3.0)
# magrittr               2.0.3      2022-03-30 [1] CRAN (R 4.3.0)
# MASS                   7.3-60     2023-05-04 [1] CRAN (R 4.3.0)
# Matrix                 1.6-0      2023-07-08 [1] CRAN (R 4.3.0)
# MatrixGenerics       * 1.13.0     2023-05-20 [1] Bioconductor
# matrixStats          * 1.0.0      2023-06-02 [1] CRAN (R 4.3.0)
# minqa                  1.2.5      2022-10-19 [1] CRAN (R 4.3.0)
# munsell                0.5.0      2018-06-12 [1] CRAN (R 4.3.0)
# mvtnorm                1.2-2      2023-06-08 [1] CRAN (R 4.3.0)
# nlme                   3.1-162    2023-01-31 [1] CRAN (R 4.3.0)
# nloptr                 2.0.3      2022-05-26 [1] CRAN (R 4.3.0)
# nnet                   7.3-19     2023-05-03 [1] CRAN (R 4.3.0)
# numDeriv               2016.8-1.1 2019-06-06 [1] CRAN (R 4.3.0)
# pbkrtest               0.5.2      2023-01-19 [1] CRAN (R 4.3.0)
# pillar                 1.9.0      2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.3.0)
# plyr                   1.8.8      2022-11-11 [1] CRAN (R 4.3.0)
# purrr                  1.0.1      2023-01-10 [1] CRAN (R 4.3.0)
# pwr                    1.3-0      2020-03-17 [1] CRAN (R 4.3.0)
# R6                     2.5.1      2021-08-19 [1] CRAN (R 4.3.0)
# rafalib              * 1.0.0      2015-08-09 [1] CRAN (R 4.3.0)
# ragg                   1.2.5      2023-01-12 [1] CRAN (R 4.3.0)
# rbibutils              2.2.13     2023-01-13 [1] CRAN (R 4.3.0)
# RColorBrewer           1.1-3      2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.11     2023-07-06 [1] CRAN (R 4.3.0)
# RCurl                  1.98-1.12  2023-03-27 [1] CRAN (R 4.3.0)
# Rdpack                 2.4        2022-07-20 [1] CRAN (R 4.3.0)
# readxl               * 1.4.3      2023-07-06 [1] CRAN (R 4.3.0)
# remaCor                0.0.16     2023-06-21 [1] CRAN (R 4.3.0)
# reshape2               1.4.4      2020-04-09 [1] CRAN (R 4.3.0)
# RhpcBLASctl            0.23-42    2023-02-11 [1] CRAN (R 4.3.0)
# rlang                * 1.1.1      2023-04-28 [1] CRAN (R 4.3.0)
# rmarkdown              2.23       2023-07-01 [1] CRAN (R 4.3.0)
# rpart                  4.1.19     2022-10-21 [1] CRAN (R 4.3.0)
# rprojroot              2.0.3      2022-04-02 [1] CRAN (R 4.3.0)
# rstatix                0.7.2      2023-02-01 [1] CRAN (R 4.3.0)
# rstudioapi             0.15.0     2023-07-07 [1] CRAN (R 4.3.0)
# S4Arrays               1.1.4      2023-06-02 [1] Bioconductor
# S4Vectors            * 0.39.1     2023-06-02 [1] Bioconductor
# scales                 1.2.1      2022-08-20 [1] CRAN (R 4.3.0)
# sdamr                  0.2.0      2022-11-16 [1] CRAN (R 4.3.0)
# segmented              1.6-4      2023-04-13 [1] CRAN (R 4.3.0)
# sessioninfo          * 1.2.2      2021-12-06 [1] CRAN (R 4.3.0)
# smplot2              * 0.1.0      2023-06-07 [1] Github (smin95/smplot2@836f909)
# stringi                1.7.12     2023-01-11 [1] CRAN (R 4.3.0)
# stringr                1.5.0      2022-12-02 [1] CRAN (R 4.3.0)
# SummarizedExperiment * 1.30.2     2023-06-06 [1] Bioconductor
# systemfonts            1.0.4      2022-02-11 [1] CRAN (R 4.3.0)
# textshaping            0.3.6      2021-10-13 [1] CRAN (R 4.3.0)
# tibble                 3.2.1      2023-03-20 [1] CRAN (R 4.3.0)
# tidyr                  1.3.0      2023-01-24 [1] CRAN (R 4.3.0)
# tidyselect             1.2.0      2022-10-10 [1] CRAN (R 4.3.0)
# utf8                   1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
# variancePartition    * 1.32.2     2023-11-14 [1] Bioconductor
# vctrs                  0.6.4      2023-10-12 [1] CRAN (R 4.3.1)
# withr                  2.5.1      2023-09-26 [1] CRAN (R 4.3.1)
# xfun                   0.39       2023-04-20 [1] CRAN (R 4.3.0)
# XVector                0.41.1     2023-06-02 [1] Bioconductor
# zlibbioc               1.47.0     2023-05-20 [1] Bioconductor
# zoo                    1.8-12     2023-04-13 [1] CRAN (R 4.3.0)
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

