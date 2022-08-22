library(here)
library(SummarizedExperiment)
library(stats)
library(edgeR)
library(limma)
library(ExploreModelMatrix)
library(ggplot2)
library(rlang)
library(cowplot)


# 1. Differential Expression Analysis 

load(here("processed-data/03_EDA/02_QC/rse_gene_blood_qc.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_adults_nicotine.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_adults_smoking.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_nicotine.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_smoking.Rdata"))


## 1.1 Modelling

## DEA for experiment mice (nicotine/smoking) vs ctrls
DEA_expt_vs_ctl<- function(RSE, formula, name, coef){
  
  pdf(file = paste("plots/04_DEA/DEA_plots_", name, ".pdf", sep="" ))
  par(mfrow=c(2,2))
  
  ## Model matrix
  model=model.matrix(formula, data=colData(RSE))
  ## Compute normalization factors to scale the raw library sizes
  RSE_scaled = calcNormFactors(RSE)
  ## Transform counts to log2(CPM) 
  ## Estimate mean-variance relationship for each gene
  vGene = voom(RSE_scaled, model, plot=TRUE)
  
  ## Fit linear model for each gene 
  fitGene = lmFit(vGene)
  
  ## Empirical Bayesian calculation to obtain our significant genes:
  ## compute F and t-statistics, p-value and logFC for DE 
  eBGene = eBayes(fitGene)

  ## Plot average log expression vs logFC
  plotMA(eBGene, coef = coef, xlab = "Mean of normalized counts", 
         ylab="logFC")
  ## Plot -log(p-value) vs logFC
  volcanoplot(eBGene, coef = coef)
  
  ## Select top-ranked genes for Group (expt vs ctrl)
  top_genes = topTable(eBGene, coef=coef, p.value = 1, number=nrow(RSE), sort.by="none")
  ## Histogram of adjusted p values
  hist(top_genes$adj.P.Val, xlab="p-value", main="")
  
  dev.off()
  
  return(list(top_genes, vGene, eBGene))
  
}  


## Plots for DE genes
plots_DE<-function(top_gene, vGene, p_value, name) {
  ## NS/Down/Upregulated genes
  DE<-vector()
  for (i in 1:dim(top_genes)[1]) {
    if (top_genes$adj.P.Val[i]>p_value) {
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
  
  ## MA plot for DE genes
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
    #geom_vline(xintercept = c(log2(0.5), log2(2)),
    #          linetype = "dashed") +
    scale_fill_manual(values = cols) + 
    scale_size_manual(values = sizes) + 
    scale_alpha_manual(values = alphas) +
    labs(x="Mean of normalized counts")
  
  
  ## Volcano plot for DE genes
  p2<-ggplot(data = top_genes, 
        aes(x = logFC,y = -log10(adj.P.Val),
               fill = DE,    
               size = DE,
               alpha = DE)) + 
    geom_point(shape = 21,    
               colour = "black") + 
    geom_hline(yintercept = -log10(p_value),
               linetype = "dashed") + 
    #geom_vline(xintercept = c(log2(0.5), log2(2)),
    #          linetype = "dashed") +
    scale_fill_manual(values = cols) + 
    scale_size_manual(values = sizes) + 
    scale_alpha_manual(values = alphas) 
  
  plot_grid(p1, p2, ncol=2)
  ggsave(paste("plots/04_DEA/DEG_plots_", name, ".pdf", sep=""), 
         width = 35, height = 15, units = "cm")
}


 
## Boxplot of lognorm counts for top 3 DE genes
DE_boxplots <- function(RSE, vGene, de_genes){
  ## Order genes by adjusted p-value
  de_genes<-de_genes[order(de_genes$adj.P.Val),]
  ## Obtain lognorm counts of DE genes
  lognorm_DE<-vGene$E[rownames(de_genes),]
  ## Samples as rows and genes as columns
  if (dim(de_genes)[1]!=1) {
    lognorm_DE<-t(lognorm_DE)
  }
  else {
    lognorm_DE<-matrix(lognorm_DE, ncol = 1)
    colnames(lognorm_DE)<-rownames(de_genes)
  }
  ## Add samples' Group information
  lognorm_DE<-data.frame(lognorm_DE, "Group"=colData(RSE)$Group)
  
  ## Boxplot for each DE gene
  if (dim(de_genes)[1]< 4) {
    n=dim(de_genes)[1]
    }
  else {
    n=3
    }
  plots<-list()
  i=1
  for (DEgene in colnames(lognorm_DE)[1:n]){
     p<-ggplot(data=as.data.frame(lognorm_DE), 
     aes(x=Group,y=eval(parse_expr(DEgene)))) + 
     ## Hide outliers
     geom_boxplot(outlier.color = "#FFFFFFFF") +
     ## Samples colored by Group + noise
     geom_jitter(aes(colour=Group),shape=16, 
                 position=position_jitter(0.2)) +
     theme_classic() +
     ggtitle(DEgene) +
     theme(legend.position="right", plot.margin=unit (c (1,1.5,1,1), 'cm'),
           plot.title = element_text(hjust=0.5, size=10, face="bold")) +
     labs(x= "Group", y = "lognorm counts")
     
     plots[[i]]<-p
     i=i+1
  }
  
  if (n==3){
    plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol=3)
    ## Save plot
    ggsave(here(paste("plots/04_DEA/DE_boxplots_",name, ".pdf", sep="")), 
        width = 35, height = 10, units = "cm")
  }
  else if(n==2) {
    plot_grid(plots[[1]], plots[[2]], ncol=2)
    ggsave(here(paste("plots/04_DEA/DE_boxplots_",name, ".pdf", sep="")), 
        width = 12, height = 10, units = "cm")    
  }
  else {
    plot_grid(plots[[1]])
    ggsave(here(paste("plots/04_DEA/DE_boxplots_",name, ".pdf", sep="")), 
        width = 12, height = 10, units = "cm")
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
top_genes<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[1]]
# vGene<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[2]]

## DE genes
table(top_genes$adj.P.Val < 0.85)
# FALSE   
# 19974   


## Model fitted by Pregnancy
formula<- ~ Group + Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate
name<-"blood_smoking_fitted"
top_genes<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[1]]
vGene<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[2]]

## DE genes
table(top_genes$adj.P.Val < 0.6)
# FALSE  TRUE 
# 19808   166 
de_genes_blood_smoking_fitted<-top_genes[top_genes$adj.P.Val < 0.6,]
## Plots for DE genes
plots_DE(top_gene, vGene, 0.6, name)
## Boxplots of lognorm counts for DEG
DE_boxplots(RSE, vGene, de_genes_blood_smoking_fitted)



## Model with interaction between Group and Pregnancy
formula<- ~ Group*Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate 
name<-"blood_smoking_interaction"
top_genes<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental:PregnancyYes")[[1]]
# vGene<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental:PregnancyYes")[[2]]

## DE genes
table(top_genes$adj.P.Val < 0.9)
# FALSE   
# 19974   







##################################
#     Nicotine adults DEA
#    Nicotine vs ctrls in adults
##################################

RSE<-rse_gene_brain_adults_nicotine

## Naive model
formula<- ~ Group + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr 
name<-"adults_nicotine_naive"
top_genes<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[1]]
# vGene<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[2]]

## DE genes
table(top_genes$adj.P.Val < 0.9)
# FALSE   
# 19974     



## Model fitted by Pregnancy
formula<- ~ Group + Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate 
name<-"adults_nicotine_fitted"
top_genes<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[1]]
# vGene<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[2]]

## DE genes
table(top_genes$adj.P.Val < 0.9)
# FALSE  
# 19974  


## Model with interaction between Group and Pregnancy
formula<- ~ Group*Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate
name<-"adults_nicotine_interaction"
top_genes<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental:PregnancyYes")[[1]]
# vGene<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental:PregnancyYes")[[2]]

## DE genes
table(top_genes$adj.P.Val < 0.9)
# FALSE 
# 19974





####################################
#       Smoking adults DEA
#    Smoking vs ctrls in adults
####################################

RSE<-rse_gene_brain_adults_smoking

## Naive model
formula<- ~ Group + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr 
name<-"adults_smoking_naive"
top_genes<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[1]]
vGene<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[2]]

## DE genes
table(top_genes$adj.P.Val < 0.65)
# FALSE  TRUE 
# 19973     1 
de_genes_adults_smoking_naive<-top_genes[top_genes$adj.P.Val < 0.65,]
## Plots for DE genes
plots_DE(top_gene, vGene, 0.65, name)
DE_boxplots(RSE, vGene, de_genes_adults_smoking_naive)



## Model fitted by Pregnancy
formula<- ~ Group + Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate 
name<-"adults_smoking_fitted"
top_genes<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[1]]
# vGene<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[2]]

## DE genes
table(top_genes$adj.P.Val < 0.9)
# FALSE  
# 19974  



## Model with interaction between Group and Pregnancy
formula<- ~ Group*Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate 
name<-"adults_smoking_interaction"
top_genes<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental:PregnancyYes")[[1]]
vGene<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental:PregnancyYes")[[2]]

## DE genes
table(top_genes$adj.P.Val < 0.30)
# FALSE  TRUE 
# 19973     1  
de_genes_adults_smoking_interaction<-top_genes[top_genes$adj.P.Val < 0.30,]
## Plots for DE genes
plots_DE(top_gene, vGene, 0.30, name)
DE_boxplots(RSE, vGene, de_genes_adults_smoking_interaction)







##################################
#     Nicotine pups DEA
#    Nicotine vs ctrls in pups
##################################

RSE<-rse_gene_brain_pups_nicotine

## Naive model
formula<- ~ Group + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr + mitoRate
name<-"pups_nicotine_naive"
top_genes<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[1]]
vGene<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[2]]

## DE genes
table(top_genes$adj.P.Val < 0.05)
# FALSE  TRUE 
# 19001   973 
de_genes_pups_nicotine_naive<-top_genes[top_genes$adj.P.Val < 0.05,]
## Plots for DE genes
plots_DE(top_gene, vGene, 0.05, name)
DE_boxplots(RSE, vGene, de_genes_pups_nicotine_naive)



## Model fitted by Sex
formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
name<-"pups_nicotine_fitted"
top_genes<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[1]]
vGene<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[2]]

## DE genes
table(top_genes$adj.P.Val < 0.05)
# FALSE  TRUE 
# 18887  1087  
de_genes_pups_nicotine_fitted<-top_genes[top_genes$adj.P.Val < 0.05,]
## Plots for DE genes
plots_DE(top_gene, vGene, 0.05, name)
DE_boxplots(RSE, vGene, de_genes_pups_nicotine_fitted)



## Model with interaction between Group and Sex
formula<- ~ Group*Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
name<-"pups_nicotine_interaction"
top_genes<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental:SexM")[[1]]
vGene<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental:SexM")[[2]]

## DE genes
table(top_genes$adj.P.Val < 0.3)
# FALSE  TRUE 
# 19958    16 
de_genes_pups_nicotine_interaction<-top_genes[top_genes$adj.P.Val < 0.3,]
## Plots for DE genes
plots_DE(top_gene, vGene, 0.3, name)
DE_boxplots(RSE, vGene, de_genes_pups_nicotine_interaction)







##################################
#     Smoking pups DEA
#    Smoking vs ctrls in pups
##################################

RSE<-rse_gene_brain_pups_smoking

## Naive model
formula<- ~ Group + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr + mitoRate
name<-"pups_smoking_naive"
top_genes<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[1]]
vGene<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[2]]

## DE genes
table(top_genes$adj.P.Val < 0.05)
# FALSE  TRUE 
# 15645  4329
de_genes_pups_smoking_naive<-top_genes[top_genes$adj.P.Val < 0.05,]
## Plots for DE genes
plots_DE(top_gene, vGene, 0.05, name)
DE_boxplots(RSE, vGene, de_genes_pups_smoking_naive)



## Model fitted by Sex
formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
name<-"pups_smoking_fitted"
top_genes<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[1]]
vGene<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental")[[2]]

## DE genes
table(top_genes$adj.P.Val < 0.05)
# FALSE  TRUE 
# 15592  4382  
de_genes_pups_smoking_fitted<-top_genes[top_genes$adj.P.Val < 0.05,]
## Plots for DE genes
plots_DE(top_gene, vGene, 0.05, name)
DE_boxplots(RSE, vGene, de_genes_pups_smoking_fitted)



## Model with interaction between Group and Sex
formula<- ~ Group*Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
name<-"pups_smoking_interaction"
top_genes<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental:SexM")[[1]]
vGene<-DEA_expt_vs_ctl(RSE, formula, name, "GroupExperimental:SexM")[[2]]

## DE genes
table(top_genes$adj.P.Val < 0.57)
# FALSE  TRUE 
# 18435  1539 
de_genes_pups_smoking_interaction<-top_genes[top_genes$adj.P.Val < 0.57,]
## Plots for DE genes
plots_DE(top_gene, vGene, 0.57, name)
DE_boxplots(RSE, vGene, de_genes_pups_smoking_interaction)
