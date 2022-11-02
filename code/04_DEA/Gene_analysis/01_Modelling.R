
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
library(sessioninfo)


# 1. Differential Expression Analysis at the gene level

load(here("raw-data/rse_gene_smoking_mouse_n208.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_gene_blood_qc.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_adults_nicotine.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_adults_smoking.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_nicotine.Rdata"))
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_smoking.Rdata"))


## 1.1 Modelling

## Extract previous output from calcNormFactors for all samples
norm_factors<-calcNormFactors(rse_gene, method = "TMM")
samples_factors<-data.frame(SAMPLE_ID=norm_factors$samples$SAMPLE_ID,
                            norm.factors=norm_factors$samples$norm.factors,
                            lib.size=norm_factors$samples$lib.size)


## DEA for experiment mice (nicotine/smoking) vs ctrls
DEA_expt_vs_ctl<- function(RSE, formula, name, coef){
  
  ## Previous lib sizes of each sample
  match_samples <- match(RSE$SAMPLE_ID, samples_factors$SAMPLE_ID)
  stopifnot(all(!is.na(match_samples)))
  factors<-samples_factors[match_samples, ]
  
  pdf(file = paste("plots/04_DEA/01_Modelling/DEA_plots_", name, ".pdf", sep="" ))
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
  
  ## Empirical Bayesian calculation to obtain our significant genes:
  ## compute F and t-statistics, p-value and logFC for DE 
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
  
  ## Gene symbols for DEG with |logFC|>1
  DEG_symbol<-vector()
  for (i in 1:dim(top_genes)[1]) {
    if (top_genes$DE[i]!="ns" & abs(top_genes$logFC[i])>1) {
      DEG_symbol<-append(DEG_symbol, top_genes$Symbol[i])
    }
    else {
      DEG_symbol<-append(DEG_symbol, NA)
    }
  }
  top_genes$DEG_symbol<- DEG_symbol
  
  
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
    scale_fill_manual(values = cols) + 
    scale_size_manual(values = sizes) + 
    scale_alpha_manual(values = alphas) +
    labs(x="Mean of normalized counts")
  
  
  ## Volcano plot for DE genes
  p2<-ggplot(data = top_genes, 
        aes(x = logFC,y = -log10(adj.P.Val),
               fill = DE,    
               size = DE,
               alpha = DE,
               label= DEG_symbol)) +
    geom_point(shape = 21) + 
    geom_hline(yintercept = -log10(FDR),
               linetype = "dashed") + 
    geom_vline(xintercept = c(-1,1),
               linetype = "dashed") +
    geom_label_repel(fill="white", size=2, max.overlaps = Inf,  
                     box.padding = 0.2, 
                     show.legend=FALSE) +
    labs(y="-log10(FDR)")+
    scale_fill_manual(values = cols) + 
    scale_size_manual(values = sizes) + 
    scale_alpha_manual(values = alphas) 
  
  plot_grid(p1, p2, ncol=2)
  ggsave(paste("plots/04_DEA/01_Modelling/DEG_plots_", name, ".pdf", sep=""), 
         width = 35, height = 15, units = "cm")
}



## Boxplot of a single gene
DE_one_boxplot <- function (de_genes, lognorm_DE, DEgene){
  
  ## q-value for the gene
  q_value<-signif(de_genes[which(de_genes$Symbol==DEgene), "adj.P.Val"], digits = 3)
  
  ## Boxplot for each DE gene
  DEgene<-chartr(":", ".",DEgene)  
  ggplot(data=as.data.frame(lognorm_DE), 
  aes(x=Group,y=eval(parse_expr(DEgene)))) + 
  ## Hide outliers
  geom_boxplot(outlier.color = "#FFFFFFFF") +
  ## Samples colored by Group + noise
  geom_jitter(aes(colour=Group),shape=16, 
             position=position_jitter(0.2)) +
  theme_classic() +
  labs(x = "Group", y = "norm counts",
      title = DEgene, 
      subtitle = paste("FDR:", q_value)) +
  theme(plot.margin=unit (c (1,1.5,1,1), 'cm'), legend.position = "none",
       plot.title = element_text(hjust=0.5, size=10, face="bold"), 
       plot.subtitle = element_text(size = 9)) 
  
}


 
## Boxplot of lognorm counts for top 3 DE genes
## Obtain lognorm counts of DE genes
DE_boxplots <- function(RSE, vGene, model, de_genes){
  ## Regress out residuals to remove batch effects
  vGene$E<-cleaningY(vGene$E, model, P=2)
  ## Order genes by q-value
  de_genes<-de_genes[order(de_genes$adj.P.Val),]
  lognorm_DE<-vGene$E[rownames(de_genes),]
  ## Samples as rows and genes as columns
  lognorm_DE<-t(lognorm_DE)
  colnames(lognorm_DE)<-rowData(RSE)[rownames(de_genes),"Symbol"]
  ## Add samples' Group information
  lognorm_DE<-data.frame(lognorm_DE, "Group"=colData(RSE)$Group)
  
 
  plots<-list()
  for (i in 1:3){
    DEgene<-colnames(lognorm_DE)[i]
    p<-DE_one_boxplot(de_genes, lognorm_DE, DEgene)
    plots[[i]]<-p
  }
  plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol = 3)
  ggsave(here(paste("plots/04_DEA/01_Modelling/DE_boxplots_",name, ".pdf", sep="")), 
         width = 25, height = 10, units = "cm")
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
save(results_blood_smoking_naive, file="processed-data/04_DEA/results_blood_smoking_naive.Rdata")
save(top_genes_blood_smoking_naive, file="processed-data/04_DEA/top_genes_blood_smoking_naive.Rdata")



## Model fitted by Pregnancy
formula<- ~ Group + Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate
name<-"blood_smoking_fitted"
coef<-"GroupExperimental"
results_blood_smoking_fitted<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes"
top_genes_blood_smoking_fitted<-results_blood_smoking_fitted[[1]]
save(results_blood_smoking_fitted, file="processed-data/04_DEA/results_blood_smoking_fitted.Rdata")
save(top_genes_blood_smoking_fitted, file="processed-data/04_DEA/top_genes_blood_smoking_fitted.Rdata")



## Model with interaction between Group and Pregnancy
formula<- ~ Group*Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate 
name<-"blood_smoking_interaction"
coef<-"GroupExperimental:PregnancyYes"
results_blood_smoking_interaction<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes" 
top_genes_blood_smoking_interaction<-results_blood_smoking_interaction[[1]]
save(results_blood_smoking_interaction, 
     file="processed-data/04_DEA/results_blood_smoking_interaction.Rdata")
save(top_genes_blood_smoking_interaction,
     file="processed-data/04_DEA/top_genes_blood_smoking_interaction.Rdata")






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
save(results_adults_nicotine_naive, file="processed-data/04_DEA/results_adults_nicotine_naive.Rdata")
save(top_genes_adults_nicotine_naive, file="processed-data/04_DEA/top_genes_adults_nicotine_naive.Rdata")
 


## Model fitted by Pregnancy
formula<- ~ Group + Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate 
name<-"adults_nicotine_fitted"
coef<-"GroupExperimental"
results_adults_nicotine_fitted<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes" 
top_genes_adults_nicotine_fitted<-results_adults_nicotine_fitted[[1]]
save(results_adults_nicotine_fitted, file="processed-data/04_DEA/results_adults_nicotine_fitted.Rdata")
save(top_genes_adults_nicotine_fitted, file="processed-data/04_DEA/top_genes_adults_nicotine_fitted.Rdata")



## Model with interaction between Group and Pregnancy
formula<- ~ Group*Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate
name<-"adults_nicotine_interaction"
coef<-"GroupExperimental:PregnancyYes"
results_adults_nicotine_interaction<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes" 
top_genes_adults_nicotine_interaction<-results_adults_nicotine_interaction[[1]]
save(results_adults_nicotine_interaction, 
     file="processed-data/04_DEA/results_adults_nicotine_interaction.Rdata")
save(top_genes_adults_nicotine_interaction,
     file="processed-data/04_DEA/top_genes_adults_nicotine_interaction.Rdata")





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
save(results_adults_smoking_naive, file="processed-data/04_DEA/results_adults_smoking_naive.Rdata")
save(top_genes_adults_smoking_naive, file="processed-data/04_DEA/top_genes_adults_smoking_naive.Rdata")



## Model fitted by Pregnancy
formula<- ~ Group + Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate 
name<-"adults_smoking_fitted"
coef<-"GroupExperimental"
results_adults_smoking_fitted<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes"
top_genes_adults_smoking_fitted<-results_adults_smoking_fitted[[1]]
save(results_adults_smoking_fitted, file="processed-data/04_DEA/results_adults_smoking_fitted.Rdata")
save(top_genes_adults_smoking_fitted, file="processed-data/04_DEA/top_genes_adults_smoking_fitted.Rdata")



## Model with interaction between Group and Pregnancy
formula<- ~ Group*Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate 
name<-"adults_smoking_interaction"
coef<-"GroupExperimental:PregnancyYes"
results_adults_smoking_interaction<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes"
top_genes_adults_smoking_interaction<-results_adults_smoking_interaction[[1]]
save(results_adults_smoking_interaction,
     file="processed-data/04_DEA/results_adults_smoking_interaction.Rdata")
save(top_genes_adults_smoking_interaction,
     file="processed-data/04_DEA/top_genes_adults_smoking_interaction.Rdata")





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
save(results_pups_nicotine_naive, file="processed-data/04_DEA/results_pups_nicotine_naive.Rdata")
save(top_genes_pups_nicotine_naive, file="processed-data/04_DEA/top_genes_pups_nicotine_naive.Rdata")
save(de_genes_pups_nicotine_naive, file="processed-data/04_DEA/de_genes_pups_nicotine_naive.Rdata")



## Model fitted by Sex
formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
name<-"pups_nicotine_fitted"
coef<-"GroupExperimental"
results_pups_nicotine_fitted<-apply_DEA(RSE, formula, name, coef)
# "1010 differentially expressed genes"
top_genes_pups_nicotine_fitted<-results_pups_nicotine_fitted[[1]][[1]]
## DE genes
de_genes_pups_nicotine_fitted<-results_pups_nicotine_fitted[[2]]
save(results_pups_nicotine_fitted, file="processed-data/04_DEA/results_pups_nicotine_fitted.Rdata")
save(top_genes_pups_nicotine_fitted, file="processed-data/04_DEA/top_genes_pups_nicotine_fitted.Rdata")
save(de_genes_pups_nicotine_fitted, file="processed-data/04_DEA/de_genes_pups_nicotine_fitted.Rdata")



## Model with interaction between Group and Sex
formula<- ~ Group*Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
name<-"pups_nicotine_interaction"
coef<-"GroupExperimental:SexM"
results_pups_nicotine_interaction<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes"
top_genes_pups_nicotine_interaction<-results_pups_nicotine_interaction[[1]]
save(results_pups_nicotine_interaction, 
     file="processed-data/04_DEA/results_pups_nicotine_interaction.Rdata")
save(top_genes_pups_nicotine_interaction,
     file="processed-data/04_DEA/top_genes_pups_nicotine_interaction.Rdata")





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
save(results_pups_smoking_naive, file="processed-data/04_DEA/results_pups_smoking_naive.Rdata")
save(top_genes_pups_smoking_naive, file="processed-data/04_DEA/top_genes_pups_smoking_naive.Rdata")
save(de_genes_pups_smoking_naive, file="processed-data/04_DEA/de_genes_pups_smoking_naive.Rdata")



## Model fitted by Sex
formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
name<-"pups_smoking_fitted"
coef<-"GroupExperimental"
results_pups_smoking_fitted<-apply_DEA(RSE, formula, name, coef)
# "4165 differentially expressed genes"
top_genes_pups_smoking_fitted<-results_pups_smoking_fitted[[1]][[1]]
## DE genes
de_genes_pups_smoking_fitted<-results_pups_smoking_fitted[[2]]
save(results_pups_smoking_fitted, file="processed-data/04_DEA/results_pups_smoking_fitted.Rdata")
save(top_genes_pups_smoking_fitted, file="processed-data/04_DEA/top_genes_pups_smoking_fitted.Rdata")
save(de_genes_pups_smoking_fitted, file="processed-data/04_DEA/de_genes_pups_smoking_fitted.Rdata")



## Model with interaction between Group and Sex
formula<- ~ Group*Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
name<-"pups_smoking_interaction"
coef<-"GroupExperimental:SexM"
results_pups_smoking_interaction<-apply_DEA(RSE, formula, name, coef)
# "No differentially expressed genes"
top_genes_pups_smoking_interaction<-results_pups_smoking_interaction[[1]]
save(results_pups_smoking_interaction, 
     file="processed-data/04_DEA/results_pups_smoking_interaction.Rdata")
save(top_genes_pups_smoking_interaction, 
     file="processed-data/04_DEA/top_genes_pups_smoking_interaction.Rdata")


