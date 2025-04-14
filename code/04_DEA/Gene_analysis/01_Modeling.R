
## All libraries used for the whole DEA

library(here)
library(SummarizedExperiment)
library(stats)
library(edgeR)
library(limma)
library(ggplot2)
library(rlang)
library(smplot2)
library(cowplot)
library(ggrepel)
library(jaffelab)
library(VennDiagram) 
library(gridExtra)
library(biomartr)
library(Hmisc)
library(readxl)
library(qvalue)
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


## DEA for experimental (nicotine/smoking) vs ctrls, or F vs M 
DEA<- function(RSE, formula, name, coef){
  
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
  
  ## Select top-ranked genes for coef
  top_genes = topTable(eBGene, coef=coef, p.value = 1, number=nrow(RSE), sort.by="none")
  
  ## Histogram of adjusted p values
  hist(top_genes$adj.P.Val, xlab="FDR", main="")
  
  dev.off()
  
  return(list(top_genes, vGene, eBGene))
  
}  



## Plots for DE genes
plots_DE<-function(top_genes, vGene, FDR=0.05, name, comparison) {
  
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
  
  ## axis title
  if(comparison == "Sex"){
    lab="Log2FC (Male vs Female)"
  }else{
    lab="Log2FC (Exposed vs Ctrl)"
  }
  
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
    labs(x="Mean of log-normalized counts", y=lab) +
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
    labs(y="-log10(FDR)", x=lab)+
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
DE_one_boxplot <- function (de_genes, lognorm_DE, DEgene, comparison){
  
  ## q-value for the gene
  q_value<-signif(de_genes[which(de_genes$Symbol==DEgene), "adj.P.Val"], digits = 2)
  
  ## FC
  FC<-signif(2**(de_genes[which(de_genes$Symbol==DEgene), "logFC"]), digits=2)
  
  ## Boxplot for each DE gene
  DEgene<-chartr(":", ".",DEgene)  
  
  ## axis title
  if(comparison == "Sex"){
    lab="Log2FC (Male vs Female)"
    colors = c("F" = "hotpink1", "M" = "dodgerblue")
    labs = c("F"="Female","M"="Male")
  }else{
    lab="Log2FC (Exposed vs Ctrl)"
    colors = c("Control" = "seashell3", "Experimental" = "orange3")
    labs = c("Control"="Ctrl","Experimental"="Expt")
  }
  
   
  ggplot(data=as.data.frame(lognorm_DE), 
    aes(x=coef, y=eval(parse_expr(DEgene)))) + 
    ## Hide outliers
    geom_boxplot(outlier.color = "#FFFFFFFF", width=0.35) +
    ## Samples colored by Group + noise
    geom_jitter(aes(colour=coef),shape=16, 
               position=position_jitter(0.2), size=2.1) +
    theme_bw() +
    scale_color_manual(values=colors) +
    scale_x_discrete(labels=labs) +
    labs(x = comparison, y = "lognorm counts",
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
DE_boxplots <- function(RSE, vGene, model, de_genes, comparison){
  
  ## Order genes by q-value
  de_genes<-de_genes[order(de_genes$adj.P.Val),]
  lognorm_DE<-vGene$E[rownames(de_genes),]
  ## Samples as rows and genes as columns
  lognorm_DE<-t(lognorm_DE)
  colnames(lognorm_DE)<-rowData(RSE)[rownames(de_genes),"Symbol"]
  ## Add samples' coef information
  lognorm_DE<-data.frame(lognorm_DE, "coef"=colData(RSE)[, comparison])
  
  ## Top most significant up and down genes
  top2up <- subset(de_genes[order(de_genes$adj.P.Val, decreasing = FALSE), ], logFC>0)[1:2, 'Symbol']
  top2down <- subset(de_genes[order(de_genes$adj.P.Val, decreasing = FALSE), ], logFC<0)[1:2, 'Symbol']
  
  plots<-list()
  for (i in 1:4){
    DEgene<-c(top2up, top2down)[i]
    p<-DE_one_boxplot(de_genes, lognorm_DE, DEgene, comparison)
    plots[[i]]<-p
  }
  plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol = 2, align = 'hv')
  ggsave(here(paste("plots/04_DEA/01_Modeling/Gene_analysis/DE_boxplots_", name, ".pdf", sep="")), 
         width = 13, height = 14, units = "cm")
}



## Perform DEA for each group of samples
apply_DEA<-function(RSE, formula, name, coef, comparison){
  ## DEA
  results<-DEA(RSE, formula, name, coef)
  top_genes<-results[[1]]
  
  ## If there are DEG
  if (length(which(top_genes$adj.P.Val<0.05))>0){
    model=model.matrix(formula, data=colData(RSE))
    vGene<-results[[2]]
    de_genes<-top_genes[top_genes$adj.P.Val < 0.05,]
    ## Plots for DE genes
    plots_DE(top_genes, vGene, 0.05, name, comparison)
    DE_boxplots(RSE, vGene, model, de_genes, comparison)
    print(paste(length(which(top_genes$adj.P.Val<0.05)), "differentially expressed genes", sep=" "))
    return(list(results, de_genes))
  }
  else {
    print("No differentially expressed genes")
    return(results)
  }
}



##########################################
#        Smoking Blood Adults DEA 
#   Smoking vs ctrls in blood in adults
##########################################
 
RSE<-rse_gene_blood_qc

## Naive model
formula<- ~ Group + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate 
name<-"blood_smoking_naive"
coef<-"GroupExperimental"
results_blood_smoking_naive<-apply_DEA(RSE, formula, name, coef, "Group")
# "No differentially expressed genes"
top_genes_blood_smoking_naive<-results_blood_smoking_naive[[1]]
save(results_blood_smoking_naive, file="processed-data/04_DEA/Gene_analysis/results_blood_smoking_naive.Rdata")
save(top_genes_blood_smoking_naive, file="processed-data/04_DEA/Gene_analysis/top_genes_blood_smoking_naive.Rdata")

## Model adjusted for Pregnancy
formula<- ~ Group + Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate
name<-"blood_smoking_fitted"
coef<-"GroupExperimental"
results_blood_smoking_fitted<-apply_DEA(RSE, formula, name, coef, "Group")
# "No differentially expressed genes"
top_genes_blood_smoking_fitted<-results_blood_smoking_fitted[[1]]
save(results_blood_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/results_blood_smoking_fitted.Rdata")
save(top_genes_blood_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_blood_smoking_fitted.Rdata")

## Model with interaction between Group and Pregnancy
formula<- ~ Group*Pregnancy + Group + Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate 
name<-"blood_smoking_interaction"
coef<-"GroupExperimental:PregnancyYes"
results_blood_smoking_interaction<-apply_DEA(RSE, formula, name, coef, "Group")
# "No differentially expressed genes" 
top_genes_blood_smoking_interaction<-results_blood_smoking_interaction[[1]]
save(results_blood_smoking_interaction, 
     file="processed-data/04_DEA/Gene_analysis/results_blood_smoking_interaction.Rdata")
save(top_genes_blood_smoking_interaction,
     file="processed-data/04_DEA/Gene_analysis/top_genes_blood_smoking_interaction.Rdata")



############################################
#         Nicotine Brain Adults DEA
#    Nicotine vs ctrls in brain in adults
############################################

RSE<-rse_gene_brain_adults_nicotine

## Naive model
formula<- ~ Group + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr 
name<-"adults_nicotine_naive"
coef<-"GroupExperimental"
results_adults_nicotine_naive<-apply_DEA(RSE, formula, name, coef, "Group")
# "No differentially expressed genes" 
top_genes_adults_nicotine_naive<-results_adults_nicotine_naive[[1]]
save(results_adults_nicotine_naive, file="processed-data/04_DEA/Gene_analysis/results_adults_nicotine_naive.Rdata")
save(top_genes_adults_nicotine_naive, file="processed-data/04_DEA/Gene_analysis/top_genes_adults_nicotine_naive.Rdata")

## Model adjusted for Pregnancy
formula<- ~ Group + Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate 
name<-"adults_nicotine_fitted"
coef<-"GroupExperimental"
results_adults_nicotine_fitted<-apply_DEA(RSE, formula, name, coef, "Group")
# "No differentially expressed genes" 
top_genes_adults_nicotine_fitted<-results_adults_nicotine_fitted[[1]]
save(results_adults_nicotine_fitted, file="processed-data/04_DEA/Gene_analysis/results_adults_nicotine_fitted.Rdata")
save(top_genes_adults_nicotine_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_adults_nicotine_fitted.Rdata")

## Model with interaction between Group and Pregnancy
formula<- ~ Group*Pregnancy +  Group + Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate
name<-"adults_nicotine_interaction"
coef<-"GroupExperimental:PregnancyYes"
results_adults_nicotine_interaction<-apply_DEA(RSE, formula, name, coef, "Group")
# "No differentially expressed genes" 
top_genes_adults_nicotine_interaction<-results_adults_nicotine_interaction[[1]]
save(results_adults_nicotine_interaction, 
     file="processed-data/04_DEA/Gene_analysis/results_adults_nicotine_interaction.Rdata")
save(top_genes_adults_nicotine_interaction,
     file="processed-data/04_DEA/Gene_analysis/top_genes_adults_nicotine_interaction.Rdata")



###########################################
#          Smoking Brain Adults DEA
#    Smoking vs ctrls in brain in adults
###########################################

RSE<-rse_gene_brain_adults_smoking

## Naive model
formula<- ~ Group + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr 
name<-"adults_smoking_naive"
coef<-"GroupExperimental"
results_adults_smoking_naive<-apply_DEA(RSE, formula, name, coef, "Group")
# "No differentially expressed genes"
top_genes_adults_smoking_naive<-results_adults_smoking_naive[[1]]
save(results_adults_smoking_naive, file="processed-data/04_DEA/Gene_analysis/results_adults_smoking_naive.Rdata")
save(top_genes_adults_smoking_naive, file="processed-data/04_DEA/Gene_analysis/top_genes_adults_smoking_naive.Rdata")

## Model adjusted for Pregnancy
formula<- ~ Group + Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate 
name<-"adults_smoking_fitted"
coef<-"GroupExperimental"
results_adults_smoking_fitted<-apply_DEA(RSE, formula, name, coef, "Group")
# "No differentially expressed genes"
top_genes_adults_smoking_fitted<-results_adults_smoking_fitted[[1]]
save(results_adults_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/results_adults_smoking_fitted.Rdata")
save(top_genes_adults_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_adults_smoking_fitted.Rdata")

## Model with interaction between Group and Pregnancy
formula<- ~ Group*Pregnancy +  Group + Pregnancy + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate 
name<-"adults_smoking_interaction"
coef<-"GroupExperimental:PregnancyYes"
results_adults_smoking_interaction<-apply_DEA(RSE, formula, name, coef, "Group")
# "No differentially expressed genes"
top_genes_adults_smoking_interaction<-results_adults_smoking_interaction[[1]]
save(results_adults_smoking_interaction,
     file="processed-data/04_DEA/Gene_analysis/results_adults_smoking_interaction.Rdata")
save(top_genes_adults_smoking_interaction,
     file="processed-data/04_DEA/Gene_analysis/top_genes_adults_smoking_interaction.Rdata")



##########################################
#         Nicotine Brain Pups DEA
#    Nicotine vs ctrls in brain in pups
##########################################

RSE<-rse_gene_brain_pups_nicotine ## 42 samples

## Naive model
formula<- ~ Group + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr + mitoRate
name<-"pups_nicotine_naive"
coef<-"GroupExperimental"
results_pups_nicotine_naive<-apply_DEA(RSE, formula, name, coef, "Group")
# "1038 differentially expressed genes"
top_genes_pups_nicotine_naive<-results_pups_nicotine_naive[[1]][[1]]
## DE genes
de_genes_pups_nicotine_naive<-results_pups_nicotine_naive[[2]]
save(results_pups_nicotine_naive, file="processed-data/04_DEA/Gene_analysis/results_pups_nicotine_naive.Rdata")
save(top_genes_pups_nicotine_naive, file="processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_naive.Rdata")
save(de_genes_pups_nicotine_naive, file="processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_naive.Rdata")


## x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
##                        Analysis of Sex differences
## x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
## Model adjusted for Sex
formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
name<-"pups_nicotine_fitted"
coef<-"GroupExperimental"
results_pups_nicotine_fitted<-apply_DEA(RSE, formula, name, coef, "Group")
# "1010 differentially expressed genes"
top_genes_pups_nicotine_fitted<-results_pups_nicotine_fitted[[1]][[1]]
## DE genes
de_genes_pups_nicotine_fitted<-results_pups_nicotine_fitted[[2]]
save(results_pups_nicotine_fitted, file="processed-data/04_DEA/Gene_analysis/results_pups_nicotine_fitted.Rdata")
save(top_genes_pups_nicotine_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata")
save(de_genes_pups_nicotine_fitted, file="processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata")
de_genes_brain_pup_nicotine <- de_genes_pups_nicotine_fitted[order(de_genes_pups_nicotine_fitted$adj.P.Val),]
write.table(de_genes_pups_nicotine_fitted, file = "processed-data/04_DEA/Gene_analysis/de_genes_brain_pup_nicotine.csv", row.names = FALSE, col.names = TRUE, sep = '\t')

## DGE analysis for Sex
name<-"pups_nicotine_sex"
coef<-"SexM"
results_pups_nicotine_sex<-apply_DEA(RSE, formula, name, coef, "Sex")
#  "743 differentially expressed genes"
top_genes_pups_nicotine_sex<-results_pups_nicotine_sex[[1]][[1]]
de_genes_pups_nicotine_sex<-results_pups_nicotine_sex[[2]]
save(results_pups_nicotine_sex, file="processed-data/04_DEA/Gene_analysis/results_pups_nicotine_sex.Rdata")
save(top_genes_pups_nicotine_sex, file="processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_sex.Rdata")
save(de_genes_pups_nicotine_sex, file="processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_sex.Rdata")
de_genes_brain_pup_nicotine_sex <- de_genes_pups_nicotine_sex[order(de_genes_pups_nicotine_sex$adj.P.Val),]
write.table(de_genes_brain_pup_nicotine_sex, file = "processed-data/04_DEA/Gene_analysis/de_genes_brain_pup_nicotine_sex.csv", row.names = FALSE, col.names = TRUE, sep = '\t')

## Intersection between Sex and Group DEGs (before and after adjusting for Sex)
length(intersect(de_genes_pups_nicotine_naive$Symbol, de_genes_pups_nicotine_sex$Symbol)) # 36/1038
length(intersect(de_genes_pups_nicotine_fitted$Symbol, de_genes_pups_nicotine_sex$Symbol)) # 81/1010

## DEGs lost after adjusting for Sex are sex DEGs?
degs_lost_byadjSex <- de_genes_pups_nicotine_naive$Symbol[which(! de_genes_pups_nicotine_naive$Symbol %in% de_genes_pups_nicotine_fitted$Symbol)]
length(intersect(degs_lost_byadjSex, de_genes_pups_nicotine_sex$Symbol)) # 14 / 304


## Model with interaction between Group and Sex
formula<- ~ Group*Sex + Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
name<-"pups_nicotine_interaction"
coef<-"GroupExperimental:SexM"
results_pups_nicotine_interaction<-apply_DEA(RSE, formula, name, coef, "Group")
# "No differentially expressed genes"
top_genes_pups_nicotine_interaction<-results_pups_nicotine_interaction[[1]]
save(results_pups_nicotine_interaction, 
     file="processed-data/04_DEA/Gene_analysis/results_pups_nicotine_interaction.Rdata")
save(top_genes_pups_nicotine_interaction,
     file="processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_interaction.Rdata")


## DGE for exposed vs ctrl in Females only
RSE_F <- RSE[, RSE$Sex == "F"] ## 22 samples (11 expt and 11 ctrl)
formula<- ~ Group + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr + mitoRate
name<-"pups_nicotine_females"
coef<-"GroupExperimental"
results_pups_nicotine_females<-apply_DEA(RSE_F, formula, name, coef, "Group")
# [1] "No differentially expressed genes"
top_genes_pups_nicotine_females<-results_pups_nicotine_females[[1]]
save(results_pups_nicotine_females, 
     file="processed-data/04_DEA/Gene_analysis/results_pups_nicotine_females.Rdata")
save(top_genes_pups_nicotine_females,
     file="processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_females.Rdata")

## DGE for exposed vs ctrl in Males only 
RSE_M <- RSE[, RSE$Sex == "M"] ## 20 samples (8 expt and 12 ctrl)
formula<- ~ Group + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr + mitoRate
name<-"pups_nicotine_males"
coef<-"GroupExperimental"
results_pups_nicotine_males<-apply_DEA(RSE_M, formula, name, coef, "Group")
# [1] "No differentially expressed genes"
top_genes_pups_nicotine_males<-results_pups_nicotine_males[[1]]
save(results_pups_nicotine_males, 
     file="processed-data/04_DEA/Gene_analysis/results_pups_nicotine_males.Rdata")
save(top_genes_pups_nicotine_males,
     file="processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_males.Rdata")



########################################
#         Smoking Brain Pups DEA
#    Smoking vs ctrls in brain in pups
########################################

RSE<-rse_gene_brain_pups_smoking ## 88 samples

## Naive model
formula<- ~ Group + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr + mitoRate
name<-"pups_smoking_naive"
coef<-"GroupExperimental"
results_pups_smoking_naive<-apply_DEA(RSE, formula, name, coef, "Group")
# "4108 differentially expressed genes"
top_genes_pups_smoking_naive<-results_pups_smoking_naive[[1]][[1]]
## DE genes
de_genes_pups_smoking_naive<-results_pups_smoking_naive[[2]]
save(results_pups_smoking_naive, file="processed-data/04_DEA/Gene_analysis/results_pups_smoking_naive.Rdata")
save(top_genes_pups_smoking_naive, file="processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_naive.Rdata")
save(de_genes_pups_smoking_naive, file="processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_naive.Rdata")


## x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
##                        Analysis of Sex differences
## x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
## Model adjusted for Sex
formula<- ~ Group + Sex + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr +  mitoRate
name<-"pups_smoking_fitted"
coef<-"GroupExperimental"
results_pups_smoking_fitted<-apply_DEA(RSE, formula, name, coef, "Group")
# "4165 differentially expressed genes"
top_genes_pups_smoking_fitted<-results_pups_smoking_fitted[[1]][[1]]
## DE genes
de_genes_pups_smoking_fitted<-results_pups_smoking_fitted[[2]]
save(results_pups_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/results_pups_smoking_fitted.Rdata")
save(top_genes_pups_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata")
save(de_genes_pups_smoking_fitted, file="processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata")
de_genes_brain_pup_smoking <- de_genes_pups_smoking_fitted[order(de_genes_pups_smoking_fitted$adj.P.Val),]
write.table(de_genes_brain_pup_smoking, file = "processed-data/04_DEA/Gene_analysis/de_genes_brain_pup_smoking.csv", row.names = FALSE, col.names = TRUE, sep = '\t')

## DGE analysis for Sex
name<-"pups_smoking_sex"
coef<-"SexM"
results_pups_smoking_sex<-apply_DEA(RSE, formula, name, coef, "Sex")
#  "55 differentially expressed genes"
top_genes_pups_smoking_sex<-results_pups_smoking_sex[[1]][[1]]
de_genes_pups_smoking_sex<-results_pups_smoking_sex[[2]]
save(results_pups_smoking_sex, file="processed-data/04_DEA/Gene_analysis/results_pups_smoking_sex.Rdata")
save(top_genes_pups_smoking_sex, file="processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_sex.Rdata")
save(de_genes_pups_smoking_sex, file="processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_sex.Rdata")
de_genes_brain_pup_smoking_sex <- de_genes_pups_smoking_sex[order(de_genes_pups_smoking_sex$adj.P.Val),]
write.table(de_genes_brain_pup_smoking_sex, file = "processed-data/04_DEA/Gene_analysis/de_genes_brain_pup_smoking_sex.csv", row.names = FALSE, col.names = TRUE, sep = '\t')


## Intersection between Sex and Group DEGs (before and after adjusting for Sex)
length(intersect(de_genes_pups_smoking_naive$Symbol, de_genes_pups_smoking_sex$Symbol)) # 7/4108
length(intersect(de_genes_pups_smoking_fitted$Symbol, de_genes_pups_smoking_sex$Symbol)) # 10/4165

## DEGs lost after adjusting for Sex are sex DEGs?
degs_lost_byadjSex <- de_genes_pups_smoking_naive$Symbol[which(! de_genes_pups_smoking_naive$Symbol %in% de_genes_pups_smoking_fitted$Symbol)]
length(intersect(degs_lost_byadjSex, de_genes_pups_smoking_sex$Symbol)) # 0


## Model with interaction between Group and Sex
formula<- ~ Group*Sex + Group + Sex + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr + mitoRate
name<-"pups_smoking_interaction"
coef<-"GroupExperimental:SexM"
results_pups_smoking_interaction<-apply_DEA(RSE, formula, name, coef, "Group")
# "No differentially expressed genes"
top_genes_pups_smoking_interaction<-results_pups_smoking_interaction[[1]]
save(results_pups_smoking_interaction, 
     file="processed-data/04_DEA/Gene_analysis/results_pups_smoking_interaction.Rdata")
save(top_genes_pups_smoking_interaction, 
     file="processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_interaction.Rdata")


## DGE for exposed vs ctrl in Females only
RSE_F <- RSE[, which(RSE$Sex == "F")] ## 44 samples (22 expt and 22 ctrl)
formula<- ~ Group + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr + mitoRate
name<-"pups_smoking_females"
coef<-"GroupExperimental"
results_pups_smoking_females<-apply_DEA(RSE_F, formula, name, coef, "Group")
# [1] "2686 differentially expressed genes"
top_genes_pups_smoking_females<-results_pups_smoking_females[[1]][[1]]
de_genes_pups_smoking_females<-results_pups_smoking_females[[2]]
save(results_pups_smoking_females, file="processed-data/04_DEA/Gene_analysis/results_pups_smoking_females.Rdata")
save(top_genes_pups_smoking_females, file="processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_females.Rdata")
save(de_genes_pups_smoking_females, file="processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_females.Rdata")
de_genes_brain_pup_smoking_females <- de_genes_pups_smoking_females[order(de_genes_pups_smoking_females$adj.P.Val),]
write.table(de_genes_brain_pup_smoking_females, file = "processed-data/04_DEA/Gene_analysis/de_genes_brain_pup_smoking_females.csv", row.names = FALSE, col.names = TRUE, sep = '\t')


## DGE for exposed vs ctrl in Males only 
RSE_M <- RSE[, which(RSE$Sex == "M")] ## 44 samples (20 expt and 24 ctrl)
formula<- ~ Group + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr + mitoRate
name<-"pups_smoking_males"
coef<-"GroupExperimental"
results_pups_smoking_males<-apply_DEA(RSE_M, formula, name, coef, "Group")
# [1] "886 differentially expressed genes"
top_genes_pups_smoking_males<-results_pups_smoking_males[[1]][[1]]
de_genes_pups_smoking_males<-results_pups_smoking_males[[2]]
save(results_pups_smoking_males, file="processed-data/04_DEA/Gene_analysis/results_pups_smoking_males.Rdata")
save(top_genes_pups_smoking_males, file="processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_males.Rdata")
save(de_genes_pups_smoking_males, file="processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_males.Rdata")
de_genes_brain_pup_smoking_males <- de_genes_pups_smoking_males[order(de_genes_pups_smoking_males$adj.P.Val),]
write.table(de_genes_brain_pup_smoking_males, file = "processed-data/04_DEA/Gene_analysis/de_genes_brain_pup_smoking_males.csv", row.names = FALSE, col.names = TRUE, sep = '\t')




## Condense results for Group DGE across all genes (don't take replication in human brain)

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
# version  R version 4.4.2 (2024-10-31)
# os       macOS Sequoia 15.3.2
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       Europe/London
# date     2025-04-02
# rstudio  2024.04.2+764 Chocolate Cosmos (desktop)
# pandoc   3.1.11 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
# quarto   1.4.555 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-8     2024-09-12 [1] CRAN (R 4.4.1)
# AnnotationDbi          1.68.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
# backports              1.5.0     2024-05-23 [1] CRAN (R 4.4.0)
# base64enc              0.1-3     2015-07-28 [1] CRAN (R 4.4.0)
# Biobase              * 2.66.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
# BiocFileCache          2.14.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
# BiocGenerics         * 0.52.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
# biomaRt                2.62.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
# biomartr             * 1.0.7     2023-12-02 [1] CRAN (R 4.4.0)
# Biostrings             2.74.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
# bit                    4.5.0.1   2024-12-03 [1] CRAN (R 4.4.1)
# bit64                  4.5.2     2024-09-22 [1] CRAN (R 4.4.1)
# bitops                 1.0-9     2024-10-03 [1] CRAN (R 4.4.1)
# blob                   1.2.4     2023-03-17 [1] CRAN (R 4.4.0)
# broom                  1.0.8     2025-03-28 [1] CRAN (R 4.4.1)
# cachem                 1.1.0     2024-05-16 [1] CRAN (R 4.4.0)
# car                    3.1-3     2024-09-27 [1] CRAN (R 4.4.1)
# carData                3.0-5     2022-01-06 [1] CRAN (R 4.4.1)
# cellranger             1.1.0     2016-07-27 [1] CRAN (R 4.4.0)
# checkmate              2.3.2     2024-07-29 [1] CRAN (R 4.4.0)
# cli                    3.6.3     2024-06-21 [1] CRAN (R 4.4.0)
# cluster                2.1.7     2024-12-08 [1] CRAN (R 4.4.1)
# colorspace             2.1-1     2024-07-26 [1] CRAN (R 4.4.0)
# cowplot              * 1.1.3     2024-01-22 [1] CRAN (R 4.4.0)
# crayon                 1.5.3     2024-06-20 [1] CRAN (R 4.4.0)
# curl                   6.0.1     2024-11-14 [1] CRAN (R 4.4.1)
# data.table             1.16.4    2024-12-06 [1] CRAN (R 4.4.1)
# DBI                    1.2.3     2024-06-02 [1] CRAN (R 4.4.0)
# dbplyr                 2.5.0     2024-03-19 [1] CRAN (R 4.4.0)
# DelayedArray           0.32.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
# digest                 0.6.37    2024-08-19 [1] CRAN (R 4.4.1)
# dplyr                  1.1.4     2023-11-17 [1] CRAN (R 4.4.0)
# edgeR                * 4.4.2     2025-01-27 [1] Bioconductor 3.20 (R 4.4.2)
# evaluate               1.0.1     2024-10-10 [1] CRAN (R 4.4.1)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.4.0)
# farver                 2.1.2     2024-05-13 [1] CRAN (R 4.4.0)
# fastmap                1.2.0     2024-05-15 [1] CRAN (R 4.4.0)
# filelock               1.0.3     2023-12-11 [1] CRAN (R 4.4.0)
# foreign                0.8-87    2024-06-26 [1] CRAN (R 4.4.2)
# formatR                1.14      2023-01-17 [1] CRAN (R 4.4.0)
# Formula                1.2-5     2023-02-24 [1] CRAN (R 4.4.1)
# fs                     1.6.5     2024-10-30 [1] CRAN (R 4.4.1)
# futile.logger        * 1.4.3     2016-07-10 [1] CRAN (R 4.4.0)
# futile.options         1.0.1     2018-04-20 [1] CRAN (R 4.4.0)
# gargle                 1.5.2     2023-07-20 [1] CRAN (R 4.4.0)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.4.0)
# GenomeInfoDb         * 1.42.1    2024-11-28 [1] Bioconductor 3.20 (R 4.4.2)
# GenomeInfoDbData       1.2.13    2024-12-05 [1] Bioconductor
# GenomicRanges        * 1.58.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
# ggplot2              * 3.5.1     2024-04-23 [1] CRAN (R 4.4.0)
# ggpubr                 0.6.0     2023-02-10 [1] CRAN (R 4.4.0)
# ggrepel              * 0.9.6     2024-09-07 [1] CRAN (R 4.4.1)
# ggsignif               0.6.4     2022-10-13 [1] CRAN (R 4.4.0)
# glue                   1.8.0     2024-09-30 [1] CRAN (R 4.4.1)
# googledrive            2.1.1     2023-06-11 [1] CRAN (R 4.4.0)
# gridExtra            * 2.3       2017-09-09 [1] CRAN (R 4.4.1)
# gtable                 0.3.6     2024-10-25 [1] CRAN (R 4.4.1)
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.4.1)
# Hmisc                * 5.2-3     2025-03-16 [1] CRAN (R 4.4.1)
# hms                    1.1.3     2023-03-21 [1] CRAN (R 4.4.0)
# htmlTable              2.4.3     2024-07-21 [1] CRAN (R 4.4.0)
# htmltools              0.5.8.1   2024-04-04 [1] CRAN (R 4.4.0)
# htmlwidgets            1.6.4     2023-12-06 [1] CRAN (R 4.4.0)
# httr                   1.4.7     2023-08-15 [1] CRAN (R 4.4.0)
# httr2                  1.0.7     2024-11-26 [1] CRAN (R 4.4.1)
# IRanges              * 2.40.1    2024-12-05 [1] Bioconductor 3.20 (R 4.4.2)
# jaffelab             * 0.99.34   2025-03-31 [1] Github (LieberInstitute/jaffelab@e7a26ca)
# jsonlite               1.8.9     2024-09-20 [1] CRAN (R 4.4.1)
# KEGGREST               1.46.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
# knitr                  1.49      2024-11-08 [1] CRAN (R 4.4.1)
# lambda.r               1.2.4     2019-09-18 [1] CRAN (R 4.4.0)
# lattice                0.22-6    2024-03-20 [1] CRAN (R 4.4.2)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.4.0)
# limma                * 3.62.2    2025-01-09 [1] Bioconductor 3.20 (R 4.4.2)
# locfit                 1.5-9.11  2025-02-03 [1] CRAN (R 4.4.1)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.4.0)
# MASS                   7.3-61    2024-06-13 [1] CRAN (R 4.4.2)
# Matrix                 1.7-1     2024-10-18 [1] CRAN (R 4.4.2)
# MatrixGenerics       * 1.18.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
# matrixStats          * 1.4.1     2024-09-08 [1] CRAN (R 4.4.1)
# memoise                2.0.1     2021-11-26 [1] CRAN (R 4.4.0)
# munsell                0.5.1     2024-04-01 [1] CRAN (R 4.4.0)
# nlme                   3.1-166   2024-08-14 [1] CRAN (R 4.4.2)
# nnet                   7.3-19    2023-05-03 [1] CRAN (R 4.4.2)
# patchwork              1.3.0     2024-09-16 [1] CRAN (R 4.4.1)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.4.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.4.0)
# png                    0.1-8     2022-11-29 [1] CRAN (R 4.4.0)
# prettyunits            1.2.0     2023-09-24 [1] CRAN (R 4.4.0)
# progress               1.2.3     2023-12-06 [1] CRAN (R 4.4.0)
# purrr                  1.0.2     2023-08-10 [1] CRAN (R 4.4.0)
# pwr                    1.3-0     2020-03-17 [1] CRAN (R 4.4.1)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.4.0)
# rafalib              * 1.0.3     2025-03-25 [1] CRAN (R 4.4.1)
# rappdirs               0.3.3     2021-01-31 [1] CRAN (R 4.4.0)
# RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.4.0)
# Rcpp                   1.0.13-1  2024-11-02 [1] CRAN (R 4.4.1)
# RCurl                  1.98-1.16 2024-07-11 [1] CRAN (R 4.4.0)
# readxl               * 1.4.5     2025-03-07 [1] CRAN (R 4.4.1)
# rlang                * 1.1.5     2025-01-17 [1] CRAN (R 4.4.1)
# rmarkdown              2.29      2024-11-04 [1] CRAN (R 4.4.1)
# rpart                  4.1.23    2023-12-05 [1] CRAN (R 4.4.2)
# rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.4.0)
# RSQLite                2.3.9     2024-12-03 [1] CRAN (R 4.4.1)
# rstatix                0.7.2     2023-02-01 [1] CRAN (R 4.4.0)
# rstudioapi             0.17.1    2024-10-22 [1] CRAN (R 4.4.1)
# S4Arrays               1.6.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
# S4Vectors            * 0.44.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.4.0)
# segmented              2.1-4     2025-02-28 [1] CRAN (R 4.4.1)
# sessioninfo          * 1.2.3     2025-02-05 [1] CRAN (R 4.4.1)
# smplot2              * 0.2.5     2025-04-02 [1] Github (smin95/smplot2@5c892f6)
# SparseArray            1.6.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
# statmod                1.5.0     2023-01-06 [1] CRAN (R 4.4.1)
# stringi                1.8.4     2024-05-06 [1] CRAN (R 4.4.0)
# stringr                1.5.1     2023-11-14 [1] CRAN (R 4.4.0)
# SummarizedExperiment * 1.36.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
# tibble                 3.2.1     2023-03-20 [1] CRAN (R 4.4.0)
# tidyr                  1.3.1     2024-01-24 [1] CRAN (R 4.4.0)
# tidyselect             1.2.1     2024-03-11 [1] CRAN (R 4.4.0)
# UCSC.utils             1.2.0     2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.4.0)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.4.0)
# VennDiagram          * 1.7.3     2022-04-12 [1] CRAN (R 4.4.0)
# withr                  3.0.2     2024-10-28 [1] CRAN (R 4.4.1)
# xfun                   0.49      2024-10-31 [1] CRAN (R 4.4.1)
# XML                    3.99-0.17 2024-06-25 [1] CRAN (R 4.4.0)
# xml2                   1.3.6     2023-12-04 [1] CRAN (R 4.4.0)
# XVector                0.46.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
# zlibbioc               1.52.0    2024-11-08 [1] Bioconductor 3.20 (R 4.4.1)
# zoo                    1.8-12    2023-04-13 [1] CRAN (R 4.4.0)
# 
# [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library
# * ── Packages attached to the search path.
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

