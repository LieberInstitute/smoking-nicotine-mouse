
# 1. Explore gene level effects


## Load data without outliers and rare samples
load(here("processed-data/03_EDA/02_QC/rse_gene_blood_qc.Rdata"))
load(here("processed-data/03_EDA/03_PCA/rse_gene_brain_adults_qc_afterPCA.Rdata"))
load(here("processed-data/03_EDA/03_PCA/rse_gene_brain_pups_qc_afterPCA.Rdata"))
load(here("processed-data/03_EDA/03_PCA/rse_exon_brain_adults_qc_afterPCA.Rdata"))
load(here("processed-data/03_EDA/03_PCA/rse_exon_brain_pups_qc_afterPCA.Rdata"))
load(here("processed-data/03_EDA/03_PCA/rse_tx_brain_adults_qc_afterPCA.Rdata"))
load(here("processed-data/03_EDA/03_PCA/rse_tx_brain_pups_qc_afterPCA.Rdata"))
load(here("processed-data/03_EDA/03_PCA/rse_jx_brain_adults_qc_afterPCA.Rdata"))
load(here("processed-data/03_EDA/03_PCA/rse_jx_brain_pups_qc_afterPCA.Rdata"))

## Reduce objects' names
rse_gene_brain_adults_qc <- rse_gene_brain_adults_qc_afterPCA
rse_gene_brain_pups_qc <- rse_gene_brain_pups_qc_afterPCA
rse_exon_brain_adults_qc <- rse_exon_brain_adults_qc_afterPCA
rse_exon_brain_pups_qc <- rse_exon_brain_pups_qc_afterPCA
rse_tx_brain_adults_qc <- rse_tx_brain_adults_qc_afterPCA
rse_tx_brain_pups_qc <- rse_tx_brain_pups_qc_afterPCA
rse_jx_brain_adults_qc <- rse_jx_brain_adults_qc_afterPCA
rse_jx_brain_pups_qc <- rse_jx_brain_pups_qc_afterPCA


## Brain data separated by type of experiment
rse_gene_brain_adults_smoking<-rse_gene_brain_adults_qc[,rse_gene_brain_adults_qc$Expt=="Smoking"]
save(rse_gene_brain_adults_smoking, 
     file="processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_adults_smoking.Rdata")
rse_gene_brain_adults_nicotine<-rse_gene_brain_adults_qc[,rse_gene_brain_adults_qc$Expt=="Nicotine"]
save(rse_gene_brain_adults_nicotine,
     file="processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_adults_nicotine.Rdata")
rse_gene_brain_pups_smoking<-rse_gene_brain_pups_qc[,rse_gene_brain_pups_qc$Expt=="Smoking"]
save(rse_gene_brain_pups_smoking, 
     file="processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_smoking.Rdata")
rse_gene_brain_pups_nicotine<-rse_gene_brain_pups_qc[,rse_gene_brain_pups_qc$Expt=="Nicotine"]
save(rse_gene_brain_pups_nicotine, 
     file="processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_nicotine.Rdata")


# ------------------------------------------------------------------------------
## 1.1 Explanatory Variables
# ------------------------------------------------------------------------------

## Plot density function for % of variance explained 
expl_var<- function(type, tissue, age, expt){
  
  ## For blood or all adults/pups data
  if (is.null(expt)){
    if (is.null(age)){ 
      RSE<-eval(parse_expr(paste("rse", type, tissue, "qc", sep="_")))
      fileName=paste("plots/03_EDA/04_Expl_Var_partition/Expl_vars_density_",type,"_",tissue,".pdf",
                     sep="")
    }
    else {
      RSE<-eval(parse_expr(paste("rse", type, tissue, age, "qc", sep="_")))
      fileName=paste("plots/03_EDA/04_Expl_Var_partition/Expl_vars_density_",type,"_",tissue,"_", 
                     age, ".pdf", sep="")
    }
  }
 
  ## Smoking/nicotine brain data for adults/pups
  else {
    RSE<-eval(parse_expr(paste("rse", type, tissue, age, expt, sep="_")))
    fileName=paste("plots/03_EDA/04_Expl_Var_partition/Expl_vars_density_",type,"_",tissue, "_", 
                   age, "_", expt, ".pdf", sep="")
  }
  
  ## % of variance in gene expression explained by each variable
  variables=c("Sex", "Pregnancy","Expt", "Group", "plate", "flowcell", 
                                 "mitoRate", "rRNA_rate","overallMapRate", "totalAssignedGene", "ERCCsumLogErr")
  exp_vars<-getVarianceExplained(RSE, variables=variables, exprs_values = "logcounts")
  ## Plot explanatory variables 
  ## Variables' colors
  colors=c("Sex"="tomato3", "Pregnancy"="lightpink", "Expt"="slateblue2","Group"="seagreen3","plate"="tan1",
           "flowcell"="mediumorchid1", "mitoRate"="royalblue1","rRNA_rate"="yellow3","overallMapRate"="wheat3", 
           "totalAssignedGene"="violetred1", "ERCCsumLogErr"="ivory4")
  values<-c()
  for (variable in variables){
    levels<-eval(parse_expr(paste("RSE$",variable, sep="")))
    ## Drop variables with fewer than 2 unique levels
    if (length(table(levels))>1){
      values<-append(values, variable)
    }
  }
  ## Plot density graph for each variable
  p<-plotExplanatoryVariables(exp_vars, theme_size = 20, nvars_to_plot = 12)
  p + scale_colour_manual(values = colors[c(values)]) + labs(color="Variables")
  ## Save plot
  ggsave(fileName, width = 25, height = 15, units = "cm")
  return(NULL)

}

## Plots
expl_var("gene", "blood", NULL, NULL)
expl_var("gene", "brain", "adults", NULL)
expl_var("gene", "brain", "adults", "smoking")
expl_var("gene", "brain", "adults", "nicotine")
expl_var("gene", "brain", "pups", NULL)
expl_var("gene", "brain", "pups", "smoking")
expl_var("gene", "brain", "pups", "nicotine")





# ------------------------------------------------------------------------------
## 1.2 Variance Partition
# ------------------------------------------------------------------------------

### 1.2.1 Canonical Correlation Analysis (CCA) 

## Plot Heatmap of CC between variables
plot_CCA<- function(tissue, age, expt){
  
  type <- 'gene'
  
  if (is.null(expt)){
    
    ## For blood 
    if (is.null(age)){
      RSE<-eval(parse_expr(paste("rse", type, tissue, "qc", sep="_")))
      formula <- ~ Pregnancy + Group + plate + flowcell + mitoRate + rRNA_rate + 
      overallMapRate + totalAssignedGene + ERCCsumLogErr
      
    }
    
    ## For adults
    else if (age=="adults"){
      RSE<-eval(parse_expr(paste("rse", type, tissue, age, "qc", sep="_")))
      formula <- ~ Pregnancy + Expt + Group + plate + flowcell + mitoRate + rRNA_rate + 
      overallMapRate + totalAssignedGene + ERCCsumLogErr      
    }
    
    ## For pups
    else if (age=="pups"){
      RSE<-eval(parse_expr(paste("rse", type, tissue, age, "qc", sep="_")))
      formula <- ~ Sex + Expt + Group + plate + flowcell + mitoRate + rRNA_rate + 
      overallMapRate + totalAssignedGene + ERCCsumLogErr     
    }
  }
  
  ## For brain adults/pups and smoking/nicotine
  else {
    RSE<-eval(parse_expr(paste("rse_gene", tissue, age, expt, sep="_")))
    if (age=="adults"){
      formula <- ~ Pregnancy + Group + plate + flowcell + mitoRate + rRNA_rate + 
      overallMapRate + totalAssignedGene + ERCCsumLogErr         
    }
    else{
      formula <- ~ Sex + Group + plate + flowcell + mitoRate + rRNA_rate + 
      overallMapRate + totalAssignedGene + ERCCsumLogErr   
    }

  }

  ## Assess correlation between all pairs of variables
  C=canCorPairs(formula, colData(RSE))
  ## Heatmap
  plotCorrMatrix(C)
  
  return(NULL)
}

## Plots
plot_CCA("blood", NULL, NULL)
plot_CCA("brain", "adults", NULL)
plot_CCA("brain", "pups", NULL)
plot_CCA("brain", "adults", "smoking")
plot_CCA("brain", "adults", "nicotine")
plot_CCA("brain", "pups", "smoking")
plot_CCA("brain", "pups", "nicotine")




### 1.2.2 Model fit

var_part<-function(tissue, age, expt) {
  if (is.null(expt)){
    ## For blood
    if (is.null(age)){
      RSE<-rse_gene_blood_qc
      ## Expression values are the response
      ## Specify categorical variables as random effects for lmer
      ## Try with all variables
      formula <- ~ (1|Pregnancy) + (1|Group) + (1|plate) + (1|flowcell) + mitoRate + rRNA_rate + 
                   overallMapRate + totalAssignedGene + ERCCsumLogErr
      fileName<-paste("plots/03_EDA/04_Expl_Var_partition/ViolinPlot_", tissue, ".pdf", sep="")
    }
    
    ## For all adults and pups
    else{
      fileName<-paste("plots/03_EDA/04_Expl_Var_partition/ViolinPlot_", tissue, "_", age, ".pdf", sep="")
      RSE<-eval(parse_expr(paste("rse", type, tissue, age, "qc", sep="_")))
      
      if (age=="adults"){
        formula <- ~ (1|Pregnancy) + (1|Expt) + (1|Group) + (1|plate) + (1|flowcell) + mitoRate + rRNA_rate + 
                      overallMapRate + totalAssignedGene + ERCCsumLogErr}
      else{
        formula <- ~ (1|Sex) + (1|Expt) + (1|Group) + (1|plate) + (1|flowcell) + mitoRate + rRNA_rate + 
                      overallMapRate + totalAssignedGene + ERCCsumLogErr}
    }
  }
  
  else {
    RSE<-eval(parse_expr(paste("rse_gene", tissue, age, expt, sep="_")))
    fileName<-paste("plots/03_EDA/04_Expl_Var_partition/ViolinPlot_", tissue, "_", 
                       age, "_", expt, ".pdf", sep="")
    ## For brain adults
    if (age=="adults"){
      formula <- ~ (1|Pregnancy) + (1|Group) + (1|plate) + (1|flowcell) + mitoRate + rRNA_rate + 
        overallMapRate + totalAssignedGene + ERCCsumLogErr}

    ## For brain pups
    if (age=="pups"){
      formula <- ~ (1|Sex) + (1|Group) + (1|plate) + (1|flowcell) + mitoRate + rRNA_rate + overallMapRate +
        totalAssignedGene + ERCCsumLogErr}
  }

  
  ## Genes with variance of 0
  genes_var_zero<-which(apply(assays(RSE)$logcounts, 1, var)==0)
  
  if (length(genes_var_zero)>0){
    ## Fit linear mixed model without those genes
    varPart <- fitExtractVarPartModel(assays(RSE)$logcounts[-genes_var_zero,],formula, colData(RSE))
  }
  else {
  ## Fit linear mixed model with all genes
    varPart <- fitExtractVarPartModel(assays(RSE)$logcounts,formula, colData(RSE))
  }
  
  ## Sort variables by median fraction of variance explained
  sort_vars <- sortCols(varPart)
  # Violin plot of contribution of each variable to total variance
  p<-plotVarPart(sort_vars, label.angle=60)
  ggsave(fileName,  p, width = 40, height = 20, units = "cm")

}

var_part("blood", NULL, NULL)
var_part("brain", "adults", NULL)
var_part("brain", "pups", NULL)
var_part("brain", "adults", "smoking")
var_part("brain", "adults", "nicotine")
var_part("brain", "pups", "smoking")
var_part("brain", "pups", "nicotine")




### 1.2.3 Variation within subsets 

var_part_subsets<-function(tissue, age, expt) {
  if (is.null(expt)){
    ## For blood
    if (is.null(age)){
      RSE<-rse_gene_blood_qc
      ## Specify formula to model Group variance separately for each Pregnancy group
      ## Drop highly correlated variables 
      formula <- ~ (Pregnancy+0|Group) + (1|Pregnancy) + (1|plate) + (1|flowcell) + rRNA_rate + 
                    totalAssignedGene + ERCCsumLogErr + overallMapRate
      fileName<-paste("plots/03_EDA/04_Expl_Var_partition/ViolinPlot_subsets_", tissue, ".pdf", sep="")
    }
    
    ## For all adults and pups
    else {
      RSE<-eval(parse_expr(paste("rse", type, tissue, age, "qc", sep="_")))
      fileName<-paste("plots/03_EDA/04_Expl_Var_partition/ViolinPlot_subsets_", tissue, "_", age, ".pdf", sep="")
     
      if (age=="adults"){
        ## Specify formula to model Group variance separately for each Expt group
        formula <- ~ (1|Pregnancy) + (Expt+0|Group) + (1|Expt) + (1|plate) + (1|flowcell) + totalAssignedGene + 
                      rRNA_rate + overallMapRate + ERCCsumLogErr}
      else{
        formula <- ~ (1|Sex) + (Expt+0|Group) + (1|Expt) +  (1|plate) + (1|flowcell) + totalAssignedGene +
                      rRNA_rate + overallMapRate + ERCCsumLogErr + mitoRate}
    }
  }
  
  else {
    RSE<-eval(parse_expr(paste("rse_gene", tissue, age, expt, sep="_")))
    fileName<-paste("plots/03_EDA/04_Expl_Var_partition/ViolinPlot_subsets_", tissue, "_", 
                       age, "_", expt, ".pdf", sep="")
    ## For brain adults
    if (age=="adults"){
      ## Group variance for each Pregnancy group
      formula <- ~ (Pregnancy+0|Group) + (1|Pregnancy) + (1|plate) + (1|flowcell) + rRNA_rate + 
                   overallMapRate + totalAssignedGene + ERCCsumLogErr
      }
      
    ## For brain pups
    if (age=="pups"){
      ## Group variance for each sex
      formula <- ~ (Sex+0|Group) + (1|Sex) + (1|plate) + (1|flowcell) + rRNA_rate + overallMapRate +
                    totalAssignedGene + ERCCsumLogErr + mitoRate
    }
  }

  
  ## Genes with variance of 0
  genes_var_zero<-which(apply(assays(RSE)$logcounts, 1, var)==0)
  
  if (length(genes_var_zero)>0){
    ## Fit linear mixed model without those genes
    varPart <- fitExtractVarPartModel(assays(RSE)$logcounts[-genes_var_zero,],formula, colData(RSE))
  }
  else {
  ## Fit linear mixed model with all genes
    varPart <- fitExtractVarPartModel(assays(RSE)$logcounts,formula, colData(RSE))
  }
  
  ## Sort variables by median fraction of variance explained
  sort_vars <- sortCols(varPart)
  # Violin plot of contribution of each variable to total variance
  p<-plotVarPart(sort_vars, label.angle=60)
  if (!is.null(age) && is.null(expt)){
    ggsave(fileName,  p, width = 35, height = 18, units = "cm")
  }
  else{
    ggsave(fileName,  p, width = 40, height = 20, units = "cm")
  }

}


var_part_subsets("blood", NULL, NULL)
var_part_subsets("brain", "adults", NULL)
var_part_subsets("brain", "pups", NULL)
var_part_subsets("brain", "adults", "smoking")
var_part_subsets("brain", "adults", "nicotine")
var_part_subsets("brain", "pups", "smoking")
var_part_subsets("brain", "pups", "nicotine")







## Reproducibility information
options(width = 120)
session_info()
# version  R version 4.2.0 (2022-04-22 ucrt)
#  os       Windows 10 x64 (build 19044)
#  system   x86_64, mingw32
#  ui       RStudio
#  language (EN)
#  collate  Spanish_Mexico.utf8
#  ctype    Spanish_Mexico.utf8
#  tz       America/Mexico_City
#  date     2022-08-06
#  rstudio  2022.02.3+492 Prairie Trillium (desktop)
#  pandoc   2.17.1.1 @ C:/Program Files/RStudio/bin/quarto/bin/ (via rmarkdown)
