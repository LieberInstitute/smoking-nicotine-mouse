# 1. Explore gene level effects

## Load data without outliers and rare samples
load(here("processed-data/03_EDA/02_QC/rse_gene_blood_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_gene_brain_adults_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_gene_brain_pups_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_exon_brain_adults_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_exon_brain_pups_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_tx_brain_adults_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_tx_brain_pups_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_jx_brain_adults_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_jx_brain_pups_qc.Rdata"))

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


## 1.1 Explanatory Variables
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
           "flowcell"="mediumorchid1", "mitoRate"="royalblue1","rRNA_rate"="yellow3","overallMapRate"="wheat4", 
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
  p<-plotExplanatoryVariables(exp_vars, theme_size = 10, nvars_to_plot = 12)
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







## 1.2 Variance Partition
### 1.2.1 Canonical Correlation Analysis (CCA) 

## Plot Heatmap of CC between variables
plot_CCA<- function(tissue, age, expt){
  ## For blood
  if (is.null(expt)){
    RSE<-rse_gene_blood_qc
    formula <- ~ Pregnancy + Group + plate + flowcell + mitoRate + rRNA_rate + 
      overallMapRate + totalAssignedGene + ERCCsumLogErr
  }
  
  ## For brain adults
  else if (age=="adults"){
    formula <- ~ Pregnancy + Group + plate + flowcell + mitoRate + rRNA_rate + overallMapRate +
      totalAssignedGene + ERCCsumLogErr
    if ( expt=="smoking") {
      RSE<-eval(parse_expr(paste("rse_gene_brain", age, "smoking", sep="_")))
    }
    else {
      RSE<-eval(parse_expr(paste("rse_gene_brain", age, "nicotine", sep="_")))
    }
  }

  
  ## For brain pups
  else if (age=="pups"){
    formula <- ~ Sex + Group + plate + flowcell + mitoRate + rRNA_rate + overallMapRate +
      totalAssignedGene + ERCCsumLogErr
    if ( expt=="smoking") {
      RSE<-eval(parse_expr(paste("rse_gene_brain", age, "smoking", sep="_")))
    }
    else {
      RSE<-eval(parse_expr(paste("rse_gene_brain", age, "nicotine", sep="_")))
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
plot_CCA("brain", "adults", "smoking")
plot_CCA("brain", "adults", "nicotine")
plot_CCA("brain", "pups", "smoking")
plot_CCA("brain", "pups", "nicotine")




