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


## 1.1 Explanatory Variables
## Plot density function for % of variance explained 
expl_var<- function(type, tissue, age){
  if (is.null(age)){
    RSE<-eval(parse_expr(paste("rse", type, tissue, "qc", sep="_")))
  }
  else {
    RSE<-eval(parse_expr(paste("rse", type, tissue, age, "qc", sep="_")))
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
  if (is.null(age)){
    fileName=paste("plots/03_EDA/04_Expl_Var_partition/Expl_vars_density_",type,"_",tissue,".pdf", sep="")
  }
  else {
    fileName=paste("plots/03_EDA/04_Expl_Var_partition/Expl_vars_density_",type,"_",tissue, "_", 
                   age, ".pdf", sep="")
  }
  ggsave(fileName, width = 25, height = 15, units = "cm")
  return(NULL)

}

## Plots
expl_var("gene", "blood", NULL)
expl_var("gene", "brain", "adults")
expl_var("gene", "brain", "pups")
