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
  exp_vars<-getVarianceExplained(RSE, variables=c("Sex", "Pregnancy", 
                    "Expt", "Group", "plate", "flowcell", "mitoRate", "rRNA_rate", 
                    "overallMapRate", "totalAssignedGene", "ERCCsumLogErr"), 
                    exprs_values = "logcounts")
  ## Plot explanatory variables ordered by % of variance explained
  plotExplanatoryVariables(exp_vars, theme_size = 10, nvars_to_plot = 12)
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







