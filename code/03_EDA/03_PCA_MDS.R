
# 1. Explore sample effects
## 1.1 Dimensionality Reduction (PCA)
### 1.1.1 Explore Samples' expression variation

## Load data
load(here("processed-data/03_EDA/02_QC/rse_gene_blood_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_gene_brain_adults_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_exon_brain_adults_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_tx_brain_adults_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_jx_brain_adults_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_gene_brain_pups_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_exon_brain_pups_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_tx_brain_pups_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_jx_brain_pups_qc.Rdata"))

## Generate PCA data
PCA<-function(tissue, type, age){
  if (is.null(age)) {
    RSE<-eval(parse_expr(paste("rse", type, tissue, "qc", sep="_"))) }
  else {
    RSE<-eval(parse_expr(paste("rse", type, tissue, age, "qc", sep="_")))}
  pca<-prcomp(t(assays(RSE)$logcounts))
  ## % of the variance explained by each PC
  pca_vars<- getPcaVars(pca)
  pca_vars_labs<- paste0(
      "PC", seq(along = pca_vars), ": ",
      pca_vars, "% Var Expl")
  
  ## Join PCs and samples' info 
  pca_data<-cbind(pca$x,colData(RSE))
  ## Add samples' phenotypes
  pca_data<-as.data.frame(pca_data)
  return(list(pca_data, pca_vars_labs))
}


## PCx vs PCy Plots 
PCx_vs_PCy <- function (PCx, PCy, pca_data, pca_vars, pheno_var) {
  plot=ggplot(data=pca_data, 
       aes(x=eval(parse_expr(PCx)),y=eval(parse_expr(PCy)),
           color=eval(parse_expr(pheno_var))) )+ 
       theme(legend.position="right", plot.margin=unit (c (1,1.5,1,1), 'cm')) +
       geom_point() + 
       labs(x= pca_vars[strtoi(gsub("PC","", PCx))], y = pca_vars[strtoi(gsub("PC","", PCy))],  
            color=pheno_var)
  return(plot)
}

## All PCA plots 
plot_PCAs<-function(type, tissue, age){
  pca_data<-PCA(tissue, type, age)[[1]]
  pca_vars<-PCA(tissue, type, age)[[2]]
  ## Plot for type and tissue
  if (is.null(age)){
    for (PCs in list(c("PC1", "PC2"), c("PC1", "PC3"))){
      plots<-list()
      i=1
        for (pheno_var in c("Age", "plate","Expt", "Sex", "Group", "Pregnancy", "medium", "flowcell")){
           p<-PCx_vs_PCy(PCs[1], PCs[2], pca_data, pca_vars, pheno_var)
           plots[[i]]=p
           i=i+1
        }
      plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]],
                nrow = 2)
      ## Save plots
      ggsave(paste("plots/03_EDA/03_PCA_MDS/",PCs[1],"_vs_",PCs[2],"_", type, "_", tissue ,".pdf", sep=""), 
             width = 65, height = 25, units = "cm")
    }
  }
  ## Plot for type and age
  else {
    for (PCs in list(c("PC1", "PC2"), c("PC1", "PC3"))){
      plots<-list()
      i=1
      for (pheno_var in c("plate","Expt", "Sex", "Group", "Pregnancy", "flowcell")){
         p<-PCx_vs_PCy(PCs[1], PCs[2], pca_data, pca_vars, pheno_var)
         plots[[i]]=p
         i=i+1
      }
      plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]],nrow = 2)
      ## Save plots
      ggsave(paste("plots/03_EDA/03_PCA_MDS/",PCs[1],"_vs_",PCs[2],"_", type, "_", tissue ,"_", age, ".pdf", sep=""), 
             width = 45, height = 25, units = "cm")
    }
  }
}

## Plots
plot_PCAs("gene", "blood", NULL)
plot_PCAs("gene", "brain", "adults")
plot_PCAs("gene", "brain", "pups")
plot_PCAs("exon", "brain", "adults")
plot_PCAs("exon", "brain", "pups")
plot_PCAs("tx", "brain", "adults")
plot_PCAs("tx", "brain", "pups")
plot_PCAs("jx", "brain", "adults")
plot_PCAs("jx", "brain", "pups")
