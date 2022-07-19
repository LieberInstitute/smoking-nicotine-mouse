
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




## Rare samples in brain PCA plots
## PC data
pca_data_gene_brain_adults<-PCA("brain", "gene", "adults")[[1]]
pca_data_exon_brain_adults<-PCA("brain", "exon", "adults")[[1]]
pca_data_tx_brain_adults<-PCA("brain", "tx", "adults")[[1]]
pca_data_jx_brain_adults<-PCA("brain", "jx", "adults")[[1]]
pca_data_gene_brain_pups<-PCA("brain", "gene", "pups")[[1]]
pca_data_exon_brain_pups<-PCA("brain", "exon", "pups")[[1]]
pca_data_tx_brain_pups<-PCA("brain", "tx", "pups")[[1]]
pca_data_jx_brain_pups<-PCA("brain", "jx", "pups")[[1]]

## Rare samples in adult plots
pca_data_gene_brain_adults[which.max(pca_data_gene_brain_adults$PC2), "SAMPLE_ID"]
# "Sample_FE3P2"
pca_data_gene_brain_adults[which.max(pca_data_gene_brain_adults$PC3), "SAMPLE_ID"]
# "Sample_FE3P2"
pca_data_exon_brain_adults[which.max(pca_data_exon_brain_adults$PC2), "SAMPLE_ID"]
# "Sample_FE3P2"
pca_data_exon_brain_adults[which.max(pca_data_exon_brain_adults$PC3), "SAMPLE_ID"]
# "Sample_4067"
pca_data_tx_brain_adults[which.max(pca_data_tx_brain_adults$PC2), "SAMPLE_ID"]
# "Sample_FE3P2"
pca_data_jx_brain_adults[which.min(pca_data_jx_brain_adults$PC3), "SAMPLE_ID"]
# "Sample_FE3P2"

## Rare samples in pup plots
pca_data_gene_brain_pups[which.min(pca_data_gene_brain_pups$PC2), "SAMPLE_ID"]
# "Sample_P2_fe2_022019"
## Male samples in females at gene level 
pca_data_gene_brain_pups[which(pca_data_gene_brain_pups$Sex=="M"& pca_data_gene_brain_pups$PC3>0),"SAMPLE_ID"]
# "Sample_P1_fe3_021819" "Sample_P2_fe2_022019"
pca_data_exon_brain_pups[which.min(pca_data_exon_brain_pups$PC1), "SAMPLE_ID"]
# "Sample_P2_fe2_022019"
pca_data_tx_brain_pups[which.min(pca_data_tx_brain_pups$PC1), "SAMPLE_ID"]
# "Sample_P2_fe2_022019"
pca_data_jx_brain_pups[which.max(pca_data_jx_brain_pups$PC1), "SAMPLE_ID"]
# "Sample_P2_fe2_022019"

## Explore Sample_FE3P2 info:
colData(rse_gene_brain_adults_qc)[which(rse_gene_brain_adults_qc$SAMPLE_ID=="Sample_FE3P2"),
                          c("mitoRate", "rRNA_rate", "totalAssignedGene", "overallMapRate", 
                            "ERCCsumLogErr", "sum", "Sex")]
## This sample has the max rRNA rate
pca_data_gene_brain_adults[which.max(pca_data_gene_brain_adults$rRNA_rate), "SAMPLE_ID"]
# "Sample_FE3P2"


## Explore Sample_4067 info:
colData(rse_gene_brain_adults_qc)[which(rse_gene_brain_adults_qc$SAMPLE_ID=="Sample_4067"),
                          c("mitoRate", "rRNA_rate", "totalAssignedGene", "overallMapRate", 
                            "ERCCsumLogErr", "sum", "Sex")]
## This sample has the max mito rate, % of mt and ribo counts and min prop of reads assigned to genes
pca_data_gene_brain_adults[which.max(pca_data_gene_brain_adults$mitoRate), "SAMPLE_ID"]
# "Sample_4067"
pca_data_gene_brain_adults[which.max(pca_data_gene_brain_adults$subsets_Mito_percent), "SAMPLE_ID"]
# "Sample_4067"
pca_data_gene_brain_adults[which.max(pca_data_gene_brain_adults$subsets_Ribo_percent), "SAMPLE_ID"]
# "Sample_4067"
pca_data_gene_brain_adults[which.min(pca_data_gene_brain_adults$totalAssignedGene), "SAMPLE_ID"]
# "Sample_4067"



## Explore Sample_P2_fe2_022019 info:
colData(rse_gene_brain_pups_qc)[which(rse_gene_brain_pups_qc$SAMPLE_ID=="Sample_P2_fe2_022019"),
                          c("mitoRate", "rRNA_rate", "totalAssignedGene", "overallMapRate", 
                            "ERCCsumLogErr", "sum", "Sex")]
## This sample has the min prop of reads assigned to genes 
pca_data_gene_brain_pups[which.min(pca_data_gene_brain_pups$totalAssignedGene), "SAMPLE_ID"]
# "Sample_P2_fe2_022019"



## Explore Sample_P1_fe3_021819 info:
colData(rse_gene_brain_pups_qc)[which(rse_gene_brain_pups_qc$SAMPLE_ID=="Sample_P1_fe3_021819"),
                          c("mitoRate", "rRNA_rate", "totalAssignedGene", "overallMapRate", 
                            "ERCCsumLogErr", "sum", "Sex")]
## This sample has the max Error value
pca_data_gene_brain_pups[which.max(abs(pca_data_gene_brain_pups$ERCCsumLogErr)), "SAMPLE_ID"]
# "Sample_P1_fe3_021819"



## Remove those samples
poorQC_samples<-c("Sample_FE3P2", "Sample_4067")
for (sample in poorQC_samples){
  rse_gene_brain_adults_qc<-rse_gene_brain_adults_qc[,-which(rse_gene_brain_adults_qc$SAMPLE_ID==sample)]
  rse_exon_brain_adults_qc<-rse_exon_brain_adults_qc[,-which(rse_exon_brain_adults_qc$SAMPLE_ID==sample)]
  rse_tx_brain_adults_qc<-rse_tx_brain_adults_qc[,-which(rse_tx_brain_adults_qc$SAMPLE_ID==sample)]
  rse_jx_brain_adults_qc<-rse_jx_brain_adults_qc[,-which(rse_jx_brain_adults_qc$SAMPLE_ID==sample)]
}

poorQC_samples<-c("Sample_P2_fe2_022019", "Sample_P1_fe3_021819")
for (sample in poorQC_samples){
  rse_gene_brain_pups_qc<-rse_gene_brain_pups_qc[,-which(rse_gene_brain_pups_qc$SAMPLE_ID==sample)]
  rse_exon_brain_pups_qc<-rse_exon_brain_pups_qc[,-which(rse_exon_brain_pups_qc$SAMPLE_ID==sample)]
  rse_tx_brain_pups_qc<-rse_tx_brain_pups_qc[,-which(rse_tx_brain_pups_qc$SAMPLE_ID==sample)]
  rse_jx_brain_pups_qc<-rse_jx_brain_pups_qc[,-which(rse_jx_brain_pups_qc$SAMPLE_ID==sample)]
}

## PCA plots without those samples
plot_PCAs("gene", "brain", "adults")
plot_PCAs("gene", "brain", "pups")
plot_PCAs("exon", "brain", "adults")
plot_PCAs("exon", "brain", "pups")
plot_PCAs("tx", "brain", "adults")
plot_PCAs("tx", "brain", "pups")
plot_PCAs("jx", "brain", "adults")
plot_PCAs("jx", "brain", "pups")




### 1.1.2 Explore samples' gene expression variation across phenotypes by Group
## PC boxplots 
PC_boxplots <- function (PCx, pheno_var1, pheno_var2, tissue, type, age) {
    pca_data<-PCA(tissue, type, age)[[1]]
    pca_vars<-PCA(tissue, type, age)[[2]]
    plot=ggplot(data=pca_data, 
         aes(x=eval(parse_expr(pheno_var1)),y=eval(parse_expr(PCx)))) + 
         ## Hide outliers
         geom_boxplot(outlier.color = "#FFFFFFFF") +
         ## PC dots colored by a group + noise
         geom_jitter(aes(colour=eval(parse_expr(pheno_var2))),shape=16, 
                     position=position_jitter(0.2)) +
         theme_classic() +
         theme(legend.position="right", plot.margin=unit (c (1,1.5,1,1), 'cm')) +
         labs(x= pheno_var1, y = pca_vars[strtoi(gsub("PC","", PCx))],  
              color=pheno_var2)
    return(plot)
}

## Boxplots of PCx data separated by Group and phenotypes
plot_PC_boxplots<-function(PCx, type, tissue, age) {
  i=1
  plots<-list()
  for (pheno_var1 in c("Age", "plate","Expt", "Sex", "Pregnancy", "flowcell")){
      p<-PC_boxplots(PCx, pheno_var1, "Group", tissue, type, age)
      plots[[i]]=p
      i=i+1
    }
  plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], 
            plots[[6]], nrow = 2)
  ## Save plot
  if (is.null(age)){
    ggsave(paste("plots/03_EDA/03_PCA_MDS/boxplot_",PCx,"_", type, "_", tissue, ".pdf", 
         sep=""), width = 50, height = 35, units = "cm")
  }
  else {
      ggsave(paste("plots/03_EDA/03_PCA_MDS/boxplot_",PCx,"_", type, "_", tissue, "_", 
                   age, ".pdf", sep=""), width = 50, height = 35, units = "cm")
  }

}

## Plots
plot_PC_boxplots("PC1", "gene", "blood", NULL)
plot_PC_boxplots("PC1", "gene", "brain", "adults")
plot_PC_boxplots("PC1", "gene", "brain", "pups")
plot_PC_boxplots("PC2", "gene", "brain", "pups")




