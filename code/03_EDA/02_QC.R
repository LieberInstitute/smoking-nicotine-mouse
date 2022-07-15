## Load all libraries for EDA 
library(SummarizedExperiment)
library(recount)
library(edgeR)
library(here)
library(ggplot2)
library(rlang)
library(scater)
library(jaffelab)
library(cowplot)
library(variancePartition)
library(gridExtra)
library(sessioninfo)


## 1. Quality Control anlaysis

## 1.1 Load data
load(here("processed-data/02_build_objects/rse_gene_brain.Rdata"))
load(here("processed-data/02_build_objects/rse_gene_blood.Rdata"))
load(here("processed-data/02_build_objects/rse_gene_brain_adults.Rdata"))
load(here("processed-data/02_build_objects/rse_gene_brain_pups.Rdata"))



## 1.2 Relationships between QC variables and mouse phenotypes
## Boxplots 
create_boxplots <- function(pheno_var, qc_var, tissue, age) {
  if (is.null(age)){
    ## Tissue data
    RSE<-eval(parse_expr(paste("rse_gene_", tissue, sep="")))
    }
  else {
    ## Tissue and Age data
    RSE<-eval(parse_expr(paste("rse_gene", tissue, age, sep="_")))
    }

  ## Quantitative QC values grouped by a qualitative phenotype variable
  plot=ggplot(as.data.frame(colData(RSE)), 
      aes(x=eval(parse_expr(pheno_var)), y=eval(parse_expr(qc_var)), 
          fill=eval(parse_expr(pheno_var)))) +
      geom_boxplot() +
      theme_classic(base_size = 10) +
      theme(legend.position="none", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
            axis.text.x = element_text(vjust = 0.45) ) +
     labs(x=pheno_var, y=qc_var) 
return(plot)
}

## Plot all boxplots for a single qc variable
plot_boxplots <- function(tissue, age){
  for (qc_var in c("rRNA_rate","mitoRate","totalAssignedGene", "ERCCsumLogErr", "overallMapRate")){
    plots<-list()
    i=1
      for (pheno_var in c("Age", "plate","Expt", "Sex", "Group", "medium", "Pregnancy", "flowcell")){
         p<-create_boxplots(pheno_var, qc_var, tissue, age)
         plots[[i]]=p
         i=i+1
      }
    plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], 
              plots[[7]], plots[[8]], nrow = 2)
    ## Save plots
    if (is.null(age)) {fileName=paste("plots/03_EDA/02_QC/boxplot_",qc_var,"_", 
                                        tissue, ".pdf", sep="")}
    else {fileName=paste("plots/03_EDA/02_QC/boxplot_",qc_var,"_", tissue, "_", 
                       age, ".pdf", sep="")}
    ggsave(fileName, width = 55, height = 25, units = "cm")
  }
}

## Plots
plot_boxplots("brain", NULL)
plot_boxplots("blood", NULL)
plot_boxplots("brain", "adults")
plot_boxplots("brain", "pups")

## Trimmed read libraries
table(rse_gene_brain$trimmed)
# FALSE 
#   184 
table(rse_gene_blood$trimmed)
# FALSE 
#   24 





## 1.3 Plot quality metrics
## Detected/mt/ribo genes VS total counts per sample 
## Samples separated by phenotypes
sum_vs_qc<- function (pheno_var, qc_stats, qc_stats_lab, tissue, age){
    if (is.null(age)){
      RSE<-eval(parse_expr(paste("rse_gene_", tissue, sep="")))
      }
    else {
      RSE<-eval(parse_expr(paste("rse_gene", tissue, age, sep="_")))
      }
    plot=ggplot(data=as.data.frame(colData(RSE)), 
           aes(x=sum, y=eval(parse_expr(qc_stats)), color=eval(parse_expr(pheno_var))))+ 
        geom_point() +
        theme(text = element_text(size = 10)) +
        theme(legend.position="right", plot.margin=unit (c (1.5,2,1,2), 'cm')) +
        labs(x="Total read counts", y=qc_stats_lab, color=pheno_var)
    return(plot)
}

## All plots for a qc variable
plot_sum_vs_qc<- function(tissue, age){
  for (qc_stats in c("detected","subsets_Mito_percent", "subsets_Ribo_percent")){
    if (qc_stats=="detected")
       {qc_stats_lab="Detected number of genes"
       qc_var="Detected"}
    if (qc_stats=="subsets_Mito_percent")
       {qc_stats_lab="Percentage of mt counts "
        qc_var="mtCounts"}
    if (qc_stats=="subsets_Ribo_percent")
       {qc_stats_lab="Percentage of ribosomal counts "
        qc_var="riboCounts"}
    plots<-list()
    i=1
      for (pheno_var in c("Age", "plate","Expt", "Sex", "Group", "medium", "Pregnancy", "flowcell")){
         p<-sum_vs_qc(pheno_var, qc_stats, qc_stats_lab, tissue, age)
         plots[[i]]=p
         i=i+1
      }
  
    plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], 
              plots[[7]], plots[[8]], nrow = 2)
    ## Save plots
    if (is.null(age)){fileName=paste("plots/03_EDA/02_QC/totalCounts_vs_",qc_var,"_",tissue,".pdf", sep="")}
    else {fileName=paste("plots/03_EDA/02_QC/totalCounts_vs_",qc_var,"_",tissue,"_", age,".pdf", sep="")}
    ggsave(fileName, width = 65, height = 25, units = "cm")
  }
}

## Plots
plot_sum_vs_qc("brain", NULL)
plot_sum_vs_qc("blood", NULL)
plot_sum_vs_qc("brain", "adults")
plot_sum_vs_qc("brain", "pups")

## Brain
## Correlation between number of genes and total read counts
cor(rse_gene_brain$sum, rse_gene_brain$detected)
# 0.7358813

## Blood
## Correlation between number of genes and total read counts
cor(rse_gene_blood$sum, rse_gene_blood$detected)
# -0.1830506

## Adult brain 
## Correlation between number of genes and total read counts
cor(rse_gene_brain_adults$sum, rse_gene_brain_adults$detected)
# 0.6474778

## Pup brain
## Correlation between number of genes and total read counts
cor(rse_gene_brain_pups$sum, rse_gene_brain_pups$detected)
# 0.7703833





## 1.4 Sample filtering by QC 
### 1.4.1 Detection-based and Mito/Ribo filtering 

## Filter brain samples
rse_gene_brain_qc<-quickPerCellQC(
  rse_gene_brain,
  ## Filter by detected number of genes and library size 
  sum.field = "sum",
  detected.field = "detected",
  ## Filter by mito/ribo content
  sub.fields = TRUE,
)
## Number of samples retained
dim(rse_gene_brain_qc)[2]
# 184


## Filter adult brain samples
rse_gene_brain_adults_qc<-quickPerCellQC(
  rse_gene_brain_adults,
  sum.field = "sum",
  detected.field = "detected",
  sub.fields = TRUE,
)
## Number of samples retained
dim(rse_gene_brain_adults_qc)[2]
# 47


# Filter pup brain samples
rse_gene_brain_pups_qc<-quickPerCellQC(
  rse_gene_brain_pups,
  sum.field = "sum",
  detected.field = "detected",
  sub.fields = TRUE,
)
## Number of samples retained
dim(rse_gene_brain_pups_qc)[2]
# 136
## Save data
save(rse_gene_brain_pups_qc, file = 'processed-data/03_EDA/02_QC/rse_gene_brain_pups_qc.Rdata')

## Filter data of exons, tx and jx by QC
rse_exon_brain_pups_qc<-rse_exon_brain_pups[,(rse_exon_brain_pups$SAMPLE_ID %in% rse_gene_brain_pups_qc$SAMPLE_ID)]
save(rse_exon_brain_pups_qc, file = 'processed-data/03_EDA/02_QC/rse_exon_brain_pups_qc.Rdata')
rse_tx_brain_pups_qc<-rse_tx_brain_pups[,(rse_tx_brain_pups$SAMPLE_ID %in% rse_gene_brain_pups_qc$SAMPLE_ID)]
save(rse_tx_brain_pups_qc, file = 'processed-data/03_EDA/02_QC/rse_tx_brain_pups_qc.Rdata')
rse_jx_brain_pups_qc<-rse_jx_brain_pups[,(rse_jx_brain_pups$SAMPLE_ID %in% rse_gene_brain_pups_qc$SAMPLE_ID)]
save(rse_jx_brain_pups_qc, file = 'processed-data/03_EDA/02_QC/rse_jx_brain_pups_qc.Rdata')


## Filter blood samples 
rse_gene_blood_qc<-quickPerCellQC(
  rse_gene_blood,
  ## Filter by detected number of genes and library size 
  sum.field = "sum",
  detected.field = "detected",
  ## Filter by mito/ribo content
  sub.fields = TRUE,
)
## Number of samples retained
dim(rse_gene_blood_qc)[2]
# 24


### 1.4.2 Plot QC metrics for samples retained and samples dropped
## Generate column with samples retained/dropped
qc_filt_col <- function(tissue, age){
  if (is.null(age)){
    RSE<-eval(parse_expr(paste("rse_gene_", tissue, sep="")))
    RSE_qc<-eval(parse_expr(paste("rse_gene", tissue, "qc", sep="_")))
  }
  else {
    RSE<-eval(parse_expr(paste("rse_gene", tissue, age, sep="_")))
    RSE_qc<-eval(parse_expr(paste("rse_gene", tissue, age, "qc", sep="_")))
  }
  qc_filt<-vector()
  for (sample in RSE$SAMPLE_ID){
    if(sample %in% RSE_qc$SAMPLE_ID){
      ## Samples retained
      qc_filt<-append(qc_filt, "Yes")
    }
    ## Samples dropped
    else {qc_filt<-append(qc_filt, "No")}
  }
  return(qc_filt)
}

## Add column for pup brain samples
rse_gene_brain_pups$Retention<-qc_filt_col("brain", "pups")

## Generate plots
qc_stats="detected"
qc_stats_lab="Detected number of genes"
p1<-sum_vs_qc("Retention", qc_stats, qc_stats_lab, "brain", "pups")

qc_stats="subsets_Mito_percent"
qc_stats_lab="% of mt genes' counts of the total counts"
p2<-sum_vs_qc("Retention", qc_stats, qc_stats_lab, "brain", "pups")

qc_stats="subsets_Ribo_percent"
qc_stats_lab="% of ribosomal genes' counts of the total counts"
p3<-sum_vs_qc("Retention", qc_stats, qc_stats_lab, "brain", "pups")

plot_grid(p1, p2, p3, nrow = 1)
## Save plots
fileName=paste("plots/03_EDA/02_QC/samples_Retained_brain_pups.pdf", sep="")
ggsave(fileName, width = 55, height = 15, units = "cm")


## Removed sample information:
colData(rse_gene_brain_pups)[which(rse_gene_brain_pups$Retention=="No"), 
        c("SAMPLE_ID", "Sex", "Expt", "Group", "Pregnancy", "mitoRate", "rRNA_rate", "totalAssignedGene", "overallMapRate", "subsets_Mito_percent", "subsets_Ribo_percent")]
