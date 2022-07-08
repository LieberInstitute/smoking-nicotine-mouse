# Exploratory Data Analysis

## Load libraries and data 
library(SummarizedExperiment)
library(recount)
library(edgeR)
library(jaffelab)
library(here)
library(ggplot2)
library(rlang)
library(scater)
library(cowplot)
library(variancePartition)
library(gridExtra)
library(sessioninfo)

load(here("processed-data/02_build_objects/rse_gene_brain.Rdata"))
load(here("processed-data/02_build_objects/rse_gene_blood.Rdata"))
load(here("processed-data/02_build_objects/rse_exon_brain.Rdata"))
load(here("processed-data/02_build_objects/rse_tx_brain.Rdata"))
load(here("processed-data/02_build_objects/rse_jx_brain.Rdata"))

## 1. Study design 

## Number of samples of each pair of phenotypes
samples_pheno <- function(pheno_var1, pheno_var2){
    ## Pheno variables levels
    levels_var1<-sort(unique(pheno[,pheno_var1]))
    levels_var2<-sort(unique(pheno[,pheno_var2]))
    for (level1 in levels_var1){
        for (level2 in levels_var2){
            print(paste("Number of samples of (", pheno_var1, ") ", level1, " and (", pheno_var2, ") ", level2, ": ", 
                        length(which(pheno[,pheno_var1]==level1 & pheno[,pheno_var2]==level2)), sep=""))
        }
    }
}

## All pairs of phenotypes
vars1<-vector()
`%nin%` = Negate(`%in%`)
for (pheno_var1 in c("Age", "Tissue", "plate", "medium","Expt", "Sex", "Group", "Pregnancy")){
    for (pheno_var2 in c("Age", "Tissue", "plate", "medium","Expt", "Sex", "Group", "Pregnancy")){
        if (pheno_var1!=pheno_var2 & pheno_var2 %nin% vars1){
            vars1<-append(vars1, pheno_var1)
            samples_pheno(pheno_var1, pheno_var2)}
    }
} 

three_pheno <- function(pheno_var1, pheno_var2, pheno_var3){
  ## Pheno variables levels
    levels_var1<-sort(unique(pheno[,pheno_var1]))
    levels_var2<-sort(unique(pheno[,pheno_var2]))
    levels_var3<-sort(unique(pheno[,pheno_var3]))
    for (level1 in levels_var1){
      for (level2 in levels_var2){
        for (level3 in levels_var3){
          print(paste("Number of samples of (", pheno_var1, ") ", level1, ", (", pheno_var2, ") ", level2, 
                " and (", pheno_var3, ") ", level3, ": ", length(which(pheno[,pheno_var1]==level1 &
                pheno[,pheno_var2]==level2 & pheno[,pheno_var3]==level3)), sep=""))
        }
      }
    }
}

## Number of samples 
three_pheno("Pregnancy", "Age", "Tissue")
three_pheno("Pregnancy", "Age", "Group")
three_pheno("Tissue", "Expt", "Group")
three_pheno("Pregnancy", "Expt", "Group")





# 2. QC analysis
## 2.1 Relationships between QC variables and mouse phenotypes

## Boxplots 
create_boxplots <- function(pheno_var, qc_var, tissue) {
    ## Tissue data
    RSE<-eval(parse_expr(paste("rse_gene_", tissue, sep="")))
    ## Quantitative QC values grouped by a qualitative phenotype variables
    plot=ggplot(as.data.frame(colData(RSE)), 
        aes(x=eval(parse_expr(pheno_var)), y=eval(parse_expr(qc_var)), fill=eval(parse_expr(pheno_var)))) +
        geom_boxplot() +
        theme_classic(base_size = 6) +
        theme(legend.position="none", plot.title = element_text (hjust = 0.5), 
              plot.margin=unit (c (1.5,2,1,2), 'cm')) +
       labs(title = paste("Boxplots of", qc_var, "separated by", pheno_var, "in", tissue), x=pheno_var, y=qc_var, ) 
  return(plot)
}


plot_boxplots <- function(tissue){

  for (qc_var in c("rRNA_rate","mitoRate","totalAssignedGene", "ERCCsumLogErr", "overallMapRate")){
    plots<-list()
    i=1
      for (pheno_var in c("Age", "plate","Expt", "Sex", "Group", "medium", "Pregnancy", "flowcell")){
         p<-create_boxplots(pheno_var, qc_var, tissue)
         plots[[i]]=p
         i=i+1
      }
    plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], 
              plots[[7]], plots[[8]], nrow = 2)
    ## Save plots
    ggsave(paste("plots/03_EDA/02_QC/boxplot_",qc_var,"_", tissue, ".pdf", sep=""), width = 35, height = 25, units = "cm")
  }
}

## Plots
plot_boxplots("brain")
plot_boxplots("blood")

## Trimmed read libraries
table(rse_gene_brain$trimmed)
# FALSE 
#   184 
table(rse_gene_blood$trimmed)
# FALSE 
#   24 



## 2.2 Calculate quality metrics

## Brain samples QC metrics
## Mt and ribosomal genes 
subsets=list(Mito=which(seqnames(rse_gene_brain)=="chrM"), 
             Ribo=grep("rRNA",rowData(rse_gene_brain)$gene_type))
## Add general qc-stats
rse_gene_brain <-addPerCellQC(rse_gene_brain, subsets)


## Blood samples QC metrics
## Mt and ribosomal genes 
subsets=list(Mito=which(seqnames(rse_gene_blood)=="chrM"), 
             Ribo=grep("rRNA",rowData(rse_gene_blood)$gene_type))
## Add general qc-stats
rse_gene_blood <-addPerCellQC(rse_gene_blood, subsets)


## Detected/mt/ribo genes VS total counts per sample 
## Samples separated by phenotypes
sum_vs_qc<- function (pheno_var, qc_stats, qc_stats_title, qc_stats_lab, tissue){
    RSE<-eval(parse_expr(paste("rse_gene_", tissue, sep="")))
    plot=ggplot(data=as.data.frame(colData(RSE)), 
           aes(x=sum, y=eval(parse_expr(qc_stats)), color=eval(parse_expr(pheno_var))))+ 
        geom_point() +
        theme(plot.title = element_text(size = 13)) +
        theme(legend.position="right", plot.title = element_text (hjust = 0.5), 
                  plot.margin=unit (c (1.5,2,1,2), 'cm')) +
        labs(title=paste(qc_stats_title, "by" ,pheno_var, "in", tissue), 
            x="Total read counts", y=qc_stats_lab, color=pheno_var)
    return(plot)
}

plot_sum_vs_qc<- function(tissue){
  for (qc_stats in c("detected","subsets_Mito_percent", "subsets_Ribo_percent")){
    if (qc_stats=="detected")
       {qc_stats_lab="Detected number of genes"
       qc_stats_title="Detected genes vs total read counts per sample,"
       qc_var="Detected"}
    if (qc_stats=="subsets_Mito_percent")
       {qc_stats_lab="% of mt genes' counts of the total counts"
        qc_stats_title="Total read counts vs mitochondrial genes' counts per sample,"
        qc_var="mtCounts"}
    if (qc_stats=="subsets_Ribo_percent")
       {qc_stats_lab="% of ribosomal genes' counts of the total counts"
        qc_stats_title="Total read counts vs ribosomal genes' counts per sample,"
        qc_var="riboCounts"}
    plots<-list()
    i=1
      for (pheno_var in c("Age", "plate","Expt", "Sex", "Group", "medium", "Pregnancy", "flowcell")){
         p<-sum_vs_qc(pheno_var, qc_stats, qc_stats_title, qc_stats_lab, tissue)
         plots[[i]]=p
         i=i+1
      }
  
    plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], 
              plots[[7]], plots[[8]], nrow = 2)
    ## Save plots
    ggsave(paste("plots/03_EDA/02_QC/totalCounts_vs_",qc_var,"_",tissue,".pdf", sep=""), 
           width = 85, height = 35, units = "cm")
  }
}

## Plots
plot_sum_vs_qc("brain")
plot_sum_vs_qc("blood")


## Filter samples by age
## Separate samples by Age in brain (all blood samples are from adult)
rse_gene_brain_adults<-rse_gene_brain[,which(rse_gene_brain$Age=="Adult")]
## 47 samples
rse_gene_brain_pups<-rse_gene_brain[,which(rse_gene_brain$Age=="Pup")]
## 137 samples




## 2.3 Samples filtering by QC 
### 2.3.1 Detection-based and Mito/Ribo filtering

## Filter brain samples
rse_gene_brain_qc<-quickPerCellQC(
  rse_gene_brain,
  ## Filter by detected number of genes and library size 
  sum.field = "sum",
  detected.field = "detected",
  ## Filter by mito/ribo content
  sub.fields = TRUE,
)
## Save data
save(rse_gene_brain_qc, file = 'processed-data/03_EDA/02_QC/rse_gene_brain_qc.Rdata')

## Filter by QC and save data of exons, tx and jx 
rse_exon_brain_qc<-rse_exon_brain[,(rse_exon_brain$SAMPLE_ID %in% rse_gene_brain_qc$SAMPLE_ID)]
save(rse_exon_brain_qc, file = 'processed-data/03_EDA/02_QC/rse_exon_brain_qc.Rdata')
rse_tx_brain_qc<-rse_tx_brain[,(rse_tx_brain$SAMPLE_ID %in% rse_gene_brain_qc$SAMPLE_ID)]
save(rse_tx_brain_qc, file = 'processed-data/03_EDA/02_QC/rse_tx_brain_qc.Rdata')
rse_jx_brain_qc<-rse_jx_brain[,(rse_jx_brain$SAMPLE_ID %in% rse_gene_brain_qc$SAMPLE_ID)]
save(rse_jx_brain_qc, file = 'processed-data/03_EDA/02_QC/rse_jx_brain_qc.Rdata')

## Number of samples retained
dim(rse_gene_brain_qc)[2]
# 161
## Number of genes per sample
table(rse_gene_brain_qc$detected)
# 8876 
#  161 


## Filter adult brain samples
rse_gene_brain_adults_qc<-quickPerCellQC(
  rse_gene_brain_adults,
  sum.field = "sum",
  detected.field = "detected",
  sub.fields = TRUE,
)
## Save data
save(rse_gene_brain_adults_qc, file = 'processed-data/03_EDA/02_QC/rse_gene_brain_adults_qc.Rdata')

## Filter by QC and save data of exons, tx and jx 
rse_exon_brain_adults_qc<-rse_exon_brain[,(rse_exon_brain$SAMPLE_ID %in% rse_gene_brain_adults_qc$SAMPLE_ID)]
save(rse_exon_brain_adults_qc, file = 'processed-data/03_EDA/02_QC/rse_exon_brain_adults_qc.Rdata')
rse_tx_brain_adults_qc<-rse_tx_brain[,(rse_tx_brain$SAMPLE_ID %in% rse_gene_brain_adults_qc$SAMPLE_ID)]
save(rse_tx_brain_adults_qc, file = 'processed-data/03_EDA/02_QC/rse_tx_brain_adults_qc.Rdata')
rse_jx_brain_adults_qc<-rse_jx_brain[,(rse_jx_brain$SAMPLE_ID %in% rse_gene_brain_adults_qc$SAMPLE_ID)]
save(rse_jx_brain_adults_qc, file = 'processed-data/03_EDA/02_QC/rse_jx_brain_adults_qc.Rdata')

## Number of samples retained
dim(rse_gene_brain_adults_qc)[2]
# 24
## Number of genes per sample
table(rse_gene_brain_adults_qc$detected)
# 8876 
#   24 


# Filter pup brain samples
rse_gene_brain_pups_qc<-quickPerCellQC(
  rse_gene_brain_pups,
  sum.field = "sum",
  detected.field = "detected",
  sub.fields = TRUE,
)
## Save data
save(rse_gene_brain_pups_qc, file = 'processed-data/03_EDA/02_QC/rse_gene_brain_pups_qc.Rdata')

## Filter by QC and save data of exons, tx and jx 
rse_exon_brain_pups_qc<-rse_exon_brain[,(rse_exon_brain$SAMPLE_ID %in% rse_gene_brain_pups_qc$SAMPLE_ID)]
save(rse_exon_brain_pups_qc, file = 'processed-data/03_EDA/02_QC/rse_exon_brain_pups_qc.Rdata')
rse_tx_brain_pups_qc<-rse_tx_brain[,(rse_tx_brain$SAMPLE_ID %in% rse_gene_brain_pups_qc$SAMPLE_ID)]
save(rse_tx_brain_pups_qc, file = 'processed-data/03_EDA/02_QC/rse_tx_brain_pups_qc.Rdata')
rse_jx_brain_pups_qc<-rse_jx_brain[,(rse_jx_brain$SAMPLE_ID %in% rse_gene_brain_pups_qc$SAMPLE_ID)]
save(rse_jx_brain_pups_qc, file = 'processed-data/03_EDA/02_QC/rse_jx_brain_pups_qc.Rdata')

## Number of samples retained
dim(rse_gene_brain_pups_qc)[2]
# 137
## Number of genes per sample
table(rse_gene_brain_pups_qc$detected)
# 8876 
#  137 

## Compare all brain samples retained with those by Age
all_byAge<-union(rse_gene_brain_adults_qc$SAMPLE_ID, rse_gene_brain_pups_qc$SAMPLE_ID)
setdiff(rse_gene_brain_qc$SAMPLE_ID, all_byAge)
# character(0)




## Filter blood samples 
rse_gene_blood_qc<-quickPerCellQC(
  rse_gene_blood,
  ## Filter by detected number of genes and library size 
  sum.field = "sum",
  detected.field = "detected",
  ## Filter by mito/ribo content
  sub.fields = TRUE,
)
## Save data
save(rse_gene_blood_qc, file = 'processed-data/03_EDA/02_QC/rse_gene_blood_qc.Rdata')

## Number of samples retained
dim(rse_gene_blood_qc)[2]
# 24

## Number of genes per sample
table(rse_gene_blood_qc$detected)
# 8035 8048 8066 8070 8078 8080 8086 8102 8113 8114 8120 8145 8164 8165 8214 8216 8223 8243 
#    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    2 
# 8262 8286 8437 8448 8702 
#    1    1    1    1    1 


## Plot QC metrics for samples retained and samples dropped
## Generate column with samples retained/dropped
qc_filt_col <- function(tissue){
  RSE<-eval(parse_expr(paste("rse_gene_", tissue, sep="")))
  RSE_qc<-eval(parse_expr(paste("rse_gene", tissue, "qc", sep="_")))
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

## Add column for brain (all blood samples were retained)
rse_gene_brain$Retention<-qc_filt_col("brain")


## Plot samples separated by two phenotypes
sum_vs_qc_byTwo<- function (pheno_var1, pheno_var2, qc_stats, qc_stats_title, qc_stats_lab, tissue){
    RSE<-eval(parse_expr(paste("rse_gene_", tissue, sep="")))
    plot=ggplot(data=as.data.frame(colData(RSE)), 
           aes(x=sum, y=eval(parse_expr(qc_stats)), color=eval(parse_expr(pheno_var1)), 
           shape=eval(parse_expr(pheno_var2))))+ 
        geom_point() +
        theme(plot.title = element_text(size = 12)) +
        theme(legend.position="right", plot.title = element_text (hjust = 0.5), 
                  plot.margin=unit (c (1.5,2,1,2), 'cm')) +
        labs(title=paste(qc_stats_title, "by" ,pheno_var1, "and", pheno_var2, "in", tissue), 
            x="Total read counts", y=qc_stats_lab, color=pheno_var1, shape=pheno_var2)
    return(plot)
}

## Generate plots
qc_stats="detected"
qc_stats_lab="Detected number of genes"
qc_stats_title="Detected genes vs total read counts per sample,"
p1<-sum_vs_qc_byTwo("Retention", "Age", qc_stats, qc_stats_title, qc_stats_lab, "brain")

qc_stats="subsets_Mito_percent"
qc_stats_lab="% of mt genes' counts of the total counts"
qc_stats_title="Total read counts vs mitochondrial genes' counts per sample,"
p2<-sum_vs_qc_byTwo("Retention", "Age", qc_stats, qc_stats_title, qc_stats_lab, "brain")

qc_stats="subsets_Ribo_percent"
qc_stats_lab="% of ribosomal genes' counts of the total counts"
qc_stats_title="Total read counts vs ribosomal genes' counts per sample,"
p3<-sum_vs_qc_byTwo("Retention", "Age", qc_stats, qc_stats_title, qc_stats_lab, "brain")

plot_grid(p1, p2, p3, nrow = 1)
## Save plots
ggsave(paste("plots/03_EDA/02_QC/totalCounts_vs_QCvars_Retained_",tissue,".pdf", sep=""), 
         width = 55, height = 15, units = "cm")



## Study design after QC filtering 
samples_pheno_afterQC <- function(pheno_var1, pheno_var2){
  RSE<-rbind(colData(rse_gene_brain_qc), colData(rse_gene_blood_qc))
  ## Pheno variables levels
  levels_var1<-sort(unique(RSE[,pheno_var1]))
  levels_var2<-sort(unique(RSE[,pheno_var2]))
  for (level1 in levels_var1){
      for (level2 in levels_var2){
          print(paste("Number of samples of (", pheno_var1, ") ", level1, " and (", 
                      pheno_var2, ") ", level2, ": ", 
                      length(which(RSE[,pheno_var1]==level1 & 
                                     RSE[,pheno_var2]==level2)), sep=""))
      }
  }
}
## Obtain number of samples for each pair of phenotypes
vars1<-vector()
`%nin%` = Negate(`%in%`)
for (pheno_var1 in c("Age", "Tissue", "plate","Expt", "Sex", "Group", "Pregnancy")){
  for (pheno_var2 in c("Age", "Tissue", "plate","Expt", "Sex", "Group", "Pregnancy")){
      if (pheno_var1!=pheno_var2 & pheno_var2 %nin% vars1){
          vars1<-append(vars1, pheno_var1)
          samples_pheno_afterQC(pheno_var1, pheno_var2)}
      }
} 
