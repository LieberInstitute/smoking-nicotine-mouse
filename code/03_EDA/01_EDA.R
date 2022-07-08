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





# 3. Explore sample effects
# 3.1 Dimensionality Reduction (PCA)
### 3.1.1 Explore Samples' expression variation

## Generate PCA data
PCA<-function(tissue, type, age){
  if (is.null(age)) {
    RSE<-eval(parse_expr(paste("rse", type, tissue, "qc", sep="_"))) }
  else {
    RSE<-eval(parse_expr(paste("rse", type, tissue, age, "qc", sep="_")))}
  pca<-prcomp(t(assay(RSE)))
  # % of the variance explained by each PC
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
PCx_vs_PCy <- function (PCx, PCy, pheno_var, tissue, type, age) {
  pca_data<-PCA(tissue, type, age)[[1]]
  pca_vars<-PCA(tissue, type, age)[[2]]
  if (is.null(age)){
    title=paste("PCA of", type, "data, colored by", pheno_var, "in", tissue)
  }
  else {
    title=paste("PCA of", type, "data, colored by", pheno_var, "in", age, "and", tissue)
  }
  plot=ggplot(data=pca_data, 
       aes(x=eval(parse_expr(PCx)),y=eval(parse_expr(PCy)),
           color=eval(parse_expr(pheno_var))) )+ 
       theme(legend.position="right", plot.title = element_text (hjust = 0.5, size = 10), 
            plot.margin=unit (c (1,1.5,1,1), 'cm')) +
       geom_point() + 
       ggtitle(title) +
       labs(x= pca_vars[strtoi(gsub("PC","", PCx))], y = pca_vars[strtoi(gsub("PC","", PCy))],  
            color=pheno_var)
  return(plot)
}

## Plots 
plot_PCAs<-function(type, tissue, age){
  ## Plot for type and tissue
  if (is.null(age)){
      for (PCs in list(c("PC1", "PC2"), c("PC1", "PC3"), c("PC2", "PC3"))){
      plots<-list()
      i=1
        for (pheno_var in c("Age", "plate","Expt", "Sex", "Group", "Pregnancy", "medium", "flowcell")){
           p<-PCx_vs_PCy(PCs[1], PCs[2], pheno_var, tissue, type, age)
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
      plots<-list()
      i=1
      for (pheno_var in c("Age", "plate","Expt", "Sex", "Group", "Pregnancy", "medium", "flowcell")){
         p<-PCx_vs_PCy("PC1", "PC2", pheno_var, tissue, type, age)
         plots[[i]]=p
         i=i+1
      }
      plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]],
                nrow = 2)
      ## Save plots
      ggsave(paste("plots/03_EDA/03_PCA_MDS/PC1_vs_PC2_", type, "_", tissue ,"_", age, ".pdf", sep=""), 
             width = 65, height = 25, units = "cm")
    }
}

## Plot for type and tissue
plot_PCAs("gene", "brain", NULL)
plot_PCAs("exon", "brain", NULL)
plot_PCAs("tx", "brain", NULL)
plot_PCAs("jx", "brain", NULL)
plot_PCAs("gene", "blood", NULL)


## Rare samples in brain plots
pca_data_gene_brain[which.min(pca_data_gene_brain$PC2), "SAMPLE_ID"]
# "Sample_P12_fc2_022019"
pca_data_exon_brain[which.max(pca_data_exon_brain$PC2), "SAMPLE_ID"]
# "Sample_P2_fe2_022019"
pca_data_tx_brain[which.min(pca_data_tx_brain$PC2), "SAMPLE_ID"]
# "Sample_P2_fe2_022019"
pca_data_jx_brain[which.max(pca_data_jx_brain$PC1), "SAMPLE_ID"]
# "Sample_P2_fe2_022019"


## Explore Sample_P12_fc2_022019 info:
colData(rse_gene_brain_qc)[which(rse_gene_brain_qc$SAMPLE_ID=="Sample_P12_fc2_022019"),
                          c("mitoRate", "rRNA_rate", "totalAssignedGene", "overallMapRate", 
                            "ERCCsumLogErr", "sum", "Age", "Sex")]
#    mitoRate rRNA_rate totalAssignedGene overallMapRate ERCCsumLogErr       sum         Age         Sex
#   <numeric> <numeric>         <numeric>      <numeric>     <numeric> <numeric> <character> <character>
#   0.0349937 0.0342992          0.764755         0.9609      -68.3041    869084         Pup           M

## This sample has the max rRNA rate of all brain samples and min overall mapping rate 
colData(rse_gene_brain)[which.max(rse_gene_brain$rRNA_rate),"SAMPLE_ID"]
# "Sample_P12_fc2_022019"
colData(rse_gene_brain_qc)[which.min(rse_gene_brain_qc$overallMapRate),"SAMPLE_ID"]
# "Sample_P12_fc2_022019"


## Explore Sample_P2_fe2_022019 info:
colData(rse_gene_brain_qc)[which(rse_gene_brain_qc$SAMPLE_ID=="Sample_P2_fe2_022019"),
                          c("mitoRate", "rRNA_rate", "totalAssignedGene", "overallMapRate", 
                            "ERCCsumLogErr", "sum", "Age", "Sex")]
#    mitoRate  rRNA_rate totalAssignedGene overallMapRate ERCCsumLogErr       sum         Age         Sex
#   <numeric>  <numeric>         <numeric>      <numeric>     <numeric> <numeric> <character> <character>
#   0.0304307 0.00612927          0.734151         0.9787      -50.9228   1002923         Pup           M

## This sample has the max number of total read counts in pups
colData(rse_gene_brain_pups_qc)[which.max(rse_gene_brain_pups_qc$sum),"SAMPLE_ID"]
# "Sample_P2_fe2_022019"



## Remove those pup samples
for (sample in c("Sample_P12_fc2_022019","Sample_P2_fe2_022019" )){
  
  rse_gene_brain_qc<-rse_gene_brain_qc[,-which(rse_gene_brain_qc$SAMPLE_ID==sample)]
  rse_gene_brain_pups_qc<-rse_gene_brain_pups_qc[,-which(rse_gene_brain_pups_qc$SAMPLE_ID==sample)]
  rse_exon_brain_qc<-rse_exon_brain_qc[,-which(rse_exon_brain_qc$SAMPLE_ID==sample)]
  rse_exon_brain_pups_qc<-rse_exon_brain_pups_qc[,-which(rse_exon_brain_pups_qc$SAMPLE_ID==sample)]
  rse_tx_brain_qc<-rse_tx_brain_qc[,-which(rse_tx_brain_qc$SAMPLE_ID==sample)]
  rse_tx_brain_pups_qc<-rse_tx_brain_pups_qc[,-which(rse_tx_brain_pups_qc$SAMPLE_ID==sample)]
  rse_jx_brain_qc<-rse_jx_brain_qc[,-which(rse_jx_brain_qc$SAMPLE_ID==sample)]
  rse_jx_brain_pups_qc<-rse_jx_brain_pups_qc[,-which(rse_jx_brain_pups_qc$SAMPLE_ID==sample)]
}

## Plot samples of brain again without those samples
plot_PCAs("gene", "brain", NULL)
plot_PCAs("exon", "brain", NULL)
plot_PCAs("tx", "brain", NULL)
plot_PCAs("jx", "brain", NULL)

## Plot samples separated by age (without those samples)
plot_PCAs("gene", "brain", "adults")
plot_PCAs("gene", "brain", "pups")
plot_PCAs("exon", "brain", "adults")
plot_PCAs("exon", "brain", "pups")
plot_PCAs("tx", "brain", "adults")
plot_PCAs("tx", "brain", "pups")
plot_PCAs("jx", "brain", "adults")
plot_PCAs("jx", "brain", "pups")



### 3.1.2 Explore samples' gene expression variation across phenotypes by Group




