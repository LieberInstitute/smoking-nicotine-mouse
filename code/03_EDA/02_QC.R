
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
library(ggrepel)
library(variancePartition)
library(gridExtra)
library(Hmisc)
library(stringr)
library(sessioninfo)


## 1. Quality Control analysis


load(here("processed-data/02_build_objects/rse_gene_brain.Rdata"))
load(here("processed-data/02_build_objects/rse_gene_blood.Rdata"))
load(here("processed-data/02_build_objects/rse_gene_brain_adults.Rdata"))
load(here("processed-data/02_build_objects/rse_exon_brain_adults.Rdata"))
load(here("processed-data/02_build_objects/rse_tx_brain_adults.Rdata"))
load(here("processed-data/02_build_objects/rse_jx_brain_adults.Rdata"))
load(here("processed-data/02_build_objects/rse_gene_brain_pups.Rdata"))
load(here("processed-data/02_build_objects/rse_exon_brain_pups.Rdata"))
load(here("processed-data/02_build_objects/rse_tx_brain_pups.Rdata"))
load(here("processed-data/02_build_objects/rse_jx_brain_pups.Rdata"))


# ------------------------------------------------------------------------------
##  1.1 Evaluate QC metrics for groups of samples
# ------------------------------------------------------------------------------

## 1.1.1 Compare QC metrics in sample groups

## Define QC metrics of interest
qc_metrics <- c("mitoRate", "overallMapRate", "totalAssignedGene", "rRNA_rate", "sum", "detected", "ERCCsumLogErr")

## Define sample variables of interest
sample_variables <- c("Group", "Expt", "Age", "Sex", "Pregnancy", "plate", "flowcell")

## Function to create boxplots of QC metrics for groups of samples
QC_boxplots <- function(qc_metric, sample_var, tissue, age) {
  
  if (is.null(age)){
    ## Tissue data
    RSE<-eval(parse_expr(paste("rse_gene_", tissue, sep="")))
  }
  else {
    ## Tissue and Age data
    RSE<-eval(parse_expr(paste("rse_gene", tissue, age, sep="_")))
  }
  
  ## Define sample colors depending on the sample variable
  if (sample_var == "Group") {
    colors <- c("Control" = "seashell3", "Experimental" = "orange3")
  } else if(sample_var=='Expt'){
    colors <- c("Nicotine" = "lightblue3", "Smoking" = "salmon")
  } else if (sample_var == "Age") {
    colors <- c("Adult" = "slateblue3", "Pup" = "yellow3")
  } else if (sample_var == "Sex") {
    colors <- c("F" = "hotpink1", "M" = "dodgerblue")
  } else if (sample_var == "Pregnancy") {
    colors <- c("Yes" = "darkorchid3", "No" = "darkolivegreen4")
  } else if (sample_var == "plate") {
    colors <- c("Plate1" = "darkorange", "Plate2" = "lightskyblue", "Plate3" = "deeppink1")
  } else if (sample_var == "flowcell") {
    colors <- c(
      "HKCG7DSXX" = "chartreuse2", "HKCMHDSXX" = "magenta", "HKCNKDSXX" = "turquoise3",
      "HKCTMDSXX" = "tomato", "HK7JHDSXX"="seagreen3", "HKCJCDSXX"="palevioletred2"
    )
  }
  
  ## Axis labels
  x_label <- capitalize(sample_var)
  y_label <- str_replace_all(qc_metric, c("_" = ""))
  
  ## x-axis text angle and position
  if (sample_var == "flowcell") {
    x_axis_angle <- 18
    x_axis_hjust <- 0.5
    x_axis_vjust <- 0.7
    x_axis_size <- 8
  } else {
    x_axis_angle <- 0
    x_axis_hjust <- 0.5
    x_axis_vjust <- 0.5
    x_axis_size <- 10
  }
  
  data <- data.frame(colData(RSE))
  
  ## Sample variable separating samples in x-axis and QC metric in y-axis
  ## (Coloring by sample variable)
  plot <- ggplot(data = data, mapping = aes(x = !!rlang::sym(sample_var), y = !!rlang::sym(qc_metric), color = !!rlang::sym(sample_var))) +
    ## Add violin plots
    geom_violin(alpha = 0, size = 0.4, color = "black", width = 0.7) +
    ## Spread dots
    geom_jitter(width = 0.05, alpha = 0.7, size = 1.3) +
    ## Add boxplots
    geom_boxplot(alpha = 0, size = 0.4, width = 0.1, color = "black") +
    ## Set colors
    scale_color_manual(values = colors) +
    ## Define axis labels
    labs(y = y_label, x = x_label) +
    ## Get rid of the background
    theme_bw() +
    ## Hide legend and define plot margins and size of axis title and text
    theme(
      legend.position = "none",
      plot.margin = unit(c(0.5, 0.4, 0.5, 0.4), "cm"),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = x_axis_size),
      axis.text.x = element_text(angle = x_axis_angle, hjust = x_axis_hjust, vjust = x_axis_vjust),
    )
  
  return(plot)
}

## Plot all boxplots for a single qc variable
plot_boxplots <- function(tissue, age){
  for (qc_var in qc_metrics){
    plots<-list()
    i=1
    for (pheno_var in sample_variables){
      p<-QC_boxplots(qc_var, pheno_var, tissue, age)
      plots[[i]]=p
      i=i+1
    }
    plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], 
              plots[[7]], nrow = 2)
    ## Save plots
    if (is.null(age)) {fileName=paste("plots/03_EDA/02_QC/boxplot_",qc_var,"_", 
                                      tissue, ".pdf", sep="")}
    else {fileName=paste("plots/03_EDA/02_QC/boxplot_",qc_var,"_", tissue, "_", 
                         age, ".pdf", sep="")}
    ggsave(fileName, width = 35, height = 17, units = "cm")
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



## 1.1.2 Examine relationships between QC metrics of samples

QC_scatterplots <- function(sample_var, qc_metric1, qc_metric2, tissue, age) {
  
  if (is.null(age)){
    ## Tissue data
    RSE<-eval(parse_expr(paste("rse_gene_", tissue, sep="")))
  }
  else {
    ## Tissue and Age data
    RSE<-eval(parse_expr(paste("rse_gene", tissue, age, sep="_")))
  }
  
  ## Define sample colors
  if (sample_var == "Group") {
    colors <- c("Control" = "seashell3", "Experimental" = "orange3")
  } else if(sample_var=='Expt'){
    colors <- c("Nicotine" = "lightblue3", "Smoking" = "salmon")
  } else if (sample_var == "Age") {
    colors <- c("Adult" = "slateblue3", "Pup" = "yellow3")
  } else if (sample_var == "Sex") {
    colors <- c("F" = "hotpink1", "M" = "dodgerblue")
  } else if (sample_var == "Pregnancy") {
    colors <- c("Yes" = "darkorchid3", "No" = "darkolivegreen4")
  } else if (sample_var == "plate") {
    colors <- c("Plate1" = "darkorange", "Plate2" = "lightskyblue", "Plate3" = "deeppink1")
  } else if (sample_var == "flowcell") {
    colors <- c(
      "HKCG7DSXX" = "chartreuse2", "HKCMHDSXX" = "magenta", "HKCNKDSXX" = "turquoise3",
      "HKCTMDSXX" = "tomato", "HK7JHDSXX"="seagreen3", "HKCJCDSXX"="palevioletred2"
    )
  } else if(sample_var=='Retention'){
    colors <- c('False'='indianred', 'True'='darkseagreen3')
  }
  
  ## Axis text size
  if (qc_metric1 == "sum") {
    x_axis_size <- 7
  } else {
    x_axis_size <- 9
  }
  
  
  ## Data
  data <- colData(RSE)
  
  ## Axis labels
  if (qc_metric1=="detected"){
    lab_x="Number of detected genes"}
  else if (qc_metric1=="subsets_Mito_percent"){
    lab_x="Percentage of mt counts"}
  else if (qc_metric1=="subsets_Ribo_percent"){
    lab_x="Percentage of ribosomal counts"}
  else if(qc_metric1=='sum'){
    lab_x='Library size'}
  else{
    lab_x = gsub("_", " ", qc_metric1)}
  
  if (qc_metric2=="detected"){
    lab_y="Number of detected genes"}
  else if (qc_metric2=="subsets_Mito_percent"){
    lab_y="Percentage of mt counts"}
  else if (qc_metric2=="subsets_Ribo_percent"){
    lab_y="Percentage of ribosomal counts"}
  else if(qc_metric2=='sum'){
    lab_y='Library size'}
  else{
    lab_y = gsub("_", " ", qc_metric2)}
  
  ## Sample labels (for the ones that were filtered out)
  if(sample_var!='Retention' | (is.null(age) & tissue=='brain')){
    data$Dropped_samples <- rep(NA, dim(data)[1])
  }
  
  ## Scatterplots: first QC metric in x-axis and second QC metric in y-axis
  plot <- ggplot(as.data.frame(data), aes(
    x = !!rlang::sym(qc_metric1),
    y = !!rlang::sym(qc_metric2),
    ## Color samples by a variable
    color = eval(parse_expr(sample_var)), 
    label = Dropped_samples
  )) +
    ## Add scatterplot
    geom_point(size = 1) +
    ## Colors
    scale_color_manual(name = capitalize(sample_var), values = colors) +
    ## Labels
    geom_text_repel(aes(label=Dropped_samples), size=2.6, color='black', min.segment.length = 0., 
                    max.overlaps=Inf, box.padding = 0.4) +
    theme_bw() +
    ## Add Pearson correlation coefficient between the metrics as subtitle
    labs(
      ## Add axis labels
      y = lab_y,
      x = lab_x
    ) +
    ## Plot margins and text size
    theme(
      plot.margin = unit(c(0.1, 1.2, 0.1, 1.2), "cm"),
      axis.title = element_text(size = (10)),
      axis.text = element_text(size = (x_axis_size)),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10)
    )
  
  if(sample_var!='Retention'){
    plot <- plot +  
      ## Add regression line
      stat_smooth(geom = "line", alpha = 0.7, size = 0.7, span = 0.25, method = lm, color = "orangered3") +
      ## Add Pearson correlation coefficient between the metrics as subtitle
      labs(subtitle = paste0("Corr: ", signif(cor(data[, qc_metric1], data[, qc_metric2], method = "pearson"), digits = 3))) +
      theme(plot.subtitle = element_text(size = 9, color = "gray40"))
  }
  
  return(plot)
}


## QC scatterplots coloring by all sample variables
multiple_QC_scatterplots <- function(qc_metric1, qc_metric2, tissue, age) {
  
  i <- 1
  plots <- list()
  for (sample_var in sample_variables) {
    plots[[i]] <- QC_scatterplots(sample_var, qc_metric1, qc_metric2, tissue, age)
    i <- i + 1
  }
  plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], nrow = 3, rel_widths = c(1, 1))
  ## Save plots
  if (is.null(age)){fileName=paste("plots/03_EDA/02_QC/",qc_metric1,"_vs_",qc_metric2, '_', tissue,".pdf", sep="")}
  else {fileName=paste("plots/03_EDA/02_QC/",qc_metric1,"_vs_", qc_metric2, '_', tissue,"_", age,".pdf", sep="")}
  ggsave(fileName, width = 35, height = 23, units = "cm")
}

## Plots of mito VS ribo percentages per sample 
multiple_QC_scatterplots('subsets_Mito_percent', 'subsets_Ribo_percent', 'brain', NULL)
multiple_QC_scatterplots('subsets_Mito_percent', 'subsets_Ribo_percent', 'blood', NULL)
multiple_QC_scatterplots('subsets_Mito_percent', 'subsets_Ribo_percent', 'brain', 'pups')
multiple_QC_scatterplots('subsets_Mito_percent', 'subsets_Ribo_percent', 'brain', 'adults')

## Plots of library size VS mito percentages per sample 
multiple_QC_scatterplots('sum', 'subsets_Mito_percent', 'brain', NULL)
multiple_QC_scatterplots('sum', 'subsets_Mito_percent', 'blood', NULL)
multiple_QC_scatterplots('sum', 'subsets_Mito_percent', 'brain', 'pups')
multiple_QC_scatterplots('sum', 'subsets_Mito_percent', 'brain', 'adults')

## Plots of library size VS ribo percentages per sample 
multiple_QC_scatterplots('sum', 'subsets_Ribo_percent', 'brain', NULL)
multiple_QC_scatterplots('sum', 'subsets_Ribo_percent', 'blood', NULL)
multiple_QC_scatterplots('sum', 'subsets_Ribo_percent', 'brain', 'pups')
multiple_QC_scatterplots('sum', 'subsets_Ribo_percent', 'brain', 'adults')

## Plots of library size VS number of detected genes per sample 
multiple_QC_scatterplots('sum', 'detected', 'brain', NULL)
multiple_QC_scatterplots('sum', 'detected', 'blood', NULL)
multiple_QC_scatterplots('sum', 'detected', 'brain', 'pups')
multiple_QC_scatterplots('sum', 'detected', 'brain', 'adults')





# ------------------------------------------------------------------------------
##  1.2  Sample filtering by QC 
# ------------------------------------------------------------------------------
## (Based on their detected number of genes, library sizes and proportions of ribo/mt counts)

# Filter brain samples
## Drop samples with lower library sizes and number of genes
outliers_sum<-isOutlier(rse_gene_brain$sum, nmads = 3, type="lower")
outliers_det<-isOutlier(rse_gene_brain$detected, nmads = 3, type="lower")
## Drop samples with higher mt and ribo percentages
outliers_mito<-isOutlier(rse_gene_brain$subsets_Mito_percent, nmads = 3, type="higher")
outliers_ribo<-isOutlier(rse_gene_brain$subsets_Ribo_percent, nmads = 3, type="higher")
not_outliers<-which(! (outliers_sum | outliers_det | outliers_mito | outliers_ribo))
rse_gene_brain_qc<-rse_gene_brain[,not_outliers]

## Number of samples retained
dim(rse_gene_brain_qc)[2]
# 135
## Save data
save(rse_gene_brain_qc, file = 'processed-data/03_EDA/02_QC/rse_gene_brain_qc.Rdata')


## Filter adult brain samples
outliers_sum<-isOutlier(rse_gene_brain_adults$sum, nmads = 3, type="lower")
outliers_det<-isOutlier(rse_gene_brain_adults$detected, nmads = 3, type="lower")
outliers_mito<-isOutlier(rse_gene_brain_adults$subsets_Mito_percent, nmads = 3, type="higher")
outliers_ribo<-isOutlier(rse_gene_brain_adults$subsets_Ribo_percent, nmads = 3, type="higher")
not_outliers<-which(! (outliers_sum | outliers_det | outliers_mito | outliers_ribo))
rse_gene_brain_adults_qc<-rse_gene_brain_adults[,not_outliers]

## Number of samples retained
dim(rse_gene_brain_adults_qc)[2]
# 42
## Save data
save(rse_gene_brain_adults_qc, file = 'processed-data/03_EDA/02_QC/rse_gene_brain_adults_qc.Rdata')

## Filter data of exons, tx and jx by QC
rse_exon_brain_adults_qc<-rse_exon_brain_adults[,(rse_exon_brain_adults$SAMPLE_ID %in% rse_gene_brain_adults_qc$SAMPLE_ID)]
save(rse_exon_brain_adults_qc, file = 'processed-data/03_EDA/02_QC/rse_exon_brain_adults_qc.Rdata')
rse_tx_brain_adults_qc<-rse_tx_brain_adults[,(rse_tx_brain_adults$SAMPLE_ID %in% rse_gene_brain_adults_qc$SAMPLE_ID)]
save(rse_tx_brain_adults_qc, file = 'processed-data/03_EDA/02_QC/rse_tx_brain_adults_qc.Rdata')
rse_jx_brain_adults_qc<-rse_jx_brain_adults[,(rse_jx_brain_adults$SAMPLE_ID %in% rse_gene_brain_adults_qc$SAMPLE_ID)]
save(rse_jx_brain_adults_qc, file = 'processed-data/03_EDA/02_QC/rse_jx_brain_adults_qc.Rdata')


# Filter pup brain samples
outliers_sum<-isOutlier(rse_gene_brain_pups$sum, nmads = 3, type="lower")
outliers_det<-isOutlier(rse_gene_brain_pups$detected, nmads = 3, type="lower")
outliers_mito<-isOutlier(rse_gene_brain_pups$subsets_Mito_percent, nmads = 3, type="higher")
outliers_ribo<-isOutlier(rse_gene_brain_pups$subsets_Ribo_percent, nmads = 3, type="higher")
not_outliers<-which(! (outliers_sum | outliers_det | outliers_mito | outliers_ribo))
rse_gene_brain_pups_qc<-rse_gene_brain_pups[,not_outliers]

## Number of samples retained
dim(rse_gene_brain_pups_qc)[2]
# 133
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
outliers_sum<-isOutlier(rse_gene_blood$sum, nmads = 3, type="lower")
outliers_det<-isOutlier(rse_gene_blood$detected, nmads = 3, type="lower")
outliers_mito<-isOutlier(rse_gene_blood$subsets_Mito_percent, nmads = 3, type="higher")
outliers_ribo<-isOutlier(rse_gene_blood$subsets_Ribo_percent, nmads = 3, type="higher")
not_outliers<-which(! (outliers_sum | outliers_det | outliers_mito | outliers_ribo))
rse_gene_blood_qc<-rse_gene_blood[,not_outliers]

## Number of samples retained
dim(rse_gene_blood_qc)[2]
# 23
## Save data
save(rse_gene_blood_qc, file = 'processed-data/03_EDA/02_QC/rse_gene_blood_qc.Rdata')



## 1.2.1 Examine QC metrics of samples retained and dropped

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
  names<-vector()
  for (sample in RSE$SAMPLE_ID){
    ## Samples retained
    if(sample %in% RSE_qc$SAMPLE_ID){
      qc_filt<-append(qc_filt, "True")
      ## Samples' names
      names<-append(names, "")
      }
    ## Samples dropped
    else {
      qc_filt<-append(qc_filt, "False")
      names<-append(names, sample)
      }
  }
  return(list(qc_filt, names))
}

## Add Retention column to all RSE
rse_gene_brain$Retention<-qc_filt_col("brain", NULL)[[1]]
rse_gene_blood$Retention<-qc_filt_col("blood", NULL)[[1]]
rse_gene_brain_adults$Retention<-qc_filt_col("brain", "adults")[[1]]
rse_gene_brain_pups$Retention<-qc_filt_col("brain", "pups")[[1]]

## Add column with names of samples dropped 
rse_gene_brain$Dropped_samples<-qc_filt_col("brain", NULL)[[2]]
rse_gene_blood$Dropped_samples<-qc_filt_col("blood", NULL)[[2]]
rse_gene_brain_adults$Dropped_samples<-qc_filt_col("brain", "adults")[[2]]
rse_gene_brain_pups$Dropped_samples<-qc_filt_col("brain", "pups")[[2]]


## Generate plots with samples separated by Retention and with labels
plots_Retained <- function(tissue, age){
  ## Plot total counts VS number of genes
  p1<-QC_scatterplots("Retention", 'sum', 'detected', tissue, age)

  ## Plot total counts VS mt percentage
  p2<-QC_scatterplots("Retention", 'sum', 'subsets_Mito_percent', tissue, age)
  
  ## Plot total counts VS ribo percentage
  p3<-QC_scatterplots("Retention", 'sum', 'subsets_Ribo_percent', tissue, age)
  
  ## Plot mt VS ribosomal percentage
  p4<-QC_scatterplots("Retention", 'subsets_Mito_percent',  'subsets_Ribo_percent', tissue, age)
  
  plot_grid(p1, p2, p3, p4, nrow = 2)
  ## Save plots
  if (is.null(age)){
    fileName=paste("plots/03_EDA/02_QC/samples_Retained_", tissue,".pdf", sep="")
  }
  else {
    fileName=paste("plots/03_EDA/02_QC/samples_Retained_", tissue,"_", age,".pdf", sep="")
  }
  ggsave(fileName, width = 23.6, height = 14, units = "cm")
}

## Plots
plots_Retained("brain", NULL)
plots_Retained("blood", NULL)
plots_Retained("brain", "adults")
plots_Retained("brain", "pups")



## 1.2.2 Map QC metrics after sample filtering 

## Boxplots of QC metrics for samples filtered and retained
boxplots_after_QC_filtering <- function(qc_metric, sample_var, tissue, age) {
  
  if (is.null(age)){
    ## Tissue data
    rse_gene<-eval(parse_expr(paste("rse_gene_", tissue, sep="")))
  }
  else {
    ## Tissue and Age data
    rse_gene<-eval(parse_expr(paste("rse_gene", tissue, age, sep="_")))
  }
  
  ## Color samples
  colors <- c("True" = "darkseagreen3", "False" = "indianred")
  
  ## Sample shape by sample variables
  if (sample_var == "Group") {
    shapes <- c("Control" = 0, "Experimental" = 15)
  } else if(sample_var == 'Expt'){
    shapes <- c("Nicotine" = 8, "Smoking" = 6)
  }else if (sample_var == "Age") {
    shapes <- c("Adult" = 16, "Pup" = 1)
  } else if (sample_var == "Sex") {
    shapes <- c("F" = 11, "M" = 17)
  } else if (sample_var == "Pregnancy") {
    shapes <- c("Yes" = 10, "No" = 1)
  } else if (sample_var == "plate") {
    shapes <- c("Plate1" = 12, "Plate2" = 5, "Plate3" = 4)
  } else if (sample_var == "flowcell") {
    shapes <- c(
      "HKCG7DSXX" = 3, "HKCMHDSXX" = 9, "HKCNKDSXX" = 14,
      "HKCTMDSXX" = 7, "HK7JHDSXX"=2, "HKCJCDSXX"=13
    )
  }
  
  y_label <- str_replace_all(qc_metric, c("_" = " "))
  
  data <- data.frame(colData(rse_gene))
  
  ## Median of the QC var values
  median <- median(eval(parse_expr(paste("rse_gene$", qc_metric, sep = ""))))
  ## Median-absolute-deviation of the QC var values
  mad <- mad(eval(parse_expr(paste("rse_gene$", qc_metric, sep = ""))))
  
  options(repr.plot.width = 4, repr.plot.height =6) 
  plot <- ggplot(data = data, mapping = aes(
    x = "", y = !!rlang::sym(qc_metric),
    color = Retention
  )) +
    geom_jitter(alpha = 0.8, size = 2.5, aes(shape = eval(parse_expr((sample_var)))), stroke = 1) +
    geom_boxplot(alpha = 0, size = 0.45, color = "black", ) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    guides(color = guide_legend(order = 1), 
           shape = guide_legend(order = 2)) +
    labs(x = "", y = y_label, color = "Retention after QC filtering", shape = capitalize(sample_var)) +
    theme_bw() +
    ## Median line
    geom_hline(yintercept = median, size = 0.5) +
    ## Line of median + 3 MADs
    geom_hline(yintercept = median + (3 * mad), size = 0.5, linetype = 2) +
    ## Line of median - 3 MADs
    geom_hline(yintercept = median - (3 * mad), size = 0.5, linetype = 2) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12),
      plot.margin = unit(c(0.1, 1.2, 0.1, 1.2), "cm"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11)
    )
  
  return(plot)
}

## Multiple boxplots
plot_QC_boxplots<- function(tissue, age){
  for (sample_var in sample_variables){
    plots<-list()
    i=1
    for (qc_var in c("sum", "detected", "subsets_Mito_percent", "subsets_Ribo_percent")){
      p<-boxplots_after_QC_filtering(qc_var, sample_var, tissue, age)
      plots[[i]]=p
      i=i+1
    }
    plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], nrow = 2, align="v")
    ## Save plots
    if (is.null(age)) {fileName=paste("plots/03_EDA/02_QC/QC_boxplot_",sample_var,"_", tissue, ".pdf", sep="")}
    else {fileName=paste("plots/03_EDA/02_QC/QC_boxplot_",sample_var,"_", tissue, "_", age, ".pdf", sep="")}
    ggsave(fileName, width = 30, height = 14, units = "cm")
  }
}

## QC boxplots
plot_QC_boxplots("brain", NULL)
plot_QC_boxplots("blood", NULL)
plot_QC_boxplots("brain", "adults")
plot_QC_boxplots("brain", "pups")







## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

#  version  R version 4.2.0 (2022-07-16 ucrt)
#  os       Windows 10 x64 (build 19044)
#  system   x86_64, mingw32
#  ui       RStudio
#  language (EN)
#  collate  Spanish_Mexico.utf8
#  ctype    Spanish_Mexico.utf8
#  tz       America/Mexico_City
