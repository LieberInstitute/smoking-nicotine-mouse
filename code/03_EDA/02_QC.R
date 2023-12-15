
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
library(ggnewscale)
library(colorblindr)
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
      "HKCTMDSXX" = "tomato", "HK7JHDSXX"="seagreen3", "HKCJCDSXX"="palevioletred2")
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
    colors <- c('False'='chocolate', 'True'='lightsteelblue3')
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
  colors <- c("True" = "lightsteelblue3", "False" = "chocolate")
  
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

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.0 (2023-04-21)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Mexico_City
# date     2023-12-12
# rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
# pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version    date (UTC) lib source
# AnnotationDbi          1.63.2     2023-07-03 [1] Bioconductor
# aod                    1.3.2      2022-04-02 [1] CRAN (R 4.3.0)
# backports              1.4.1      2021-12-13 [1] CRAN (R 4.3.0)
# base64enc              0.1-3      2015-07-28 [1] CRAN (R 4.3.0)
# beachmat               2.16.0     2023-05-08 [1] Bioconductor
# beeswarm               0.4.0      2021-06-01 [1] CRAN (R 4.3.0)
# Biobase              * 2.61.0     2023-06-02 [1] Bioconductor
# BiocFileCache          2.9.1      2023-07-14 [1] Bioconductor
# BiocGenerics         * 0.47.0     2023-06-02 [1] Bioconductor
# BiocIO                 1.11.0     2023-06-02 [1] Bioconductor
# BiocNeighbors          1.18.0     2023-05-08 [1] Bioconductor
# BiocParallel         * 1.35.3     2023-07-07 [1] Bioconductor
# BiocSingular           1.16.0     2023-05-08 [1] Bioconductor
# biomaRt                2.57.1     2023-06-14 [1] Bioconductor
# Biostrings             2.69.2     2023-07-05 [1] Bioconductor
# bit                    4.0.5      2022-11-15 [1] CRAN (R 4.3.0)
# bit64                  4.0.5      2020-08-30 [1] CRAN (R 4.3.0)
# bitops                 1.0-7      2021-04-24 [1] CRAN (R 4.3.0)
# blob                   1.2.4      2023-03-17 [1] CRAN (R 4.3.0)
# boot                   1.3-28.1   2022-11-22 [1] CRAN (R 4.3.0)
# broom                  1.0.5      2023-06-09 [1] CRAN (R 4.3.0)
# BSgenome               1.69.0     2023-07-07 [1] Bioconductor
# bumphunter             1.43.0     2023-06-14 [1] Bioconductor
# cachem                 1.0.8      2023-05-01 [1] CRAN (R 4.3.0)
# caTools                1.18.2     2021-03-28 [1] CRAN (R 4.3.0)
# checkmate              2.2.0      2023-04-27 [1] CRAN (R 4.3.0)
# cli                    3.6.1      2023-03-23 [1] CRAN (R 4.3.0)
# cluster                2.1.4      2022-08-22 [1] CRAN (R 4.3.0)
# codetools              0.2-19     2023-02-01 [1] CRAN (R 4.3.0)
# colorblindr          * 0.1.0      2023-12-13 [1] Github (clauswilke/colorblindr@90d64f8)
# colorspace             2.1-0      2023-01-23 [1] CRAN (R 4.3.0)
# corpcor                1.6.10     2021-09-16 [1] CRAN (R 4.3.0)
# cowplot              * 1.1.1      2020-12-30 [1] CRAN (R 4.3.0)
# crayon                 1.5.2      2022-09-29 [1] CRAN (R 4.3.0)
# curl                   5.0.1      2023-06-07 [1] CRAN (R 4.3.0)
# data.table             1.14.8     2023-02-17 [1] CRAN (R 4.3.0)
# DBI                    1.1.3      2022-06-18 [1] CRAN (R 4.3.0)
# dbplyr                 2.3.3      2023-07-07 [1] CRAN (R 4.3.0)
# DelayedArray           0.26.6     2023-07-02 [1] Bioconductor
# DelayedMatrixStats     1.23.0     2023-04-25 [1] Bioconductor
# derfinder              1.35.0     2023-07-07 [1] Bioconductor
# derfinderHelper        1.35.0     2023-06-02 [1] Bioconductor
# digest                 0.6.33     2023-07-07 [1] CRAN (R 4.3.0)
# doParallel             1.0.17     2022-02-07 [1] CRAN (R 4.3.0)
# doRNG                  1.8.6      2023-01-16 [1] CRAN (R 4.3.0)
# downloader             0.4        2015-07-09 [1] CRAN (R 4.3.0)
# dplyr                  1.1.2      2023-04-20 [1] CRAN (R 4.3.0)
# dynamicTreeCut       * 1.63-1     2016-03-11 [1] CRAN (R 4.3.0)
# edgeR                * 3.43.7     2023-06-21 [1] Bioconductor
# EnvStats               2.8.0      2023-07-08 [1] CRAN (R 4.3.0)
# evaluate               0.21       2023-05-05 [1] CRAN (R 4.3.0)
# fANCOVA                0.6-1      2020-11-13 [1] CRAN (R 4.3.0)
# fansi                  1.0.5      2023-10-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1      2022-07-06 [1] CRAN (R 4.3.0)
# fastcluster          * 1.2.3      2021-05-24 [1] CRAN (R 4.3.0)
# fastmap                1.1.1      2023-02-24 [1] CRAN (R 4.3.0)
# filelock               1.0.2      2018-10-05 [1] CRAN (R 4.3.0)
# foreach                1.5.2      2022-02-02 [1] CRAN (R 4.3.0)
# foreign                0.8-84     2022-12-06 [1] CRAN (R 4.3.0)
# Formula                1.2-5      2023-02-24 [1] CRAN (R 4.3.0)
# fs                     1.6.3      2023-07-20 [1] CRAN (R 4.3.0)
# gargle                 1.5.2      2023-07-20 [1] CRAN (R 4.3.0)
# generics               0.1.3      2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.37.2     2023-06-21 [1] Bioconductor
# GenomeInfoDbData       1.2.10     2023-05-28 [1] Bioconductor
# GenomicAlignments      1.37.0     2023-07-07 [1] Bioconductor
# GenomicFeatures        1.53.1     2023-06-22 [1] Bioconductor
# GenomicFiles           1.37.0     2023-07-07 [1] Bioconductor
# GenomicRanges        * 1.53.1     2023-06-02 [1] Bioconductor
# GEOquery               2.69.0     2023-06-16 [1] Bioconductor
# ggbeeswarm             0.7.2      2023-04-29 [1] CRAN (R 4.3.0)
# ggnewscale           * 0.4.9      2023-05-25 [1] CRAN (R 4.3.0)
# ggplot2              * 3.4.4      2023-10-12 [1] CRAN (R 4.3.1)
# ggrepel              * 0.9.3      2023-02-03 [1] CRAN (R 4.3.0)
# glue                   1.6.2      2022-02-24 [1] CRAN (R 4.3.0)
# GO.db                  3.17.0     2023-05-28 [1] Bioconductor
# googledrive            2.1.1      2023-06-11 [1] CRAN (R 4.3.0)
# gplots                 3.1.3      2022-04-25 [1] CRAN (R 4.3.0)
# gridExtra            * 2.3        2017-09-09 [1] CRAN (R 4.3.0)
# gtable                 0.3.4      2023-08-21 [1] CRAN (R 4.3.0)
# gtools                 3.9.4      2022-11-27 [1] CRAN (R 4.3.0)
# here                 * 1.0.1      2020-12-13 [1] CRAN (R 4.3.0)
# Hmisc                * 5.1-0      2023-05-08 [1] CRAN (R 4.3.0)
# hms                    1.1.3      2023-03-21 [1] CRAN (R 4.3.0)
# htmlTable              2.4.1      2022-07-07 [1] CRAN (R 4.3.0)
# htmltools              0.5.5      2023-03-23 [1] CRAN (R 4.3.0)
# htmlwidgets            1.6.2      2023-03-17 [1] CRAN (R 4.3.0)
# httr                   1.4.6      2023-05-08 [1] CRAN (R 4.3.0)
# IRanges              * 2.35.2     2023-06-23 [1] Bioconductor
# irlba                  2.3.5.1    2022-10-03 [1] CRAN (R 4.3.0)
# iterators              1.0.14     2022-02-05 [1] CRAN (R 4.3.0)
# jaffelab             * 0.99.32    2023-05-28 [1] Github (LieberInstitute/jaffelab@21e6574)
# jsonlite               1.8.7      2023-06-29 [1] CRAN (R 4.3.0)
# KEGGREST               1.41.0     2023-07-07 [1] Bioconductor
# KernSmooth             2.23-22    2023-07-10 [1] CRAN (R 4.3.0)
# knitr                  1.43       2023-05-25 [1] CRAN (R 4.3.0)
# labeling               0.4.3      2023-08-29 [1] CRAN (R 4.3.0)
# lattice                0.21-8     2023-04-05 [1] CRAN (R 4.3.0)
# lifecycle              1.0.3      2022-10-07 [1] CRAN (R 4.3.0)
# limma                * 3.57.6     2023-06-21 [1] Bioconductor
# lme4                   1.1-34     2023-07-04 [1] CRAN (R 4.3.0)
# lmerTest               3.1-3      2020-10-23 [1] CRAN (R 4.3.0)
# locfit                 1.5-9.8    2023-06-11 [1] CRAN (R 4.3.0)
# magrittr               2.0.3      2022-03-30 [1] CRAN (R 4.3.0)
# MASS                   7.3-60     2023-05-04 [1] CRAN (R 4.3.0)
# Matrix                 1.6-0      2023-07-08 [1] CRAN (R 4.3.0)
# MatrixGenerics       * 1.13.0     2023-05-20 [1] Bioconductor
# matrixStats          * 1.0.0      2023-06-02 [1] CRAN (R 4.3.0)
# memoise                2.0.1      2021-11-26 [1] CRAN (R 4.3.0)
# minqa                  1.2.5      2022-10-19 [1] CRAN (R 4.3.0)
# munsell                0.5.0      2018-06-12 [1] CRAN (R 4.3.0)
# mvtnorm                1.2-2      2023-06-08 [1] CRAN (R 4.3.0)
# nlme                   3.1-162    2023-01-31 [1] CRAN (R 4.3.0)
# nloptr                 2.0.3      2022-05-26 [1] CRAN (R 4.3.0)
# nnet                   7.3-19     2023-05-03 [1] CRAN (R 4.3.0)
# numDeriv               2016.8-1.1 2019-06-06 [1] CRAN (R 4.3.0)
# pbkrtest               0.5.2      2023-01-19 [1] CRAN (R 4.3.0)
# pillar                 1.9.0      2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.3.0)
# plyr                   1.8.8      2022-11-11 [1] CRAN (R 4.3.0)
# png                    0.1-8      2022-11-29 [1] CRAN (R 4.3.0)
# prettyunits            1.1.1      2020-01-24 [1] CRAN (R 4.3.0)
# progress               1.2.2      2019-05-16 [1] CRAN (R 4.3.0)
# purrr                  1.0.1      2023-01-10 [1] CRAN (R 4.3.0)
# qvalue                 2.33.0     2023-05-11 [1] Bioconductor
# R6                     2.5.1      2021-08-19 [1] CRAN (R 4.3.0)
# rafalib              * 1.0.0      2015-08-09 [1] CRAN (R 4.3.0)
# ragg                   1.2.5      2023-01-12 [1] CRAN (R 4.3.0)
# rappdirs               0.3.3      2021-01-31 [1] CRAN (R 4.3.0)
# rbibutils              2.2.13     2023-01-13 [1] CRAN (R 4.3.0)
# RColorBrewer           1.1-3      2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.11     2023-07-06 [1] CRAN (R 4.3.0)
# RCurl                  1.98-1.12  2023-03-27 [1] CRAN (R 4.3.0)
# Rdpack                 2.4        2022-07-20 [1] CRAN (R 4.3.0)
# readr                  2.1.4      2023-02-10 [1] CRAN (R 4.3.0)
# recount              * 1.27.0     2023-04-25 [1] Bioconductor
# remaCor                0.0.16     2023-06-21 [1] CRAN (R 4.3.0)
# rentrez                1.2.3      2020-11-10 [1] CRAN (R 4.3.0)
# reshape2               1.4.4      2020-04-09 [1] CRAN (R 4.3.0)
# restfulr               0.0.15     2022-06-16 [1] CRAN (R 4.3.0)
# RhpcBLASctl            0.23-42    2023-02-11 [1] CRAN (R 4.3.0)
# rjson                  0.2.21     2022-01-09 [1] CRAN (R 4.3.0)
# rlang                * 1.1.1      2023-04-28 [1] CRAN (R 4.3.0)
# rmarkdown              2.23       2023-07-01 [1] CRAN (R 4.3.0)
# rngtools               1.5.2      2021-09-20 [1] CRAN (R 4.3.0)
# rpart                  4.1.19     2022-10-21 [1] CRAN (R 4.3.0)
# rprojroot              2.0.3      2022-04-02 [1] CRAN (R 4.3.0)
# Rsamtools              2.17.0     2023-07-07 [1] Bioconductor
# RSQLite                2.3.1      2023-04-03 [1] CRAN (R 4.3.0)
# rstudioapi             0.15.0     2023-07-07 [1] CRAN (R 4.3.0)
# rsvd                   1.0.5      2021-04-16 [1] CRAN (R 4.3.0)
# rtracklayer            1.61.0     2023-07-07 [1] Bioconductor
# S4Arrays               1.1.4      2023-06-02 [1] Bioconductor
# S4Vectors            * 0.39.1     2023-06-02 [1] Bioconductor
# ScaledMatrix           1.9.1      2023-05-03 [1] Bioconductor
# scales                 1.2.1      2022-08-20 [1] CRAN (R 4.3.0)
# scater               * 1.29.1     2023-06-15 [1] Github (davismcc/scater@eb6b801)
# scuttle              * 1.9.4      2023-01-23 [1] Bioconductor
# segmented              1.6-4      2023-04-13 [1] CRAN (R 4.3.0)
# sessioninfo          * 1.2.2      2021-12-06 [1] CRAN (R 4.3.0)
# SingleCellExperiment * 1.23.0     2023-04-25 [1] Bioconductor
# sparseMatrixStats      1.13.0     2023-05-20 [1] Bioconductor
# stringi                1.7.12     2023-01-11 [1] CRAN (R 4.3.0)
# stringr              * 1.5.0      2022-12-02 [1] CRAN (R 4.3.0)
# SummarizedExperiment * 1.30.2     2023-06-06 [1] Bioconductor
# survival               3.5-5      2023-03-12 [1] CRAN (R 4.3.0)
# systemfonts            1.0.4      2022-02-11 [1] CRAN (R 4.3.0)
# textshaping            0.3.6      2021-10-13 [1] CRAN (R 4.3.0)
# tibble                 3.2.1      2023-03-20 [1] CRAN (R 4.3.0)
# tidyr                  1.3.0      2023-01-24 [1] CRAN (R 4.3.0)
# tidyselect             1.2.0      2022-10-10 [1] CRAN (R 4.3.0)
# tzdb                   0.4.0      2023-05-12 [1] CRAN (R 4.3.0)
# utf8                   1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
# variancePartition    * 1.32.2     2023-11-14 [1] Bioconductor
# VariantAnnotation      1.47.1     2023-07-07 [1] Bioconductor
# vctrs                  0.6.4      2023-10-12 [1] CRAN (R 4.3.1)
# vipor                  0.4.5      2017-03-22 [1] CRAN (R 4.3.0)
# viridis                0.6.3      2023-05-03 [1] CRAN (R 4.3.0)
# viridisLite            0.4.2      2023-05-02 [1] CRAN (R 4.3.0)
# withr                  2.5.1      2023-09-26 [1] CRAN (R 4.3.1)
# xfun                   0.39       2023-04-20 [1] CRAN (R 4.3.0)
# XML                    3.99-0.14  2023-03-19 [1] CRAN (R 4.3.0)
# xml2                   1.3.5      2023-07-06 [1] CRAN (R 4.3.0)
# XVector                0.41.1     2023-06-02 [1] Bioconductor
# yaml                   2.3.7      2023-01-23 [1] CRAN (R 4.3.0)
# zlibbioc               1.47.0     2023-05-20 [1] Bioconductor
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

