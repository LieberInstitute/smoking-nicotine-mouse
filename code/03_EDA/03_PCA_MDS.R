
# 1. Explore sample effects
##  (Dimensionality Reduction) 

## Load data
load(here("processed-data/02_build_objects/rse_gene_brain.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_gene_blood_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_gene_brain_adults_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_exon_brain_adults_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_tx_brain_adults_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_jx_brain_adults_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_gene_brain_pups_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_exon_brain_pups_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_tx_brain_pups_qc.Rdata"))
load(here("processed-data/03_EDA/02_QC/rse_jx_brain_pups_qc.Rdata"))
## Use not filtered brain data
rse_gene_brain_qc<-rse_gene_brain


# ------------------------------------------------------------------------------
##  1.1  Principal Component Analysis (PCA)
# ------------------------------------------------------------------------------
### Explore samples' expression variation

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

## Colors 
colors = list("Group"=c("Control" = "seashell3", "Experimental" = "orange3"),
              "Expt"=c("Nicotine" = "lightblue3", "Smoking" = "salmon"),
              "Age"=c("Adult" = "slateblue3", "Pup" = "yellow3"),
              "Sex"=c("F" = "hotpink1", "M" = "dodgerblue"),
              "Pregnancy"=c("Yes" = "darkorchid3", "No" = "darkolivegreen4"),
              "plate"=c("Plate1" = "darkorange", "Plate2" = "lightskyblue", "Plate3" = "deeppink1"),
              "flowcell"=c("HKCG7DSXX" = "chartreuse2", "HKCMHDSXX" = "magenta", "HKCNKDSXX" = "turquoise3",
                           "HKCTMDSXX" = "tomato", "HK7JHDSXX"="seagreen3", "HKCJCDSXX"="palevioletred2")
)

## Colors to highlight poor QC samples
poorQC_samples_colors <- c("Sample_FE3P2"="magenta2", "Sample_4067"="darkorange2", "Sample_FC41"="blue", 
                           "Sample_P2_fe2_022019"='cyan2', "Sample_P1_fe3_021819"='purple2', "Sample_P7_fe3_021719"='green1')

PCx_vs_PCy <- function (PCx, PCy, pca_data, pca_vars, sample_var, level) {
  
  if(unique(pca_data$Tissue)=='Brain' & unique(pca_data$Age)=='Pup' & level=='jx'){
    axis_text_size = 8
  }
  else{
    axis_text_size = 10
  }
  
  plot=ggplot(data=pca_data, 
       aes(x=eval(parse_expr(PCx)),y=eval(parse_expr(PCy))))+ 
        geom_point(aes(color=eval(parse_expr(sample_var))), size=2) + 
        scale_color_manual(name = capitalize(sample_var), values = colors[[sample_var]]) +
        theme_bw() + 
        labs(x= pca_vars[strtoi(gsub("PC","", PCx))], y = pca_vars[strtoi(gsub("PC","", PCy))]) +
        ## Add square around segregated samples
        new_scale_color() +
        geom_point(data=subset(pca_data, !is.na(label)), aes(color=label), pch = 0, size=4, stroke = 1) +
        scale_color_manual(name=label, values = poorQC_samples_colors) +
        guides(color = 'none') + 
        new_scale_color() +
        geom_point(data=subset(pca_data, !is.na(label)), aes(color=eval(parse_expr(sample_var))), size=2) +
        scale_color_manual(values = colors[[sample_var]]) +
        guides(color = 'none') + 
        theme(plot.margin=unit (c (1,1.5,1,1), 'cm'), 
              legend.position = "right",
              legend.text = element_text(size = 11),
              legend.title = element_text(size = 12),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = axis_text_size)) 
  return(plot)
}

## All PCA plots 
plot_PCAs<-function(type, tissue, age, afterFiltering){
  
  ## PC data
  pca_results<-PCA(tissue, type, age)
  pca_data<-pca_results[[1]]
  pca_vars<-pca_results[[2]]
  
  ## Add label of rare samples
  if(!is.null(age)){
    pca_data$label <- sapply(pca_data$SAMPLE_ID, function(x){if(x %in% c("Sample_FE3P2", "Sample_4067", "Sample_FC41", 
                                                                       "Sample_P2_fe2_022019", "Sample_P1_fe3_021819", "Sample_P7_fe3_021719")){x} else{NA}})
  }
  else{
    pca_data$label <- rep(NA, dim(pca_data)[1])}
  
  ## File name termination for plots after manual sample filtering  
  if(afterFiltering==TRUE){
    term <- '_afterManualFilt.pdf'
  }
  else{
    term <- '.pdf'
  }
  
  
  ## PC pairs to plot
  PC_pairs <- list(c("PC1", "PC2"), c("PC3", "PC4"), c("PC5", "PC6"))
  
  ## Plots for blood and brain 
  if (is.null(age)){
    if (tissue=="blood"){
      for (PCs in PC_pairs){
        plots<-list()
        i=1
        for (sample_var in c("Group", "Pregnancy", "plate", "flowcell")){
           p<-PCx_vs_PCy(PCs[1], PCs[2], pca_data, pca_vars, sample_var, type)
           plots[[i]]=p
           i=i+1
        }
        plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], nrow = 2, align="v")
        ## Save plots
        ggsave(paste("plots/03_EDA/03_PCA_MDS/",PCs[1],"_vs_",PCs[2],"_", type, "_", tissue , term, sep=""), 
                 width = 26, height = 16, units = "cm")
      }
    }  
    ## For brain
    else {
       for (PCs in PC_pairs){
        plots<-list()
        i=1
        for (sample_var in c("Group", "Expt", "Age", "Sex", "Pregnancy", "plate", "flowcell")){
           p<-PCx_vs_PCy(PCs[1], PCs[2], pca_data, pca_vars, sample_var, type)
           plots[[i]]=p
           i=i+1
        }
        plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], nrow = 2, align="v")
        ## Save plots
        ggsave(paste("plots/03_EDA/03_PCA_MDS/",PCs[1],"_vs_",PCs[2],"_", type, "_", tissue , term, sep=""), 
                 width = 52, height = 16, units = "cm")
      }
    }
  }
  
  ## Plots for adult brain 
  else if (age=="adults") {
    for (PCs in PC_pairs){
      plots<-list()
      i=1
      for (pheno_var in c("Group", "Expt", "Pregnancy", "plate", "flowcell")){
         p<-PCx_vs_PCy(PCs[1], PCs[2], pca_data, pca_vars, pheno_var, type)
         plots[[i]]=p
         i=i+1
      }
      plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], nrow = 2, align="v")
      ## Save plots
      ggsave(paste("plots/03_EDA/03_PCA_MDS/",PCs[1],"_vs_",PCs[2],"_", type, "_", tissue ,"_", age, term, sep=""), 
             width = 39, height = 16, units = "cm")
    }
  }
  
  ## Plots for brain and pups
  else if  (age=="pups") {
    for (PCs in PC_pairs){
      plots<-list()
      i=1
      for (pheno_var in c("Group", "Expt", "Sex", "plate", "flowcell")){
         p<-PCx_vs_PCy(PCs[1], PCs[2], pca_data, pca_vars, pheno_var, type)
         plots[[i]]=p
         i=i+1
      }
      plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], nrow = 2, align="v")
      ## Save plots
      ggsave(paste("plots/03_EDA/03_PCA_MDS/",PCs[1],"_vs_",PCs[2],"_", type, "_", tissue ,"_", age, term, sep=""), 
             width = 39, height = 16, units = "cm")

  }
 }
}

## Plots
## Gene level
plot_PCAs("gene", "blood", NULL, FALSE)
plot_PCAs("gene", "brain", NULL, FALSE)
plot_PCAs("gene", "brain", "adults", FALSE)
plot_PCAs("gene", "brain", "pups", FALSE)
## Exon level 
plot_PCAs("exon", "brain", "adults", FALSE)
plot_PCAs("exon", "brain", "pups", FALSE)
## Tx level
plot_PCAs("tx", "brain", "adults", FALSE)
plot_PCAs("tx", "brain", "pups", FALSE)
## Jxn level
plot_PCAs("jx", "brain", "adults", FALSE)
plot_PCAs("jx", "brain", "pups", FALSE)



### 1.1.1 Manual sample filtering 
## Removal of rare samples in brain PCA plots

## PC data
pca_data_gene_brain_adults<-PCA("brain", "gene", "adults")[[1]]
pca_data_exon_brain_adults<-PCA("brain", "exon", "adults")[[1]]
pca_data_tx_brain_adults<-PCA("brain", "tx", "adults")[[1]]
pca_data_jx_brain_adults<-PCA("brain", "jx", "adults")[[1]]
pca_data_gene_brain_pups<-PCA("brain", "gene", "pups")[[1]]
pca_data_exon_brain_pups<-PCA("brain", "exon", "pups")[[1]]
pca_data_tx_brain_pups<-PCA("brain", "tx", "pups")[[1]]
pca_data_jx_brain_pups<-PCA("brain", "jx", "pups")[[1]]


## In adult plots
pca_data_gene_brain_adults[which.max(pca_data_gene_brain_adults$PC2), "SAMPLE_ID"]
# "Sample_FE3P2"
pca_data_gene_brain_adults[which.max(pca_data_gene_brain_adults$PC3), "SAMPLE_ID"]
# "Sample_FE3P2"
pca_data_exon_brain_adults[which.max(pca_data_exon_brain_adults$PC2), "SAMPLE_ID"]
# "Sample_FE3P2"
pca_data_exon_brain_adults[which.min(pca_data_exon_brain_adults$PC3), "SAMPLE_ID"]
# "Sample_4067"
pca_data_tx_brain_adults[which.max(pca_data_tx_brain_adults$PC2), "SAMPLE_ID"]
# "Sample_FE3P2"
pca_data_jx_brain_adults[which.min(pca_data_jx_brain_adults$PC3), "SAMPLE_ID"]
# "Sample_FE3P2"
pca_data_jx_brain_adults[which.min(pca_data_jx_brain_adults$PC5), "SAMPLE_ID"]
# "Sample_4057"
pca_data_jx_brain_adults[which.min(pca_data_jx_brain_adults$PC6), "SAMPLE_ID"]
# "Sample_FC41"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## Explore Sample_FE3P2 info:
colData(rse_gene_brain_adults_qc)[which(rse_gene_brain_adults_qc$SAMPLE_ID=="Sample_FE3P2"),
                          c("mitoRate", "rRNA_rate", "totalAssignedGene", "overallMapRate", 
                            "ERCCsumLogErr", "sum", "Sex")]
## This sample has the max rRNA rate
pca_data_gene_brain_adults[which.max(pca_data_gene_brain_adults$rRNA_rate), "SAMPLE_ID"]
# "Sample_FE3P2"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## Explore Sample_4067 info:
colData(rse_gene_brain_adults_qc)[which(rse_gene_brain_adults_qc$SAMPLE_ID=="Sample_4067"),
                          c("mitoRate", "rRNA_rate", "totalAssignedGene", "overallMapRate", 
                            "ERCCsumLogErr", "sum", "Sex")]
## This sample has the max mito rate, % of mt and ribo counts and min prop of reads assigned to genes 
## and number of genes
pca_data_gene_brain_adults[which.max(pca_data_gene_brain_adults$mitoRate), "SAMPLE_ID"]
# "Sample_4067"
pca_data_gene_brain_adults[which.max(pca_data_gene_brain_adults$subsets_Mito_percent), "SAMPLE_ID"]
# "Sample_4067"
pca_data_gene_brain_adults[which.max(pca_data_gene_brain_adults$subsets_Ribo_percent), "SAMPLE_ID"]
# "Sample_4067"
pca_data_gene_brain_adults[which.min(pca_data_gene_brain_adults$totalAssignedGene), "SAMPLE_ID"]
# "Sample_4067"
pca_data_gene_brain_adults[which.min(pca_data_gene_brain_adults$detected), "SAMPLE_ID"]
# "Sample_4067"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## Explore Sample_FC41 info:
colData(rse_gene_brain_adults_qc)[which(rse_gene_brain_adults_qc$SAMPLE_ID=="Sample_FC41"),
                          c("mitoRate", "rRNA_rate", "totalAssignedGene", "overallMapRate", 
                            "ERCCsumLogErr", "sum", "Sex")]
## This sample has the min overall mapping rate and library size
pca_data_gene_brain_adults[which.min(pca_data_gene_brain_adults$overallMapRate), "SAMPLE_ID"]
# "Sample_FC41"
pca_data_gene_brain_adults[which.min(pca_data_gene_brain_adults$sum), "SAMPLE_ID"]
# "Sample_FC41"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## Explore Sample_4057 info:
colData(rse_gene_brain_adults_qc)[which(rse_gene_brain_adults_qc$SAMPLE_ID=="Sample_4057"),
                          c("mitoRate", "rRNA_rate", "totalAssignedGene", "overallMapRate", 
                            "ERCCsumLogErr", "sum", "Sex")]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


## In pup plots
pca_data_gene_brain_pups[which.min(pca_data_gene_brain_pups$PC2), "SAMPLE_ID"]
# "Sample_P2_fe2_022019"
## Male samples in females at gene level 
pca_data_gene_brain_pups[which(pca_data_gene_brain_pups$Sex=="M" & pca_data_gene_brain_pups$PC3>0),"SAMPLE_ID"]
# "Sample_P1_fe3_021819" "Sample_P2_fe2_022019"
pca_data_exon_brain_pups[which.min(pca_data_exon_brain_pups$PC1), "SAMPLE_ID"]
# "Sample_P2_fe2_022019"
pca_data_tx_brain_pups[which.min(pca_data_tx_brain_pups$PC1), "SAMPLE_ID"]
# "Sample_P2_fe2_022019"
pca_data_tx_brain_pups[which.min(pca_data_tx_brain_pups$PC6), "SAMPLE_ID"]
# "Sample_P2_fe2_022019"
pca_data_jx_brain_pups[which.max(pca_data_jx_brain_pups$PC1), "SAMPLE_ID"]
# "Sample_P2_fe2_022019"
pca_data_jx_brain_pups[which.max(pca_data_jx_brain_pups$PC3), "SAMPLE_ID"]
# "Sample_P5_fe1_021919"
pca_data_jx_brain_pups[which.min(pca_data_jx_brain_pups$PC5), "SAMPLE_ID"]
# "Sample_P7_fe3_021719"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## Explore Sample_P2_fe2_022019 info:
colData(rse_gene_brain_pups_qc)[which(rse_gene_brain_pups_qc$SAMPLE_ID=="Sample_P2_fe2_022019"),
                          c("mitoRate", "rRNA_rate", "totalAssignedGene", "overallMapRate", 
                            "ERCCsumLogErr", "sum", "Sex")]
## This sample has the min prop of reads assigned to genes 
pca_data_gene_brain_pups[which.min(pca_data_gene_brain_pups$totalAssignedGene), "SAMPLE_ID"]
# "Sample_P2_fe2_022019"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## Explore Sample_P1_fe3_021819 info:
colData(rse_gene_brain_pups_qc)[which(rse_gene_brain_pups_qc$SAMPLE_ID=="Sample_P1_fe3_021819"),
                          c("mitoRate", "rRNA_rate", "totalAssignedGene", "overallMapRate", 
                            "ERCCsumLogErr", "sum", "Sex")]
## This sample has the max Error value
pca_data_gene_brain_pups[which.max(abs(pca_data_gene_brain_pups$ERCCsumLogErr)), "SAMPLE_ID"]
# "Sample_P1_fe3_021819"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## Explore Sample_P5_fe1_021919 info:
colData(rse_gene_brain_pups_qc)[which(rse_gene_brain_pups_qc$SAMPLE_ID=="Sample_P5_fe1_021919"),
                          c("mitoRate", "rRNA_rate", "totalAssignedGene", "overallMapRate", 
                            "ERCCsumLogErr", "sum", "Sex")]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## Explore Sample_P7_fe3_021719 info:
colData(rse_gene_brain_pups_qc)[which(rse_gene_brain_pups_qc$SAMPLE_ID=="Sample_P7_fe3_021719"),
                          c("mitoRate", "rRNA_rate", "totalAssignedGene", "overallMapRate", 
                            "ERCCsumLogErr", "sum", "Sex")]
## This sample has the min overall mapping rate and library size
pca_data_gene_brain_pups[which.min(pca_data_gene_brain_pups$overallMapRate), "SAMPLE_ID"]
# "Sample_P7_fe3_021719"
pca_data_gene_brain_pups[which.min(pca_data_gene_brain_pups$sum), "SAMPLE_ID"]
# "Sample_P7_fe3_021719"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


## Remove those samples
## In adults
poorQC_samples<-c("Sample_FE3P2", "Sample_4067", "Sample_FC41")
rse_gene_brain_adults_qc <- rse_gene_brain_adults_qc_afterPCA<-rse_gene_brain_adults_qc[,-which(rse_gene_brain_adults_qc$SAMPLE_ID %in% poorQC_samples)]
rse_exon_brain_adults_qc <- rse_exon_brain_adults_qc_afterPCA<-rse_exon_brain_adults_qc[,-which(rse_exon_brain_adults_qc$SAMPLE_ID %in% poorQC_samples)]
rse_tx_brain_adults_qc <- rse_tx_brain_adults_qc_afterPCA<-rse_tx_brain_adults_qc[,-which(rse_tx_brain_adults_qc$SAMPLE_ID %in% poorQC_samples)]
rse_jx_brain_adults_qc <- rse_jx_brain_adults_qc_afterPCA<-rse_jx_brain_adults_qc[,-which(rse_jx_brain_adults_qc$SAMPLE_ID %in% poorQC_samples)]

## Save RSE without those samples
save(rse_gene_brain_adults_qc_afterPCA, file="processed-data/03_EDA/03_PCA/rse_gene_brain_adults_qc_afterPCA.Rdata")
save(rse_exon_brain_adults_qc_afterPCA, file="processed-data/03_EDA/03_PCA/rse_exon_brain_adults_qc_afterPCA.Rdata")
save(rse_tx_brain_adults_qc_afterPCA, file="processed-data/03_EDA/03_PCA/rse_tx_brain_adults_qc_afterPCA.Rdata")
save(rse_jx_brain_adults_qc_afterPCA, file="processed-data/03_EDA/03_PCA/rse_jx_brain_adults_qc_afterPCA.Rdata")


## In pups
poorQC_samples<-c("Sample_P2_fe2_022019", "Sample_P1_fe3_021819", "Sample_P7_fe3_021719")
rse_gene_brain_pups_qc <- rse_gene_brain_pups_qc_afterPCA<-rse_gene_brain_pups_qc[,-which(rse_gene_brain_pups_qc$SAMPLE_ID %in% poorQC_samples)]
rse_exon_brain_pups_qc <- rse_exon_brain_pups_qc_afterPCA<-rse_exon_brain_pups_qc[,-which(rse_exon_brain_pups_qc$SAMPLE_ID %in% poorQC_samples)]
rse_tx_brain_pups_qc <- rse_tx_brain_pups_qc_afterPCA<-rse_tx_brain_pups_qc[,-which(rse_tx_brain_pups_qc$SAMPLE_ID %in% poorQC_samples)]
rse_jx_brain_pups_qc <- rse_jx_brain_pups_qc_afterPCA<-rse_jx_brain_pups_qc[,-which(rse_jx_brain_pups_qc$SAMPLE_ID %in% poorQC_samples)]

## Save RSE without those samples
save(rse_gene_brain_pups_qc_afterPCA, file="processed-data/03_EDA/03_PCA/rse_gene_brain_pups_qc_afterPCA.Rdata")
save(rse_exon_brain_pups_qc_afterPCA, file="processed-data/03_EDA/03_PCA/rse_exon_brain_pups_qc_afterPCA.Rdata")
save(rse_tx_brain_pups_qc_afterPCA, file="processed-data/03_EDA/03_PCA/rse_tx_brain_pups_qc_afterPCA.Rdata")
save(rse_jx_brain_pups_qc_afterPCA, file="processed-data/03_EDA/03_PCA/rse_jx_brain_pups_qc_afterPCA.Rdata")


## PCA plots without those samples
plot_PCAs("gene", "brain", "adults", TRUE)
plot_PCAs("gene", "brain", "pups", TRUE)
plot_PCAs("exon", "brain", "adults", TRUE)
plot_PCAs("exon", "brain", "pups", TRUE)
plot_PCAs("tx", "brain", "adults", TRUE)
plot_PCAs("tx", "brain", "pups", TRUE)
plot_PCAs("jx", "brain", "adults", TRUE)
plot_PCAs("jx", "brain", "pups", TRUE)




### 1.1.2 Explore the separation of nicotine and smoking samples by Group

## PCA plots for Group and Expt
PCA_Expt_Group<- function(type, tissue, age){
  
  RSE<-eval(parse_expr(paste("rse", type, tissue, age, "qc", sep="_")))

  ## PCA plots for smoking samples separated by Group
  RSE_smoking<-RSE[,RSE$Expt=="Smoking"]
  pca<-prcomp(t(assays(RSE_smoking)$logcounts))
  # % of the variance explained by each PC
  pca_vars<- getPcaVars(pca)
  pca_vars_labs<- paste0(
      "PC", seq(along = pca_vars), ": ",
      pca_vars, "% Var Expl")
  
  ## Join PCs and samples' info 
  pca_data<-cbind(pca$x,colData(RSE_smoking))
  pca_data<-as.data.frame(pca_data)
  pca_data$label <- rep(NA, dim(pca_data)[1])
  
  ## Plot
  p1<-PCx_vs_PCy("PC1", "PC2", pca_data, pca_vars_labs, "Group", type) + ggtitle("Smoking")
  p2<-PCx_vs_PCy("PC3", "PC4", pca_data, pca_vars_labs, "Group", type) + ggtitle("Smoking")
  p3<-PCx_vs_PCy("PC5", "PC6", pca_data, pca_vars_labs, "Group", type) + ggtitle("Smoking")
  
  
  ## PCA plots for nicotine samples separated by Group
  RSE_nicotine<-RSE[,RSE$Expt=="Nicotine"]
  pca<-prcomp(t(assays(RSE_nicotine)$logcounts))
  pca_vars<- getPcaVars(pca)
  pca_vars_labs<- paste0(
      "PC", seq(along = pca_vars), ": ",
      pca_vars, "% Var Expl")
  pca_data<-cbind(pca$x,colData(RSE_nicotine))
  pca_data<-as.data.frame(pca_data)
  pca_data$label <- rep(NA, dim(pca_data)[1])
  
  ## Plot
  p4<-PCx_vs_PCy("PC1", "PC2", pca_data, pca_vars_labs, "Group", type) + ggtitle("Nicotine")
  p5<-PCx_vs_PCy("PC3", "PC4", pca_data, pca_vars_labs, "Group", type) + ggtitle("Nicotine")
  p6<-PCx_vs_PCy("PC5", "PC6", pca_data, pca_vars_labs, "Group", type) + ggtitle("Nicotine")
  
  plot_grid(p1,p2,p3,p4,p5,p6, nrow=2)
  ggsave(paste("plots/03_EDA/03_PCA_MDS/Expt_samples_", type, "_", tissue, "_", age, ".pdf", sep=""), 
           width = 39, height = 18, units = "cm")
  return(NULL)
}

## Plots
PCA_Expt_Group("gene", "brain", "adults")
PCA_Expt_Group("gene", "brain", "pups")
PCA_Expt_Group("exon", "brain", "adults")
PCA_Expt_Group("exon", "brain", "pups")
PCA_Expt_Group("tx", "brain", "adults")
PCA_Expt_Group("tx", "brain", "pups")
PCA_Expt_Group("jx", "brain", "adults")
PCA_Expt_Group("jx", "brain", "pups")




### 1.1.3 Explore separation of samples by sample variable and Group in a PC 

## PC boxplots 
PC_boxplots <- function (PCx, sample_var1, sample_var2, pca_data, pca_vars) {
  
    if(sample_var1=='flowcell'){
      axis_text_size=7
    }
    else{
      axis_text_size=10
    }
  
    plot=ggplot(data=pca_data, 
         aes(x=eval(parse_expr(sample_var1)),y=eval(parse_expr(PCx)))) + 
         ## Hide outliers
         geom_boxplot(outlier.color = "#FFFFFFFF") +
         ## PC dots colored by a group + noise
         geom_jitter(aes(colour=eval(parse_expr(sample_var2))),shape=16, size=2, 
                     position=position_jitter(0.2)) +
         scale_color_manual(values = colors[[sample_var2]]) +
         theme_bw() +
         labs(x=capitalize(sample_var1), y = pca_vars[strtoi(gsub("PC","", PCx))], color=capitalize(sample_var2)) +
         theme(legend.position="right", 
               plot.margin=unit (c (0.1,0.1,0.1,0.1), 'cm'),
               legend.text = element_text(size = 11),
               legend.title = element_text(size = 12),
               axis.title = element_text(size = 12),
               axis.text.x = element_text(size = axis_text_size)) 
         
    return(plot)
}

## Boxplots of PC1 data of samples separated by Group and phenotypes
plot_PC_boxplots<-function(PCx, type, tissue, age) {
  
  pca_data<-PCA(tissue, type, age)[[1]]
  pca_vars<-PCA(tissue, type, age)[[2]]

  if (is.null(age)){
    i=1
    plots<-list()
    for (sample_var1 in c("Pregnancy", "plate", "flowcell")){
        p<-PC_boxplots(PCx, sample_var1, "Group", pca_data, pca_vars)
        plots[[i]]=p
        i=i+1
      }
    plot_grid(plots[[1]], plots[[2]], plots[[3]], nrow = 1)
    ## Save plot
    ggsave(paste("plots/03_EDA/03_PCA_MDS/boxplot_",PCx,"_", type, "_", tissue, ".pdf", 
         sep=""), width = 30, height = 7, units = "cm")
  }

  else if (age=="adults") {
    i=1
    plots<-list()
    for (sample_var1 in c("Expt", "Pregnancy", "plate", "flowcell")){
        p<-PC_boxplots(PCx, sample_var1, "Group", pca_data, pca_vars)
        plots[[i]]=p
        i=i+1
      }
    plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], nrow = 2)
    ## Save plot
    ggsave(paste("plots/03_EDA/03_PCA_MDS/boxplot_",PCx,"_", type, "_", tissue, "_", age, ".pdf", 
         sep=""), width = 24, height = 17, units = "cm")
  }

  else if (age=="pups") {
    i=1
    plots<-list()
    for (sample_var1 in c("Expt", "Sex", "plate", "flowcell")){
        p<-PC_boxplots(PCx, sample_var1, "Group", pca_data, pca_vars)
        plots[[i]]=p
        i=i+1
      }
    plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], nrow = 2)
    ## Save plot
    ggsave(paste("plots/03_EDA/03_PCA_MDS/boxplot_",PCx,"_", type, "_", tissue, "_", age, ".pdf", 
         sep=""), width = 24, height = 17, units = "cm")
  }
}


## Plots
plot_PC_boxplots("PC1", "gene", "blood", NULL)
plot_PC_boxplots("PC1", "gene", "brain", "adults")
plot_PC_boxplots("PC1", "gene", "brain", "pups")
plot_PC_boxplots("PC2", "gene", "brain", "pups")





# ------------------------------------------------------------------------------
##  1.2  Multidimensional scaling (MDS)
# ------------------------------------------------------------------------------

## MDS plot
MDS<- function(sample_var, MDS){
  plot=ggplot(data=as.data.frame(MDS), 
         aes(x=V1,y=V2, color=eval(parse_expr(sample_var)))) + 
         geom_point(size=2) + 
         scale_color_manual(values = colors[[sample_var]]) +
         theme_bw() +
         labs(x= "Component 1", y = "Component 2",  color=capitalize(sample_var)) +
         theme(plot.margin=unit (c (1,1.5,1,1), 'cm'), 
              legend.position = "right",
              legend.text = element_text(size = 11),
              legend.title = element_text(size = 12),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 10)) 
         
    return(plot)
} 

## All MDS plots
plot_MDS<- function(tissue, type, age) {
  
  if (is.null(age)){
    RSE <-eval(parse_expr(paste("rse", type, tissue, "qc", sep="_")))
    ## MDS in 2 components
    MDS <-calculateMDS(assays(RSE)$logcounts, ncomponents=2)
    MDS <- cbind(MDS, colData(RSE)[49:61])
    
    i=1
    plots<-list()
    for (pheno_var in c("Group", "Pregnancy", "plate", "flowcell")){
        p<-MDS(pheno_var, MDS)
        plots[[i]]=p
        i=i+1
      }
    plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], nrow = 2, align='v')
    ## Save plot
    ggsave(paste("plots/03_EDA/03_PCA_MDS/MDS_", type, "_", tissue, ".pdf", 
           sep=""), width = 26, height = 16, units = "cm")
  }
  
  else if (age=="adults") {
    RSE <-eval(parse_expr(paste("rse", type, tissue, age, "qc", sep="_")))
    MDS <-calculateMDS(assays(RSE)$logcounts, ncomponents=2)
    MDS <- cbind(MDS, colData(RSE)[49:61])
    
    i=1
    plots<-list()
    for (pheno_var in c("Expt", "Group", "Pregnancy", "plate", "flowcell")){
        p<-MDS(pheno_var, MDS)
        plots[[i]]=p
        i=i+1
      }
    plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], nrow = 2, align='v')
    ## Save plot
    ggsave(paste("plots/03_EDA/03_PCA_MDS/MDS_", type, "_", tissue, "_", age, ".pdf", 
           sep=""), width = 39, height = 16, units = "cm")
  }
  
  else if (age=="pups") {
    RSE <-eval(parse_expr(paste("rse", type, tissue, age, "qc", sep="_")))
    MDS <-calculateMDS(assays(RSE)$logcounts, ncomponents=2)
    MDS <- cbind(MDS, colData(RSE)[49:61])
    
    i=1
    plots<-list()
    for (pheno_var in c("Expt", "Group", "Sex", "plate", "flowcell")){
        p<-MDS(pheno_var, MDS)
        plots[[i]]=p
        i=i+1
      }
    plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], nrow = 2, align='v')
    ## Save plot
    ggsave(paste("plots/03_EDA/03_PCA_MDS/MDS_", type, "_", tissue, "_", age, ".pdf", 
           sep=""), width = 39, height = 16, units = "cm")
  }
}

## Plots
plot_MDS("blood", "gene", NULL)
plot_MDS("brain", "gene", "adults")
plot_MDS("brain", "gene", "pups")
plot_MDS("brain", "exon", "adults")
plot_MDS("brain", "exon", "pups")
plot_MDS("brain", "tx", "adults")
plot_MDS("brain", "tx", "pups")
plot_MDS("brain", "jx", "adults")
plot_MDS("brain", "jx", "pups")







## Reproducibility information
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
