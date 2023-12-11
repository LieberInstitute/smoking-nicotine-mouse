
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

## Colors to highlight poor QC plots
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
# version  R version 4.2.0 (2022-04-22 ucrt)
# os       Windows 10 x64 (build 19044)
# system   x86_64, mingw32
# ui       RStudio
# language (EN)
# collate  Spanish_Mexico.utf8
# ctype    Spanish_Mexico.utf8
# tz       America/Mexico_City
# date     2022-07-19

