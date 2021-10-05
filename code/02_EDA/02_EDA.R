#with R 4.0.x

## libraries
library("here")
library("sessioninfo")
library("ggplot2")
library("recount")
library("jaffelab")
library("stats")
library("ggfortify")
library("scater")
library("limma")
library("edgeR")

# load data
load(here("processed-data","build_objects","rse_gene.Rdata"), verbose = TRUE)
load(here("processed-data","build_objects","rse_exon.Rdata"), verbose = TRUE)
load(here("processed-data","build_objects","rse_jx.Rdata"), verbose = TRUE)
load(here("processed-data","build_objects","rse_tx.Rdata"), verbose = TRUE)

# Obtain PCa
geneExprs <- log2(recount::getRPKM(rse_gene, "Length") + 1)
set.seed(20201006)
pca <- stats::prcomp(t(geneExprs))
pca_vars <- jaffelab::getPcaVars(pca)
pca_vars_lab <- paste0(
  "PC", seq(along = pca_vars), ": ",
  pca_vars, "% Var Expl"
)

# Create data frame
coldata <- as.data.frame(colData(rse_gene))
pca_data <- cbind(coldata, pca$x)


# Explore gene expression variability: PC vs PC
## Age
### PC1 vs PC2
pdf(here("plots", "02_EDA", "pc-Age-plot.pdf"))
ggplot(data =pca_data, aes(x=PC1,y=PC2,color= Age.x)) + 
  geom_point() + 
  ggtitle("PCA of the smoking mouse data, colored by age") +
  labs(x= pca_vars_lab[1], y = pca_vars_lab[2])

### PC1 vs PC3
ggplot(data =pca_data, aes(x=PC1,y=PC3,color= Age.x)) + 
  geom_point() + 
  ggtitle("PCA of the smoking mouse data, colored by age") +
  labs(x= pca_vars_lab[1], y = pca_vars_lab[3])

### PC2 vs PC3
ggplot(data =pca_data, aes(x=PC2,y=PC3,color= Age.x)) + 
  geom_point() +
  ggtitle("PCA of the smoking mouse data, colored by age") +
  labs(x= pca_vars_lab[2], y = pca_vars_lab[3])
dev.off()


## Exposition
### PC1 vs PC2
pdf(here("plots", "02_EDA", "pc-exposition-plot.pdf"))
ggplot(data =pca_data, aes(x=PC1,y=PC2,color=Expt.x)) + 
  geom_point(aes(shape=Age.x))  +
  ggtitle("PCA of the smoking mouse data, colored by exposition") +
  labs(x= pca_vars_lab[1], y = pca_vars_lab[2])

### PC1 vs PC3
ggplot(data =pca_data, aes(x=PC1,y=PC3,color=Expt.x)) + 
  geom_point(aes(shape=Age.x))  +
  ggtitle("PCA of the smoking mouse data, colored by exposition") +
  labs(x= pca_vars_lab[1], y = pca_vars_lab[3])

### PC2 vs PC3
ggplot(data =pca_data, aes(x=PC2,y=PC3,color=Expt.x)) + 
  geom_point(aes(shape=Age.x))  +
  ggtitle("PCA of the smoking mouse data, colored by exposition") +
  labs(x= pca_vars_lab[2], y = pca_vars_lab[3])
dev.off()


##Control vs cases
### PC1 vs PC2
pdf(here("plots", "02_EDA", "pc-casecontrol-plot.pdf"))
ggplot(data =pca_data, aes(x=PC1,y=PC2,color=Group.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by group") +
  labs(x= pca_vars_lab[1], y = pca_vars_lab[2])

### PC1 vs PC3
ggplot(data =pca_data, aes(x=PC1,y=PC3,color=Group.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by group") +
  labs(x= pca_vars_lab[1], y = pca_vars_lab[3])

### PC2 vs PC3
ggplot(data =pca_data, aes(x=PC2,y=PC3,color=Group.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by group") +
  labs(x= pca_vars_lab[2], y = pca_vars_lab[3])
dev.off()


## Location
### PC1 vs PC2
pdf(here("plots", "02_EDA", "pc-location-plot.pdf"))
ggplot(data =pca_data, aes(x=PC1,y=PC2,color=location.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by location") +
  labs(x= pca_vars_lab[1], y = pca_vars_lab[2])

### PC1 vs PC3
ggplot(data =pca_data, aes(x=PC1,y=PC3,color=location.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by location") +
  labs(x= pca_vars_lab[1], y = pca_vars_lab[3])

### PC2 vs PC3
ggplot(data =pca_data, aes(x=PC2,y=PC3,color=location.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by location") +
  labs(x= pca_vars_lab[2], y = pca_vars_lab[3])
dev.off()


## Plate
### PC1 vs PC2
pdf(here("plots", "02_EDA", "pc-plate-plot.pdf"))
ggplot(data =pca_data, aes(x=PC1,y=PC2,color=plate.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by plate") +
  labs(x= pca_vars_lab[1], y = pca_vars_lab[2])

### PC1 vs PC3
ggplot(data =pca_data, aes(x=PC1,y=PC3,color=plate.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by plate") +
  labs(x= pca_vars_lab[1], y = pca_vars_lab[3])

### PC2 vs PC3
ggplot(data =pca_data, aes(x=PC2,y=PC3,color=plate.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by plate") +
  labs(x= pca_vars_lab[2], y = pca_vars_lab[3])
dev.off()


## Tissue
### PC1 vs PC2
pdf(here("plots", "02_EDA", "pc-tissue-plot.pdf"))
ggplot(data =pca_data, aes(x=PC1,y=PC2,color=Tissue.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(x= pca_vars_lab[1], y = pca_vars_lab[2])

### PC1 vs PC3
ggplot(data =pca_data, aes(x=PC1,y=PC3,color=Tissue.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(x= pca_vars_lab[1], y = pca_vars_lab[3])

### PC2 vs PC3
ggplot(data =pca_data, aes(x=PC2,y=PC3,color=Tissue.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(x= pca_vars_lab[2], y = pca_vars_lab[3])
dev.off()


## PCs vs SPEAQeasy output summary metrics
pdf(here("plots", "02_EDA", "pc-metricsSPEAQeasy-plot.pdf"))
### PC1
ggplot(data =pca_data, aes(x=PC1,y=mitoMapped,color=Tissue.x)) + 
  geom_point() +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(x= pca_vars_lab[1])

ggplot(data =pca_data, aes(x=PC1, y=mitoRate,color=Tissue.x)) + 
  geom_point() +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(x= pca_vars_lab[1])

ggplot(data =pca_data, aes(x=PC1, y=totalAssignedGene,color=Tissue.x)) + 
  geom_point() +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(x= pca_vars_lab[1])

ggplot(data =pca_data, aes(x=PC1, y=rRNA_rate, color=Tissue.x)) + 
  geom_point() +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(x= pca_vars_lab[1])

ggplot(data =pca_data, aes(x=PC1, y=ERCCsumLogErr, color=Tissue.x)) + 
  geom_point() +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(x= pca_vars_lab[1])

### PC2
ggplot(data =pca_data, aes(x=PC2,y=mitoMapped,color=Tissue.x)) + 
  geom_point() +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(x= pca_vars_lab[2])

ggplot(data =pca_data, aes(x=PC2, y=mitoRate,color=Tissue.x)) + 
  geom_point() +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(x= pca_vars_lab[2])

ggplot(data =pca_data, aes(x=PC2, y=totalAssignedGene,color=Tissue.x)) + 
  geom_point() +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(x= pca_vars_lab[2])

ggplot(data =pca_data, aes(x=PC2, y=rRNA_rate, color=Tissue.x)) + 
  geom_point() +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(x= pca_vars_lab[2])

ggplot(data =pca_data, aes(x=PC2, y=ERCCsumLogErr, color=Tissue.x)) + 
  geom_point() +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(x= pca_vars_lab[2])


### PC3
ggplot(data =pca_data, aes(x=PC3,y=mitoMapped,color=Tissue.x)) + 
  geom_point() +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(x= pca_vars_lab[3])

ggplot(data =pca_data, aes(x=PC3, y=mitoRate,color=Tissue.x)) + 
  geom_point() +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(x= pca_vars_lab[3])

ggplot(data =pca_data, aes(x=PC3, y=totalAssignedGene,color=Tissue.x)) + 
  geom_point() +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(x= pca_vars_lab[3])

ggplot(data =pca_data, aes(x=PC3, y=rRNA_rate, color=Tissue.x)) + 
  geom_point() +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(x= pca_vars_lab[3])

ggplot(data =pca_data, aes(x=PC3, y=ERCCsumLogErr, color=Tissue.x)) + 
  geom_point() +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(x= pca_vars_lab[3])
dev.off()



# Explore gene expression variability: PCs boxplots

pc_list = pca_data[, which(colnames(pca_data)=="PC1" | colnames(pca_data)=="PC2"
                           | colnames(pca_data)=="PC3" | colnames(pca_data)=="PC4"
                           | colnames(pca_data)=="PC5") ]
var_list = pca_data[, which(colnames(pca_data)=="Tissue.x" | colnames(pca_data)=="Age.x"
                           | colnames(pca_data)=="Sex.x" | colnames(pca_data)=="Expt.x"
                           | colnames(pca_data)=="Group.x" | colnames(pca_data)=="plate.x"
                           | colnames(pca_data)=="location.x")]

colnum_var= length(colnames(var_list))
colnum_pc= length(colnames(pc_list))

plot_list=list()
n=1
 for (j in 1:colnum_pc){
   for (i in 1:colnum_var){
 # j=2
 # i=1
    
    pca_subdata= cbind(pc_list[j], var_list[i], var_list$Group.x)
    colnames(pca_subdata)=c("PC", "var", "group")
    
    plot_list[[n]]= ggplot(data =pca_subdata, aes(x=var, y=PC)) +
      geom_boxplot(colour = "grey", fill = "light grey") +
      geom_point(aes(color=group), position = "jitter") +
      ggtitle("Boxplot, PCA of the smoking mouse data") +
      labs(x=colnames(var_list[i]) ,y = pca_vars_lab[j])
    n=n+1
   }
 }
pdf(here("plots", "02_EDA", "boxplots.pdf"))
lapply(plot_list, print)
dev.off()



# Explore gene expression variability: model for each gene

## Diagnostic plots for quality control
mito= subsets=list(Mito=grep("Mt_", rowData(rse_gene)$gene_type))
rse_gene <- addPerCellQC(rse_gene, subsets=mito)

### detected genes VS total count. 
pdf(here("plots", "02_EDA", "diagnostic_qual_plots"))
ggplot(data =coldata, aes(x=sum, y=detected, color= Tissue.x))+ 
  geom_point() 

### Mitochondrial content VS total count. 
ggplot(data =coldata, aes(x=sum, y=subsets_Mito_percent, color= Tissue.x)) + 
  geom_point() 

### Explanatory Variables
rse_gene <- logNormCounts(rse_gene)
vars <- getVarianceExplained(rse_gene, 
                             variables=c("Tissue.x", "Age.x", "Sex.x", "Age.x", 
                                         "Expt.x", "Group.x", "plate.x", "location.x"))
plotExplanatoryVariables(vars)
dev.off()



## Applying variancePartition 
## UNFINISHED code taken from http://bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/variancePartition.pdf

### identify genes that pass expression cutoff
geneCounts <- assays(rse_gene)$counts
isexpr <- rowSums(cpm(geneCounts)>1) >= 0.5 * ncol(geneCounts)
### create data structure with only expressed genes
gExpr <- DGEList(counts=geneCounts[isexpr,])
### Perform TMM normalization
gExpr <- calcNormFactors(gExpr)
### Specify variables to be included in the voom() estimates of uncertainty.
design <- model.matrix(~ Tissue.x + Age.x, data=coldata)
### Estimate precision weights for each gene and sample
vobjGenes <- voom(gExpr, design)




# Explore covariates
## UNFINISHED code taken from https://github.com/LieberInstitute/brainseq_phase2/blob/master/expr_cutoff/explore_metrics.R

pdf(here("plots", "02_EDA", "explore_metrics.pdf"))

ggplot(data =coldata, aes(x=colData(rse_gene)$mitoRate, y=colData(rse_gene)$rRNA_rate, color= Tissue.x)) + 
  geom_point() 
ggplot(data =coldata, aes(x=colData(rse_gene)$mitoRate, y=colData(rse_gene)$ERCCsumLogErr, color= Tissue.x)) + 
  geom_point() 
ggplot(data =coldata, aes(x=colData(rse_gene)$rRNA_rate, y=colData(rse_gene)$ERCCsumLogErr, color= Tissue.x)) + 
  geom_point() 

dev.off()
