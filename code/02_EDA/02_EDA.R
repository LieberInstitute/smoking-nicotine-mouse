#with R 4.0.x

## libraries
library("here")
library("sessioninfo")
library("ggplot2")
library("recount")
library("jaffelab")
library("stats")
library("ggfortify")

#load data
load(here("processed-data","build_objects","rse_gene.Rdata"), verbose = TRUE)
load(here("processed-data","build_objects","rse_exon.Rdata"), verbose = TRUE)
load(here("processed-data","build_objects","rse_jx.Rdata"), verbose = TRUE)
load(here("processed-data","build_objects","rse_tx.Rdata"), verbose = TRUE)

#Obtain PC
geneExprs <- log2(recount::getRPKM(rse_gene, "Length") + 1) #check why add 1
set.seed(20201006)
pca <- stats::prcomp(t(geneExprs))
pca_vars <- jaffelab::getPcaVars(pca)
pca_var_expl <- paste0(
  "PC", seq(along = pca_vars), ": ",
  pca_vars, "% Var Expl"
)

#Create data frame
coldata <- as.data.frame(colData(rse_gene))
pca_data <- cbind(coldata, pca$x)

#Plot PCs

#Age
#PC1 vs PC2
pdf(here("plots", "02_EDA", "pc-Age-plot.pdf"))
ggplot(data =pca_data, aes(x=PC1,y=PC2,color= Age.x)) + 
  geom_point() + 
  ggtitle("PCA of the smoking mouse data, colored by age") +
  labs(y= pca_vars_lab[1], x = pca_vars_lab[2])

#PC1 vs PC3
ggplot(data =pca_data, aes(x=PC1,y=PC3,color= Age.x)) + 
  geom_point() + 
  ggtitle("PCA of the smoking mouse data, colored by age") +
  labs(y= pca_vars_lab[1], x = pca_vars_lab[3])

#PC2 vs PC3
ggplot(data =pca_data, aes(x=PC2,y=PC3,color= Age.x)) + 
  geom_point() +
  ggtitle("PCA of the smoking mouse data, colored by age") +
  labs(y= pca_vars_lab[2], x = pca_vars_lab[3])
dev.off()


#Exposition
#PC1 vs PC2
pdf(here("plots", "02_EDA", "pc-exposition-plot.pdf"))
ggplot(data =pca_data, aes(x=PC1,y=PC2,color=Expt.x)) + 
  geom_point(aes(shape=Age.x))  +
  ggtitle("PCA of the smoking mouse data, colored by exposition") +
  labs(y= pca_vars_lab[1], x = pca_vars_lab[2])

#PC1 vs PC3
ggplot(data =pca_data, aes(x=PC1,y=PC3,color=Expt.x)) + 
  geom_point(aes(shape=Age.x))  +
  ggtitle("PCA of the smoking mouse data, colored by exposition") +
  labs(y= pca_vars_lab[1], x = pca_vars_lab[3])

#PC2 vs PC3
ggplot(data =pca_data, aes(x=PC2,y=PC3,color=Expt.x)) + 
  geom_point(aes(shape=Age.x))  +
  ggtitle("PCA of the smoking mouse data, colored by exposition") +
  labs(y= pca_vars_lab[2], x = pca_vars_lab[3])
dev.off()


#Control vs cases
#PC1 vs PC2
pdf(here("plots", "02_EDA", "pc-casecontrol-plot.pdf"))
ggplot(data =pca_data, aes(x=PC1,y=PC2,color=Group.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by group") +
  labs(y= pca_vars_lab[1], x = pca_vars_lab[2])

#PC1 vs PC3
ggplot(data =pca_data, aes(x=PC1,y=PC3,color=Group.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by group") +
  labs(y= pca_vars_lab[1], x = pca_vars_lab[3])

#PC2 vs PC3
ggplot(data =pca_data, aes(x=PC2,y=PC3,color=Group.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by group") +
  labs(y= pca_vars_lab[2], x = pca_vars_lab[3])
dev.off()


#Location
#PC1 vs PC2
pdf(here("plots", "02_EDA", "pc-location-plot.pdf"))
ggplot(data =pca_data, aes(x=PC1,y=PC2,color=location.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by location") +
  labs(y= pca_vars_lab[1], x = pca_vars_lab[2])

#PC1 vs PC3
ggplot(data =pca_data, aes(x=PC1,y=PC3,color=location.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by location") +
  labs(y= pca_vars_lab[1], x = pca_vars_lab[3])

#PC2 vs PC3
ggplot(data =pca_data, aes(x=PC2,y=PC3,color=location.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by location") +
  labs(y= pca_vars_lab[2], x = pca_vars_lab[3])
dev.off()


#Plate
#PC1 vs PC2
pdf(here("plots", "02_EDA", "pc-plate-plot.pdf"))
ggplot(data =pca_data, aes(x=PC1,y=PC2,color=plate.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by plate") +
  labs(y= pca_vars_lab[1], x = pca_vars_lab[2])

#PC1 vs PC3
ggplot(data =pca_data, aes(x=PC1,y=PC3,color=plate.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by plate") +
  labs(y= pca_vars_lab[1], x = pca_vars_lab[3])

#PC2 vs PC3
ggplot(data =pca_data, aes(x=PC2,y=PC3,color=plate.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by plate") +
  labs(y= pca_vars_lab[2], x = pca_vars_lab[3])
dev.off()


#Tissue
#PC1 vs PC2
pdf(here("plots", "02_EDA", "pc-tissue-plot.pdf"))
ggplot(data =pca_data, aes(x=PC1,y=PC2,color=Tissue.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(y= pca_vars_lab[1], x = pca_vars_lab[2])

#PC1 vs PC3
ggplot(data =pca_data, aes(x=PC1,y=PC3,color=Tissue.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(y= pca_vars_lab[1], x = pca_vars_lab[3])

#PC2 vs PC3
ggplot(data =pca_data, aes(x=PC2,y=PC3,color=Tissue.x)) + 
  geom_point(aes(shape=Age.x)) +
  ggtitle("PCA of the smoking mouse data, colored by tissue") +
  labs(y= pca_vars_lab[2], x = pca_vars_lab[3])
dev.off()


#Boxplots

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
      labs(y = pca_var_expl[j])
    n=n+1
   }
   # print("done")
 }
pdf(here("plots", "02_EDA", "boxplots.pdf"))
lapply(plot_list, print)
dev.off()

