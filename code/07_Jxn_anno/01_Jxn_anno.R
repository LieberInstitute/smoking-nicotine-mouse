
## 1. Junction annotation 

library(ggplot2)
library(here)
library(GenomicRanges)

## Only for DE jxns from pup and brain samples

load(here("raw-data/rse_jx_smoking_mouse_n208.Rdata"))
load(here("processed-data/04_DEA/Jx_analysis/top_jxns_nic.Rdata"))
load(here("processed-data/04_DEA/Jx_analysis/de_jxns_nic.Rdata"))
load(here("processed-data/04_DEA/Jx_analysis/top_jxns_smo.Rdata"))
load(here("processed-data/04_DEA/Jx_analysis/de_jxns_smo.Rdata"))



## Explore jxn classes

## "Novel" jxns have unknown (not in GENCODE) start and end sites
table(rowData(rse_jx)[which(rowData(rse_jx)$Class=="Novel"), c("inGencodeStart", "inGencodeEnd")])
#                 inGencodeEnd
#   inGencodeStart FALSE
#            FALSE 847606

## Note: jxns without associated gene are all Novel
table(rowData(rse_jx)[which(is.na(rowData(rse_jx)$newGeneID)), "Class"])
# Novel 
# 700460 

## "AltStartEnd" jxns have only one known site 
table(rowData(rse_jx)[which(rowData(rse_jx)$Class=="AltStartEnd"), c("inGencodeStart", "inGencodeEnd")])
#                   inGencodeEnd
#   inGencodeStart  FALSE   TRUE
#       FALSE        0    579436
#       TRUE        8342     0

## "ExonSkip" jxns have sites from non-successive exons, both known individually but not together 
table(rowData(rse_jx)[which(rowData(rse_jx)$Class=="ExonSkip"), c("inGencodeStart", "inGencodeEnd")])
#                  inGencodeEnd
#   inGencodeStart TRUE
#             TRUE  142

## "InGen" jxns are in GENCODE (both sites are known)
table(rowData(rse_jx)[which(rowData(rse_jx)$Class=="InGen"), c("inGencodeStart", "inGencodeEnd")])
#                  inGencodeEnd
#   inGencodeStart TRUE
#             TRUE  542

## "Fusion" jxns span  multiple genes (start and end in jxns that can be known or unknown)
table(rowData(rse_jx)[which(rowData(rse_jx)$isFusion=="TRUE"), c("inGencodeStart", "inGencodeEnd")])
#                 inGencodeEnd
#   inGencodeStart FALSE TRUE
#            FALSE   4   74
#             TRUE   11   15



## Create table with the following info for each jxn class:
##        - Class: jxn class
##        - number: total number of jxns from that class
##        - number_DEjxns_nic: number of nicotine DE jxns from that class
##        - percentage_DEjxns_nic: % of nicotine DE jxns that were from that class
##        - number_DEjxns_smo: number of smoking DE jxns from that class
##        - percentage_DEjxns_smo: % of smoking DE jxns that were from that class

jxn_classes <- lapply(unique(rowData(rse_jx)$Class), function(x){cbind(Class=x, 
                                                                       number_jxs=length(which(rowData(rse_jx)$Class==x)), 
                                                                       number_DEjxns_nic=length(which(de_jxns_nic$Class==x)),
                                                                       percentage_DEjxns_nic=signif(length(which(de_jxns_nic$Class==x))/dim(de_jxns_nic)[1]*100, 3),
                                                                       number_DEjxns_smo=length(which(de_jxns_smo$Class==x)),
                                                                       percentage_DEjxns_smo=signif(length(which(de_jxns_smo$Class==x))/dim(de_jxns_smo)[1]*100, 3))})
jxn_classes <- rbind(jxn_classes[[1]], jxn_classes[[2]], jxn_classes[[3]], jxn_classes[[4]])

## Add "isFusion" class information 
jxn_classes <- as.data.frame(rbind(jxn_classes, cbind("isFusion", 
                                                       length(which(rowData(rse_jx)$isFusion=="TRUE")), 
                                                       length(which(de_jxns_nic$isFusion=="TRUE")),
                                                       signif(length(which(de_jxns_nic$isFusion=="TRUE"))/dim(de_jxns_nic)[1]*100, 3),
                                                       length(which(de_jxns_smo$isFusion=="TRUE")),
                                                       signif(length(which(de_jxns_smo$isFusion=="TRUE"))/dim(de_jxns_smo)[1]*100, 3))))

jxn_classes$number_jxs <- as.numeric(jxn_classes$number_jxs)
jxn_classes$number_DEjxns_nic <- as.numeric(jxn_classes$number_DEjxns_nic)
jxn_classes$percentage_DEjxns_nic<- as.numeric(jxn_classes$percentage_DEjxns_nic)
jxn_classes$number_DEjxns_smo <- as.numeric(jxn_classes$number_DEjxns_smo)
jxn_classes$percentage_DEjxns_smo<- as.numeric(jxn_classes$percentage_DEjxns_smo)

save(jxn_classes, file="processed-data/07_Jxn_anno/jxn_classes.Rdata")



## Create table with the information of DE jxns' genes (the known ones):
##        - Gene: ID of gene with DE jxns
##        - number_DEjxns_nic: number of nicotine DE jxns the gene has
##        - number_DEjxns_nic_Novel: number of nicotine DE novel jxns the gene has
##        - number_DEjxns_nic_AltStartEnd: number of nicotine DE jxns with alternative start/end, the gene has
##        - number_DEjxns_nic_InGen: number of nicotine DE known (in GENCODE) jxns the gene has
##        - number_DEjxns_nic_ExonSkip: number of nicotine DE jxns from non-successive exons, the gene has
##        - number_DEjxns_nic_isFusion: number of nicotine DE jxns that span multiple genes, the gene has
##        - number_DEjxns_smo: number of smoking DE jxns the gene has
##        - number_DEjxns_smo_Novel: number of smoking DE novel jxns the gene has
##        - number_DEjxns_smo_AltStartEnd: number of smoking DE jxns with alternative start/end, the gene has
##        - number_DEjxns_smo_InGen: number of smoking DE known (in GENCODE) jxns the gene has
##        - number_DEjxns_smo_ExonSkip: number of smoking DE jxns from non-successive exons, the gene has
##        - number_DEjxns_smo_isFusion: number of smoking DE jxns that span multiple genes, the gene has

## Genes with DE jxns
DEjxns_genes <- union(unique(de_jxns_nic$newGeneID), unique(de_jxns_smo$newGeneID))
DEjxns_genes <- DEjxns_genes[which(! is.na(DEjxns_genes))]
## Info
DEjxns_genes_info <- vector()
for (gene in DEjxns_genes){
          DEjxns_genes_info <- rbind(DEjxns_genes_info, 
                                    c(Gene=gene,
                                      ## newGeneID: gene name(s) associated with the exons that each junction spans
                                      number_DEjxns_nic=length(which(de_jxns_nic$newGeneID==gene)),
                                      number_DEjxns_nic_Novel=length(which(de_jxns_nic$newGeneID==gene & de_jxns_nic$Class=="Novel")),
                                      number_DEjxns_nic_AltStartEnd=length(which(de_jxns_nic$newGeneID==gene & de_jxns_nic$Class=="AltStartEnd")),
                                      number_DEjxns_nic_InGen=length(which(de_jxns_nic$newGeneID==gene & de_jxns_nic$Class=="InGen")),
                                      number_DEjxns_nic_ExonSkip=length(which(de_jxns_nic$newGeneID==gene & de_jxns_nic$Class=="ExonSkip")),
                                      number_DEjxns_nic_isFusion=length(which(de_jxns_nic$newGeneID==gene & de_jxns_nic$isFusion==TRUE)),
                                      
                                      number_DEjxns_smo=length(which(de_jxns_smo$newGeneID==gene)),
                                      number_DEjxns_smo_Novel=length(which(de_jxns_smo$newGeneID==gene & de_jxns_smo$Class=="Novel")),
                                      number_DEjxns_smo_AltStartEnd=length(which(de_jxns_smo$newGeneID==gene & de_jxns_smo$Class=="AltStartEnd")),
                                      number_DEjxns_smo_InGen=length(which(de_jxns_smo$newGeneID==gene & de_jxns_smo$Class=="InGen")),
                                      number_DEjxns_smo_ExonSkip=length(which(de_jxns_smo$newGeneID==gene & de_jxns_smo$Class=="ExonSkip")),
                                      number_DEjxns_smo_isFusion=length(which(de_jxns_smo$newGeneID==gene & de_jxns_smo$isFusion==TRUE))))
}
## Char to numeric
DEjxns_genes_info <- cbind(Gene=DEjxns_genes_info[,1], as.data.frame(apply(DEjxns_genes_info[,2:13], 2, function(x){as.numeric(as.character(x))})))
save(DEjxns_genes_info, file="processed-data/07_Jxn_anno/DEjxns_genes_info.Rdata")



## Histograms of the number of total, novel, AltStartEnd, InGen, ExonSkip and Fusion DE jxns of the genes in nic and smo

## Data frame
nic=cbind(DEjxns_genes_info[,2:7], expt=rep("Nicotine", dim(DEjxns_genes_info)[1]))
total_DEjxns <- as.data.frame(rbind(nic, 
                                    setNames(cbind(DEjxns_genes_info[,8:13], rep("Smoking", dim(DEjxns_genes_info)[1])), colnames(nic))))
colnames(total_DEjxns) <- c("number_DEjxns", "number_DEjxns_Novel", "number_DEjxns_AltStartEnd", 
                            "number_DEjxns_InGen", "number_DEjxns_ExonSkip", "number_DEjxns_isFusion", "expt")
total_DEjxns_expt <- as.data.frame(apply(total_DEjxns[,1:6], 2, function(x){as.numeric(x)}))
total_DEjxns_expt$expt <-total_DEjxns$expt

## Histograms ignoring zeros

data=total_DEjxns_expt[which(total_DEjxns_expt$number_DEjxns!=0),]
h1 <- ggplot(data, aes(x=number_DEjxns, fill=expt)) +
      geom_histogram(color="black", alpha=0.9, position="dodge") +
      scale_x_continuous(breaks = seq(1, max(data$number_DEjxns), 10)) +
      xlab("Number of DE jxns per gene") +
      ylab("Frecuency") +
      theme(axis.title=element_text(size=10,face="bold"), legend.position = "None")+
      facet_wrap(~expt)

data=total_DEjxns_expt[which(total_DEjxns_expt$number_DEjxns_Novel!=0),]
h2 <- ggplot(data, aes(x=number_DEjxns_Novel, fill=expt)) +
      geom_histogram(color="black", alpha=0.9, position="dodge")+
      xlab("Number of Novel DE jxns each gene has")+
      ylab("Frecuency") +
      scale_x_continuous(breaks=seq(1,max(data$number_DEjxns_Novel),1)) + 
      scale_y_continuous(breaks=seq(1, max(table(data$number_DEjxns_Novel)),1)) +
      theme(axis.title=element_text(size=10,face="bold"), legend.position = "None", plot.margin = margin(0.2, 5, 0.2, 5, "cm"))+
      facet_wrap(~expt)

data=total_DEjxns_expt[which(total_DEjxns_expt$number_DEjxns_AltStartEnd!=0),]
h3 <- ggplot(data, aes(x=number_DEjxns_AltStartEnd, color=expt, fill=expt)) +
      geom_histogram(color="black", alpha=0.9, position="dodge")+
      scale_x_continuous(breaks = seq(1, max(data$number_DEjxns_AltStartEnd), 10)) +    
      xlab("Number of DE jxns with alternative start/end, each gene has")+
      ylab("Frecuency")+
      theme(axis.title=element_text(size=10,face="bold"), legend.position = "None")+
      facet_wrap(~expt)

data=total_DEjxns_expt[which(total_DEjxns_expt$number_DEjxns_InGen!=0),]
h4 <- ggplot(data, aes(x=number_DEjxns_InGen, color=expt, fill=expt)) +
      geom_histogram(color="black", alpha=0.9, position="dodge")+
      xlab("Number of DE known (in GENCODE) jxns each gene has")+
      ylab("Frecuency")+
      scale_x_continuous(breaks=seq(1,max(data$number_DEjxns_InGen),1)) + 
      theme(axis.title=element_text(size=10,face="bold"), legend.position = "None", plot.margin = margin(0.2, 6.5, 0.2, 6.5, "cm"))+
      facet_wrap(~expt)

data=total_DEjxns_expt[which(total_DEjxns_expt$number_DEjxns_ExonSkip!=0),]
h5 <- ggplot(data, aes(x=number_DEjxns_ExonSkip, color=expt, fill=expt)) +
      geom_histogram(color="black", alpha=0.9, position="dodge")+
      xlab("Number of DE jxns from non-successive exons, each gene has")+
      ylab("Frecuency")+
      scale_x_continuous(breaks=seq(1,max(data$number_DEjxns_ExonSkip),1)) + 
      theme(axis.title=element_text(size=10,face="bold"), legend.position = "None", plot.margin = margin(0.2, 5, 0.2, 5, "cm"))+
      facet_wrap(~expt)

## There were not DE fusion jxns 
length(which(total_DEjxns_expt$number_DEjxns_isFusion!=0))
## [1] 0


plot_grid(h1, h2, h3, h4, h5, ncol=3)
ggsave(here("plots/07_Jxn_anno/histograms_numberDEjxns.pdf"), width = 50, height = 25, units = "cm")





## 1. Find the closest upstream gene of DE Novel jxns (without gene)

## Obtain novel DE introns/jxns without associated gene

## Nicotine
novel_jxns_nic <- de_jxns_nic[which(! (de_jxns_nic$Class=="InGen")),]
## Smoking
novel_jxns_smo <- de_jxns_smo[which(! (de_jxns_smo$Class=="InGen")),]







## 3. Find nearest genes of novel jxns without assigned gene






