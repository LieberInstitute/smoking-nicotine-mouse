
## 1. Junction annotation 

library(ggplot2)
library(here)
library(SummarizedExperiment)
library(GenomicRanges)
library(rtracklayer)
library(AcidGenomes)

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





## 1. Find the nearest, preceding and following genes of DE Novel jxns (without associated gene)

## Obtain novel DE jxns without associated gene
## Nicotine
novel_jxns_nic <- de_jxns_nic[which(de_jxns_nic$Class=="Novel" & is.na(de_jxns_nic$newGeneID)),]
## Smoking
novel_jxns_smo <- de_jxns_smo[which(de_jxns_smo$Class=="Novel" & is.na(de_jxns_smo$newGeneID)),]

## GRanges with novel nic and smo DE jxns = querys
GRanges_novel_jxns_nic <- rowRanges(rse_jx_brain_pups_nicotine)[rownames(novel_jxns_nic)]
GRanges_novel_jxns_smo <- rowRanges(rse_jx_brain_pups_smoking)[rownames(novel_jxns_smo)]
## Same level style as subject to enable comparison between them
seqlevelsStyle(GRanges_novel_jxns_nic) <- "NCBI"
seqlevelsStyle(GRanges_novel_jxns_smo) <- "NCBI"

## GRanges of all mouse genes = subject (GRCm38 genome assembly)
GRanges_mouse_genes <- makeGRangesFromEnsembl(
                              organism="Mus musculus",
                              level = "genes",
                              genomeBuild = "GRCm38",
                              release = NULL)

############## Nearest genes ##############
nearest_genes_nic <- GRanges_mouse_genes[IRanges::nearest(GRanges_novel_jxns_nic, GRanges_mouse_genes)]
nearest_genes_smo <- GRanges_mouse_genes[IRanges::nearest(GRanges_novel_jxns_smo, GRanges_mouse_genes)]


############## Following genes ##############
## Note that precede() returns the preceded gene by the jxn, i.e. the following gene in + strands
following_genes_nic <- GRanges_mouse_genes[IRanges::precede(GRanges_novel_jxns_nic, GRanges_mouse_genes)]
following_genes_smo <- GRanges_mouse_genes[IRanges::precede(GRanges_novel_jxns_smo, GRanges_mouse_genes)]


############## Preceding genes ##############
## Note that follow() returns the gene directly followed by the jxn, i.e. the preceding gene in + strands
preceding_genes_nic <- GRanges_mouse_genes[IRanges::follow(GRanges_novel_jxns_nic, GRanges_mouse_genes)]
preceding_genes_smo <- GRanges_mouse_genes[IRanges::follow(GRanges_novel_jxns_smo, GRanges_mouse_genes)]


## List of all GRanges for nic and smo DE Novel jxns' genes

novel_jxns_foundGenes <- list (
  DE_Novel_jxns_nic = GRanges_novel_jxns_nic,
  nearest_genes_nic = nearest_genes_nic,
  following_genes_nic = following_genes_nic, 
  preceding_genes_nic = preceding_genes_nic,
  DE_Novel_jxns_smo = GRanges_novel_jxns_smo,
  nearest_genes_smo = nearest_genes_smo,
  following_genes_smo = following_genes_smo, 
  preceding_genes_smo = preceding_genes_smo
)

save(novel_jxns_foundGenes, file="processed-data/07_Jxn_anno/novel_jxns_foundGenes.Rdata")


## Data frame with DE novel jxns in nic and smo, their nearest, following and preceding genes IDs:
##        - DE_Novel_jxn: ID of the Novel DE jxn
##        - nearest_gene: ID of the nearest gene to the jxn
##        - following_gene: ID of the following gene to the jxn
##        - preceding_gene: ID of the preceding gene to the jxn
##        - expt: the experiment (smo/nic) in which the jxn is DE

NovelDEJxns_foundGenesIDs<- data.frame(
  DE_Novel_jxn=c(names(novel_jxns_foundGenes[["DE_Novel_jxns_nic"]]), names(novel_jxns_foundGenes[["DE_Novel_jxns_smo"]])),
  nearest_gene=c(names(novel_jxns_foundGenes[["nearest_genes_nic"]]), names(novel_jxns_foundGenes[["nearest_genes_smo"]])),
  following_gene=c(names(novel_jxns_foundGenes[["following_genes_nic"]]), names(novel_jxns_foundGenes[["following_genes_smo"]])),
  preceding_gene=c(names(novel_jxns_foundGenes[["preceding_genes_nic"]]), names(novel_jxns_foundGenes[["preceding_genes_smo"]])),
  expt=c(rep("Nicotine", length(names(novel_jxns_foundGenes[["DE_Novel_jxns_nic"]]))), rep("Smoking", length(names(novel_jxns_foundGenes[["DE_Novel_jxns_smo"]]))))
)
save(NovelDEJxns_foundGenesIDs, file="processed-data/07_Jxn_anno/NovelDEJxns_foundGenesIDs.Rdata")







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
# date     2023-12-29
# rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
# pandoc   NA
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                * version   date (UTC) lib source
# AcidBase                 0.7.3     2023-12-15 [1] Bioconductor
# AcidCLI                  0.3.0     2023-10-03 [1] Bioconductor
# AcidGenerics             0.7.6     2023-12-15 [1] Bioconductor
# AcidGenomes            * 0.7.2     2023-12-06 [1] Bioconductor
# AcidPlyr                 0.5.3     2023-12-13 [1] local
# AnnotationDbi          * 1.63.2    2023-07-03 [1] Bioconductor
# AnnotationFilter       * 1.26.0    2023-10-26 [1] Bioconductor
# AnnotationHub            3.9.1     2023-06-14 [1] Bioconductor
# Biobase                * 2.61.0    2023-06-02 [1] Bioconductor
# BiocFileCache            2.9.1     2023-07-14 [1] Bioconductor
# BiocGenerics           * 0.48.1    2023-11-02 [1] Bioconductor
# BiocIO                   1.11.0    2023-06-02 [1] Bioconductor
# BiocManager              1.30.21.1 2023-07-18 [1] CRAN (R 4.3.0)
# BiocParallel             1.35.3    2023-07-07 [1] Bioconductor
# BiocVersion              3.18.0    2023-05-11 [1] Bioconductor
# biomaRt                  2.57.1    2023-06-14 [1] Bioconductor
# biomartr               * 1.0.7     2023-12-02 [1] CRAN (R 4.3.1)
# Biostrings               2.69.2    2023-07-05 [1] Bioconductor
# bit                      4.0.5     2022-11-15 [1] CRAN (R 4.3.0)
# bit64                    4.0.5     2020-08-30 [1] CRAN (R 4.3.0)
# bitops                   1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
# blob                     1.2.4     2023-03-17 [1] CRAN (R 4.3.0)
# cachem                   1.0.8     2023-05-01 [1] CRAN (R 4.3.0)
# cli                      3.6.1     2023-03-23 [1] CRAN (R 4.3.0)
# codetools                0.2-19    2023-02-01 [1] CRAN (R 4.3.0)
# colorspace               2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# cowplot                * 1.1.1     2020-12-30 [1] CRAN (R 4.3.0)
# crayon                   1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# curl                     5.0.1     2023-06-07 [1] CRAN (R 4.3.0)
# data.table               1.14.8    2023-02-17 [1] CRAN (R 4.3.0)
# DBI                      1.1.3     2022-06-18 [1] CRAN (R 4.3.0)
# dbplyr                   2.3.3     2023-07-07 [1] CRAN (R 4.3.0)
# DelayedArray             0.26.6    2023-07-02 [1] Bioconductor
# digest                   0.6.33    2023-07-07 [1] CRAN (R 4.3.0)
# dplyr                    1.1.2     2023-04-20 [1] CRAN (R 4.3.0)
# edgeR                  * 3.43.7    2023-06-21 [1] Bioconductor
# ellipsis                 0.3.2     2021-04-29 [1] CRAN (R 4.3.0)
# ensembldb              * 2.26.0    2023-10-26 [1] Bioconductor
# fansi                    1.0.5     2023-10-08 [1] CRAN (R 4.3.1)
# farver                   2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
# fastmap                  1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
# filelock                 1.0.2     2018-10-05 [1] CRAN (R 4.3.0)
# formatR                  1.14      2023-01-17 [1] CRAN (R 4.3.0)
# fs                       1.6.3     2023-07-20 [1] CRAN (R 4.3.0)
# futile.logger          * 1.4.3     2016-07-10 [1] CRAN (R 4.3.0)
# futile.options           1.0.1     2018-04-20 [1] CRAN (R 4.3.0)
# gargle                   1.5.2     2023-07-20 [1] CRAN (R 4.3.0)
# generics                 0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb           * 1.37.2    2023-06-21 [1] Bioconductor
# GenomeInfoDbData         1.2.10    2023-05-28 [1] Bioconductor
# GenomicAlignments        1.37.0    2023-07-07 [1] Bioconductor
# GenomicFeatures        * 1.53.1    2023-06-22 [1] Bioconductor
# GenomicRanges          * 1.54.1    2023-10-30 [1] Bioconductor
# ggplot2                * 3.4.4     2023-10-12 [1] CRAN (R 4.3.1)
# ggrepel                * 0.9.3     2023-02-03 [1] CRAN (R 4.3.0)
# glue                     1.6.2     2022-02-24 [1] CRAN (R 4.3.0)
# goalie                   0.7.7     2023-12-04 [1] Bioconductor
# googledrive              2.1.1     2023-06-11 [1] CRAN (R 4.3.0)
# gridExtra              * 2.3       2017-09-09 [1] CRAN (R 4.3.0)
# gtable                   0.3.4     2023-08-21 [1] CRAN (R 4.3.0)
# here                   * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
# hms                      1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
# htmltools                0.5.5     2023-03-23 [1] CRAN (R 4.3.0)
# httpuv                   1.6.11    2023-05-11 [1] CRAN (R 4.3.0)
# httr                     1.4.6     2023-05-08 [1] CRAN (R 4.3.0)
# interactiveDisplayBase   1.39.0    2023-06-02 [1] Bioconductor
# IRanges                * 2.36.0    2023-10-26 [1] Bioconductor
# jaffelab               * 0.99.32   2023-05-28 [1] Github (LieberInstitute/jaffelab@21e6574)
# jsonlite                 1.8.8     2023-12-04 [1] CRAN (R 4.3.1)
# KEGGREST                 1.41.0    2023-07-07 [1] Bioconductor
# labeling                 0.4.3     2023-08-29 [1] CRAN (R 4.3.0)
# lambda.r                 1.2.4     2019-09-18 [1] CRAN (R 4.3.0)
# later                    1.3.1     2023-05-02 [1] CRAN (R 4.3.0)
# lattice                  0.21-8    2023-04-05 [1] CRAN (R 4.3.0)
# lazyeval                 0.2.2     2019-03-15 [1] CRAN (R 4.3.0)
# lifecycle                1.0.3     2022-10-07 [1] CRAN (R 4.3.0)
# limma                  * 3.57.6    2023-06-21 [1] Bioconductor
# locfit                   1.5-9.8   2023-06-11 [1] CRAN (R 4.3.0)
# magrittr                 2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# MASS                     7.3-60    2023-05-04 [1] CRAN (R 4.3.0)
# Matrix                   1.6-4     2023-11-30 [1] CRAN (R 4.3.1)
# MatrixGenerics         * 1.13.0    2023-05-20 [1] Bioconductor
# matrixStats            * 1.0.0     2023-06-02 [1] CRAN (R 4.3.0)
# memoise                  2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
# mime                     0.12      2021-09-28 [1] CRAN (R 4.3.0)
# munsell                  0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# nlme                     3.1-162   2023-01-31 [1] CRAN (R 4.3.0)
# pillar                   1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pipette                  0.15.2    2023-12-15 [1] Bioconductor
# pkgconfig                2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# png                      0.1-8     2022-11-29 [1] CRAN (R 4.3.0)
# prettyunits              1.1.1     2020-01-24 [1] CRAN (R 4.3.0)
# progress                 1.2.2     2019-05-16 [1] CRAN (R 4.3.0)
# promises                 1.2.0.1   2021-02-11 [1] CRAN (R 4.3.0)
# ProtGenerics             1.34.0    2023-10-26 [1] Bioconductor
# purrr                    1.0.1     2023-01-10 [1] CRAN (R 4.3.0)
# R.methodsS3            * 1.8.2     2022-06-13 [1] CRAN (R 4.3.0)
# R.oo                   * 1.25.0    2022-06-12 [1] CRAN (R 4.3.0)
# R.utils                * 2.12.2    2022-11-11 [1] CRAN (R 4.3.0)
# R6                       2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# rafalib                * 1.0.0     2015-08-09 [1] CRAN (R 4.3.0)
# ragg                     1.2.5     2023-01-12 [1] CRAN (R 4.3.0)
# rappdirs                 0.3.3     2021-01-31 [1] CRAN (R 4.3.0)
# RColorBrewer             1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                     1.0.11    2023-07-06 [1] CRAN (R 4.3.0)
# RCurl                    1.98-1.12 2023-03-27 [1] CRAN (R 4.3.0)
# restfulr                 0.0.15    2022-06-16 [1] CRAN (R 4.3.0)
# rjson                    0.2.21    2022-01-09 [1] CRAN (R 4.3.0)
# rlang                  * 1.1.1     2023-04-28 [1] CRAN (R 4.3.0)
# rprojroot                2.0.3     2022-04-02 [1] CRAN (R 4.3.0)
# Rsamtools                2.17.0    2023-07-07 [1] Bioconductor
# RSQLite                  2.3.1     2023-04-03 [1] CRAN (R 4.3.0)
# rstudioapi               0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
# rtracklayer            * 1.61.0    2023-07-07 [1] Bioconductor
# S4Arrays                 1.1.4     2023-06-02 [1] Bioconductor
# S4Vectors              * 0.40.2    2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
# scales                   1.2.1     2022-08-20 [1] CRAN (R 4.3.0)
# segmented                1.6-4     2023-04-13 [1] CRAN (R 4.3.0)
# sessioninfo            * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
# shiny                    1.7.4.1   2023-07-06 [1] CRAN (R 4.3.0)
# stringi                  1.7.12    2023-01-11 [1] CRAN (R 4.3.0)
# stringr                  1.5.0     2022-12-02 [1] CRAN (R 4.3.0)
# SummarizedExperiment   * 1.30.2    2023-06-06 [1] Bioconductor
# syntactic                0.7.1     2023-10-27 [1] Bioconductor
# systemfonts              1.0.4     2022-02-11 [1] CRAN (R 4.3.0)
# textshaping              0.3.6     2021-10-13 [1] CRAN (R 4.3.0)
# tibble                   3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyselect               1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# utf8                     1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                    0.6.4     2023-10-12 [1] CRAN (R 4.3.1)
# VennDiagram            * 1.7.3     2022-04-12 [1] CRAN (R 4.3.0)
# withr                    2.5.2     2023-10-30 [1] CRAN (R 4.3.1)
# XML                      3.99-0.14 2023-03-19 [1] CRAN (R 4.3.0)
# xml2                     1.3.5     2023-07-06 [1] CRAN (R 4.3.0)
# xtable                   1.8-4     2019-04-21 [1] CRAN (R 4.3.0)
# XVector                  0.41.1    2023-06-02 [1] Bioconductor
# yaml                     2.3.8     2023-12-11 [1] CRAN (R 4.3.1)
# zlibbioc                 1.47.0    2023-05-20 [1] Bioconductor
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
