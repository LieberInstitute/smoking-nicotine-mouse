
## 1. Junction annotation 

## Only for DE jxns from pup and brain samples

load(here("raw-data/rse_jx_smoking_mouse_n208.Rdata"))
load(here("processed-data/04_DEA/Jx_analysis/top_jxns_nic.Rdata"))
load(here("processed-data/04_DEA/Jx_analysis/de_jxns_nic.Rdata"))
load(here("processed-data/04_DEA/Jx_analysis/top_jxns_smo.Rdata"))
load(here("processed-data/04_DEA/Jx_analysis/de_jxns_smo.Rdata"))



## Explore jxns classes

## "Novel" jxns have unknown (not in Gencode) start and end sites
table(rowData(rse_jx)[which(rowData(rse_jx)$Class=="Novel"), c("inGencodeStart", "inGencodeEnd")])
#                 inGencodeEnd
#   inGencodeStart FALSE
#            FALSE 847606

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

## "InGen" jxns are in Gencode (both sites are known)
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


## Create table with the information of DE jxns' genes:
##        - Gene: ID of gene with DE jxns
##        - number_DEjxns_nic: number of nicotine DE jxns the gene has
##        - number_DEjxns_nic_Novel: number of nicotine DE novel jxns the gene has
##        - number_DEjxns_nic_AltStartEnd: number of nicotine DE jxns with alternative start/end, the gene has
##        - number_DEjxns_nic_InGen: number of nicotine DE known (in Gencode) jxns the gene has
##        - number_DEjxns_nic_ExonSkip: number of nicotine DE jxns from non-successive exons, the gene has
##        - number_DEjxns_smo: number of smoking DE jxns the gene has
##        - number_DEjxns_smo_Novel: number of smoking DE novel jxns the gene has
##        - number_DEjxns_smo_AltStartEnd: number of smoking DE jxns with alternative start/end, the gene has
##        - number_DEjxns_smo_InGen: number of smoking DE known (in Gencode) jxns the gene has
##        - number_DEjxns_smo_ExonSkip: number of smoking DE jxns from non-successive exons, the gene has

## Genes with DE jxns
DEjxns_genes <- union(unique(de_jxns_nic$newGeneID), unique(de_jxns_smo$newGeneID))
DEjxns_genes <- DEjxns_genes[which(! is.na(DEjxns_genes))]
## Info
DEjxns_genes_info <- vector()
for (gene in DEjxns_genes){
          DEjxns_genes_info <- rbind(DEjxns_genes_info, 
                                    c(Gene=gene,
                                      number_DEjxns_nic=length(which(de_jxns_nic$newGeneID==gene)),
                                      number_DEjxns_nic_Novel=length(which(de_jxns_nic$newGeneID==gene & de_jxns_nic$Class=="Novel")),
                                      number_DEjxns_nic_AltStartEnd=length(which(de_jxns_nic$newGeneID==gene & de_jxns_nic$Class=="AltStartEnd")),
                                      number_DEjxns_nic_InGen=length(which(de_jxns_nic$newGeneID==gene & de_jxns_nic$Class=="InGen")),
                                      number_DEjxns_nic_ExonSkip=length(which(de_jxns_nic$newGeneID==gene & de_jxns_nic$Class=="ExonSkip")),
                                      number_DEjxns_smo=length(which(de_jxns_smo$newGeneID==gene)),
                                      number_DEjxns_smo_Novel=length(which(de_jxns_smo$newGeneID==gene & de_jxns_smo$Class=="Novel")),
                                      number_DEjxns_smo_AltStartEnd=length(which(de_jxns_smo$newGeneID==gene & de_jxns_smo$Class=="AltStartEnd")),
                                      number_DEjxns_smo_InGen=length(which(de_jxns_smo$newGeneID==gene & de_jxns_smo$Class=="InGen")),
                                      number_DEjxns_smo_ExonSkip=length(which(de_jxns_smo$newGeneID==gene & de_jxns_smo$Class=="ExonSkip"))))
}
## Char to numeric
DEjxns_genes_info <- cbind(Gene=DEjxns_genes_info[,1], as.data.frame(apply(DEjxns_genes_info[,2:11], 2, function(x){as.numeric(as.character(x))})))
save(DEjxns_genes_info, file="processed-data/07_Jxn_anno/DEjxns_genes_info.Rdata")


## Boxplot of the number of DE jxns per gene
pdf(file = paste("plots/07_Jxn_anno/Number_DEjxns_per_gene.pdf", sep="" ))
par(mfrow=c(2,2))

ggplot(data=DEjxns_genes_info, aes(x=0,y=number_DEjxns_nic)) + 
      geom_boxplot() +
      theme_classic() +
      labs(x = 0, y = "Number of nicotine DE jxns per gene") 

ggplot(data=DEjxns_genes_info, aes(x=0,y=number_DEjxns_smo)) + 
      geom_boxplot() +
      theme_classic() +
      labs(x = 0, y = "Number of smoking DE jxns per gene") 

dev.off()




## Obtain novel DE introns/jxns (with at least one unknown end or if the combination of ends is unknown)

## Nicotine
novel_jxns_nic <- de_jxns_nic[which(! (de_jxns_nic$Class=="InGen")),]
## Smoking
novel_jxns_smo <- de_jxns_smo[which(! (de_jxns_smo$Class=="InGen")),]

## Find nearest genes of novel jxns




