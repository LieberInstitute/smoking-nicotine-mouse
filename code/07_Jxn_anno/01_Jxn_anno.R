
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

## "ExonSkip" jxns have sites from non-successive exons, both known
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




## 1.1 Extract exons and genes from novel jxns??

## Novel DE introns/jxns
## Nicotine
novel_jxns_nic <- de_jxns_nic[which(! (de_jxns_nic$Class=="InGen" | de_jxns_nic$Class=="ExonSkip")),]
## Smoking
novel_jxns_smo <- de_jxns_smo[which(! (de_jxns_smo$Class=="InGen" | de_jxns_smo$Class=="ExonSkip")),]





