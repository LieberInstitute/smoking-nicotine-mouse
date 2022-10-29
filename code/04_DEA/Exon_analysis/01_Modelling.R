
library(here)


## 1. Differential Expression Analysis at the exon level



## Only for brain and pups samples

load(here("processed-data/03_EDA/02_QC/rse_exon_brain_pups_qc.Rdata"))



## 1.1 Modelling

## Extract previous output from calcNormFactors for all samples
norm_factors<-calcNormFactors(rse_exon_brain_pups_qc, method = "TMM")
samples_factors<-data.frame(SAMPLE_ID=norm_factors$samples$SAMPLE_ID,
                            norm.factors=norm_factors$samples$norm.factors,
                            lib.size=norm_factors$samples$lib.size)

