## 1. Study design
library(sessioninfo)

## 1.1 Load data
load(here("processed-data/02_build_objects/pheno.tsv"))

## 1.2 Number of samples for each pair of phenotypes
samples_pheno <- function(pheno_var1, pheno_var2){
    ## Pheno variables levels
    levels_var1<-sort(unique(pheno[,pheno_var1]))
    levels_var2<-sort(unique(pheno[,pheno_var2]))
    for (level1 in levels_var1){
        for (level2 in levels_var2){
            print(paste("Number of samples of (", pheno_var1, ") ", level1, " and (", pheno_var2, ") ", level2, ": ", 
                        length(which(pheno[,pheno_var1]==level1 & pheno[,pheno_var2]==level2)), sep=""))
        }
    }
}

## All pairs of phenotypes
vars1<-vector()
`%nin%` = Negate(`%in%`)
for (pheno_var1 in c("Age", "Tissue", "plate", "medium","Expt", "Sex", "Group", "Pregnancy")){
    for (pheno_var2 in c("Age", "Tissue", "plate", "medium","Expt", "Sex", "Group", "Pregnancy")){
        if (pheno_var1!=pheno_var2 & pheno_var2 %nin% vars1){
            vars1<-append(vars1, pheno_var1)
            samples_pheno(pheno_var1, pheno_var2)}
    }
} 

## 1.3 Number of samples for each triad of phenotypes
three_pheno <- function(pheno_var1, pheno_var2, pheno_var3){
  ## Pheno variables levels
    levels_var1<-sort(unique(pheno[,pheno_var1]))
    levels_var2<-sort(unique(pheno[,pheno_var2]))
    levels_var3<-sort(unique(pheno[,pheno_var3]))
    for (level1 in levels_var1){
      for (level2 in levels_var2){
        for (level3 in levels_var3){
          print(paste("Number of samples of (", pheno_var1, ") ", level1, ", (", pheno_var2, ") ", level2, 
                " and (", pheno_var3, ") ", level3, ": ", length(which(pheno[,pheno_var1]==level1 &
                pheno[,pheno_var2]==level2 & pheno[,pheno_var3]==level3)), sep=""))
        }
      }
    }
}

## Number of samples 
three_pheno("Pregnancy", "Age", "Tissue")
three_pheno("Pregnancy", "Age", "Group")
three_pheno("Tissue", "Expt", "Group")
three_pheno("Pregnancy", "Expt", "Group")



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

#  version  R version 4.2.0 (2022-04-22 ucrt)
#  os       Windows 10 x64 (build 19044)
#  system   x86_64, mingw32
#  ui       RStudio
#  language (EN)
#  collate  Spanish_Mexico.utf8
#  ctype    Spanish_Mexico.utf8
#  tz       America/Mexico_City
