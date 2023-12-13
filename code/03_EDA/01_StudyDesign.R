
## 1. Study design


library(sessioninfo)
load(here("processed-data/02_build_objects/pheno.tsv"))


## 1.1 Number of samples for each pair of phenotypes

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



## 1.2 Number of samples for each triad of phenotypes

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

