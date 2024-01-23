
# Script summary

## 01. SPEAQeasy 

Run [_SPEAQeasy_](https://doi.org/10.1186/s12859-021-04142-3) pipeline for quantification of expression features. 

* [01_make_manifest.R](01_SPEAQeasy/01_make_manifest.R): construct a manifest listing sample IDs and associated raw FASTQ files.

* [02_run_pipeline.sh](01_SPEAQeasy/02_run_pipeline.sh): invoke _SPEAQeasy_ to align data, quantify features, and produce [`RangedSummarizedExperiment`](https://www.bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html) objects for downstream analysis.


## 02. Build objects
This initial part of the code builds the necessary objects to analyze in downstream steps.

* [01_add_sequencing_lane.R](02_build_objects/01_add_sequencing_lane.R): Extract flowcell information of the samples.

* [02_build_objects.R](02_build_objects/02_build_objects.R): Explore, clean, correct and format data. Here raw counts are log-normalized, lowly-expressed features filtered out and datasets separated by tissue and age. 


## 03. Exploratory Data Analysis (EDA)

* [01_StudyDesign.R](03_EDA/01_StudyDesign.R): Compute the number of samples belonging to each pair/triad of sample phenotypes.

* [02_QC.R](03_EDA/02_QC.R): Evaluate and compare quality control (QC) metrics of sample groups; explore the relationships between QC metrics, and sample filtering by QC, followed by examination of QC metrics of samples that were kept and removed.

* [03_PCA_MDS.R](03_EDA/03_PCA_MDS.R): Perform Principal Component Analysis (PCA) and Multidimensional Scaling Analysis (MDS) to explore sample-level variability in gene expression and identify drivers of such variation. QC metrics of segregated samples in PCA plots were further analyzed and those samples were removed if turned out to be poor-quality. 

* [04_Expl_Var_partition.R](03_EDA/04_Expl_Var_partition.R): Explore gene-level effects: Explanatory variables and Variance partition.


## 04. Differential Expression Analysis (DEA)
Separated in gene, tx, exon and jx level analyses
* *01_Modeling.R*: Perform DEA 
* *02_Comparisons.R*: Create plots comparing t-stats of different groups of features and Venn diagrams to quantify the number of common DE features between different experiments and/or regulation directions. Analysis of blood vs brain biomarkers and comparisons of human vs mouse genes are within the gene analysis. 

## 05. Functionl Enrichment Analysis (GO & KEGG)
Separated in gene, tx and exon level analyses
* *01_GO_KEGG_analyses.R*: Gene Ontology and KEGG enrichment analyses at the gene level

## 06. DGE visualization
Only at the gene level
* *01_Heatmap_DEG.R*: Create heatmaps of DEG 

## 07. Novel DE junction gene annotation
* *01_Jxn_anno.R*: Obtain novel DE jxns and explore potential novel isoforms



# Code organization

Files are organized following the structure from [LieberInstitute/template_project](https://github.com/LieberInstitute/template_project). Scripts include the  R session information with details about version numbers of the packages we used.
